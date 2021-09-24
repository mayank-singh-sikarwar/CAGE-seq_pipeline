###############################################################################
#creation date: 31-August-2021
#Description: This workflow consumes cage fastq files and delineates it into a
#dataframe.
#Argument: 
###############################################################################
arg <- commandArgs(trailingOnly = TRUE)
#default Argument
for( i in 1:length(arg)){
  if(length(arg) == 0 ){
    stop("No argument provided, refer Readme regarding argument and its format.")
  }
}
suppressMessages(library(ORFik))
parent_folder <- file.path(getwd())
wd <- file.path(parent_folder,"configur_details/conf.csv")
conf<- read.csv(wd )
conf<- setNames(as.character(conf$x), conf$X)
organism <- "Homo Sapiens"
annotation <- getGenomeAndAnnotation(organism,
                                     output.dir = conf["ref"],
                                     db = "ensembl",
                                     GTF = FALSE,
                                     genome = FALSE)
annotation["genome"] <- "/home/mayank/ref/GRCh38.primary_assembly.genome.fa"
annotation["gtf"] <- "/home/mayank/ref/gencode.v35.annotation.gtf"


gtf<- annotation["gtf"] 
genome <- annotation["genome"]
message("GTF file and Genome file path instantiated")
bam_dir <- file.path(conf["bam CAGE"],"aligned/")
bam_file_name <- list.files(bam_dir, pattern = ".bam")
processed_bam <- file.path(bam_dir, bam_file_name)

#directory instances
ofst_name <- unlist(strsplit(bam_file_name,"[.]"))[1]
ofst_directory <- paste0(file.path(bam_dir,ofst_name),".ofst")
message("ofst file location: ",ofst_directory)

#check if ofst exist; if TRUE load ofst; If False read bam file,
#export ofst and import ofst file.
if(file.exists(ofst_directory)){
  cage_ofst <- import.ofst(ofst_directory)
}else{
  cage_bam <- readBam(processed_bam)
  export.ofst(cage_bam ,ofst_directory)
  cage_ofst <- import.ofst(ofst_directory)
}
message("ofst file loaded.")
#load gtf as TXdb object
gtf_db_as_txdb <- loadTxdb(gtf , chrStyle = "UCSC")
#load 5' prime leaders from annotation db
annotated_5utrs <- loadRegion(gtf_db_as_txdb, "leader", by = "tx")
#instantiating grangeslist object to retrive gene id and gene name
suppressMessages(library(rtracklayer))
gene_annotations<- import(gtf)
gene_info_df <-mcols(gene_annotations)[,c("gene_id","gene_name","transcript_id")]

###############################################################################
#creation date: 5-AUgust-2021
#Description: Function returns grangelist object of transcripts with reassigned 
#TSS and a holistic plot of all the TSS reassigned.
###############################################################################
reassignTSS_wrapper <- function(gtf_db_as_txdb , cage_ofst,annotated_5utrs){
  
 if(arg[1]=="yes"){
   ####################################################################################################################
   #creating experiment
   message("Creating experiment ",conf["exp CAGE"])
   create.experiment(
     dir = bam_dir,
     exper = conf["exp CAGE"],
     saveDir = file.path(parent_folder,conf["exp CAGE"]),
     txdb = gtf,
     fa = genome,
     organism = organism,
     pairedEndBam = FALSE,
     viewTemplate = FALSE,
     types = c("bam", "bed", "wig"),
     libtype = "CAGE",
     stage = "Prefrontal_cortex",
     rep = "1",
     condition = "Normal",
     fraction = "auto"
   )
   exp.df <- read.experiment(conf["exp CAGE"],file.path(parent_folder,conf["exp CAGE"]))
   #TSS heatmap
   heatMapRegion(
     exp.df,
     region = "TSS",
     outdir = file.path(conf["bam CAGE"],"Plots"),
     scores = c( "transcriptNormalized", "sum"),
     type = "ofst",
     cage = cage_ofst,
     #   plot.ext = c(".png"),
     acceptedLengths = 21:75,
     upstream = 500,
     downstream = 500,
     shifting = c("5prime"),
     #   longestPerGene = FALSE
   )
   #####################################################################################################################
   
   
   #load 5' prime leaders from annotation db and use them to reassignTSS
   median_read_per_leader <- median(countOverlaps(annotated_5utrs,cage_ofst))
   message("Reassiging TSS of the leaders using cage data")
   cage_reannotated_5utrs <- reassignTSSbyCage(annotated_5utrs, 
                                               cage_ofst,
                                               extension = 1000,#no. of bps upstream from TSS 
                                               #to look for cage peaks
                                               filterValue = median_read_per_leader, 
                                               restrictUpstreamToTx = FALSE,#restrict upstream TSS 
                                               #overlapping of one transcript to another
                                               removeUnused = TRUE,#remove leader with no cage support
                                               preCleanup = TRUE,#input CAGE reads, not on the leaders
                                               #avoid changing leaders, when the change 
                                               #is less than 5 nt away from original
                                               cageMcol = TRUE
   )
 }else if(arg[1] == "no"){
   #load 5' prime leaders from annotation db and use them to reassignTSS
   median_read_per_leader <- median(countOverlaps(annotated_5utrs,cage_ofst))
   message("Reassiging TSS of the leaders using cage data")
   cage_reannotated_5utrs <- reassignTSSbyCage(annotated_5utrs, 
                                               cage_ofst,
                                               extension = 1000,#no. of bps upstream from TSS 
                                               #to look for cage peaks
                                               filterValue = median_read_per_leader, 
                                               restrictUpstreamToTx = FALSE,#restrict upstream TSS 
                                               #overlapping of one transcript to another
                                               removeUnused = TRUE,#remove leader with no cage support
                                               preCleanup = TRUE,#input CAGE reads, not on the leaders
                                               #avoid changing leaders, when the change 
                                               #is less than 5 nt away from original
                                               cageMcol = TRUE)
   
 }
  
  
  return(cage_reannotated_5utrs)
}
#calling function 
cage_reannotated_5utrs <- reassignTSS_wrapper(gtf_db_as_txdb , cage_ofst,annotated_5utrs)

###############################################################################
#creation date: 28-August-2021
#Description: Function returns stable id for a gene and gene acronym using a tr
#anscript stable id
###############################################################################
get_gene_info <- function(gene_info_df, transcript_id){
  tx_from_gtf <- gene_info_df$transcript_id
  row_ind <- match(transcript_id, tx_from_gtf)
  gene_id <- gene_info_df[row_ind,]$gene_id
  gene_name <- gene_info_df[row_ind,]$gene_name
  
  gene_info <- list(gene_id, gene_name)
  
  return(gene_info)
}

#################################################################################
#creation date: 30-August-2021
#Description: Function returns indexs of exons with changed tss in a transcript 
#list object 
#################################################################################
index_by_exon <-function(ann_start_or_end,cage_start_or_end){
  
  index <- which(!(cage_start_or_end %in% ann_start_or_end))
  
  if(length(index)==0){
    index <- as.numeric(1)
  }else{
    index <- index
  }
  return(index)
}
#######################################################################################
#creation date: 30-August-2021
#Description: Function returns start site and end sites from cage and annoateted  
#file while handling exon mismatch case observed for transcript like:"ENST00000371714.5"
#######################################################################################
index_pos_strand <- function(annotated_gr,cage_gr){
  #storing annotated & cage exon names
  ann_exon_name <- mcols(annotated_gr)$exon_name
  cage_exon_name <- mcols(cage_gr)$exon_name
  if(length(ann_exon_name) == length(cage_exon_name)){
    #storing annotated & cage start site
    ann_start <- start(annotated_gr)
    cage_start <- start(cage_gr)
    index <- index_by_exon(ann_start,cage_start)
    index_vec <- c(index, index)
  }else if(length(ann_exon_name) != length(cage_exon_name)){
    ann_exon_index <- match(cage_exon_name, ann_exon_name)
    ann_exon_index <- ann_exon_index[!is.na(ann_exon_index)]
    ann_exon_index <- ann_exon_index[1]
    #         print(ann_exon_index)
    cage_exon_index <- match(ann_exon_name, cage_exon_name)
    cage_exon_index <- cage_exon_index[!is.na(cage_exon_index)]
    cage_exon_index <- cage_exon_index[1]
    
    index_vec <- c(ann_exon_index, cage_exon_index)
    
  }
  return(index_vec)
}
#######################################################################################
#creation date: 30-August-2021
#Description: Function returns start site and end sites from cage and annoateted  
#file while handling exon mismatch case observed for transcript like:"ENST00000371714.5"
#######################################################################################
index_neg_strand <- function(annotated_gr,cage_gr){
  #storing annotated & cage exon names
  ann_exon_name <- mcols(annotated_gr)$exon_name
  cage_exon_name <- mcols(cage_gr)$exon_name
  if(length(ann_exon_name) == length(cage_exon_name)){
    #storing annotated & cage start site
    ann_end <- end(annotated_gr)
    cage_end <- end(cage_gr)
    index <- index_by_exon(ann_end,cage_end)
    index_vec <- c(index, index)
  }else if(length(ann_exon_name) != length(cage_exon_name)){
    ann_exon_index <- match(cage_exon_name, ann_exon_name)
    ann_exon_index <- ann_exon_index[!is.na(ann_exon_index)]
    ann_exon_index <- ann_exon_index[1]
    #         print(ann_exon_index)
    cage_exon_index <- match(ann_exon_name, cage_exon_name)
    cage_exon_index <- cage_exon_index[!is.na(cage_exon_index)]
    cage_exon_index <- cage_exon_index[1]
    
    index_vec <- c(ann_exon_index, cage_exon_index)
    
  }
  return(index_vec)
}
#################################################################################
#creation date: 30-August-2021
#Description: Function returns returns reads count for a tranbscript by checking
#index in a vector containing read counts for all the transcripts present in data
#################################################################################
read_counts <- function(transcript_id, read_count, diff_index){
  transcript_name_fullset <- names(read_count)
  name_match_index <- match(transcript_id, transcript_name_fullset)
  if(diff_index > 1){
    read_index = (name_match_index + diff_index - 1)
    read_count_numeric = as.numeric(read_count[read_index])
  }else{
    read_count_numeric = as.numeric(read_count[name_match_index])
  }
  return(read_count_numeric)
}
###############################################################################
#creation date: 4-August-2021
#Description: Function returns shift flag for changed TSS, on the basis of diff
#erence observed between TSS_by_cage and Original_TSS
###############################################################################
shift_flag <- function(diff, strand){
  if(runValue(strand) =="+"){
    if(diff < 0){
      flag <- "upstream"
    }else if (diff > 0){
      flag <- "downstream"
    }else {
      flag <- "No change"
    }
  }else if(runValue(strand) == "-"){
    if(diff < 0){
      flag <- "downstream"
    }else if (diff > 0){
      flag <- "upstream"
    }else {
      flag <- "No change"
    }
  }
  return(flag)
}
###############################################################################
#creation date: 4-August-2021
#Description: Function returns shift flag for changed TSS, on the basis of diff
#erence observed between TSS_by_cage and Original_TSS
###############################################################################
getTSSchangeinDF <- function(cage_reannotated_5utrs,gtf_db_as_txdb, gene_info_df){
  #selecting 5' UTRS with changed TSS from original set of 5'UTRs
  newleader_names <- names(cage_reannotated_5utrs)
  message("Number of Transcript with reassigned TSS:",length(newleader_names))
  selected_original_5UTRs <- loadRegion(gtf_db_as_txdb, "leader",
                                        names.keep = newleader_names ,by = "tx")
  message("Number of Transcripts selected on the basis of reassigned TSS from annotated Transcripts set: "
          ,length(selected_original_5UTRs))
  
  #sorting GRangeslist objects according to genomic ranges
  cage_reannotated_5utrs <- sort(cage_reannotated_5utrs) 
  selected_original_5UTRs <- sort(selected_original_5UTRs)
  
  
  #storing all genomic overlaps between cage Galigment object and granges object
  reads_per_cage_reannotatedTSS_start <- countOverlaps(unlistGrl(cage_reannotated_5utrs), 
                                                       cage_ofst, type=c("start"))
  reads_per_annotatedTSS_start <- countOverlaps(unlistGrl(selected_original_5UTRs), cage_ofst
                                                , type=c("start"))
  message("Overlapping start sites stored of annotated and reassigned TSS on positive strand.")
  reads_per_cage_reannotatedTSS_end <- countOverlaps(unlistGrl(cage_reannotated_5utrs), 
                                                     cage_ofst, type=c("end"))
  reads_per_annotatedTSS_end <- countOverlaps(unlistGrl(selected_original_5UTRs), cage_ofst
                                              , type=c("end"))
  message("Overlapping start sites stored of annotated and reassigned TSS on negative strand.")
  
  #creating empty df to store values
  df <- data.frame(Seqname=factor(),
                   Strand = factor(),
                   Gene_id = character(),
                   Gene_acronym = character(),
                   Transcript_name=character(),
                   Exon_name=character(),
                   Original_TSS=integer(),
                   Original_TSS_read_count = integer(),
                   TSS_by_cage=integer(),
                   Reannotated_TSS_read_count = integer(),
                   Nucleaotide_Difference=integer(),
                   Flag = character(),
                   stringsAsFactors=FALSE)
  total_transcripts <- length(selected_original_5UTRs)
  pb <- txtProgressBar(min = 0, max = total_transcripts, style = 3)
  for(i in 1:total_transcripts){
    #         print(i)
    # update progress bar
    setTxtProgressBar(pb, i)
    
    #taking one transcript and unlisting it
    annotated <- selected_original_5UTRs[i]
    annotated_gr <- unlistGrl(annotated)
    cage <- cage_reannotated_5utrs[i]
    cage_gr <- unlistGrl(cage)
    
    if(runValue(strand(annotated_gr)) == "+"){
      
      index_list  <- index_pos_strand(annotated_gr,cage_gr)
      #             print(index_list)
      ann_site_index <- index_list[1]
      cage_site_index <- index_list[2]
      annotated_gr_single <- annotated_gr[ann_site_index, "exon_name"]
      cage_gr_single <- cage_gr[cage_site_index, "exon_name"]
      Seqname <- seqnames(annotated_gr_single)
      Strand <- strand(annotated_gr_single)
      transcript_id <- names(annotated_gr_single)
      #             print(transcript_id)
      gene_info <- get_gene_info(gene_info_df, transcript_id)
      gene_id <- as.character(gene_info[1])
      gene_name <- as.character(gene_info[2])
      Exon_name <- mcols(annotated_gr_single)$exon_name
      Original_TSS <- start(annotated_gr_single)
      Original_TSS_read_count <- read_counts(transcript_id, 
                                             reads_per_annotatedTSS_start, 
                                             ann_site_index)
      TSS_by_cage <- start(cage_gr_single)
      Reannotated_TSS_read_count <- read_counts(transcript_id, 
                                                reads_per_cage_reannotatedTSS_start, 
                                                cage_site_index)
      Diff <- (TSS_by_cage - Original_TSS)
      Flag <- shift_flag(Diff, Strand)
      
      de <- list(Seqname = Seqname,
                 Strand = Strand,
                 Gene_id = gene_id,
                 Gene_acronym = gene_name,
                 Transcript_name = transcript_id,
                 Exon_name = Exon_name,
                 Original_TSS = Original_TSS,
                 Original_TSS_read_count = Original_TSS_read_count,
                 TSS_by_cage = TSS_by_cage,
                 Reannotated_TSS_read_count = Reannotated_TSS_read_count,
                 Nucleaotide_Difference = Diff,
                 Flag = Flag
      )
      df <- rbind(df, de, stringsAsFactors=FALSE)
      
      
    }else if(runValue(strand(annotated_gr)) == "-"){
      index_list <- index_neg_strand(annotated_gr,cage_gr)
      #             print(index_list)
      ann_site_index <- index_list[1]
      cage_site_index <- index_list[2]
      annotated_gr_single <- annotated_gr[ann_site_index, "exon_name"]
      cage_gr_single <- cage_gr[cage_site_index, "exon_name"]
      Seqname <- seqnames(annotated_gr_single)
      Strand <- strand(annotated_gr_single)
      transcript_id <- names(annotated_gr_single)
      gene_info <- get_gene_info(gene_info_df, transcript_id)
      gene_id <- as.character(gene_info[1])
      gene_name <- as.character(gene_info[2])
      Exon_name <- mcols(annotated_gr_single)$exon_name
      Original_TSS <- end(annotated_gr_single)
      Original_TSS_read_count <- read_counts(transcript_id, 
                                             reads_per_annotatedTSS_end, 
                                             ann_site_index)
      TSS_by_cage <- end(cage_gr_single)
      Reannotated_TSS_read_count <- read_counts(transcript_id, 
                                                reads_per_cage_reannotatedTSS_end, 
                                                cage_site_index)
      Diff <- (TSS_by_cage - Original_TSS)
      Flag <- shift_flag(Diff, Strand) 
      
      de <- list(Seqname = Seqname,
                 Strand = Strand,
                 Gene_id = gene_id,
                 Gene_acronym = gene_name,
                 Transcript_name = transcript_id,
                 Exon_name = Exon_name,
                 Original_TSS = Original_TSS,
                 Original_TSS_read_count = Original_TSS_read_count,
                 TSS_by_cage = TSS_by_cage,
                 Reannotated_TSS_read_count = Reannotated_TSS_read_count,
                 Nucleaotide_Difference = Diff,
                 Flag = Flag
      )
      df <- rbind(df, de, stringsAsFactors=FALSE)
      
    }
    
  }
  close(pb)
  return(df)
}
#calling getTSSchangeinDF to write dataframe containing all the transcripts
#with reassigned TSS, its supporting cage read count and annotated TSS, its
#supporting cage reads and all related information.
DF <- getTSSchangeinDF(cage_reannotated_5utrs,gtf_db_as_txdb,gene_info_df)
write.csv(DF,file.path(arg[2],"CAGE_data.csv"))





