###############################################################################
#creation date: 6-July-2021
#Description: This function generates Txdb object using either local genome and 
#annotation file or download from a reputed repository like ensembl.
#Argument: <local/online>
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
Genome_annotation<- function(arg){
  organism <- "Homo Sapiens"
#control flow for the use of local file or to retrieve it from 
#online open repo.
  if(arg[1] == "local"){
    annotation <- getGenomeAndAnnotation(organism,
                                         output.dir = conf["ref"],
                                         db = "ensembl",
                                         GTF = FALSE,
                                         genome = FALSE)
    annotation["genome"] <- "/home/mayank/ref/GRCh38.primary_assembly.genome.fa"
    annotation["gtf"] <- "/home/mayank/ref/gencode.v35.annotation.gtf"
    
  }else{
    annotation <- getGenomeAndAnnotation(organism,
                                         output.dir = conf["ref"],
                                         db = "ensembl",
                                         GTF = TRUE,
                                         genome = TRUE,
                                         merge_contaminants = FALSE,
                                         phix = FALSE,ncRNA = FALSE,
                                         tRNA = FALSE,rRNA = FALSE,
                                         gunzip = TRUE,remake = FALSE,
                                         assembly_type = "primary_assembly"
    )
  }
  return(annotation)
}

annotation <- Genome_annotation(arg)

#index the genome
index <- STAR.index(annotation,
                    star.path = "/mnt/data/miniconda_python2/miniconda2/envs/STAR_2.7.7a_env/bin/STAR",
                    max.cpus = 6,
                    max.ram = 30,
                    #SAsparse = ,
                    wait = TRUE,
                    remake = FALSE,
                    #script = 
                      )

#aligning the data
alignment <- STAR.align.folder(conf["fastq CAGE"],
                               conf["bam CAGE"],
                               index,
                               star.path = "/mnt/data/miniconda_python2/miniconda2/envs/STAR_2.7.7a_env/bin/STAR",
                               fastp = "/home/mayank/scripts/fastp",
                               paired.end = FALSE,
                               steps = "tr-ge",
                               adapter.sequence = c("auto"),
                               #min.length = ,
                               mismatches = 2,
                               #trim.front =,
                               max.multimap = 10,
                               alignment.type = "EndToEnd",
                               #EndToEnd#Local,
                               max.cpus = 4,
                               #wait=,
                               #include.subfolders=,
                               #resume=,
                               #multiQC=,
                               #script.folder =,
                               #script.single =,
                               )


