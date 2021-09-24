###############################################################################
#creation date: 6-July-2021
#Description: this script either transfers local fastq files from single folder
#or download it open repo. using primary or secondary experiment accession.
#Argument: <local/Online> <"Path to local file"/"Experiment Accession(primary or 
#secondary)">
###############################################################################
arg <- commandArgs(trailingOnly = TRUE)
#default Argument
for( i in 1:length(arg)){
  if(length(arg) == 0 | nchar(arg[i]) == 0  ){
    stop("No argument provided, Refer Readme regarding argument and its format.")
  }
}
#reading configuration paths and creating a indexed character vector
parent_folder <- file.path(getwd())
wd <- file.path(parent_folder,"configur_details/conf.csv")
conf<- read.csv(wd )
conf<- setNames(as.character(conf$x), conf$X)
print(conf)
suppressMessages(library(ORFik))
#control flow for the use of local file or to retrieve it from 
#online open repo.

get_fastq <-function(arg){
  if(arg[1] == "local"){
    
    message("Transfer local file from specified path")
    message("-----------------------------------------------------")
    
    from_file_path = arg[2]
    to_file_path <- conf["fastq CAGE"]
    to_path_stat <- dir.exists(conf["fastq CAGE"])
    if(to_path_stat == FALSE){
      dir.create(to_file_path, recursive = TRUE)
    }else{
      message(to_path_stat,": Destination folders already exist")
    }
    message("-----------------------------------------------------")
    if(dir.exists(from_file_path)){
      file_list <- list.files(from_file_path, pattern = "\\.fastq.gz$",full.names = TRUE)
    }else{
      file_list <- from_file_path
    }
    
    copy_stat <-file.copy(file_list,to_file_path, overwrite = TRUE)
    false_num <- table(copy_stat)["FALSE"]
    #check if all files are copied
    if(is.na(false_num)){
      message("All files copied successfully.")
    }else{
      message("Few file transfer failed, need to check manually")
    }
    #downloading fastq files using 
  }else if (arg[1] == "online"){
    message("Using download.sra command")
    message("-----------------------------------------------------")
    #requesting input for study accession 
    accession = as.character(arg[2])
    message("Run Available accession")
    exp_meta <- download.SRA.metadata(accession, outdir = conf["fastq CAGE"])
    exp_meta <- exp_meta[,1]
    download.SRA(exp_meta, 
                 conf["fastq CAGE"],
                 fastq.dump.path = "/mnt/data/miniconda_python2/miniconda2/envs/SRA_tools_env/bin/fastq-dump",
                 rename = FALSE)
  }
}
get_fastq(arg)