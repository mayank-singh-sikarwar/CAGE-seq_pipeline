###############################################################################
#creation date: 6-July-2021
#Description: script creates folder structure directory and stores it in a con
#-fig file. It also configures the experiment in base folder structure.
#Argument; <config file name> <Experiemnt name> <Assembly name>
###############################################################################

#Argument
arg <- commandArgs(trailingOnly = TRUE)
def <- c("config_file","Experiment","Assembly")

#default Argument
for( i in 1:length(arg)){
  if(nchar(arg[i]) == 0){
    arg[i] = def[i]
  }
}
configure_base <- function(arg){
  #folder structure setup script
  suppressMessages(library(ORFik))
  #set working directory
  for_conf <- file.path
  #store current working directory
  parent_folder <- file.path(getwd())
  cat("parent_folder:",parent_folder)
  cat("\n")
  #raw data directory 
  fastq.dir<- file.path(parent_folder, "base/raw_data")
  
  #processed(trimmed and aligned) data directory 
  bam.dir <- file.path(parent_folder, "base/processed_data")
  
  #reference genome directory
  reference.dir <- file.path(parent_folder, "base/reference_genome")
  
  #config file name prompt
  config_name <- arg[1]
  
  config.filepath <- file.path(parent_folder, config_name)
  message(config.filepath)
  #create config file to save folder structure path
  config.save(config.filepath, 
              fastq.dir, 
              bam.dir, 
              reference.dir
  )
  #check if config file is created
  config_exist <- file.exists(config.filepath)
  
  if(config_exist == TRUE){
    message("config creation status: ",config_exist)
    experiment_name <- arg[2]
    assembly_name <- arg[3]

    conf <- config.exper(experiment = experiment_name,
                         assembly = assembly_name,
                         type = c("RNA-seq", "Ribo-seq", "CAGE"),
                         config = config(config.filepath))
    
  }else{
    message("config creation status: ",config_exist," Rexecute script giving correct config name.") 
  }
  dir.create("configur_details")
  wd <- file.path(parent_folder,"configur_details/conf.csv")
  write.csv(conf, wd)
}
#calling function to return dataframe
configure_base(arg)

