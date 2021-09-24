# mb6303
The development and implementation of our CAGE-seq pipeline was done in R version 4.1., bioconductor(latest) and  using ORFik as one of the base packages.
To run CAGE-pipeline, run below command in terminal.

1. To setup  folder structure, use below mentioned scripts with the arguments in same sequence.
Rscript construct_folder_struct.R <config file name> <Experiemnt name> <Assembly name>

2. To move raw files from local directory or to download from online repo, use below mentioned script.
Rscript get_fastq.R <local/Online> <"Path to local file"/"Experiment Accession(primary or secondary)">

3. To map raw reads to a genome, use below mwntioned script.
Rscript mapping.R <local/Online>

4. To reassign TSS of each transcript and generate a simple dataset:
Rscript Reassign_TSS.R <Yes/No> <"Path to save transcript data set">

