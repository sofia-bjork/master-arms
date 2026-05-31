#!/usr/bin/env Rscript

# List of packages you will need
required_packages <- c('argparse')

# Determine already installed packages
installed_packages <- installed.packages()[, 'Package']

# Loop through required packages
for (package in required_packages) {
  
  if (!(package %in% installed_packages)) {  # If package not installed...
    options(repos = "https://cran.rstudio.com/")  # ...set the CRAN mirror...
    install.packages(package)                     # ...and download the package.
  }
  
  suppressPackageStartupMessages(library(package, character.only = TRUE))  # Load package silently.
}

# create the output directory
if(!dir.exists("output")) dir.create("output")

# create the output output/miseq_fastq directory
if(!dir.exists("output/miseq_fastq")) dir.create("output/miseq_fastq")

########################################
## PARSING THE COMMAND LINE ARGUMENTS ##
########################################

# Initialize command line argument parser
parser <- ArgumentParser(description = 'DOWNLOAD ENA ACCESSIONS')

# All of your command line arguments
parser$add_argument('-f', '--file', metavar = 'fileName', type = 'character', required = TRUE, help = 'Specify the file that contains the ENA accessions.')
parser$add_argument('-d', '--directory', metavar = 'directory', type = 'character', required = TRUE, help = 'Specify the folder that should be downloaded into.')

# Parse the arguments
args <- parser$parse_args()

# Access the arguments
ENAFile <- args$file
directory <- args$directory


##########################
## FILE EXISTENCE CHECK ##
##########################

# Check if your file exists
if(!file.exists(ENAFile)){ # If the file does not exist...
  stop(paste("Error:", ENAFile, "not found.")) # ... stop the script.
}


############################
## FILE CONTENTS HANDLING ##
############################

# Read the contents of the ENA file
lines <- readLines(ENAFile, warn = F)

# Make an empty variable that, later on, will contain information about the sequencing runs and their samples
runs_samples <- list()

# Make an empty variable that, later on, will help us remember what sequencing run we are currently handling
current_run <- NULL

# Loop over all of the ENA file lines
for(line in lines){
  
  # Remove all of the leading and trailing white space characters (spaces, tab and newline)
  line <- trimws(line)
  
  # Remove all of the other white space characters (spaces)
  line <- gsub(pattern = ' ', replacement = '', line) 
  
  # Determine what the first character of the line is
  firstCharacter = substr(line,1,1)
  
  if(firstCharacter == '>'){                                # If the first character is a >...
    current_run <- sub(pattern = '>', replacement = '', line) # ... this line is the name of the current sequencing run...
    runs_samples[[current_run]] <- c()                       # ... this line is the key of a vector that will, later on, contain its samples.
  }
  
  else if(line == ''){ # If the line contains no information... 
    next              # ... skip this line.
  }
  
  else{                                                                 # If none of the conditions is valid...
    runs_samples[[current_run]] = c(runs_samples[[current_run]], line)   # ... the line is a sample name that is added to its corresponding sequencing run
  }
}


#####################################
## FETCHING ONLINE ENA INFORMATION ##
#####################################

# Loop over all sequencing runs
for(run in names(runs_samples)){
  
  # Create the path into which the samples of the current sequencing run should be downloaded
  path.download <- file.path(directory, run)
  
  if(!dir.exists(path.download)){ # If the directory does not exist yet...
    dir.create(path.download)     # ... create the directory.
  }
  
  # Loop over all the samples belonging to the current sequencing run.
  for(sample in runs_samples[[run]]){
    
    # Imagine that our current sample is ERR4914118
    # The link for this ENA sample is ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR491/008/ERR4914118/ERR4914118_1.fastq.gz
    # You can see that there is a 6 letter code in this link (ERR491)
    
    # We extract this 6 character code from the sample name
    six_letter_code <- substr(sample, start = 1, stop = 6)
    
    # You can also see that there is a 1 character code in this link, preceded by '00' (008)
    
    # We extact this 1 number code
    one_letter_code <- substr(sample, start = nchar(sample), stop = nchar(sample))
    
    # We can see the forward read filename in the link (ERR4914118_1.fastq.gz), so we construct this one
    fwd_file_name = paste0(sample, '_1.fastq.gz')
    
    # There is also a reverse read filename in the link for the reverse read file, so we construct this one
    rev_file_name = paste0(sample, '_2.fastq.gz')
    
    # Construct the entire link for the forward read
    url_fwd = paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', six_letter_code, '/', '00', one_letter_code, '/', sample, '/', fwd_file_name)
    
    # Construct the entire link for the reverse read
    url_rev = paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', six_letter_code, '/', '00', one_letter_code, '/', sample, '/', rev_file_name)
    
    # Define the location and filename for the files that will be downloaded from ENA
    dest_fwd <- file.path(path.download, fwd_file_name)
    dest_rev <- file.path(path.download, rev_file_name)
    
    # Define how many times we can attempt to download an ENA accession
    max_retries <- 3
    
    # Define a variable that defines the number of current retries
    retry_count <- 0
    
    # Define a variable that will specify if the download of the ENA file was succesful or not
    download_success <- FALSE
    
    while (!download_success && retry_count < max_retries) { # Keep on trying as long as downloading did not succeed (download_succes == FALSE) and the maximum amount of retries is lower than 3
      retry_count <- retry_count + 1                         # Add 1 to the amount of current retries
      tryCatch({                                  
        download.file(url_fwd, dest_fwd)                     # Download the information from the url to the defined destination
        download_success <- TRUE                             # Mark the download as a succes
      }, error = function(e) {                               # If something within the tryCatch statement gave an error...
        message(paste("Error downloading file, retrying (", retry_count, "/", max_retries, ")")) # ... print an error message...
        Sys.sleep(10)                                                                            # ... and wait 10 seconds.
      })
    }
    
    if (!download_success) { # If the download was not a succes (download_succes == FALSE)...
      stop("Failed to download file after", max_retries, "attempts") # stop the script and mention what accession did not work.
    }
    
    
    # The exact same thing as above, but for the reverse reads
    max_retries <- 3
    retry_count <- 0
    download_success <- FALSE
    
    while (!download_success && retry_count < max_retries) {
      retry_count <- retry_count + 1
      tryCatch({
        download.file(url_rev, dest_rev)
        download_success <- TRUE
      }, error = function(e) {
        message(paste("Error downloading file, retrying (", retry_count, "/", max_retries, ")"))
        Sys.sleep(10)
      })
    }
    
    if (!download_success) {
      stop("Failed to download file after", max_retries, "attempts")
    }
  }
}