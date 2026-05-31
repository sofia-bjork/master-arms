#!/usr/bin/env Rscript

# script to summarize all files that have been symlinked to the fastq_files dir

number <- c(1:5)

batch_list_output <- data.frame(
    Batch <- c(),
    File_name <- c()
)

for (nr in number){
    batch <- paste0("Batch_", nr)
    path    <- file.path("output", "novaseq_fastq", batch)
    batch_list <- list.files(path, pattern = ".fastq", full.names = FALSE)
    batch_list <- sort(batch_list[grep("_1_H", batch_list)])

    length <- length(batch_list)

    loop_dataframe <- data.frame(
        Batch <- rep(nr, length),
        File_name <- batch_list
    )
    colnames(loop_dataframe) <- c("Batch", "File_name") 

    batch_list_output <- rbind(batch_list_output, loop_dataframe)
    colnames(batch_list_output) <- c("Batch", "File_name")    
}

write.csv(batch_list_output, "output/extracted_sample_names.csv", 
          row.names = FALSE, quote = FALSE)
        
