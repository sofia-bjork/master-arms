#!/usr/bin/env Rscript

### dada2 COI workflow with cutadapt primer removal ###

COI_trimming_filtering_func <- function(seq_batch, cutadapt) {

run_id  <- paste("Run_", seq_batch, sep = "")

# directory containing the fastq.gz files
path    <- file.path("output", "miseq_fastq", run_id)

# COI filtering directory
output_coi_dir <- file.path("output", "COI")
if(!dir.exists(output_coi_dir)) dir.create(output_coi_dir)

coi_dir <- file.path("output", "COI", "miseq_dada2")
if(!dir.exists(coi_dir)) dir.create(coi_dir)


img_id  <- gsub("_", "", run_id)


message("The specified directory contains the following files:")
print(list.files(path))


# Use cutadapt for primer removal (prior installation of cutadapt on your machine via python, anaconda, etc. required)
# directory containing cutadapt (if using a conda env, find by running > which cutadapt in terminal)

system2(cutadapt, args = "--version") # see if R recognizes cutadapt and shows its version


# load / install necessary packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("dada2", version = 1.30.26) # if this does not work, try to install via devtools (requires prior installation of devtools)
# BiocManager::install("ShortRead", version = 1.30.26)
# BiocManager::install("Biostrings", version = 1.30.26)

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)


# generate matched lists of the forward and reverse read files, as well as parsing out the sample name

fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))


# Designate sequences [including ambiguous nucleotides (base = N, Y, W, etc.) if present) of the primers used
# The reverse COI primer jgHCO2198 contains Inosine nucleotides.
# These "I" bases are not part of IUPAC convention and are not recognized by the packages used here. Change "I"s to "N"s.

FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"  ## forward primer sequence
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"  ## reverse primer sequence


# Verify the presence and orientation of these primers in the data

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients


# Calculate number of reads containing forward and reverse primer sequences (considering all possible primer orientations. Only exact matches are found.).
# Only one set of paired end fastq.gz files will be checked (second sample in this case).
# This is is sufficient, assuming all the files were created using the same library preparation.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  read_fastq <- readFastq(fn)
  nhits <- vcountPattern(primer, sread(read_fastq), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs))

# Output:
# FWD primer should mainly be found in the forward reads in its forward orientation.
# REV primer should mainly be found in the reverse reads in its forward orientation.


# Create output filenames for the cutadapt-ed files.
# Define the parameters for the cutadapt command.
# See here for a detailed explanation of paramter settings: https://cutadapt.readthedocs.io/en/stable/guide.html#

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 
# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.1 --discard-untrimmed", R1.flags, R2.flags,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

# see here for a detailed explanation of the output:
# https://cutadapt.readthedocs.io/en/stable/guide.html#cutadapt-s-output
# Sometimes, you will see this: "WARNING: One or more of your adapter sequences may be incomplete. Please see the detailed output above."
# This usually refers to: "WARNING: The adapter is preceded by "T" (or any other base) extremely often. The provided adapter sequence could be incomplete at its 3' end."
# The amplified regions and primer binding sites are usually highly conserved, so primer sequences are often preceded by the same base.
# Cutadapt just warns us that this is the case and tells us to check if the preceding base is indeed not part of the primer. 


# Count the presence of primers in the first cutadapt-ed sample to check if cutadapt worked as intended:

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut))

# The primer-free sequence read files are now ready to be analyzed.
# Similar to the earlier steps of reading in FASTQ files, read in the names of the cutadapt-ed FASTQ files. 
# Apply some string manipulation to get the matched lists of forward and reverse fastq files.


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))


# Check if forward and reverse files match:

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutFs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

print("The extracted sample names are:")
print(head(sample.names))


# Inspect read quality profiles. 
# If there are more than 20 samples, grab 20 randomly

set.seed(1)

if(length(cutFs) <= 20) {
  fwd_qual_plots<-plotQualityProfile(cutFs) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
  rev_qual_plots<-plotQualityProfile(cutRs) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
}
# fwd_qual_plots
# rev_qual_plots


# Print out the forward quality plot

ggsave(paste0("COI_", img_id, "_quality_forward.jpg"),
       plot = fwd_qual_plots, path = coi_dir, width = 15,
       height = 8, units = "in", dpi = 300)


# Print out the reverse quality plot

ggsave(paste0("COI_", img_id, "_quality_reverse.jpg"),
       plot = rev_qual_plots, path = coi_dir, width = 15,
       height = 8, units = "in", dpi = 300)


## Filter and trim ##

# Assign filenames to the fastq.gz files of filtered and trimmed reads.

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


# Set filter and trim parameters.

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen =c(200,130),maxN = 0, maxEE = c(2,4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = T) 


# Save this output as RDS file for the read tracking table created downstream:
saveRDS(out, file = file.path(coi_dir, paste("filter_and_trim_out_", img_id, ".rds", sep = "")))


# check how many reads remain after filtering

out


# Check if file names match

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtFs) <- sample.names
names(filtRs) <- sample.namesR


# Estimate error models of the amplicon dataset. 

set.seed(100) # set seed to ensure that randomized steps are replicatable
errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)


# save error calculation as RDS files:

saveRDS(errF, file = file.path(coi_dir, paste("errF_", img_id, ".rds", sep = "")))
saveRDS(errR, file = file.path(coi_dir, paste("errR_", img_id, ".rds", sep = "")))


# As a sanity check, visualize the estimated error rates and write to file:

plot_err_F <- plotErrors(errF, nominalQ = TRUE)
plot_err_R <- plotErrors(errR, nominalQ = TRUE)

ggsave(paste0("COI_", img_id, "_error_forward.jpg"),
       plot = plot_err_F, path = coi_dir, width = 15,
       height = 8, units = "in", dpi = 300)

ggsave(paste0("COI_", img_id, "_error_reverse.jpg"),
       plot = plot_err_R, path = coi_dir, width = 15,
       height = 8, units = "in", dpi = 300)


### The dada2 tutorial implements a dereplication step at this point. 
### This does not seem to be necessary any more with the newer dada2 versions, according to what the developers stated in the dada2 github forum.

# Apply the dada2's core sequence-variant inference algorithm:

# Set pool = pseudo", see https://benjjneb.github.io/dada2/pool.html

dadaFs <- dada(filtFs, err = errF, multithread = T, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = T, pool = "pseudo")


# Apply the sample names extracted earlier (see above) to remove the long fastq.gz file names
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names


# Save sequence-variant inference output as RDS files: 

saveRDS(dadaFs, file = file.path(coi_dir, paste("dadaFs_", img_id, ".rds", sep = "")))
saveRDS(dadaRs, file = file.path(coi_dir, paste("dadaRs_", img_id, ".rds", sep = "")))


# Merge the forward and reverse reads.
# Adjust the minimum overlap (default = 12) and maximum mismatch allowed if necessary.

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                      minOverlap = 10, maxMismatch = 1, verbose = TRUE)

saveRDS(mergers, file = file.path(coi_dir, paste("mergers_", img_id, ".rds", sep = "")))


# Construct an amplicon sequence variant table (ASV) table
# If maxMismatch > 0 has been allowed in the mergePairs step,
# "Duplicate sequences detected and merged" may appear as output during the sequence table creation
# This is not a problem, just ignore it.


seqtab <- makeSequenceTable(mergers)

# if only one sample is running, makeSequenceTable will not automatically add a row name
# this becomes problematic further downstream, so here we add row names for runs with only one sample
if(length(sample.names) == 1) {
    rownames(seqtab) <- sample.names
}

# How many sequence variants were inferred? 
# The first column displays the number of samples. The second column displays the number of ASVs. 
dim(seqtab)


# Save sequence table

saveRDS(seqtab, file = file.path(coi_dir, paste("seqtab_", img_id, ".rds", sep = "")))


## Track reads throughout the pipeline ##

# Get number of reads in files prior to cutadapt application

input <- countFastq(path, pattern = ".gz")          # get statistics from input files
input$Sample <- rownames(input)
input$Sample <- gsub("_.*", "", input$Sample)       # Remove all characters after _ (incl. _) in file names
input <- aggregate(.~Sample, input, FUN = "mean")   # Aggregate forward and reverse read files 


# Get number of reads from each step of dada2 pipeline

getN <- function(x) sum(getUniques(x))

if(length(sample.names) == 1) {
    track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers))
} else {
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,getN)) }

colnames(track) <- c("cutadapt", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names


# Combine with read numbers from input files

input<-input[order(match(input[,1],rownames(track))),]
track<-cbind(input$records,track)
colnames(track)[1] <- "input"


# Save to file

write.table(track, 
            file = file.path(coi_dir, paste("track_", img_id, ".txt", sep = "")),
            sep = "\t", col.names = NA)

}