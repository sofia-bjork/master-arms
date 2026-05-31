#!/usr/bin/env bash

export PATH=/envs/git_env/bin:$PATH

set -e

cd /cfs/klemming/home/b/bjorkso/Private/temporal-arms-main

# symlink the files from genoscope in the project space (novaseq)
Rscript code/sample_download/novaseq_symlinks.R
echo "-----------------symlinks (novaseq) done-----------------"

# create a list of files symlinked (novaseq)
Rscript code/metadata_processing/novaseq_file_list.R
echo "-----------------file list (novaseq) done-----------------"

# download metadata files and specify parameters of interest for miseq
# -o observatory -p preservative -u UnitID -f fraction
Rscript code/metadata_processing/miseq_sample_summary.R -o Koster -p DMSO
echo "-----------------metadata download (miseq) done-----------------"

# create run number summary file (miseq)
Rscript code/metadata_processing/miseq_sequencing_run.R
echo "-----------------run number summary (miseq) done-----------------"

# create miseq fastq files directory 
mkdir -p output/miseq_fastq

# download the specified ENA fastq files (miseq)
Rscript code/sample_download/ENADownload.R -f metadata/generated_meta/COI_ENA_accessions.txt -d output/miseq_fastq
Rscript code/sample_download/ENADownload_8digits.R -f metadata/generated_meta/COI_ENA_accessions_8digits.txt -d output/miseq_fastq
echo "-----------------ENA download (miseq) done-----------------"

echo "-----------------script finished successfully-----------------"
