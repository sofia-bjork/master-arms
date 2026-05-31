#!/usr/bin/env bash

export PATH=/envs/git_env/bin:$PATH

set -e

cd /cfs/klemming/home/b/bjorkso/Private/temporal-arms-main

# filter and trim using cutadapt and dada2 (miseq)
Rscript code/reads_to_OTUs/miseq_run_filter_trim.R -d /envs/git_env/bin/cutadapt
echo "-----------------filter and trim (miseq) done-----------------"

# filter and trim using cutadapt and dada2 (novaseq)
Rscript code/reads_to_OTUs/novaseq_filter_trim.R -d /envs/git_env/bin/cutadapt
echo "-----------------filter and trim (novaseq) done-----------------"

echo "-----------------script finished successfully-----------------"
