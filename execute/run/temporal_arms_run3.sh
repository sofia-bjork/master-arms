#!/usr/bin/env bash

export PATH=/envs/git_env/bin:$PATH

set -e

cd /cfs/klemming/home/b/bjorkso/Private/temporal-arms-main

# remove chimeras and singletons, merge ASV data sets
Rscript code/reads_to_OTUs/COI_chimera.R
echo "-----------------chimera removal done-----------------"

# identify and remove nuclear mitochondrial DNA pseudogenes (nuMTs)
Rscript code/reads_to_OTUs/MACSE_align_pseudo.R -d /envs/git_env/share/macse-2.07-0/macse_v2.07.jar
echo "-----------------nuMT removal done-----------------"

echo "-----------------script finished successfully-----------------"
