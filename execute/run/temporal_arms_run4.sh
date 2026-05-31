#!/usr/bin/env bash

export PATH=/envs/git_env/bin:$PATH

set -e

cd /cfs/klemming/home/b/bjorkso/Private/temporal-arms-main

# subset non-nuMT ASVs from the COI_nochim_nosingle.fa file
grep -w -A 1 -f output/COI/pseudo/nonpseudo.combined.names.txt output/COI/COI_nochim_nosingle_ASVs.fa > output/COI/pseudo/COI_nochim_nosingle_nopseudo.fa --no-group-separator

# negative control correction, read filtering
Rscript code/reads_to_OTUs/COI_blank_corr.R
echo "-----------------blank correction, reads filtering done-----------------"

# subset the ASVs remaining after negative control correction
grep -w -A 1 -f output/COI/blank_corr/no_contam_headers_COI.txt output/COI/pseudo/COI_nochim_nosingle_nopseudo.fa --no-group-separator > output/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa 

# generate headers with the abundance of each ASV included
Rscript code/reads_to_OTUs/dereplication_headers.R
echo "-----------------dereplication headers done-----------------"

# replace headers with dereplicated header names
grep -v "^--" output/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa | awk 'NR%2==0' | paste -d'\n' output/COI/ASV_dereplicated.txt - > output/COI/COI_dereplicated_ASVs.fa

# change directory for swarm to run 
pushd output/COI/

# cluster ASVs into MOTUs using swarm
/envs/git_env/bin/swarm -d 13 -i swarm/internal.txt -o swarm/output.txt -s swarm/statistics.txt -u swarm/uclust.txt -w swarm/COI_cluster_reps.fa COI_dereplicated_ASVs.fa
echo "-----------------OTU clustering done-----------------"

# move back to root
popd

# generate MOTU tables from swarm output 
Rscript code/reads_to_OTUs/MOTU_tables.R

# remove read abundance line from sequence header
awk -F'_' '{print $1}' output/COI/swarm/COI_cluster_reps.fa > output/COI/MOTU/COI_cluster_reps_lulu_ready.fa

# generate match lists using BLASTn
makeblastdb -in output/COI/MOTU/COI_cluster_reps_lulu_ready.fa -parse_seqids -dbtype nucl
blastn -db output/COI/MOTU/COI_cluster_reps_lulu_ready.fa -outfmt '6 qseqid sseqid pident' -out output/COI/MOTU/match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query output/COI/MOTU/COI_cluster_reps_lulu_ready.fa

# LULU curation
Rscript code/reads_to_OTUs/LULU_curation.R
echo "-----------------LULU curation done-----------------"

# generate fasta files with remaining MOTUs after LULU curation
grep -w -A 1 -f output/COI/MOTU/lulu_curated_headers.txt output/COI/MOTU/COI_cluster_reps_lulu_ready.fa --no-group-separator > output/COI/MOTU/COI_cluster_reps_lulu_curated.fa


# taxonomic assignment using BOLDigger3 on public animal library (--db 1) on exhaustive search mode (--mode 3)
pip install boldigger3==2.1.4
pip install lxml-html-clean==0.4.3
# this loop was written with chatgpt
for f in output/COI/MOTU/chunks/*.fasta; do
    echo "Processing $f"
    
    boldigger3 identify  "$f" --db 1 --mode 3
    
    echo "Sleeping 15 minutes to avoid rate limit..."
    sleep 900
done
# boldigger3 identify output/COI/MOTU/COI_cluster_reps_lulu_curated.fa --db 1 --mode 3
echo "-----------------BOLDigger done-----------------"

# create taxonomic table from boldigger output 
Rscript code/reads_to_OTUs/BOLDigger_tax_table.R

echo "-----------------script finished successfully-----------------"
