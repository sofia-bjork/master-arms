
# FASTA files from ARMS-MBON to taxonomic assignments 

## How to run
This repository contains the code and commands to process MiSeq and NovaSeq data from ARMS-MBON from quality controlled FASTA files to ASVs/OTUs/taxonomic assignments. The code specified in `execute/run/temporal_arms_run(1-4).sh` was run on an HPC using an Apptainer sandbox container based on [docker://continuumio/miniconda3:25.3.1-1](docker://continuumio/miniconda3:25.3.1-1). The environment used for building the sandbox can be found as `environment/git_env.yml`. The code was executed using `execute/run/temporal_arms_slurm.sh` for each of `execute/run/temporal_arms_run(1-4).sh` (in order) to yield the final taxonomic assignments. 


## Important notes
In `execute/run/temporal_arms_run1.sh`, MiSeq FASTA files and metadata is downloaded from ENA and NovaSeq files are symliked from an internal directory. After download/symlinking, input FASTA files are placed in `output/miseq_fastq/Run_X` and `output/novaseq_fastq/Batch_X` respectively. X represents the batch in where samples have been sequenced. 

The NovaSeq files were at the time of completion not made publically available. The NovaSeq processing code is thus adapted for a specific file structure in the first execute script `execute/run/temporal_arms_run1.sh`. Subsequent scripts are in some parts adapted for a naming convention that might not be applied in the final files when made available. 

## Additional scripts
This repository also contains the scripts used for statistical analysis and visualization in my thesis titled "Investigating Temporal Ecosystem Changes on Marine Hard Bottoms in Kosterhavet National Park Using Long-Term Metabarcoding Data". These scripts are found in `code/post-processing`. These scripts were run on MacOS 12.7.6, using the conda environment specified in `environment/post_processing_export.yml`. 


