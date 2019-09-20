# SPORK

This repository contains the JSON script written in Common Workflow Language (CWL) for Small Powerful Orthogonal Read Mapper from KNIFE (SPORK), a new, statistically driven algorithm for detecting of criptic splicing, a splicing not necessary at exon-exon boundaries and which could be within a gene or between genes in a fusion.

# Software Requirements

- Anaconda
- SciPy 1.1.0
- Bowtie2 2.2.9
- Python 3.4 with Biopython/1.70 installed
- java
- R 3.5.1
- Trim Galore! 0.4.4

# SPORK main script

All computational steps including any alignment step required for running SPORK has been packaged in a single JSON file "KNIFE_and_SPORK_pipeline.json". This script can be run on any local cluster using an input JSON file that provides the paths for reference index files and input RNA-seq file. To run the script, Rabix toolkit should be pre-installed on the local cluster

In addition, if one has KNIFE output files, SPORK could be run with KNIFE outputs as input files. THe workflow for running SPORK using KNIFE outputs 
has been packaged in a single JSON file "SPORK_pipeline.json". 

# Input file

All input parameters required for running SPORK JSON script should be provided via an input JSON file "SPORK_input.json". The following parameters should be set in the input JSON file:

- Bowtie2 index files for the reference genome
- Bowtie2 index files for the regular junctions
- Bowtie2 index files for the scrambled junctions
- Bowtie2 index files for ribosome
- Bowtie2 index files for the transcriptome
- Bowtie2 index files for indel junctions (for junctions with up to 5 symmetric indels at the splice site)
- Fastq files for the input RNA-Seq data
- Pickle file for gene/exon annotation (gtfs_info)
- Pickle file for gene/exon sequences (gtfs_info_seq)

Note: Since index files are too large and cannot be uploaded to githiub, we provide original reference fasta files (genome, transcriptome, regular junctions, scrambled junctions, ribosome) along with needed scripts/instructions for building index files here: 


# Toolkit for executing SPORK JSON script

For running SPORK JSON script on a local cluster, Rabix should be installed first. Rabix is an open-source tool developed by Seven Bridges, that can be used to run a computational workflow written in Common Workflow Language (CWL) on a locul cluster. More information on how to install Rabix can be found in this GitHub repositiory: https://github.com/rabix/bunny 


# An example Batch script for submitting SPORK jobs on a local cluster

An example batch script "KNIFE_and_SPORK_submit_job.sbatch", based on job scheduler Slurm has been provided. In the batch script file the path to where Rabix has been installed, SPORK pipeline JSON file (KNIFE_and_SPORK_pipeline.json), and SPORK input JSON file (KNIFE_and_SPORK_input.json) should be provided. 

In addition, if one has KNIFE output files, SPORK could be run with KNIFE outputs as input files. An example batch script "SPORK_submit_job.sbatch", based on job scheduler Slurm has been provided. In the batch script file the path to where Rabix has been installed, SPORK pipeline JSON file (SPORK_pipeline.json), and SPORK input JSON file (SPORK_input.json) should be provided. 

# Output files

A main report file containing reported criptic junctions with their corresponding statsitical scores and number of various types of aligned reads. High confidence junctions are called via filtering based on the statistical scores (all needed filters and their thresholds are described in the paper).
