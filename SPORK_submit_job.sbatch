#!/bin/bash
#SBATCH --job-name=SPORK
#SBATCH --output=SPORK.%j.out
#SBATCH --error=SPORK.%j.err
#SBATCH --time=24:00:00
#SBATCH -p bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB
#SBATCH --export=RABIX.CORES=32,RABIX.MEMORY=531520,ALL
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jorda@stanford.edu
source /home/users/jorda/KNIFE_ONLY_getUnalignedReadsOverlap_python2/bin/activate
module load java
module load R
module load python/3.6.1
export PATH="$PATH:/scratch/PI/horence/jorda/bowtie2/bowtie2-2.2.9/"
export PATH="$PATH:/scratch/PI/horence/jorda/TrimGalore-0.4.4/"
Rscript $HOME/packages.R
#module load python/3.4
#pip3 install --user docopt
#pip3 install --user typing
java -jar /scratch/PI/horence/jorda/bunny-1.0.5-1/rabix-cli/local/rabix-cli-1.0.5.jar --no-container --verbose /scratch/PI/horence/jorda/SPORK/SPORK_pipeline.json /scratch/PI/horence/jorda/SPORK/SPORK_input.json
