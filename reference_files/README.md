# Building reference files

The following reference files are needed in order to run SPORK:

- Bowtie2 index files for the reference genome
- Bowtie2 index files for the regular junctions
- Bowtie2 index files for the scrambled junctions
- Bowtie2 index files for ribosome
- Bowtie2 index files for the transcriptome
- Bowtie2 index files for indel junctions (for junctions with up to 5 symmetric indels at the splice site)
- Pickle files for genes/exons annotation


# Software Requirements

- Bowtie2 2.2.9
- Python 3.4 with Biopython/1.70 installed


# Creating regular/scrambled junction fasta files

In order to create regular/scrambled junction fasta files, please run the following commands:
    
    1. mkdir index && python makeExonDB.py -f /path/to/fasta_file -a path to /path/to/gtf_file -o /path/to/makeExonsDB_output_dir 
    2. sh createJunctionIndex.sh path/to/output_dir  /path/to/makeExonsDB_output_dir prefix_string

# Creating pickle file for gene/exon annotation (gtfs_info)

In order to create pickle for gene/exon annotation, please run the following command:

    python extract_exons.py -gtf /path/to/gtf_file -rn reference_name
    
# Creating pickle file for gene/exon sequences (gtfs_info_seq)

In order to create pickle for gene/exon sequences, please run the following command:
    
    python create_gtf_info_with_seq_new_new.py -gtf /path/to/gtf_file -fa /path/to/fasta_file 


# Creating transcriptome file

In order to create transcriptome file, please use Tophat2 gtf_to_fasta tool:

    gtf_to_fasta /path/to/gtf_file /path/to/fasta_file gtf.fasta && awk '{if ($1 ~ /^>/) print ">"$2; else print $0}' gtf.fasta > transcripts.fa

# Creating indel junctions

In order to create indel junctions files, please run the following command:
    
    python subsample_reg_indels.py --input_file /path/to/reg_junctions.fa --probability p --reference_string reference_name
    
Parameter p represents the probability that regular junctions will be kept in indel junctions files. Since regular junctions fasta file for human genome is large, we have used p = 0.1

# Build bowtie2 indices

In order to build Bowtie2 index files for the reference genome, regular junctions, scrambled junctions, ribosome, transcriptome, indel junctions, known fusions, run the following command:

    bowtie2-build -f --threads numTreads ref.fa ./ref
 

