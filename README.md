Custom analysis software for ChIP-seq of DNA repair proteins after CRISPR/Cas9-mediated DSBs
====

## Software requirements
- [python 2.7](https://www.anaconda.com/distribution/) (Anaconda's python distribution comes with the required numpy and scipy libraries)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- Ensure that both `samtools` and `bowtie2` are added to path and can be called directly from bash

## Data requirements
- Raw ChIP-seq FASTQ files and processed BAM files can be downloaded from SRA.

## Generate processed BAM files (if starting from raw FASTQ files)
1. Download sequencing reads in FASTQ format from SRA
2. Download human prebuilt bowtie2 indices
    - [Human hg38](http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)
    - move to the corresponding folder named `hg38_bowtie2/`
3. Download the human (hg38) genome assembly in FASTA format
    - [hg38.fa](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
    - extract from gz files
    - move to the corresponding folder named `hg38_bowtie2/`
4. Generate FASTA file indices
    - `samtools faidx hg38_bowtie2/hg38.fa`
5. Modify the first lines of bash script `process_reads.sh`
    - Fill the list variable `filelist` with paths to each **sample name** for processing. Note that the actual file paths include the **sample name** followed by "\_1.fastq" or "\_2.fastq", denoting read1 or read2 of paired-end reads, respectively.
    - Fill in the path to the indexed genome denoted by variable `genomepath`.
6. Run the following code snippet, where `-p` denotes the number of samples to process in parallel; modify accordingly. This performs genome alignment, filtering, sorting, removal of PCR duplicates, indexing, and sample statistics output. This step takes less than one day to complete on our Intel i7-8700K, 32GB RAM desktop, though speed appears to be bottlenecked by read/writes to disk.
    ```
    bash process_reads.sh -p 6
    ```
7. *(optional)* For fair comparison between different time points in a timeseries, we subset reads from all relevant samples to the sample with the fewest reads of the set. Run the following code snippet, where `-s` inputs the number of mapped reads to subset for each sample; modify accordingly. This step takes less than 10 minutes.
    ```
    bash subset_reads.sh -p 6 -s 24400000
    ```
    - for 53BP1 DNA-PKcs inhibitor experiments, subsetted to 9,045,179 for both replicates.
    - for γH2AX DNA-PKcs inhibitor experiments, subsetted to 13,474,505 for replicate 1 and 7,626,728 for replicate 2.
    - for MRE11 DNA-PKcs inhibitor experiments, subsetted to 5,507,748 for replicate 1 and 8,057,929 for replicate 2.
    - for 53BP1 timeseries experiments, subsetted to 8,481,406 for replicate 1 and 6,504,200 for replicate 2.
    - for MRE11 timeseries experiments, subsetted to 2,686,969 for replicate 1 and 6,207,727 for replicate 2.


## Start from pre-processed BAM files
In addition to raw paired-end reads in FASTQ format, we have also uploaded pre-processed sequencing reads in BAM format to SRA. These are the output of the previous section. It is highly recommended to start from these BAM files.
1. Download pre-processed paired-end reads in BAM format from SRA.
2. Move the downloaded data for MRE11 and γH2AX ChIP-seq to the desired folder.
3. If not already done, index the downloaded BAM files:
    ```
    samtools index /path/to/output.bam
    ```

## ChIP-seq analysis scripts in Python
#### `chipseq_M_pki.py`: analyze MRE11 ChIP-seq, with/without DNA-PKcs inhibitor KU60648
#### `chipseq_BH_pki.py`: analyze 53BP1 and γh2AX ChIP-seq, with/without DNA-PKcs inhibitor KU60648
#### `chipseq_M_ts.py`: analyze MRE11 ChIP-seq, timeseries after Cas9/pcRNA deactivation
#### `chipseq_B_ts.py`: analyze 53BP1 ChIP-seq, timeseries after Cas9/pcRNA deactivation
1. Set `base` variable to be the path to the directory that holds the BAM files.
2. Create a new folder to hold the output of 53BP1 and γH2AX scripts, set its path to `bs_a`.
3. Create a new folder to hold the output of MRE11 scripts, set its path to `bs_s`.
4. Ensure that all file names are correct (if the BAM files were directly downloaded from SRA, they should be), then run script.
