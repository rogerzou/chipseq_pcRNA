#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
"/mnt/d/pcRNA_SRA1/mre11-wt" \
"/mnt/d/pcRNA_SRA1/gh2ax-wt" \
"/mnt/d/pcRNA_SRA1/53bp1-wt" \
"/mnt/d/pcRNA_SRA1/53bp1-1h-nL-nD-rep1" \
"/mnt/d/pcRNA_SRA1/53bp1-4h-nL-nD-rep1" \
"/mnt/d/pcRNA_SRA1/53bp1-4h-L-nD-rep1" \
"/mnt/d/pcRNA_SRA1/53bp1-4h-L-PKi-rep1" \
"/mnt/d/pcRNA_SRA1/gh2ax-1h-nL-nD-rep1" \
"/mnt/d/pcRNA_SRA1/gh2ax-4h-nL-nD-rep1" \
"/mnt/d/pcRNA_SRA1/gh2ax-4h-L-nD-rep1" \
"/mnt/d/pcRNA_SRA1/gh2ax-4h-L-PKi-rep1" \
"/mnt/d/pcRNA_SRA1/mre11-1h-nL-nD-rep1" \
"/mnt/d/pcRNA_SRA1/mre11-4h-nL-nD-rep1" \
"/mnt/d/pcRNA_SRA1/mre11-4h-L-nD-rep1" \
"/mnt/d/pcRNA_SRA1/mre11-4h-L-PKi-rep1" \
"/mnt/d/pcRNA_SRA1/53bp1-1h-nL-nD-rep2" \
"/mnt/d/pcRNA_SRA1/53bp1-4h-nL-nD-rep2" \
"/mnt/d/pcRNA_SRA1/53bp1-4h-L-nD-rep2" \
"/mnt/d/pcRNA_SRA1/53bp1-4h-L-PKi-rep2" \
"/mnt/d/pcRNA_SRA1/gh2ax-1h-nL-nD-rep2" \
"/mnt/d/pcRNA_SRA1/gh2ax-4h-nL-nD-rep2" \
"/mnt/d/pcRNA_SRA1/gh2ax-4h-L-nD-rep2" \
"/mnt/d/pcRNA_SRA1/gh2ax-4h-L-PKi-rep2" \
"/mnt/d/pcRNA_SRA1/mre11-1h-nL-nD-rep2" \
"/mnt/d/pcRNA_SRA1/mre11-4h-nL-nD-rep2" \
"/mnt/d/pcRNA_SRA1/mre11-4h-L-nD-rep2" \
"/mnt/d/pcRNA_SRA1/mre11-4h-L-PKi-rep2" \
"/mnt/d/pcRNA_SRA2/53bp1-00m-rep1" \
"/mnt/d/pcRNA_SRA2/53bp1-30m-rep1" \
"/mnt/d/pcRNA_SRA2/53bp1-60m-rep1" \
"/mnt/d/pcRNA_SRA2/53bp1-2h-rep1" \
"/mnt/d/pcRNA_SRA2/53bp1-00m-rep2" \
"/mnt/d/pcRNA_SRA2/53bp1-30m-rep2" \
"/mnt/d/pcRNA_SRA2/53bp1-60m-rep2" \
"/mnt/d/pcRNA_SRA2/53bp1-2h-rep2" \
"/mnt/d/pcRNA_SRA2/mre11-00m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-15m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-30m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-45m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-60m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-2h-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-00m-rep2" \
"/mnt/d/pcRNA_SRA2/mre11-15m-rep2" \
"/mnt/d/pcRNA_SRA2/mre11-30m-rep2" \
"/mnt/d/pcRNA_SRA2/mre11-45m-rep2" \
"/mnt/d/pcRNA_SRA2/mre11-60m-rep2" \
"/mnt/d/pcRNA_SRA2/mre11-2h-rep2" \
)
# Enter path to indexed genome
genomepath="/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38"
##########################################


# processing arguments, proceed with bioinformatics pipeline
main() {
  numthreads=1                      # default number of parallel processes
  while getopts 'p:' opt; do      # pass number of reads to subset via -s option
    case "$opt" in
      p) numthreads="$OPTARG"
        ;;
      \?) echo "Usage: $(basename $0) [-p number of threads] [-s number of reads to subset]"
      exit 1
        ;;
    esac
  done
  echo "Number of parallel processes: $numthreads"
  for file in "${filelist[@]}"; do  # process each file
    ((i=i%numthreads)); ((i++==0)) && wait
    align2bam ${genomepath} ${file} &
  done
}


# main bioinformatics pipeline (alignment to indexing and read statistics)
align2bam() {

  # Align reads to either hg38 or mm10 using bowtie2
  bowtie2 -p 6 -q --local -X 1000 -x $1 \
  -1 "$2_1.fastq" -2 "$2_2.fastq" -S "$2.sam" ;

  # Convert from SAM to BAM, filter for mapping quality >=25 and singleton reads
  samtools view -h -S -b -F 0x08 -q 25 \
  "$2.sam" > "$2_unsorted.bam" ;

  # Add mate score tags to ensure that paired-end reads contain correct information about the mate reads
  samtools fixmate -m \
  "$2_unsorted.bam" "$2_fixmate.bam" ;

  # Sort BAM file entries by genomic position
  samtools sort \
  "$2_fixmate.bam" > "$2_sorted.bam" ;

  # Remove potential PCR duplicates
  samtools markdup -r \
  "$2_sorted.bam" "$2_rmdup.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$2_rmdup.bam" ;

  # Retrieve read count statistics
  samtools flagstat \
  "$2_rmdup.bam" > "$2_flagstats.txt"

}

main "$@"; exit
