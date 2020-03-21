#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
"/mnt/d/pcRNA_SRA2/mre11-00m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-15m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-30m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-45m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-60m-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-2h-rep1" \
"/mnt/d/pcRNA_SRA2/mre11-wt-ts-rep1" \
)
##########################################


# processing arguments, proceed with reads subsetting
main() {
  numthreads=1                      # default number of parallel processes
  numsubset=1000000000              # upper limit for number of reads
  while getopts 'p:s:' opt; do      # pass number of reads to subset via -s option
    case "$opt" in
      p) numthreads="$OPTARG"
        ;;
      s) numsubset="$OPTARG"
        ;;
      \?) echo "Usage: $(basename $0) [-p number of threads] [-s number of reads to subset]"
      exit 1
        ;;
    esac
  done
  echo "Number of parallel processes: $numthreads"
  for file in "${filelist[@]}"; do  # process each file
    ((i=i%numthreads)); ((i++==0)) && wait
    subset ${numsubset} ${file} &
  done
}


# normalize by mapped read counts, only if script has been run before without [-s] flag so that "*_rmdup.bam" files have been generated in align2bam()
subset() {

  if test -f "$2_rmdup.bam"; then       # only run program if *_rmdup.bam exists

    # Subset BAM file by number of mapped reads
    total=$( samtools view -F 0x04 -c "$2_rmdup.bam" )
    frac=$(echo "$1/$total" | bc -l)
    if (( $(echo "$frac >= 1" | bc -l) )) ; then
      cp "$2_rmdup.bam" "$2.bam"
    else
      samtools view -bs ${frac} \
      "$2_rmdup.bam" > "$2.bam"
    fi

    # Index the subsetted BAM file
    samtools index \
    "$2.bam" ;

    # Delete index output for rmdup
#    rm "$2_rmdup.bam.bai"

    # Retrieve read count statistics
    samtools flagstat \
    "$2.bam" > "$2_flagstats.txt"

  else
    echo "ERROR: *_rmdup.bam does not exist. First run process_reads.sh without the -s flag"
  fi
}


main "$@"; exit
