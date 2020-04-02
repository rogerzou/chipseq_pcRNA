# -*- coding: utf-8 -*-
""" Useful functions for ChIP-seq analysis after Cas9 cleavage
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

from collections import defaultdict
import pysam
import os
import re
import numpy as np
import csv


def status_statement(current, final, count, chr=None):
    """ Print progress statements for long processes

    :param current: current iteration number
    :param final: total number of iterations
    :param count: number of print statements total
    :param chr: print chromosome information of current progress

    """
    if current % int(final/count) == 0:
        if chr is None:
            print("Processed %i out of %i" % (current, final))
        else:
            print("Processed %i out of %i in %s" % (current, final, chr))


def hg38_generator():
    """ Generate the number of base pairs for each chromosome in hg38 (human)

    :return generator that outputs the number of base pairs for each chromosome in hg38
            in the format: [chr7, 159345973]

    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/hg38.sizes", 'r') as f:
        for row in csv.reader(f, dialect='excel', delimiter='\t'):
            yield row


def mm10_generator():
    """ Generate the number of base pairs for each chromosome in mm10 (mouse)

    :return generator that outputs the number of base pairs for each chromosome in mm10
            in the format: [chr7, 159345973]

    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/mm10.sizes", 'r') as f:
        for row in csv.reader(f, dialect='excel', delimiter='\t'):
            yield row


def read_pair_generator(bam, region_string=None):
    """ Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.

    :param bam: pysam AlignmentFile loaded with BAM file that contains paired-end reads
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :return generator of read pairs with the following format: [read1, read2]

    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_pair_align(read1, read2):
    """ Extract read pair locations as a fragment oriented in increasing chromosome coordinates

    :param read1: read #1 of pair in pysam AlignedSegment format
    :param read2: read #2 of pair in pysam AlignedSegment format
    :return 4-item array in the following format: [fragA-start, fragA-end, fragB-start, fragB-end]
            with monotonically increasing chromosome coordinates

    """
    r1pos = [x+1 for x in read1.positions]
    r2pos = [x+1 for x in read2.positions]
    if read1.mate_is_reverse and r1pos[0] < r2pos[0]:  # read1 is earlier
        read = [r1pos[0], r1pos[-1], r2pos[0], r2pos[-1]]
    elif read2.mate_is_reverse and r2pos[0] < r1pos[0]:  # read2 is earlier
        read = [r2pos[0], r2pos[-1], r1pos[0], r1pos[-1]]
    else:
        read = []
        # print("Skipping read pair from error in alignment.")
    # print("%s--%s>  <%s--%s" % tuple(read))
    return read


def get_read_subsets(filein, fileout, region_string, target, bool_array=None):
    """ Categorize fragments (paired reads) by its relationship to the target (cut) site

    For a target (cut) site, subset the fragments within region_range to
    - fragments that span the cut site
    - fragments that do not -> fragments that start from the cut site and span left vs right

    :param filein: BAM file containing paired-end reads
    :param fileout: base output file name with extension (.bam) omitted - function generates
                    multiple files where each file only contains reads from a specific category
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :param target: coordinate within region_string that denotes the 4th nucleotide from PAM
    :param bool_array: boolean array indicating whether to output the subsetted reads in the
                       following order: [M, N, L, R, S1, S2, I]

    """

    if bool_array is None:
        bool_array = [True, True, False, False]
    bam = pysam.AlignmentFile(filein, 'rb')
    countM = 0      # fragments that span
    countN = 0      # fragments that don't span
    countL = 0      # fragments that don't span, on left
    countR = 0      # fragments that don't span, on right
    counter = 0     # count total number of reads
    fileM = fileout + "_M.bam"
    fileN = fileout + "_N.bam"
    fileL = fileout + "_L.bam"
    fileR = fileout + "_R.bam"
    readsM = pysam.AlignmentFile(fileM, "wb", template=bam)
    readsN = pysam.AlignmentFile(fileN, "wb", template=bam)
    readsL = pysam.AlignmentFile(fileL, "wb", template=bam)
    readsR = pysam.AlignmentFile(fileR, "wb", template=bam)
    for read1, read2 in read_pair_generator(bam, region_string):
        read = read_pair_align(read1, read2)
        if not read:
            continue
        counter += 1
        if read[0] < target < read[3]:          # fragments that span
            countM += 1
            readsM.write(read1)
            readsM.write(read2)
        elif target + 5 >= read[0] >= target:   # fragments that begin 5bp of cleavage site
            countR += 1
            readsR.write(read1)
            readsR.write(read2)
            countN += 1
            readsN.write(read1)
            readsN.write(read2)
        elif target - 5 <= read[-1] <= target:   # fragments that end 5bp of cleavage site
            countL += 1
            readsL.write(read1)
            readsL.write(read2)
            countN += 1
            readsN.write(read1)
            readsN.write(read2)
    readsM.close()
    readsN.close()
    readsL.close()
    readsR.close()
    file_array = [fileM, fileN, fileL, fileR]
    for i, boo in enumerate(bool_array):
        file_i = file_array[i]
        if boo:
            pysam.sort("-o", file_i, file_i)
            os.system("samtools index " + file_i)
        else:
            os.system("rm " + file_i)
    print("%i span | %i end/start | %i window total | %i mapped total"
          % (countM, countN, counter, bam.mapped))
    bam.close()


def to_wiggle_pairs(filein, fileout, region_string, endcrop=False):
    """ Constructs fragment pile-ups in wiggle format using paired-end information

    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.wig) omitted
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :param endcrop: if True, don't count the first and last read of a fragment - the output for each
                    position is effectively a safe underestimate for the number of spanning reads

    """
    [chr, sta, end] = re.split('[:-]', region_string)
    sta = int(sta)
    end = int(end)
    wlist = [0] * (end-sta+1)
    wig = open(fileout + ".wig", "w")
    wig.write("variableStep\tchrom=%s\n" % chr)
    bam = pysam.AlignmentFile(filein, 'rb')
    for read1, read2 in read_pair_generator(bam, region_string):
        read = read_pair_align(read1, read2)
        if not read:
            continue
        if endcrop:
            wlist = [x+1 if read[0]-sta < i < read[-1]-sta else x for i, x in enumerate(wlist)]
        else:
            wlist = [x+1 if read[0]-sta <= i <= read[-1]-sta else x for i, x in enumerate(wlist)]
    for i, x in enumerate(wlist):
        wig.write("%i\t%i\n" % (sta + i, x))
    wig.close()
    bam.close()


def to_wiggle_windows(filein, fileout, window, chr=None, generator=None):
    """ Outputs wiggle file that counts reads per million in each window

    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.wig) omitted
    :param window: size of window (i.e. if 500, genome is divided into 500bp windows, function
                outputs the number of reads in each window)
    :param chr: array of chromosome strings to limit analysis to particular chromosomes
                (i.e. ['chr7', 'chr8']). If not set, then all chromosomes are processed.
    :param generator: species-specific generator that outputs the number of base pairs for each
                      chromosome in the format: [chr7, 159345973]. If not set, then hg38 is used.

    """
    if not generator:
        generator = hg38_generator()
    bam = pysam.AlignmentFile(filein, 'rb')
    treads = float(bam.mapped)
    wig = open(fileout + ".wig", "w")
    for row in generator:  # iterate over each chromosome
        if chr is None or (chr is not None and row[0] in chr):      # checks for chr #
            chr_i = row[0]
            wig.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
            numwins = int(int(row[1]) / window)  # number of windows
            for i in range(numwins):              # iterate over each bin (total - 1)
                start1 = i * window
                finish1 = (i + 1) * window - 1
                rpm = bam.count(chr_i, start1, finish1) / treads * 1E6
                wig.write("%0.6f\n" % rpm)
                status_statement(i, numwins, 20, chr_i)
    wig.close()
    bam.close()


def to_bins(filein, fileout, window, numbins, chr=None, generator=None):
    """ Outputs CSV file that counts number of reads in each bin (column) for each window (row)

    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.csv) omitted; each row corresponds
                    to each window, each column corresponds to the number of reads for each bin
    :param window: size of window (i.e. if 500, genome is divided into 500bp windows)
    :param numbins: number of bins for each window
    :param chr: array of chromosome strings to limit analysis to particular chromosomes
                (i.e. ['chr7', 'chr8']). If not set, then all chromosomes are processed.
    :param generator: species-specific generator that outputs the number of base pairs for each
                      chromosome in the format: [chr7, 159345973]. If not set, then hg38 is used.

    """
    if not generator:
        generator = hg38_generator()
    bam = pysam.AlignmentFile(filein, 'rb')
    cm = []
    for row in generator:                                           # iterate over each chromosome
        if chr is None or (chr is not None and row[0] in chr):      # checks for chr #
            count = int(int(row[1]) / window)                       # number of windows
            res = int(window / numbins)
            chr_i = row[0]
            for i in range(count):                                  # iterate over each window
                win_start = i * window
                win_finish = (i + 1) * window - 1
                cm_i = np.zeros(3 + numbins, dtype=object)          # array to hold bin counts info
                cm_i[0] = chr_i
                cm_i[1] = win_start
                cm_i[2] = win_finish
                for j in range(numbins):                            # iterate over each bin
                    bin_start = win_start + j * res
                    bin_finish = win_start + (j + 1) * res - 1
                    cm_i[j + 3] = bam.count(chr_i, bin_start, bin_finish)
                cm.append(cm_i)
                status_statement(i, count, 20, chr_i)
    np.savetxt(fileout + ".csv", np.asarray(cm), fmt='%s', delimiter=',')
    bam.close()


def avgwig(file1, file2, fileout):
    """
    Given two wig files of same structure, output the average of both as a wig file
    :param file1: first wig file
    :param file2: second wig file
    :param fileout: wig file (no extension) that outputs the average of the two input files

    """
    with open(file1, "r") as f1, open(file2, "r") as f2, open(fileout + ".wig", "w") as outwig:
        for line1, line2 in zip(f1, f2):
            if line1.split("\t")[0] == "fixedStep" or line2.split("\t")[0] == "fixedStep":
                if line1 == line2:
                    outwig.write(line1)
                else:
                    raise ValueError("The files are too different to merge!")
            else:
                val1 = float(line1.split("\n")[0])
                val2 = float(line2.split("\n")[0])
                outwig.write("%0.5f\n" % np.mean([val2, val1]))


def percentchange(file1, file2, fileout):
    """
    Given two wig files, calculate the percent change from file1 to file2
    :param file1: first wig file
    :param file2: second wig file
    :param fileout: output wig file (extension omitted)

    """
    with open(file1, "r") as f1, open(file2, "r") as f2, open(fileout + ".wig", "w") as outwig:
        for l1, l2 in zip(f1, f2):
            if l1.split("\t")[0] == "fixedStep" or l2.split("\t")[0] == "fixedStep":
                if l1 == l2:
                    outwig.write(l1)
                else:
                    raise ValueError("The files are too different to compare!")
            else:
                val1 = float(l1.split("\n")[0])
                val2 = float(l2.split("\n")[0])
                if val1 == 0:
                    outwig.write("%0.5f\n" % 0)
                else:
                    outwig.write("%0.5f\n" % ((val2 - val1) / val1))
