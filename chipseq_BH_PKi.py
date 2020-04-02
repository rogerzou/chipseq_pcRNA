# -*- coding: utf-8 -*-
""" 53BP1 and gH2AX ChIP-seq analysis for Cas9/pcRNA - DNA-PKcs inhibitor effect """

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
bs = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA1/"
bs_a = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA1/analysis/"

chrs = ['chr7', 'chr8']


""" Convert BAM file to WIG file that counts the number of reads in each window span. """
win = 5000
c.to_wiggle_windows(bs + "53bp1-1h-nL-nD-rep1.bam", bs_a + "53bp1-1h-nL-nD-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-4h-nL-nD-rep1.bam", bs_a + "53bp1-4h-nL-nD-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-4h-L-nD-rep1.bam", bs_a + "53bp1-4h-L-nD-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-4h-L-PKi-rep1.bam", bs_a + "53bp1-4h-L-PKi-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-1h-nL-nD-rep2.bam", bs_a + "53bp1-1h-nL-nD-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-4h-nL-nD-rep2.bam", bs_a + "53bp1-4h-nL-nD-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-4h-L-nD-rep2.bam", bs_a + "53bp1-4h-L-nD-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-4h-L-PKi-rep2.bam", bs_a + "53bp1-4h-L-PKi-rep2", win, chrs)

c.to_wiggle_windows(bs + "gh2ax-1h-nL-nD-rep1.bam", bs_a + "gh2ax-1h-nL-nD-rep1", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-4h-nL-nD-rep1.bam", bs_a + "gh2ax-4h-nL-nD-rep1", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-4h-L-nD-rep1.bam", bs_a + "gh2ax-4h-L-nD-rep1", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-4h-L-PKi-rep1.bam", bs_a + "gh2ax-4h-L-PKi-rep1", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-1h-nL-nD-rep2.bam", bs_a + "gh2ax-1h-nL-nD-rep2", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-4h-nL-nD-rep2.bam", bs_a + "gh2ax-4h-nL-nD-rep2", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-4h-L-nD-rep2.bam", bs_a + "gh2ax-4h-L-nD-rep2", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-4h-L-PKi-rep2.bam", bs_a + "gh2ax-4h-L-PKi-rep2", win, chrs)

c.to_wiggle_windows(bs + "53bp1-wt.bam", bs_a + "53bp1-wt", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-wt-sub1.bam", bs_a + "gh2ax-wt-sub1", win, chrs)
c.to_wiggle_windows(bs + "gh2ax-wt-sub2.bam", bs_a + "gh2ax-wt-sub2", win, chrs)


""" For each window span, count number of reads in each bin. """
win = 50000
numbins = 50
c.to_bins(bs + "53bp1-1h-nL-nD-rep1.bam", bs_a + "53bp1-1h-nL-nD-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-4h-nL-nD-rep1.bam", bs_a + "53bp1-4h-nL-nD-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-4h-L-nD-rep1.bam", bs_a + "53bp1-4h-L-nD-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-4h-L-PKi-rep1.bam", bs_a + "53bp1-4h-L-PKi-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-1h-nL-nD-rep2.bam", bs_a + "53bp1-1h-nL-nD-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-4h-nL-nD-rep2.bam", bs_a + "53bp1-4h-nL-nD-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-4h-L-nD-rep2.bam", bs_a + "53bp1-4h-L-nD-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-4h-L-PKi-rep2.bam", bs_a + "53bp1-4h-L-PKi-rep2", win, numbins, chrs)

c.to_bins(bs + "gh2ax-1h-nL-nD-rep1.bam", bs_a + "gh2ax-1h-nL-nD-rep1", win, numbins, chrs)
c.to_bins(bs + "gh2ax-4h-nL-nD-rep1.bam", bs_a + "gh2ax-4h-nL-nD-rep1", win, numbins, chrs)
c.to_bins(bs + "gh2ax-4h-L-nD-rep1.bam", bs_a + "gh2ax-4h-L-nD-rep1", win, numbins, chrs)
c.to_bins(bs + "gh2ax-4h-L-PKi-rep1.bam", bs_a + "gh2ax-4h-L-PKi-rep1", win, numbins, chrs)
c.to_bins(bs + "gh2ax-1h-nL-nD-rep2.bam", bs_a + "gh2ax-1h-nL-nD-rep2", win, numbins, chrs)
c.to_bins(bs + "gh2ax-4h-nL-nD-rep2.bam", bs_a + "gh2ax-4h-nL-nD-rep2", win, numbins, chrs)
c.to_bins(bs + "gh2ax-4h-L-nD-rep2.bam", bs_a + "gh2ax-4h-L-nD-rep2", win, numbins, chrs)
c.to_bins(bs + "gh2ax-4h-L-PKi-rep2.bam", bs_a + "gh2ax-4h-L-PKi-rep2", win, numbins, chrs)

c.to_bins(bs + "53bp1-wt.bam", bs_a + "53bp1-wt", win, numbins, chrs)
c.to_bins(bs + "gh2ax-wt-sub1.bam", bs_a + "gh2ax-wt-sub1", win, numbins, chrs)
c.to_bins(bs + "gh2ax-wt-sub2.bam", bs_a + "gh2ax-wt-sub2", win, numbins, chrs)


""" Get the average 53BP1 and gH2AX peaks in RPM """
c.avgwig(bs_a + "53bp1-1h-nL-nD-rep1.wig", bs_a + "53bp1-1h-nL-nD-rep2.wig", bs_a + "53bp1-1h-nL-nD-avg")
c.avgwig(bs_a + "53bp1-4h-nL-nD-rep1.wig", bs_a + "53bp1-4h-nL-nD-rep2.wig", bs_a + "53bp1-4h-nL-nD-avg")
c.avgwig(bs_a + "53bp1-4h-L-nD-rep1.wig", bs_a + "53bp1-4h-L-nD-rep2.wig", bs_a + "53bp1-4h-L-nD-avg")
c.avgwig(bs_a + "53bp1-4h-L-PKi-rep1.wig", bs_a + "53bp1-4h-L-PKi-rep2.wig", bs_a + "53bp1-4h-L-PKi-avg")
c.avgwig(bs_a + "gh2ax-1h-nL-nD-rep1.wig", bs_a + "gh2ax-1h-nL-nD-rep2.wig", bs_a + "gh2ax-1h-nL-nD-avg")
c.avgwig(bs_a + "gh2ax-4h-nL-nD-rep1.wig", bs_a + "gh2ax-4h-nL-nD-rep2.wig", bs_a + "gh2ax-4h-nL-nD-avg")
c.avgwig(bs_a + "gh2ax-4h-L-nD-rep1.wig", bs_a + "gh2ax-4h-L-nD-rep2.wig", bs_a + "gh2ax-4h-L-nD-avg")
c.avgwig(bs_a + "gh2ax-4h-L-PKi-rep1.wig", bs_a + "gh2ax-4h-L-PKi-rep2.wig", bs_a + "gh2ax-4h-L-PKi-avg")
c.avgwig(bs_a + "gh2ax-wt-sub1.wig", bs_a + "gh2ax-wt-sub2.wig", bs_a + "gh2ax-wt-avg")


""" Calculate percent change from averaged 53BP1 and gH2AX peaks in RPM """
c.percentchange(bs_a + "53bp1-1h-nL-nD-avg.wig", bs_a + "53bp1-4h-nL-nD-avg.wig", bs_a + "53bp1-1hto4h-delta")
c.percentchange(bs_a + "53bp1-4h-nL-nD-avg.wig", bs_a + "53bp1-4h-L-nD-avg.wig", bs_a + "53bp1-4hnLtoL-delta")
c.percentchange(bs_a + "53bp1-4h-L-nD-avg.wig", bs_a + "53bp1-4h-L-PKi-avg.wig", bs_a + "53bp1-4hnDtoD-delta")
c.percentchange(bs_a + "gh2ax-1h-nL-nD-avg.wig", bs_a + "gh2ax-4h-nL-nD-avg.wig", bs_a + "gh2ax-1hto4h-delta")
c.percentchange(bs_a + "gh2ax-4h-nL-nD-avg.wig", bs_a + "gh2ax-4h-L-nD-avg.wig", bs_a + "gh2ax-4hnLtoL-delta")
c.percentchange(bs_a + "gh2ax-4h-L-nD-avg.wig", bs_a + "gh2ax-4h-L-PKi-avg.wig", bs_a + "gh2ax-4hnDtoD-delta")
