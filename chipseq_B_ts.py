# -*- coding: utf-8 -*-
""" 53BP1 and gH2AX ChIP-seq analysis for Cas9/pcRNA - timeseries after deactivation """

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
bs = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA2/"
bs_a = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA2/analysis/"

chrs = ['chr7', 'chr8']


""" Convert BAM file to WIG file that counts the number of reads in each window span. """
win = 5000
c.to_wiggle_windows(bs + "53bp1-00m-rep1.bam", bs_a + "53bp1-00m-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-30m-rep1.bam", bs_a + "53bp1-30m-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-60m-rep1.bam", bs_a + "53bp1-60m-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-2h-rep1.bam", bs_a + "53bp1-2h-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-wt-ts-rep1.bam", bs_a + "53bp1-wt-ts-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-00m-rep2.bam", bs_a + "53bp1-00m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-30m-rep2.bam", bs_a + "53bp1-30m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-60m-rep2.bam", bs_a + "53bp1-60m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-2h-rep2.bam", bs_a + "53bp1-2h-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-wt-ts-rep2.bam", bs_a + "53bp1-wt-ts-rep2", win, chrs)


""" For each window span, count number of reads in each bin. """
win = 50000
numbins = 50
c.to_bins(bs + "53bp1-00m-rep1.bam", bs_a + "53bp1-00m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-30m-rep1.bam", bs_a + "53bp1-30m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-60m-rep1.bam", bs_a + "53bp1-60m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-2h-rep1.bam", bs_a + "53bp1-2h-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-wt-ts-rep1.bam", bs_a + "53bp1-wt-ts-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-00m-rep2.bam", bs_a + "53bp1-00m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-30m-rep2.bam", bs_a + "53bp1-30m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-60m-rep2.bam", bs_a + "53bp1-60m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-2h-rep2.bam", bs_a + "53bp1-2h-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-wt-ts-rep2.bam", bs_a + "53bp1-wt-ts-rep2", win, numbins, chrs)


""" Get the average 53BP1 peak in RPM """
c.avgwig(bs_a + "53bp1-00m-rep1.wig", bs_a + "53bp1-00m-rep2.wig", bs_a + "53bp1-00m-avg")
c.avgwig(bs_a + "53bp1-30m-rep1.wig", bs_a + "53bp1-30m-rep2.wig", bs_a + "53bp1-30m-avg")
c.avgwig(bs_a + "53bp1-60m-rep1.wig", bs_a + "53bp1-60m-rep2.wig", bs_a + "53bp1-60m-avg")
c.avgwig(bs_a + "53bp1-2h-rep1.wig", bs_a + "53bp1-2h-rep2.wig", bs_a + "53bp1-2h-avg")
c.avgwig(bs_a + "53bp1-wt-ts-rep1.wig", bs_a + "53bp1-wt-ts-rep2.wig", bs_a + "53bp1-wt-ts-avg")


""" Calculate percent change from averaged 53BP1 peak in RPM """
c.percentchange(bs_a + "53bp1-00m-avg.wig", bs_a + "53bp1-30m-avg.wig", bs_a + "53bp1-00m-delta")
c.percentchange(bs_a + "53bp1-30m-avg.wig", bs_a + "53bp1-60m-avg.wig", bs_a + "53bp1-30m-delta")
c.percentchange(bs_a + "53bp1-60m-avg.wig", bs_a + "53bp1-2h-avg.wig", bs_a + "53bp1-60m-delta")
c.percentchange(bs_a + "53bp1-2h-avg.wig", bs_a + "53bp1-wt-ts-avg.wig", bs_a + "53bp1-2h-delta")
