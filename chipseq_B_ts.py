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
c.to_wiggle_windows(bs + "53bp1-15m-rep1.bam", bs_a + "53bp1-15m-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-30m-rep1.bam", bs_a + "53bp1-30m-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-45m-rep1.bam", bs_a + "53bp1-45m-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-60m-rep1.bam", bs_a + "53bp1-60m-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-2h-rep1.bam", bs_a + "53bp1-2h-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-wt-ts-rep1.bam", bs_a + "53bp1-wt-ts-rep1", win, chrs)
c.to_wiggle_windows(bs + "53bp1-00m-rep2.bam", bs_a + "53bp1-00m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-15m-rep2.bam", bs_a + "53bp1-15m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-30m-rep2.bam", bs_a + "53bp1-30m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-45m-rep2.bam", bs_a + "53bp1-45m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-60m-rep2.bam", bs_a + "53bp1-60m-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-2h-rep2.bam", bs_a + "53bp1-2h-rep2", win, chrs)
c.to_wiggle_windows(bs + "53bp1-wt-ts-rep2.bam", bs_a + "53bp1-wt-ts-rep2", win, chrs)


""" For each window span, count number of reads in each bin. """
win = 50000
numbins = 50
c.to_bins(bs + "53bp1-00m-rep1.bam", bs_a + "53bp1-00m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-15m-rep1.bam", bs_a + "53bp1-15m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-30m-rep1.bam", bs_a + "53bp1-30m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-45m-rep1.bam", bs_a + "53bp1-45m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-60m-rep1.bam", bs_a + "53bp1-60m-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-2h-rep1.bam", bs_a + "53bp1-2h-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-wt-ts-rep1.bam", bs_a + "53bp1-wt-ts-rep1", win, numbins, chrs)
c.to_bins(bs + "53bp1-00m-rep2.bam", bs_a + "53bp1-00m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-15m-rep2.bam", bs_a + "53bp1-15m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-30m-rep2.bam", bs_a + "53bp1-30m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-45m-rep2.bam", bs_a + "53bp1-45m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-60m-rep2.bam", bs_a + "53bp1-60m-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-2h-rep2.bam", bs_a + "53bp1-2h-rep2", win, numbins, chrs)
c.to_bins(bs + "53bp1-wt-ts-rep2.bam", bs_a + "53bp1-wt-ts-rep2", win, numbins, chrs)


""" Perform T-test on bins by comparing each time point to wild type """
c.ttest_two(bs_a + "53bp1-00m-rep1.csv", bs_a + "53bp1-wt-ts-rep1.csv", bs_a + "53bp1-00m-rep1", p=0.05)
c.ttest_two(bs_a + "53bp1-15m-rep1.csv", bs_a + "53bp1-wt-ts-rep1.csv", bs_a + "53bp1-15m-rep1", p=0.05)
c.ttest_two(bs_a + "53bp1-30m-rep1.csv", bs_a + "53bp1-wt-ts-rep1.csv", bs_a + "53bp1-30m-rep1", p=0.05)
c.ttest_two(bs_a + "53bp1-45m-rep1.csv", bs_a + "53bp1-wt-ts-rep1.csv", bs_a + "53bp1-45m-rep1", p=0.05)
c.ttest_two(bs_a + "53bp1-60m-rep1.csv", bs_a + "53bp1-wt-ts-rep1.csv", bs_a + "53bp1-60m-rep1", p=0.05)
c.ttest_two(bs_a + "53bp1-2h-rep1.csv", bs_a + "53bp1-wt-ts-rep1.csv", bs_a + "53bp1-2h-rep1", p=0.05)
c.ttest_two(bs_a + "53bp1-00m-rep2.csv", bs_a + "53bp1-wt-ts-rep2.csv", bs_a + "53bp1-00m-rep2", p=0.05)
c.ttest_two(bs_a + "53bp1-15m-rep2.csv", bs_a + "53bp1-wt-ts-rep2.csv", bs_a + "53bp1-15m-rep2", p=0.05)
c.ttest_two(bs_a + "53bp1-30m-rep2.csv", bs_a + "53bp1-wt-ts-rep2.csv", bs_a + "53bp1-30m-rep2", p=0.05)
c.ttest_two(bs_a + "53bp1-45m-rep2.csv", bs_a + "53bp1-wt-ts-rep2.csv", bs_a + "53bp1-45m-rep2", p=0.05)
c.ttest_two(bs_a + "53bp1-60m-rep2.csv", bs_a + "53bp1-wt-ts-rep2.csv", bs_a + "53bp1-60m-rep2", p=0.05)
c.ttest_two(bs_a + "53bp1-2h-rep2.csv", bs_a + "53bp1-wt-ts-rep2.csv", bs_a + "53bp1-2h-rep2", p=0.05)


""" Calculate the width of 53BP1 peaks in broadPeak format. """
names = ["ACTB", "MYC"]
cuts = [5529660, 127736258]
span = 500000
c.ttest_span(bs_a + "53bp1-00m-rep1_ttest.csv", bs_a + "53bp1-00m-rep1_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-15m-rep1_ttest.csv", bs_a + "53bp1-15m-rep1_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-30m-rep1_ttest.csv", bs_a + "53bp1-30m-rep1_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-45m-rep1_ttest.csv", bs_a + "53bp1-45m-rep1_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-60m-rep1_ttest.csv", bs_a + "53bp1-60m-rep1_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-2h-rep1_ttest.csv", bs_a + "53bp1-2h-rep1_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-00m-rep2_ttest.csv", bs_a + "53bp1-00m-rep2_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-15m-rep2_ttest.csv", bs_a + "53bp1-15m-rep2_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-30m-rep2_ttest.csv", bs_a + "53bp1-30m-rep2_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-45m-rep2_ttest.csv", bs_a + "53bp1-45m-rep2_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-60m-rep2_ttest.csv", bs_a + "53bp1-60m-rep2_span", chrs, cuts, names, span)
c.ttest_span(bs_a + "53bp1-2h-rep2_ttest.csv", bs_a + "53bp1-2h-rep2_span", chrs, cuts, names, span)


""" Get span average """
c.avgspan(bs_a + "53bp1-00m-rep1_span.broadPeak", bs_a + "53bp1-00m-rep2_span.broadPeak", bs_a + "53bp1-00m-span")
c.avgspan(bs_a + "53bp1-30m-rep1_span.broadPeak", bs_a + "53bp1-30m-rep2_span.broadPeak", bs_a + "53bp1-30m-span")
c.avgspan(bs_a + "53bp1-60m-rep1_span.broadPeak", bs_a + "53bp1-60m-rep2_span.broadPeak", bs_a + "53bp1-60m-span")
c.avgspan(bs_a + "53bp1-2h-rep1_span.broadPeak", bs_a + "53bp1-2h-rep2_span.broadPeak", bs_a + "53bp1-2h-span")


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
