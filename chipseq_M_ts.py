# -*- coding: utf-8 -*-
""" MRE11 ChIP-seq analysis for Cas9/pcRNA - timeseries after deactivation """

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
bs = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA2/"                # directory with input files
bs_s = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA2/subsets/"      # output directory

rr_A = "chr7:5527160-5532160"           # 5kb window centered at ACTB cut site
tar_A = 5529660                         # ACTB cleavage site
rr_M = "chr8:127733758-127738758"       # 5kb window centered at MYC cut site
tar_M = 127736258                       # MYC cleavage site


""" Subset BAM files (in 5kb window centered at cut site) by characteristics of each paired-end read
    in relation to the cut site. """
""" ACTB """
c.get_read_subsets(bs + "mre11-00m-rep1.bam", bs_s + "mre11-00m-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-15m-rep1.bam", bs_s + "mre11-15m-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-30m-rep1.bam", bs_s + "mre11-30m-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-45m-rep1.bam", bs_s + "mre11-45m-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-60m-rep1.bam", bs_s + "mre11-60m-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-2h-rep1.bam", bs_s + "mre11-2h-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-wt-ts-rep1.bam", bs_s + "mre11-wt-ts-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-00m-rep2.bam", bs_s + "mre11-00m-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-15m-rep2.bam", bs_s + "mre11-15m-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-30m-rep2.bam", bs_s + "mre11-30m-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-45m-rep2.bam", bs_s + "mre11-45m-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-60m-rep2.bam", bs_s + "mre11-60m-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-2h-rep2.bam", bs_s + "mre11-2h-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-wt-ts-rep2.bam", bs_s + "mre11-wt-ts-rep2_ACTB", rr_A, tar_A)
""" MYC """
c.get_read_subsets(bs + "mre11-00m-rep1.bam", bs_s + "mre11-00m-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-15m-rep1.bam", bs_s + "mre11-15m-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-30m-rep1.bam", bs_s + "mre11-30m-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-45m-rep1.bam", bs_s + "mre11-45m-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-60m-rep1.bam", bs_s + "mre11-60m-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-2h-rep1.bam", bs_s + "mre11-2h-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-wt-ts-rep1.bam", bs_s + "mre11-wt-ts-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-00m-rep2.bam", bs_s + "mre11-00m-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-15m-rep2.bam", bs_s + "mre11-15m-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-30m-rep2.bam", bs_s + "mre11-30m-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-45m-rep2.bam", bs_s + "mre11-45m-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-60m-rep2.bam", bs_s + "mre11-60m-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-2h-rep2.bam", bs_s + "mre11-2h-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-wt-ts-rep2.bam", bs_s + "mre11-wt-ts-rep2_MYC", rr_M, tar_M)


""" Convert all fragments (in 5kb window centered at cut site) to wiggle format. """
""" ACTB """
c.to_wiggle_pairs(bs + "mre11-00m-rep1.bam", bs_s + "mre11-00m-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-15m-rep1.bam", bs_s + "mre11-15m-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-30m-rep1.bam", bs_s + "mre11-30m-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-45m-rep1.bam", bs_s + "mre11-45m-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-60m-rep1.bam", bs_s + "mre11-60m-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-2h-rep1.bam", bs_s + "mre11-2h-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-wt-ts-rep1.bam", bs_s + "mre11-wt-ts-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-00m-rep2.bam", bs_s + "mre11-00m-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-15m-rep2.bam", bs_s + "mre11-15m-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-30m-rep2.bam", bs_s + "mre11-30m-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-45m-rep2.bam", bs_s + "mre11-45m-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-60m-rep2.bam", bs_s + "mre11-60m-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-2h-rep2.bam", bs_s + "mre11-2h-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-wt-ts-rep2.bam", bs_s + "mre11-wt-ts-rep2_ACTB", rr_A)
""" MYC """
c.to_wiggle_pairs(bs + "mre11-00m-rep1.bam", bs_s + "mre11-00m-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-15m-rep1.bam", bs_s + "mre11-15m-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-30m-rep1.bam", bs_s + "mre11-30m-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-45m-rep1.bam", bs_s + "mre11-45m-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-60m-rep1.bam", bs_s + "mre11-60m-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-2h-rep1.bam", bs_s + "mre11-2h-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-wt-ts-rep1.bam", bs_s + "mre11-wt-ts-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-00m-rep2.bam", bs_s + "mre11-00m-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-15m-rep2.bam", bs_s + "mre11-15m-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-30m-rep2.bam", bs_s + "mre11-30m-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-45m-rep2.bam", bs_s + "mre11-45m-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-60m-rep2.bam", bs_s + "mre11-60m-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-2h-rep2.bam", bs_s + "mre11-2h-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-wt-ts-rep2.bam", bs_s + "mre11-wt-ts-rep2_MYC", rr_M)


""" Convert spanning fragments to wiggle format. """
""" ACTB """
c.to_wiggle_pairs(bs_s + "mre11-00m-rep1_ACTB_M.bam", bs_s + "mre11-00m-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep1_ACTB_M.bam", bs_s + "mre11-15m-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep1_ACTB_M.bam", bs_s + "mre11-30m-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep1_ACTB_M.bam", bs_s + "mre11-45m-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep1_ACTB_M.bam", bs_s + "mre11-60m-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep1_ACTB_M.bam", bs_s + "mre11-2h-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep1_ACTB_M.bam", bs_s + "mre11-wt-ts-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-00m-rep2_ACTB_M.bam", bs_s + "mre11-00m-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep2_ACTB_M.bam", bs_s + "mre11-15m-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep2_ACTB_M.bam", bs_s + "mre11-30m-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep2_ACTB_M.bam", bs_s + "mre11-45m-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep2_ACTB_M.bam", bs_s + "mre11-60m-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep2_ACTB_M.bam", bs_s + "mre11-2h-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep2_ACTB_M.bam", bs_s + "mre11-wt-ts-rep2_ACTB_M", rr_A)
""" MYC """
c.to_wiggle_pairs(bs_s + "mre11-00m-rep1_MYC_M.bam", bs_s + "mre11-00m-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep1_MYC_M.bam", bs_s + "mre11-15m-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep1_MYC_M.bam", bs_s + "mre11-30m-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep1_MYC_M.bam", bs_s + "mre11-45m-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep1_MYC_M.bam", bs_s + "mre11-60m-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep1_MYC_M.bam", bs_s + "mre11-2h-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep1_MYC_M.bam", bs_s + "mre11-wt-ts-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-00m-rep2_MYC_M.bam", bs_s + "mre11-00m-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep2_MYC_M.bam", bs_s + "mre11-15m-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep2_MYC_M.bam", bs_s + "mre11-30m-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep2_MYC_M.bam", bs_s + "mre11-45m-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep2_MYC_M.bam", bs_s + "mre11-60m-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep2_MYC_M.bam", bs_s + "mre11-2h-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep2_MYC_M.bam", bs_s + "mre11-wt-ts-rep2_MYC_M", rr_M)


""" Convert fragments that start/end 5bp away from cut site to wiggle format. """
""" ACTB """
c.to_wiggle_pairs(bs_s + "mre11-00m-rep1_ACTB_N.bam", bs_s + "mre11-00m-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep1_ACTB_N.bam", bs_s + "mre11-15m-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep1_ACTB_N.bam", bs_s + "mre11-30m-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep1_ACTB_N.bam", bs_s + "mre11-45m-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep1_ACTB_N.bam", bs_s + "mre11-60m-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep1_ACTB_N.bam", bs_s + "mre11-2h-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep1_ACTB_N.bam", bs_s + "mre11-wt-ts-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-00m-rep2_ACTB_N.bam", bs_s + "mre11-00m-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep2_ACTB_N.bam", bs_s + "mre11-15m-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep2_ACTB_N.bam", bs_s + "mre11-30m-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep2_ACTB_N.bam", bs_s + "mre11-45m-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep2_ACTB_N.bam", bs_s + "mre11-60m-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep2_ACTB_N.bam", bs_s + "mre11-2h-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep2_ACTB_N.bam", bs_s + "mre11-wt-ts-rep2_ACTB_N", rr_A)
""" MYC """
c.to_wiggle_pairs(bs_s + "mre11-00m-rep1_MYC_N.bam", bs_s + "mre11-00m-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep1_MYC_N.bam", bs_s + "mre11-15m-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep1_MYC_N.bam", bs_s + "mre11-30m-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep1_MYC_N.bam", bs_s + "mre11-45m-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep1_MYC_N.bam", bs_s + "mre11-60m-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep1_MYC_N.bam", bs_s + "mre11-2h-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep1_MYC_N.bam", bs_s + "mre11-wt-ts-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-00m-rep2_MYC_N.bam", bs_s + "mre11-00m-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-15m-rep2_MYC_N.bam", bs_s + "mre11-15m-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-30m-rep2_MYC_N.bam", bs_s + "mre11-30m-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-45m-rep2_MYC_N.bam", bs_s + "mre11-45m-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-60m-rep2_MYC_N.bam", bs_s + "mre11-60m-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-2h-rep2_MYC_N.bam", bs_s + "mre11-2h-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-wt-ts-rep2_MYC_N.bam", bs_s + "mre11-wt-ts-rep2_MYC_N", rr_M)
