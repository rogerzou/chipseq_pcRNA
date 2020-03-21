# -*- coding: utf-8 -*-
""" MRE11 ChIP-seq analysis for Cas9/pcRNA - DNA-PKcs inhibitor effect """

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
bs = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA1/"                # directory with input files
bs_s = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/pcRNA_SRA1/subsets/"      # output directory

rr_A = "chr7:5527160-5532160"           # 5kb window centered at ACTB cut site
tar_A = 5529660                         # ACTB cleavage site
rr_M = "chr8:127733758-127738758"       # 5kb window centered at MYC cut site
tar_M = 127736258                       # MYC cleavage site


""" Subset BAM files (in 5kb window centered at cut site) by characteristics of each paired-end read
    in relation to the cut site. """
""" ACTB """
c.get_read_subsets(bs + "mre11-1h-nL-nD-rep1.bam", bs_s + "mre11-1h-nL-nD-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-4h-nL-nD-rep1.bam", bs_s + "mre11-4h-nL-nD-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-4h-L-nD-rep1.bam", bs_s + "mre11-4h-L-nD-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-4h-L-PKi-rep1.bam", bs_s + "mre11-4h-L-PKi-rep1_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-1h-nL-nD-rep2.bam", bs_s + "mre11-1h-nL-nD-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-4h-nL-nD-rep2.bam", bs_s + "mre11-4h-nL-nD-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-4h-L-nD-rep2.bam", bs_s + "mre11-4h-L-nD-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-4h-L-PKi-rep2.bam", bs_s + "mre11-4h-L-PKi-rep2_ACTB", rr_A, tar_A)
c.get_read_subsets(bs + "mre11-wt.bam", bs_s + "mre11-wt_ACTB", rr_A, tar_A)
""" MYC """
c.get_read_subsets(bs + "mre11-1h-nL-nD-rep1.bam", bs_s + "mre11-1h-nL-nD-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-4h-nL-nD-rep1.bam", bs_s + "mre11-4h-nL-nD-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-4h-L-nD-rep1.bam", bs_s + "mre11-4h-L-nD-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-4h-L-PKi-rep1.bam", bs_s + "mre11-4h-L-PKi-rep1_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-1h-nL-nD-rep2.bam", bs_s + "mre11-1h-nL-nD-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-4h-nL-nD-rep2.bam", bs_s + "mre11-4h-nL-nD-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-4h-L-nD-rep2.bam", bs_s + "mre11-4h-L-nD-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-4h-L-PKi-rep2.bam", bs_s + "mre11-4h-L-PKi-rep2_MYC", rr_M, tar_M)
c.get_read_subsets(bs + "mre11-wt.bam", bs_s + "mre11-wt_MYC", rr_M, tar_M)


""" Convert all fragments (in 5kb window centered at cut site) to wiggle format. """
""" ACTB """
c.to_wiggle_pairs(bs + "mre11-1h-nL-nD-rep1.bam", bs_s + "mre11-1h-nL-nD-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-4h-nL-nD-rep1.bam", bs_s + "mre11-4h-nL-nD-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-4h-L-nD-rep1.bam", bs_s + "mre11-4h-L-nD-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-4h-L-PKi-rep1.bam", bs_s + "mre11-4h-L-PKi-rep1_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-1h-nL-nD-rep2.bam", bs_s + "mre11-1h-nL-nD-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-4h-nL-nD-rep2.bam", bs_s + "mre11-4h-nL-nD-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-4h-L-nD-rep2.bam", bs_s + "mre11-4h-L-nD-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-4h-L-PKi-rep2.bam", bs_s + "mre11-4h-L-PKi-rep2_ACTB", rr_A)
c.to_wiggle_pairs(bs + "mre11-wt.bam", bs_s + "mre11-wt_ACTB", rr_A)
""" MYC """
c.to_wiggle_pairs(bs + "mre11-1h-nL-nD-rep1.bam", bs_s + "mre11-1h-nL-nD-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-4h-nL-nD-rep1.bam", bs_s + "mre11-4h-nL-nD-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-4h-L-nD-rep1.bam", bs_s + "mre11-4h-L-nD-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-4h-L-PKi-rep1.bam", bs_s + "mre11-4h-L-PKi-rep1_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-1h-nL-nD-rep2.bam", bs_s + "mre11-1h-nL-nD-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-4h-nL-nD-rep2.bam", bs_s + "mre11-4h-nL-nD-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-4h-L-nD-rep2.bam", bs_s + "mre11-4h-L-nD-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-4h-L-PKi-rep2.bam", bs_s + "mre11-4h-L-PKi-rep2_MYC", rr_M)
c.to_wiggle_pairs(bs + "mre11-wt.bam", bs_s + "mre11-wt_MYC", rr_M)


""" Convert spanning fragments to wiggle format. """
""" ACTB """
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep1_ACTB_M.bam", bs_s + "mre11-1h-nL-nD-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep1_ACTB_M.bam", bs_s + "mre11-4h-nL-nD-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep1_ACTB_M.bam", bs_s + "mre11-4h-L-nD-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep1_ACTB_M.bam", bs_s + "mre11-4h-L-PKi-rep1_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep2_ACTB_M.bam", bs_s + "mre11-1h-nL-nD-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep2_ACTB_M.bam", bs_s + "mre11-4h-nL-nD-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep2_ACTB_M.bam", bs_s + "mre11-4h-L-nD-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep2_ACTB_M.bam", bs_s + "mre11-4h-L-PKi-rep2_ACTB_M", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-wt_ACTB_M.bam", bs_s + "mre11-wt_ACTB_M", rr_A)
""" MYC """
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep1_MYC_M.bam", bs_s + "mre11-1h-nL-nD-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep1_MYC_M.bam", bs_s + "mre11-4h-nL-nD-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep1_MYC_M.bam", bs_s + "mre11-4h-L-nD-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep1_MYC_M.bam", bs_s + "mre11-4h-L-PKi-rep1_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep2_MYC_M.bam", bs_s + "mre11-1h-nL-nD-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep2_MYC_M.bam", bs_s + "mre11-4h-nL-nD-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep2_MYC_M.bam", bs_s + "mre11-4h-L-nD-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep2_MYC_M.bam", bs_s + "mre11-4h-L-PKi-rep2_MYC_M", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-wt_MYC_M.bam", bs_s + "mre11-wt_MYC_M", rr_M)


""" Convert fragments that start/end 5bp away from cut site to wiggle format. """
""" ACTB """
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep1_ACTB_N.bam", bs_s + "mre11-1h-nL-nD-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep1_ACTB_N.bam", bs_s + "mre11-4h-nL-nD-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep1_ACTB_N.bam", bs_s + "mre11-4h-L-nD-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep1_ACTB_N.bam", bs_s + "mre11-4h-L-PKi-rep1_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep2_ACTB_N.bam", bs_s + "mre11-1h-nL-nD-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep2_ACTB_N.bam", bs_s + "mre11-4h-nL-nD-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep2_ACTB_N.bam", bs_s + "mre11-4h-L-nD-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep2_ACTB_N.bam", bs_s + "mre11-4h-L-PKi-rep2_ACTB_N", rr_A)
c.to_wiggle_pairs(bs_s + "mre11-wt_ACTB_N.bam", bs_s + "mre11-wt_ACTB_N", rr_A)
""" MYC """
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep1_MYC_N.bam", bs_s + "mre11-1h-nL-nD-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep1_MYC_N.bam", bs_s + "mre11-4h-nL-nD-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep1_MYC_N.bam", bs_s + "mre11-4h-L-nD-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep1_MYC_N.bam", bs_s + "mre11-4h-L-PKi-rep1_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-1h-nL-nD-rep2_MYC_N.bam", bs_s + "mre11-1h-nL-nD-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-nL-nD-rep2_MYC_N.bam", bs_s + "mre11-4h-nL-nD-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-nD-rep2_MYC_N.bam", bs_s + "mre11-4h-L-nD-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-4h-L-PKi-rep2_MYC_N.bam", bs_s + "mre11-4h-L-PKi-rep2_MYC_N", rr_M)
c.to_wiggle_pairs(bs_s + "mre11-wt_MYC_N.bam", bs_s + "mre11-wt_MYC_N", rr_M)
