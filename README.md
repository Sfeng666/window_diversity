# Calculate and plot nucleotide diversity in windows across each chromosome

----
This program calculates amd plots nucleotide diversity in windows across each chromosome, in order to indirectly reflect regions of different recombination rates.

## Notes
1. The window is simply defined to be 100kb from the first bp of each mappable config.
2. We cannot confirm the order of contigs along each chromosome arms, but it seems ok to simply concatenate the windows from contigs that could be mapped to each chromosome arm, other wise many of those contigs will have only several windows, which is less informative.
3. The position of each window nucleotide diversity is not shown in the plot, because a continuous axis of genomic coordinates cannot be applied to concatenated contigs without knowing the order. Instead, the windows are only shown in the following order: (contig_1: window_1, window_2 ... window_n) (contig_2: window_1, window_2 ... window_n) ... (contig_n: window_1, window_2 ... window_n). 

----
###  Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)
###  Mail: siyuanfeng.bioinfo@gmail.com