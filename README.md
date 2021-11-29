# Calculate and plot nucleotide diversity in windows across each chromosome

----
This program calculates amd plots nucleotide diversity in windows across each chromosome, in order to indirectly reflect regions of different recombination rates.

## steps
1. Devide each contig into windows using mpileup file of raw allele count;
2. Within each window, calculate nucleotide diversity;
3. plot the distribution of window diversity of each chromosome arm for selected samples;
4. label the window of the lowest diversity within each contig;
5. gernerate lists of gene symbols of genes in those windows.
6. plot genes within those windows along with window diversity.

## Notes
1. The window is defined as a continuous genomic region that include 250000 filtered sites.
2. We cannot confirm the order of contigs along each chromosome arms, but it seems ok to simply concatenate the windows from contigs that could be mapped to each chromosome arm, other wise many of those contigs will have only several windows, which is less informative. Contigs of a chromosome arm are ordered by their total length.
3. The position of each window nucleotide diversity is not shown in the plot, because a continuous axis of genomic coordinates cannot be applied to concatenated contigs without knowing the order. Instead, the windows are only shown in the following order: (contig_1: window_1, window_2 ... window_n) (contig_2: window_1, window_2 ... window_n) ... (contig_n: window_1, window_2 ... window_n). 

----
###  Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)
###  Mail: siyuanfeng.bioinfo@gmail.com