# plot genes within low-nucleotide-diversity windows

----
This program plot genes within low-nucleotide-diversity windows along with the distribution of nucleotide diversity. The purpose is to screen for specific genes that might be under selection. 

## steps
1. Generate a file containing windows of interest (low diversity) with coordinates;
2. Within each window of interest, calculate nucleotide diversity at the resolution of shorter windows (each including 1000 analyzed sites);
3. plot the distribution of window (small) diversity within each window of interest;
4. gernerate lists of gene symbols of genes that intersect with or are within those windows;
6. plot genes along with window diversity.

## Notes
1. The within-window window is defined as a continuous genomic region that include 1000 analyzed sites.
2. Since we are looking within contigs, genomic coordinates could be directly applied to the x axis.
3. The UCSC genome browser is not usable to track those genes, because they only support unannotated gene ID for Drosophila suzukii. Instead, the genomic coordinates of genes are used to plot the track of genes using the R package [KaryoploteR](https://bernatgel.github.io/karyoploter_tutorial). 
4. Some of the gene labels are missing from KaryoploteR (genes that intersect with the window border), and some overlap with other labels. These labels have been mannually added or adjusted, then saved with a file suffix '_adj'.

----
###  Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)
###  Mail: siyuanfeng.bioinfo@gmail.com