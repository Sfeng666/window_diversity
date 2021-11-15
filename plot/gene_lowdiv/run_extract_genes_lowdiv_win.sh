#!/bin/bash

wd=plot/gene_lowdiv
script=scripts/extract_genes_lowdiv_win.py
in_win_stats=results/winlen_125000_mincov_12_mincount_1/window_stats.txt
in_win_lowdiv=plot/window_lowest_diversity.txt
in_anno=/home/siyuan/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.gff.gz
# in_anno="/Users/siyuansmac/Desktop/PhD research/Suzukii WGS analysis/dataset.nosync/annotation/Refseq/reannotate/in_parallel/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.gff.gz"
in_ortho=/home/siyuan/reference/WT3-2.0/refseq/ortholog/gene_ortholog_suzukii_melanogaster.txt
# in_ortho="/Users/siyuansmac/bioinfo/project/suzukii_WGS/reference/WT3-2.0/refseq/ortholog/gene_ortholog_suzukii_melanogaster.txt"
out_gene=$wd/list_gene_lowdiv_win.txt
conda_env=WGS_analysis

python $script \
--wd $wd \
--in_win_stats $in_win_stats \
--in_win_lowdiv $in_win_lowdiv \
--in_anno $in_anno \
--in_ortho $in_ortho \
--out_gene $out_gene \
> $wd/run_extract_genes_lowdiv_win.log 2>&1