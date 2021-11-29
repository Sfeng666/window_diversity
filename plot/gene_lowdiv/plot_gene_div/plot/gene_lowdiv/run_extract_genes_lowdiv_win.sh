#!/bin/bash

wd=plot/gene_lowdiv/plot_gene_div/plot/gene_lowdiv
script=plot/gene_lowdiv/plot_gene_div/scripts/extract_genes_coord_lowdiv_win.py
in_win_stats=plot/gene_lowdiv/plot_gene_div/selected_windows.txt
in_anno=/home/siyuan/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.gff.gz
# in_anno="/Users/siyuansmac/Desktop/PhD research/Suzukii WGS analysis/dataset.nosync/annotation/Refseq/reannotate/in_parallel/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.gff.gz"
in_ortho=/home/siyuan/reference/WT3-2.0/refseq/ortholog/gene_ortholog_suzukii_melanogaster.txt
# in_ortho="/Users/siyuansmac/bioinfo/project/suzukii_WGS/reference/WT3-2.0/refseq/ortholog/gene_ortholog_suzukii_melanogaster.txt"
out_gene=$wd/list_gene_lowdiv_win.txt
conda_env=WGS_analysis

python $script \
--wd $wd \
--in_win_stats $in_win_stats \
--in_anno $in_anno \
--in_ortho $in_ortho \
--out_gene $out_gene \
> $wd/run_extract_genes_lowdiv_win.log 2>&1