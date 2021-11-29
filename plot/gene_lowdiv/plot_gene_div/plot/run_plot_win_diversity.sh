#!/bin/bash

### 1. configure the conda enviroment ###
source /opt/miniconda3/bin/activate WGS_analysis

### 2. set parameters ###
rscript=plot/gene_lowdiv/plot_gene_div/plot/plot_gene_div.R
wd=plot/gene_lowdiv/plot_gene_div/plot
dd=plot/gene_lowdiv/plot_gene_div/results

### 3. run the plotting R script ###
Rscript $rscript $wd $dd winlen_1000_mincov_12_mincount_1 >> $wd/run_plot_win_diversity.log 2>&1