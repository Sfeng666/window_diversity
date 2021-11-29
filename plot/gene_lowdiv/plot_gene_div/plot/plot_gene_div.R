library(scales)
library(karyoploteR)

# receive cmd arguments
Args <- commandArgs(T)
wd <- Args[1]
dd <- Args[2]
para <- Args[3]

# # parameters for test only
# wd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/window_diversity/plot/gene_lowdiv/plot_gene_div/plot"
# dd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/window_diversity/plot/gene_lowdiv/plot_gene_div/results"
# para <- "winlen_1000_mincov_12_mincount_1"
# chr <- "2R"
# contig <- 'NW_023496807.1'

# parameter specified inside script
# selected.range <- names(vc_color)
selected.range <- c("Chinese", "European")

# def a function to plot dots for window nucleotide diversity
plot_dot <- function(wd, dd, para, chr, contig, selected.range) {
  
  ## read the table including Fst
  df <- read.table(file=paste(dd, "/", para, "/", "diversity_", chr, '_', contig, ".txt", sep = ""), 
                   header=T, sep="\t", check.names=F)
  df_color <- data.frame(sample=colnames(df[,-1:-3]), 
                         range=factor(c("FR-Run", "European", "European", "Chinese",
                                        "Chinese", "European", "European", "European",
                                        "European", "Japanese", "Japanese", "Chinese", 
                                        "American", "American", "Chinese", "American", 
                                        "American", "American", "US-Haw", "American",
                                        "BR-Pal", "European", "American", "American",
                                        "Japanese", "European", "Japanese", "European", "European"),
                                      levels = c("US-Haw", "American", "BR-Pal", "FR-Run", "European", "Chinese", "Japanese")))
  df_color <- df_color[order(df_color$sample),]
  df_color <- df_color[order(df_color$range),]
  df_color$range <- as.character(df_color$range)
  
  # update 10.09.2021: separate island populations from the mainland population ranges, and assign separate color to them
  vc_color <- c(rgb(232, 125,	114, maxColorValue = 255), rgb(83, 182,76, maxColorValue = 255), rgb(109, 157, 248, maxColorValue = 255), rgb(109 - 35, 157 + 35, 248, maxColorValue = 255),
                rgb(232, 125,	114 + 70, maxColorValue = 255), rgb(232 -50, 125,	114 + 20, maxColorValue = 255), rgb(83 + 70, 182,76, maxColorValue = 255))
  names(vc_color) <- c("American", "European", "Chinese", "Japanese", "US-Haw", "BR-Pal","FR-Run")
  vc_spcolor <- vc_color[df_color$range]
  df_color$color <- vc_spcolor
  
  # determine coordinates of the shaded background
  coord <- c(0, row.names(df[df$window == 0,][-1,]))
  if (length(coord) %% 2 == 1) {
    coord <- c(coord, 1e10)
  }
  xleft <- coord[seq(1, length(coord), 2) + 1]
  xright <- coord[seq(2, length(coord), 2) + 1]
  
  # show the window at the middle of each contig (longer than 50 windows) that has the lowest average diversity across samples
  # mt.out <- matrix(, nrow = 0, ncol = dim(df)[2])
  win.lowest.div <- c()
  df$meandiv <- apply(df[,-1:-3], 1, function(x){mean(x)})
  for (contig in unique(df$contig)) {
    df.contig <- df[df$contig == contig, ]
    if (dim(df.contig)[1] >= 20) {
      df.contig.mid <- df.contig[df.contig$window > .1*length(df.contig$window) & df.contig$window < .9*length(df.contig$window), ]
      meandiv.lowest <- df.contig.mid[order(df.contig.mid$meandiv), ][1,]$meandiv
      win.lowest.div <- c(win.lowest.div, 
                          df.contig.mid[df.contig.mid$meandiv <= 1.2*meandiv.lowest, "start_coord"])
      # mt.out <- rbind(mt.out, df.contig.mid[df.contig.mid$meandiv <= 1.2*meandiv.lowest, ])
    }
  }

  # # output the lowest-diversity windows to a text file
  # write.table(mt.out, file = paste(wd, '/', "window_lowest_diversity.txt", sep = ""),
  #             quote = F, sep = "\t", row.names = F, append = T, col.names = F)
  
  ## initialize the multi-lane setting
  pdf(file=paste(wd, '/', "gene_lowdivwin_", chr, "_", contig, ".pdf", sep = ""))
  par(ps = 9, cex = 1, cex.main = 1)
  par(mfrow=c(2, 1), mai = c(1, 0.8, 0.5, 0.2))
  # par(mfrow=c(2, 1))

  # plot the dot plot
  par(fig=c(0,1,0,0.65))
  selected.pop <- df_color[df_color$range %in% selected.range, 1]
  for (sample in selected.pop) {
    if (sample == selected.pop[1]){
      plot(x = df$start_coord, y = df[, sample], 
           col = alpha(df_color[df_color$sample == sample, 3], alpha = 0.3), 
           pch = 19, cex = 0.2, cex.axis = 0.9, ylim = c(0, 0.055), bty = 'l',
           xlab = "Genomic Region (bp)", ylab = "Nucleotide Diversity",
           panel.first = rect(xleft, 0 , xright, 1e10, col='lightgrey', border=NA))
    }
    else {
      points(x = df$start_coord,  y = df[, sample], col = alpha(df_color[df_color$sample == sample, 3], alpha = 0.3), pch = 19, cex = 0.5)
    }
  }
  # win.lowest.div <- c(win.lowest.div, 
  #                     8869702, df.contig[1, "start_coord"], 8789981, 8724867)
  # abline(v = win.lowest.div, lty = "dashed", lwd = 0.2, xpd = T)
  legend("topright", legend=selected.range,title="Range", cex = 0.8,
         fill=vc_color[selected.range], bty="n", xpd = T)
  
  # plot genes
  par(new = T)
  genomic.region <- paste(wd, "/gene_lowdiv/", "genomic_region_", chr, '_', contig, ".txt", sep = "")
  df.genomic.region <- read.table(file = genomic.region, 
                                  header = T, sep = '\t', check.names = F)
  custom.genome <- toGRanges(genomic.region)
  pp <- getDefaultPlotParams(plot.type=1)
  pp$bottommargin <- 650
  pp$topmargin <- 140
  pp$data1height <- 250
  pp$data1inmargin <- 10
  # pp$dataloutmargin <- 5
  pp$ideogramheight <- 10
  pp$rightmargin <- 0.0509
  # pp$rightmargin <- 0.0865
  pp$leftmargin <- 0.1439
  kp <- plotKaryotype(genome = custom.genome, labels.plotter=NULL, 
                      plot.params = pp, plot.type = 1, lwd = 0.3)
  kpAddChromosomeNames(kp, chr.names = paste(df.genomic.region$start, df.genomic.region$end, sep = " - "), cex = 0.8)
  title(main = paste(chr, contig, sep = ":"), line = -2, cex.main = 0.8)
  gene.coords <- paste(wd, "/gene_lowdiv/", "gene_coords_", chr, '_', contig, ".txt", sep = "")
  df.genes <- read.table(file=gene.coords, 
                         header=F, sep="\t", check.names=F)
  data.genes <- toGRanges(gene.coords)
  kpPlotRegions(kp, data=paste(wd, "/gene_lowdiv/", "gene_coords_", chr, '_', contig, ".txt", sep = ""), 
                col="#999999", border="black", r0=0.02, r1=0.15, lwd = 0.2)
  kpPlotMarkers(kp, data.genes, labels = df.genes$V5, ignore.chromosome.ends = T,
                adjust.label.position = T, label.dist = 0.001, 
                text.orientation = "vertical", pos = 4,
                cex = 0.8, r0=0.09, r1=0.6, max.iter = 1000, lwd = 0.2)
  dev.off()
}

# plot
plot_dot(wd, dd, para, "2L", 'NW_023496800.1', selected.range)
plot_dot(wd, dd, para, "2R", 'NW_023496807.1', selected.range)
plot_dot(wd, dd, para, "2R", 'NW_023496808.1', selected.range)
plot_dot(wd, dd, para, "3R", 'NW_023496837.1',selected.range)
plot_dot(wd, dd, para, "X", 'NW_023496846.1',selected.range)
