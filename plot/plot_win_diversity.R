library(scales)

# receive cmd arguments
Args <- commandArgs(T)
wd <- Args[1]
dd <- Args[2]
para <- Args[3]

# # parameters for test only
# wd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/window_diversity/plot"
# dd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/window_diversity/results"
# para <- "winlen_125000_mincov_12_mincount_1"
# chr <- "3R"

# parameter specified inside script
# selected.range <- names(vc_color)
selected.range <- c("Chinese", "European")

# def a function to plot dots for window nucleotide diversity
plot_dot <- function(wd, dd, para, chr, selected.range) {

## read the table including Fst
setwd(wd)
df <- read.table(file=paste(dd, "/", para, "/", "diversity_", chr, "_sorted.txt", sep = ""), 
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
mt.out <- matrix(, nrow = 0, ncol = dim(df)[2])
win.lowest.div <- c()
df$meandiv <- apply(df[,-1:-3], 1, function(x){mean(x)})
for (contig in unique(df$contig)) {
  df.contig <- df[df$contig == contig, ]
  if (dim(df.contig)[1] >= 20) {
    df.contig.mid <- df.contig[df.contig$window > .1*length(df.contig$window) & df.contig$window < .9*length(df.contig$window), ]
    meandiv.lowest <- df.contig.mid[order(df.contig.mid$meandiv), ][1,]$meandiv
    win.lowest.div <- c(win.lowest.div, 
                        rownames(df.contig.mid[df.contig.mid$meandiv <= 1.05*meandiv.lowest, ]))
    mt.out <- rbind(mt.out, df.contig.mid[df.contig.mid$meandiv <= 1.05*meandiv.lowest, ])
  }
}
# rownames(df.contig.mid[order(df.contig.mid$meandiv), ][1:(0.05*dim(df.contig.mid)[1]), ]))

# output the lowest-diversity windows to a text file
write.table(mt.out, file = paste(wd, '/', "window_lowest_diversity.txt", sep = ""),
            quote = F, sep = "\t", row.names = F, append = T, col.names = F)

## plot the dot plot
pdf(file=paste(wd, '/', "diversity_", chr, "_", para, ".pdf", sep = ""))
par(ps = 12, cex = 1, cex.main = 1)
selected.pop <- df_color[df_color$range %in% selected.range, 1]
for (sample in selected.pop) {
  if (sample == selected.pop[1]){
    plot(x = 1:nrow(df), y = df[, sample], 
         col = alpha(df_color[df_color$sample == sample, 3], alpha = 0.3), pch = 19, cex = 0.5,
         ylim = c(0, 0.05), xlab = "Windows", ylab = "Nucleotide Diversity", main = chr,
         panel.first = rect(xleft, 0 , xright, 1e10, col='lightgrey', border=NA))
  }
  else {
    points(x = 1:nrow(df),  y = df[, sample], col = alpha(df_color[df_color$sample == sample, 3], alpha = 0.3), pch = 19, cex = 0.5)
  }
}
abline(v = win.lowest.div, lty = "dotted", lwd = 0.8)
legend("topright", legend=selected.range,title="Range", cex = 0.8,
       fill=vc_color[selected.range], bty="n", xpd = T)
dev.off()
}

# plot
plot_dot(wd, dd, para, "2L", selected.range)
plot_dot(wd, dd, para, "2R", selected.range)
plot_dot(wd, dd, para, "3L", selected.range)
plot_dot(wd, dd, para, "3R", selected.range)
plot_dot(wd, dd, para, "4", selected.range)
plot_dot(wd, dd, para, "X", selected.range)
