# receive cmd arguments
Args <- commandArgs(T)
wd <- Args[1]
dd <- Args[2]
para <- Args[3]

# # parameters for test only
# wd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/window_diversity/plot"
# dd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/window_diversity/results"
# para <- "winheter_800_mincov_12_mincount_1"
# chr <- "2L"

# def a function to plot dots for window nucleotide diversity
plot_dot <- function(wd, dd, para, chr) {

## read the table including Fst
setwd(wd)
df <- read.table(file=paste(dd, "/", para, "/", "diversity_", chr, ".txt", sep = ""), 
                 header=T, sep="\t", check.names=FALSE)
df_color <- data.frame(sample=colnames(df), 
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

## plot the dot plot
pdf(file=paste(wd, '/', "diversity_", chr, "_", para, ".pdf", sep = ""))
par(ps = 12, cex = 1, cex.main = 1)
for (sample in colnames(df)) {
  if (sample == colnames(df)[1]){
    plot(x = 1:nrow(df), y = df[, which(colnames(df) == sample)], 
         col = "black", bg = df_color[df_color$sample == sample, 3], pch = 21, 
         ylim = c(0, 0.05), xlab = "Windows", ylab = "Nucleotide Diversity", main = chr)
  }
  else {
    points(x = 1:nrow(df),  y = df[, which(colnames(df) == sample)], col = "black", bg = df_color[df_color$sample == sample, 3], pch = 21)
  }
}
legend("topright", legend=names(vc_color),title="Range", cex = 0.8,
       fill=vc_color, bty="n", xpd = T)
dev.off()
}

# plot
plot_dot(wd, dd, para, "2L")
plot_dot(wd, dd, para, "2R")
plot_dot(wd, dd, para, "3L")
plot_dot(wd, dd, para, "3R")
plot_dot(wd, dd, para, "4")
plot_dot(wd, dd, para, "X")
