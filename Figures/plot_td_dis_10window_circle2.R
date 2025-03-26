library(ggplot2)
library(dplyr)
library(circlize)

input_file <- '~/Desktop/graph/MTD/distribution/fen1-24k_MTD_trf.tsv'
output_file <- '~/Desktop/graph/MTD/distribution/fen1-24k_MTD_trf_circle3.png'

# input_file <- '~/Desktop/graph/MTD/distribution/wt-24k_MTD_trf.tsv'
# output_file <- '~/Desktop/graph/MTD/distribution/wt-24k_MTD_trf_circle3.png'

# 染色体名称和长度
chromosomes <- data.frame(
  chr = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10", 
          "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6", 
          "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5", 
          "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", 
          "NC_001224.1"),
  length = c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 
             562643, 439888, 745751, 666816, 1078177, 924431, 784333, 
             1091291, 948066, 85779)
)

# Load TD data
td_data <- read.table(input_file, sep = "\t", header = FALSE)
colnames(td_data) <- c("chr", "start_before", "ref", "alt", "length", "seq_ref", "seq_alt", "score")

# Calculate the starting position of TD
td_data <- td_data %>%
  mutate(start = start_before + 1)

# Count the number of times each site is detected
td_counts <- td_data %>%
  group_by(chr, start) %>%
  summarise(count = n(), .groups = 'drop')

# Merge chromosome length information into TD data
td_counts <- td_counts %>%
  left_join(chromosomes, by = "chr")

# Calculate the cumulative starting position of each chromosome
chromosomes <- chromosomes %>%
  arrange(chr) %>%
  mutate(cumulative_length = cumsum(length),
         start_pos = lag(cumulative_length, default = 0))

# Merge the cumulative start position into TD data
td_counts <- td_counts %>%
  left_join(chromosomes %>% select(chr, start_pos), by = "chr") %>%
  mutate(genome_pos = start_pos + start)

# Set the window size
window_size <- 10

# Calculate the number of TDs per window
td_windows <- td_counts %>%
  mutate(window = floor(genome_pos / window_size)) %>%
  group_by(window) %>%
  summarise(td_count = sum(count), .groups = 'drop')

# Calculate the genome location of each window
td_windows <- td_windows %>%
  mutate(genome_pos = window * window_size + window_size / 2)

chromosome <- data.frame(chr = chromosomes$chr, start = chromosomes$start_pos, end = chromosomes$cumulative_length)

# plot
circos.clear()
png(output_file, width = 9, height = 9, units = "in", res = 300)
circos.par(track.margin=c(0,0),cell.padding = c(0, 0, 0, 0), gap.after = c(rep(1, nrow(chromosomes)-1),4))

circos.genomicInitialize(chromosome,plotType = c("axis"),axis.labels.cex = 0.6)

circos.track(ylim=c(0,1),panel.fun=function(x,y) {
    Genome=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
  
    if (Genome == "NC_001224.1") {
        x_position = mean(xlim) - 30000  # Move left
    } else if (Genome == "NC_001133.9") {
        x_position = mean(xlim) + 0  # Move to the right but will be blocked
        x_position = mean(xlim)  # Others keep default
    }
    circos.text(x_position,mean(ylim),Genome,cex=0.5,col='grey20',
                facing="bending.inside",niceFacing=TRUE)
},bg.col="grey90",bg.border=F,track.height=0.06)
circos.track(factors = chromosome$chr, ylim = c(7, 0), 
            panel.fun = function(x, y) {
            Genome = CELL_META$sector.index
            xlim = CELL_META$xlim
            ylim = CELL_META$ylim
            for (i in 0:7) {
                circos.lines(CELL_META$cell.xlim, c(i,i),lty = 1, lwd = 0.5, col = 'grey')
            }
            },
            bg.col = NA,
            bg.border = T,
            track.height = 0.6)
for (i in 1:7) {
  circos.text(12230000, i,labels = i, cex = 0.7, col = 'black')
}

td_counts_1 <- td_counts[td_counts$count == 1,]
td_counts_n <- td_counts[td_counts$count > 1,]
# can set different sizes for points with different counts for observation
circos.trackPoints(td_counts_1$chr, x = td_counts_1$genome_pos, y = td_counts_1$count,col = rgb(1,0,0,0.6), pch = 16, cex = 0.7)
circos.trackPoints(td_counts_n$chr, x = td_counts_n$genome_pos, y = td_counts_n$count,col = rgb(1,0,0,0.6), pch = 16, cex = 0.7)

title("Distribution of TDs in Yeast Genome (10bp Windows)", cex.main = 1.5, col.main = "black")

circos.clear()
