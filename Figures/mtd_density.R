library(ggplot2)
library(tidyr)
setwd("~/Desktop/graph/MTD/")

# load data
raw_data <- openxlsx::read.xlsx("mtd-figure4.xlsx", sheet = 2, colNames = FALSE)
df <- raw_data[2:3,2:ncol(raw_data) ]
df <- as.data.frame(t(df))
colnames(df) <- c("species", "nMTD")
df$nMTD <- as.numeric(df$nMTD)
df$species <- factor(df$species, levels = unique(df$species)) 

colors <- c(rep(rgb(42,54,130, maxColorValue = 255),8),rgb(16,124,54, maxColorValue = 255),rep(rgb(230,0,29, maxColorValue = 255),8))

p <- ggplot(df, aes(x = species, y = nMTD)) +
    geom_bar(stat = "identity", fill = colors, width = 0.8) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 200), breaks = seq(0, 200, 20)) +
    theme_classic() +
    coord_fixed(0.035) + 
    theme(
        legend.position = "none",
        axis.text = element_text(size = 21),
        axis.title = element_text(size = 21),
        legend.text = element_text(size = 18)
    ) +
    coord_cartesian(clip = "off") + 
    theme(
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        axis.line = element_line(size = 2),
        axis.text.x = element_text(size = 24, face = "bold.italic", color = "black", angle = 90, hjust = 1),
        axis.text.y = element_text(size = 24, face = "bold", color = "black"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 30, face = "bold", color = "black"),
        strip.text.x = element_text(size = 17),
        plot.margin = unit(c(2, 2, 2, 2), "cm")
    ) +
    ylab("MTD per million bases") 
    p

ggsave(p, filename = "Figure4/mtd_density.png", width = 11, height = 12)
