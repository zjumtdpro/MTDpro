# Proportion of unique MTD loci across population sequencing replicates for 12 species.
library(ggplot2)
library(tidyr)
library(ggsci)
setwd("~/Desktop/graph/MTD/")

raw_data <-  openxlsx::read.xlsx("mtd-figure4.xlsx", sheet = 3, colNames = FALSE)

data <- raw_data[3:(nrow(raw_data)-2),2:ncol(raw_data)]
x_axis <- raw_data[2,2:ncol(raw_data)]
colnames(data) <- x_axis
df <- data %>% gather(key = "species", value = "percentage")
df$percentage <- as.numeric(df$percentage)
df$species <- factor(df$species, levels = unique(df$species)) 

colors <- c(rep(rgb(42,54,130, maxColorValue = 255),6),rgb(16,124,54, maxColorValue = 255),rep(rgb(230,0,29, maxColorValue = 255),4),'grey70')


p <- ggplot(data = df, aes(x = species, y = percentage, color = species))+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.5, 
                                            jitter.height = 0, 
                                            dodge.width = 1),
                                            size = 3.5, alpha = 1)+
    geom_boxplot(size = 1, outlier.shape = NA, alpha = 0)+
    scale_color_manual(values=colors)+
    scale_fill_manual(values = colors)+
    scale_y_continuous(expand = c(0,0),limits = c(0, 100),breaks = seq(0, 100, 10))+
    scale_x_discrete(labels=x_axis)+
    scale_color_manual(values = colors)+
    theme_classic()+
    theme(axis.title.x = element_blank(),legend.position = "none",
    axis.text = element_text(size = 21), 
    axis.title = element_text(size = 21), 
    legend.text = element_text(size = 15),
    legend.key.height = unit(1, "cm"),
    axis.ticks = element_line(size = 1.5), 
    axis.ticks.length = unit(0.5, "cm"), 
    axis.line = element_line(size = 1.5), 
    axis.text.x = element_text(size = 22, face = "bold.italic",color = "black", angle = 50, vjust = 0.85, hjust = 0.9), 
    axis.title.y = element_text(size = 30, face = "bold",color = "black"), 
    plot.margin = unit(c(2, 2, 2, 2), "cm")
    )+
    coord_cartesian(clip = 'off')+
    ylab("% of unique MTD")

p
ggsave(p, filename = "Figure4/unique.svg", width = 10, height = 11)
