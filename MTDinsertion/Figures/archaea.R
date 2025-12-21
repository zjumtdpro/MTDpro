# Proportions of TDs with MHA at the inter-species level
library(ggplot2)
library(tidyverse)
setwd("~/Desktop/graph/MTD/")

# load data
raw_data <-  openxlsx::read.xlsx("mtd-figure2-3level.xlsx", sheet = 8, colNames = FALSE)

line <- 423 
data <- raw_data[2:line,2:5]
pvalue <- raw_data[1,8:11]
random <- as.numeric(unlist(raw_data[2,8:11]))*100
x_axis <- raw_data[1,2:5]

colnames(data) <- x_axis
value <- data %>%
    gather(key = "length",value = "percentage")
value$length <- factor(value$length, levels = x_axis)

x_axis <- ifelse(grepl("^[0-9]", x_axis), paste0("=", x_axis), gsub("^>=", "≥", x_axis))

colors <- c('blue', 'red', '#01c101', 'purple')

p <- 
ggplot(data = value, aes(x = length, y = as.numeric(percentage),color=length)) +
    geom_boxplot(size = 1, outlier.shape = NA)+
    geom_jitter(position=position_jitterdodge(jitter.width = 1, 
                                            jitter.height = 0, 
                                            dodge.width = 1),
                                            size = 3.5, alpha = 0.35)+

    scale_color_manual(values=colors)+
    scale_fill_manual(values = colors)+
    scale_y_continuous(expand = c(0,0),limits = c(0, 100),breaks = seq(0, 100, 10))+
    theme_classic()+
    coord_fixed(0.035)+ 
    geom_segment(aes(x = as.numeric(length)-0.4, xend = as.numeric(length)+0.4,
     y = random[as.numeric(length)], yend = random[as.numeric(length)], color = length), 
                 linetype = "dashed", size = 1) + 
    theme(axis.title.x = element_blank(),legend.position="none",
    axis.text = element_text(size = 21), 
    axis.title = element_text(size = 21), 
    legend.text = element_text(size = 18))+
    coord_cartesian(clip = 'off')+ 
    theme(
        axis.ticks = element_line(size = 2.5),
        axis.ticks.length = unit(0.5, "cm"),
        axis.line = element_line(size = 2.5), 
        axis.text.x = element_text(size = 30, face = "bold",color = "black", angle = 0), 
        axis.text.y = element_text(size = 34, face = "bold",color = "black"), 
        axis.title.y = element_text(size = 40, face = "bold",color = "black"), 
        plot.margin = unit(c(2, 2, 2, 2), "cm")
    )+
    ylab("% of TDs with MHA")+
    scale_x_discrete(labels=x_axis)+
    annotation_custom(grob = grid::textGrob("(bp)", x = 0.95, y = -0.043, 
        just = "left", gp = grid::gpar(fontsize = 34, fontface = "bold"))) +
    annotation_custom(grob = grid::textGrob(sprintf("%.4f" ,pvalue[1]), x = 0.15, y = 1.05, 
         gp = grid::gpar(fontsize = 23, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(sprintf("%.4f" ,pvalue[2]), x = 0.38, y = 1.05, 
         gp = grid::gpar(fontsize = 23, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(sprintf("%.4f" ,pvalue[3]), x = 0.615, y = 1.05, 
         gp = grid::gpar(fontsize = 23, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(sprintf("%.4f" ,pvalue[4]), x = 0.845, y = 1.05, 
         gp = grid::gpar(fontsize = 23, fontface = "bold",color = "black")))
p

ggsave(p,filename = 'Figure2/archaea_plot.png',width = 11,height = 9)
