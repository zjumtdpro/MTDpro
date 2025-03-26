# Proportions of TDs across evolutionary scales.
librar2y(ggplot2)
library(tidyverse)
setwd("~/Desktop/graph/MTD/")

# load data
raw_data <-  openxlsx::read.xlsx("mtd-figure2-3level.xlsx", sheet = 2, colNames = FALSE)

line <- 2247
data <- raw_data[2:line,2:4]
pvalue <- raw_data[1,6:8]
random <- as.numeric(unlist(raw_data[2,6:8]))*100
x_axis <- raw_data[1,2:4]

colnames(data) <- x_axis
value <- data %>%
    gather(key = "group",value = "percentage")%>%
    drop_na()

value$group <- factor(value$group, levels = x_axis)
colors <- c('blue', 'red', '#01c101')

p <- 
ggplot(data = value, aes(x = group, y = as.numeric(percentage),color=group, fill = group)) +
    geom_jitter(position=position_jitterdodge(jitter.width = 1, 
                                            jitter.height = 0, 
                                            dodge.width = 1),
                                            size = 3.5, alpha = 0.35)+
    geom_boxplot(size = 1, outlier.shape = NA, alpha = 0.5)+
    scale_color_manual(values=colors)+
    scale_fill_manual(values = colors)+
    scale_y_continuous(expand = c(0,0),limits = c(0, 100),breaks = seq(0, 100, 10))+
    theme_classic()+
    coord_fixed(0.035)+ 

    geom_segment(aes(x = -Inf, xend = Inf, y = random[1], yend = random[1]),
                    linetype = "dashed", size = 1,color="black")+


    theme(axis.title.x = element_blank(),legend.position="none",
    axis.text = element_text(size = 21), 
    axis.title = element_text(size = 21), 
    legend.text = element_text(size = 18))+
    coord_cartesian(clip = 'off')+ 
    theme(
        axis.ticks = element_line(size = 2.5), 
        axis.ticks.length = unit(0.5, "cm"),    
        axis.line = element_line(size = 2.5),   
        axis.text.x = element_text(size = 30, face = "bold",color = "black", angle = 45, vjust = 0.85, hjust = 0.9),
        axis.text.y = element_text(size = 34, face = "bold",color = "black"), 
        axis.title.y = element_text(size = 40, face = "bold",color = "black"), 
        plot.margin = unit(c(2, 2, 2, 2), "cm")
    )+
    ylab("% of TD with MHA ≥1 bp")+
    scale_x_discrete(labels=x_axis)+

    annotation_custom(grob = grid::textGrob(pvalue[1], x = 0.22, y = 1.05, 
         gp = grid::gpar(fontsize = 23, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(pvalue[2], x = 0.51, y = 1.05, 
         gp = grid::gpar(fontsize = 23, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(sprintf("%.4f" ,pvalue[3]), x = 0.8, y = 1.05, 
         gp = grid::gpar(fontsize = 23, fontface = "bold",color = "black")))
p
p <- p+ggsci::scale_color_npg()+ggsci::scale_fill_npg()
p
ggsave(p,filename = 'Figure2/3level-M1_plot.png',width = 11,height = 13)
