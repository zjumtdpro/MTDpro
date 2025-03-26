library(ggplot2)
library(tidyr)
setwd("~/Desktop/graph/MTD/")

raw_data <-  openxlsx::read.xlsx("mtd-figure3-human.xlsx", sheet = 4, colNames = FALSE)

data <- as.data.frame(t(raw_data))
x_axis <- data[1,2:5]
x_axis <- ifelse(grepl("^[0-9]", x_axis), paste0("=", x_axis), gsub("^>=", "≥", x_axis))
x_axis <- c("name" , x_axis)
value <- data[2:10,1:5]
colnames(value) <- x_axis
random <- c(75,25,6.25,1.5625)

#colors <- pal_npg("nrc")(9)
colors <- ggsci::pal_d3("category10")(9)
value <- value %>%
    pivot_longer(cols = !name,
    names_to = "length", values_to = "percentage")
value$length <- as.factor(value$length)
value$percentage <- as.numeric(value$percentage)
p <- ggplot(data = value, aes(x = length, y = percentage, color = name))+
    geom_point(size = 4)+
    scale_y_continuous(expand = c(0,0),limits = c(0, 100),breaks = seq(0, 100, 20))+
    scale_color_manual(values = colors)+
    theme_classic()+
    theme(axis.title.x = element_blank(),
    axis.text = element_text(size = 21), 
    axis.title = element_text(size = 21), 
    legend.text = element_text(size = 15),
    legend.key.height = unit(1, "cm"),
    axis.ticks = element_line(size = 1.5), 
    axis.ticks.length = unit(0.5, "cm"), 
    axis.line = element_line(size = 1.5), 
    axis.text.x = element_text(size = 20, face = "bold",color = "black", angle = 0),
    axis.text.y = element_text(size = 20, face = "bold",color = "black"),
    axis.title.y = element_text(size = 20, face = "bold",color = "black"), 
    plot.margin = unit(c(2, 2, 2, 2), "cm")
    )+
    coord_cartesian(clip = 'off')+
    guides(color = guide_legend(title = NULL)) + 
    ylab("percentage (%)")+
    annotation_custom(grob = grid::textGrob("(bp)", x = 0.95, y = -0.07, 
                just = "left", gp = grid::gpar(fontsize =21, fontface = "bold")))+
    annotation_custom(grob = grid::textGrob(pvalue[1], x = 0.15, y = 0.97, 
         gp = grid::gpar(fontsize = 15, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(pvalue[2], x = 0.385, y = 0.97, 
         gp = grid::gpar(fontsize = 15, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(pvalue[3], x = 0.625, y = 0.97, 
         gp = grid::gpar(fontsize = 15, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(pvalue[4], x = 0.865, y = 0.97, 
         gp = grid::gpar(fontsize = 15, fontface = "bold",color = "black")))+
    geom_segment(aes(x = as.numeric(length)-0.4, xend = as.numeric(length)+0.4, y = random[as.numeric(length)], yend = random[as.numeric(length)]), color = "black", 
                 linetype = "dashed", size = 0.8) 

p
ggsave(p,filename = 'Figure3/disease_classify.png',width = 13,height = 6)
