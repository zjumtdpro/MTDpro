library(ggplot2)
library(ggnewscale)
library(tidyr)
library(cowplot)
setwd("~/Desktop/graph/MTD/")
# load data
raw_data <-  openxlsx::read.xlsx("mtd-figure3-2.xlsx", sheet = 4, colNames = FALSE)
x_axis <- raw_data[2,2:5]
df1 <- raw_data[3:12,2:5]
colnames(df1) <- x_axis
df1$group <- raw_data[1,2]
df1$group <- "COSMIC\r\n(n = 17,180)"
df2 <- raw_data[3:12,7:10]
colnames(df2) <- x_axis
df2$group <- raw_data[1,7]
df2$group <- "ClinVar\r\n(n = 1,223)"
df <- rbind(df1,df2)
pvalue <- c(rep("<0.0001",4))
random <- c(75,25,6.25,1.5625)

df_long <- df %>% 
    pivot_longer(cols = !group,
                names_to = "length",
                values_to = "percentage")
df_long$percentage <- as.numeric(df_long$percentage)*100
df_long$length <- factor(df_long$length, levels = x_axis)
colors <- c('blue', 'red', '#01c101', 'purple')
x_axis <- ifelse(grepl("^[0-9]", x_axis), paste0("=", x_axis), gsub("^>=", "≥", x_axis))

p <- 
ggplot(data = df_long, aes(x = length, y = percentage,color=length)) +
    geom_boxplot(width = 0.5, size = 0.2, outlier.shape = NA, alpha = 1)+
    geom_jitter(position=position_jitterdodge(jitter.width = 0.5, 
                                            jitter.height = 0, 
                                            dodge.width = 1),
                                            size = 3.5, alpha = 0.35)+

    scale_color_manual(values=colors)+
    scale_fill_manual(values = colors)+
    scale_y_continuous(expand = c(0,0),limits = c(0, 100),breaks = seq(0, 100, 10))+
    theme_classic()+
    coord_fixed(0.035)+
    geom_segment(aes(x = as.numeric(length)-0.4, xend = as.numeric(length)+0.45,
     y = random[as.numeric(length)], yend = random[as.numeric(length)], color = length), 
                 linetype = "dashed", size = 1) + 
     theme(axis.title.x = element_blank(),legend.position="none",
    axis.text = element_text(size = 21), 
    axis.title = element_text(size = 21), 
    legend.text = element_text(size = 18))+
    coord_cartesian(clip = 'off')+ 
    theme(
        axis.ticks = element_line(size = 1.5), 
        axis.ticks.length = unit(0.5, "cm"),  
        axis.line = element_line(size = 1.5), 
        axis.text.x = element_text(size = 20, face = "bold",color = "black", angle = 0), 
        axis.text.y = element_text(size = 20, face = "bold",color = "black"),
        axis.title.y = element_text(size = 30, face = "bold",color = "black"), 
        strip.text.x = element_text(size = 20),
        strip.background = element_blank(),
        plot.margin = unit(c(2, 2, 2, 2), "cm")
    )+
    ylab("% of TDs with MHA")+
    scale_x_discrete(labels=x_axis)+
    annotation_custom(grob = grid::textGrob(pvalue[1], x = 0.15, y = 0.97, 
         gp = grid::gpar(fontsize = 14, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(pvalue[1], x = 0.38, y = 0.97, 
         gp = grid::gpar(fontsize = 14, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(pvalue[1], x = 0.615, y = 0.97, 
         gp = grid::gpar(fontsize = 14, fontface = "bold",color = "black")))+
    annotation_custom(grob = grid::textGrob(pvalue[1], x = 0.845, y = 0.97, 
         gp = grid::gpar(fontsize = 14, fontface = "bold",color = "black")))
p <- p + facet_wrap(~group)+
    geom_segment(aes(x = 0, xend = 4.5, y = 100, yend = 100), color = "black", size = 1)

p

p_annotated <- ggdraw(p) +
  draw_label("(bp)", x = 0.93, y = 0.12, hjust = 0.5, fontface = "bold", size = 20)
p_annotated
ggsave(p_annotated,filename = 'Figure3-2/disease.png',width = 10,height = 8.5)







