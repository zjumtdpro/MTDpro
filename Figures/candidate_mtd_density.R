library(ggplot2)
# library(ggpmisc) # Can be used to display fitting formulas and R^2 values
library(tidyr)
setwd("~/Desktop/graph/MTD/")
eq_label <- "y = 1.64x + 2.2583"
r2_label <- "R² = 0.9718"
# load data
raw_data <- openxlsx::read.xlsx("mtd-figure4.xlsx", sheet = 1, colNames = FALSE)
df <- raw_data[3:4,3:ncol(raw_data) ]
df <- as.data.frame(t(df))
colnames(df) <- c("x", "y")
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)

# fungi (red), bacteria (blue), archaea (green) and virus (grey)
colors <- c(rep(rgb(42,54,130, maxColorValue = 255),8),rgb(16,124,54, maxColorValue = 255),rep(rgb(230,0,29, maxColorValue = 255),8),'grey70')

p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = colors, size = 4.5,alpha = 0.9) + 
    geom_smooth(method = "lm", se = FALSE, color = "black") + 
    #   stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    #                formula = y ~ x,
    #                parse = TRUE,
    #                size = 5, fontface = "bold")+ 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0, 100, 10)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 50), breaks = seq(0, 50, 10)) +
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
        axis.text.x = element_text(size = 24, face = "bold", color = "black", angle = 0),
        axis.text.y = element_text(size = 24, face = "bold", color = "black"),
        axis.title.y = element_text(size = 30, face = "bold", color = "black"),
        axis.title.x = element_text(size = 30, face = "bold", color = "black"),
        strip.text.x = element_text(size = 17),
        plot.margin = unit(c(2, 2, 2, 2), "cm")
    ) +
    xlab("Genome size (Mb)") +
    ylab("No. of candidate MTD (M)") +
    annotate("text", x = 15, y = 100, label = eq_label, parse = FALSE, size = 8, fontface = "bold", color = "black") + 
    annotate("text", x = 14, y = 93, label = r2_label, parse = FALSE, size = 8, fontface = "bold", color = "black")

    p

ggsave(p, filename = "Figure4/candidate_mtd_density.png", width = 10, height = 9)
