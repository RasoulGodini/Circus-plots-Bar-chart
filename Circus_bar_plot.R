library(truncnorm)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

#-------- Generate random data ---------
set.seed(123)  # For reproducibility

n_participants <- 2
n_genes <- 30

sex <- sample(c("male", "female"), n_participants, replace = FALSE)

df <- data.frame(Gene_ID = paste0("Gene_", 1:n_genes),
                 male = rtruncnorm(n_genes, a = 0, b = 20, mean = 5, sd = 3),
                 male_pvalue = rtruncnorm(n_genes, a = 0, b = 0.9, mean = 0.08, sd = 0.07),
                 female = rtruncnorm(n_genes, a = 0, b = 20, mean = 12, sd = 8),
                 female_pvalue = rtruncnorm(n_genes, a = 0, b = 0.9, mean = 0.07, sd = 0.06),
                 Gene_family = sample(c("Family_A", "Family_B", "Family_C"), n_genes, replace = TRUE)
)


#Add coordinate of x_location on the plot for each sex
df_loc <- df %>% group_by(Gene_family) %>%
  mutate(x_location = row_number())


# Sort the genes based on the family
df_loc_sort <- df_loc %>% arrange(Gene_family)
df_loc_sort$Gene_ID <- factor(df_loc_sort$Gene_ID, 
                              levels = unique(df_loc_sort$Gene_ID))


# ----------------- Make the circus plot -----------------
circos.clear()

#To make the circus plot later
pdf("circus_bar_plot_example.pdf", width = 3.7, height = 3.7, bg = "transparent")

# Define spaces between each two gene family
costum_degree <- c("Family_A" = 5, 
                   "Family_B" = 5, 
                   "Family_C" = 20)

#Control the space around the objects of the plot
circos.par(gap.after = costum_degree, 
           cell.padding = c(0, 3, 0, 3), #Space around the bars 
           start.degree = 90) #Where the axis will be

circos.initialize(sectors = factor(df_loc_sort$Gene_family),
                  x = as.numeric(df_loc_sort$Gene_ID))

# Add the annotation layer
circos.trackPlotRegion(
  sectors = factor(df_loc_sort$Gene_family),
  x = as.numeric(df_loc_sort$Gene_ID), 
  ylim = c(21, 5), #Change based on your Y axis (here expression) values
  
  panel.fun = function(x, y) {
    # Get the current sector index (class)
    sector_index <- get.cell.meta.data("sector.index")
    
    # Define a color palette for the classes
    class_colors <- c("Family_A" = "purple", 
                      "Family_B" = "springgreen", 
                      "Family_C" = "orange")
    
    #The class annotation
    circos.rect(
      xleft = get.cell.meta.data("xlim")[1], 
      xright = get.cell.meta.data("xlim")[2], 
      ytop = 21, #Change based on your Y axis (here expression) values
      ybottom = 19, #Change based on your Y axis (here expression) values
      lwd = 0.4,
      col = class_colors[sector_index],
    )
  },bg.border = NA 
)
set_track_gap(mm_h(0.001))


# Add the outer track
circos.trackPlotRegion(
  sectors = factor(df_loc_sort$Gene_family), bg.lwd = 0.8,
  x = as.numeric(df_loc_sort$Gene_ID), bg.col = "azure",
  ylim = c(0, 21),
  panel.fun = function(x, y) {
    sector_index <- get.cell.meta.data("sector.index")
    sector_data <- df_loc_sort[df_loc_sort$Gene_family == sector_index, ]
    #To add a line 
    circos.lines(x = as.numeric(sector_data$Gene_ID),
                 y = rep(10, length(sector_data$Gene_ID)),
                 col = "black", lty = 2, lwd = 0.3, straight = TRUE)
    
    #Add expression values
    circos.barplot(value = sector_data$female, bar_width = 0.3, lwd = 0.2, 
                   pos = as.numeric(sector_data$Gene_ID),
                   col = case_when((sector_data$female >= 10 & 
                                      sector_data$female_pvalue <= 0.05) ~ "orangered",
                                   (sector_data$female < 10 & 
                                      sector_data$female_pvalue <= 0.05) ~ "dodgerblue",
                                   (sector_data$female_pvalue >= 0.05) ~ "grey")
    )
    
    # Add the labels (here gene names)
    circos.text(
      x = as.numeric(sector_data$Gene_ID), 
      y = sector_data$female - sector_data$female + 26, #To set the position of the labels
      labels = sector_data$Gene_ID, cex = 0.4, col = "black", 
      facing = "reverse.clockwise", adj = c(1, 0.5)
    )
  }
)
circos.yaxis(labels.cex = 0.37,side = "left", sector.index = "Family_A")
set_track_gap(mm_h(0.5)) 



# Add the inner track
circos.trackPlotRegion(
  sectors = factor(df_loc_sort$Gene_family), bg.lwd = 0.8,
  x = as.numeric(df_loc_sort$Gene_ID), bg.col = "ivory",
  ylim = c(0,21),
  panel.fun = function(x, y) {
    sector_index <- get.cell.meta.data("sector.index")
    sector_data <- df_loc_sort[df_loc_sort$Gene_family == sector_index, ]
    #To add a line 
    circos.lines(x = as.numeric(sector_data$Gene_ID),
                 y = rep(10, length(sector_data$Gene_ID)),
                 col = "black", lty = 2, lwd = 0.3, straight = TRUE)
    
    #Add expression values
    circos.barplot(value = sector_data$male, bar_width = 0.3, lwd = 0.2, 
                   pos = as.numeric(sector_data$Gene_ID),
                   col = case_when((sector_data$male >= 10 & 
                                      sector_data$male_pvalue <= 0.05) ~ "orangered",
                                   (sector_data$male < 10 & 
                                      sector_data$male_pvalue <= 0.05) ~ "dodgerblue",
                                   (sector_data$male_pvalue >= 0.05) ~ "grey")
    )
  }
)
circos.yaxis(labels.cex = 0.37,side = "left", sector.index = "Family_A")

# Prepare the legends
lgd_points = Legend(at = c("High expression \nPvalue <= 0.05)",  
                           "Low expression \nPvalue <= 0.05",
                           "Pvalue > 0.05"), 
                    type = "points", background = NA,
                    legend_gp = gpar(col = c("orangered", "dodgerblue", "grey")),
                    grid_height = unit(1, "mm"), gap = unit(0.2, "mm"),
                    grid_width = unit(2, "mm"),labels_gp = gpar(cex = 0.4),
                    title = expression("Gene expression"),
                    title_gp = gpar(cex = 0.4, fontface = "bold"),
                    title_gap = unit(0.2, "mm"),
                    title_position = "topcenter")


lgd_gene_family = Legend(at = c("Family_A", "Family_B", "Family_C"), 
                        type = "points", background = NA,
                        legend_gp = gpar(col = c("purple","springgreen","orange")),
                         grid_height = unit(1, "mm"), , gap = unit(0.2, "mm"),
                         grid_width = unit(2, "mm"),labels_gp = gpar(cex = 0.4),
                         title = expression("Gene family"), nrow = 2,
                         title_gp = gpar(cex = 0.4, fontface = "bold"),
                         title_gap = unit(0.2, "mm"),
                         title_position = "topcenter")

#Pack both legends
lgd_combined <- packLegend(lgd_points, lgd_gene_family, gap = unit(1, "mm"))

# Draw the legends
draw(lgd_combined, just = "center")

dev.off()

