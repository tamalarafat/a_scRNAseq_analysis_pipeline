generate_data_summary_plot <- function(seuratObject,
                                  covariate,
                                  threshold = NULL,
                                  draw_threshold = FALSE,
                                  figure_name_suffix = NULL,
                                  binwidth = 1000){
  
  # Extract the metadata information from the seurat object
  CVR_mat <- as.data.frame(seuratObject@meta.data)
  
  # Calculate the summary statistics of the covariates and create annotation label for the plot
  summ <- data.frame(min = min(CVR_mat[[covariate]]), 
                     max = max(CVR_mat[[covariate]]),
                     mean = mean(CVR_mat[[covariate]]), 
                     q1= quantile(CVR_mat[[covariate]], probs = 0.25), 
                     median = median(CVR_mat[[covariate]]), 
                     q3= quantile(CVR_mat[[covariate]], probs = 0.75),
                     sd = sd(CVR_mat[[covariate]]),
                     lab = paste("Number of cells = ", ncol(seuratObject), 
                                 "\nmin = ", round(min(CVR_mat[[covariate]]), 2), 
                                 "\nmax = ", round(max(CVR_mat[[covariate]]), 2), 
                                 "\nmean = ", round(mean(CVR_mat[[covariate]]), 2), 
                                 "\nq1 = ",round(quantile(CVR_mat[[covariate]], probs = 0.25), 2), 
                                 "\nmedian = ", round(median(CVR_mat[[covariate]]), 2),
                                 "\nq3 = ",round(quantile(CVR_mat[[covariate]], probs = 0.75), 2), 
                                 "\nsd = ",round(sd(CVR_mat[[covariate]]), 2), 
                                 "\n",names(CVR_mat[covariate]), ">", threshold, "=" , 
                                 nrow(CVR_mat[CVR_mat[[covariate]] > threshold, ]), " cells"))
  
  p <- ggplot(CVR_mat, aes(x = .data[[covariate]])) + geom_histogram(binwidth = binwidth, alpha = 1, col = "white", fill = "#00A08A") + 
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) + 
    {if (draw_threshold) geom_vline(data = CVR_mat, aes(xintercept = threshold))} + 
    geom_text(data = summ, aes(label = lab), x = Inf, y = Inf, hjust = 1, vjust = 1.2, size = 8) +
    xlab(names(CVR_mat[covariate])) + 
    ylab("Number of cells") + 
    ggtitle(str_c("Distribution of ", names(CVR_mat[covariate]))) + 
    theme_classic() + 
    theme(axis.title.x = element_text(size = 18, face = "bold", colour = "black"), 
          axis.title.y = element_text(size = 18, face = "bold", colour = "black"), 
          axis.ticks.length = unit(.30, "cm"), 
          axis.text = element_text(size = 18, face = "bold", colour = "black"),
          title =  element_text(size = 18, face = "bold", colour = "black"))
  
  ggsave(filename = str_c(covariate, "_distribution", figure_name_suffix, ".png"), plot = p, width = 16, height = 12, dpi = 300)
}
