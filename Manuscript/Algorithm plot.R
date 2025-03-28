library(ggplot2)
library(dplyr)
library(ggnewscale)
library(greekLetters)

#everything in mmol/L
mu <- 8.4 
sig <- 0.45
cut <- 8.4 


make_example <- function(data=F, siglevel_outl = 0.999, siglevel_mean = 0, sexM=T, n=20, seed=11, saveplot=F, outlier_lab =F) {
  
  set.seed(seed)
  
  x = data
  
  if (sexM){
    sig <- 10.03/10  #based on finnish
    cut <- 135/10
  } else {
    sig <- 9.7/10 #based on finnish
    cut <- 125/10
  }
  
  if (sexM) {
    dt <- round(runif((n-1), min=56, max = 85))
  } else {
    dt <- round(runif((n-1), min=122, max = 150))
  }
  
  t <- c(0,cumsum(dt))
  print(t)
  
  nsig_outl <- qnorm((1+siglevel_outl)/2)
  nsig_mean <- qnorm((1+siglevel_mean)/2)
  
  mean <- cumsum(x) / seq(1, n)
  mean[1]<- NA
  historicmean <- mean[1:n]
  std <- sig / sqrt(seq(1, n))
  low <- mean - nsig_mean * std
  high <- mean + nsig_mean * std
  high[2]<- NA
  low[2]<- NA
  
  low2 <- historicmean - nsig_outl * sig
  high2 <- historicmean + nsig_outl * sig
  
  deferred <- x < cut
  deferred_outlier <- x + nsig_outl * sig < historicmean
  deferred_new <- ifelse(deferred_outlier == TRUE & !is.na(deferred_outlier), deferred_outlier, high < cut)
  deferred_new[2]<- NA
  
  print(x)
  
  df <- data.frame(t=t, x=x, mean=mean, low=low, high=high, low2=low2, high2=high2, cut=cut,
                   deferred=deferred, deferred_new=deferred_new, deferred_outlier=deferred_outlier)
  
  x_end <- max(df$t, na.rm = TRUE)
  
  # Get corresponding high2 and low2 values at max(df$t)
  y_high_end <- df$high2[df$t == x_end]
  y_low_end <- df$low2[df$t == x_end]
  
  # Ensure there are no NA values
  if (is.na(y_high_end) | is.na(y_low_end)) {
    y_high_end <- max(df$high2, na.rm = TRUE)
    y_low_end <- min(df$low2, na.rm = TRUE)
  }
  
  #x_blue <- df[3, ]$t
  
  # Get corresponding high and low values for the blue ribbon at max(df$t)
  #y_high_end_blue <- df$high[df$t == x_blue]
  #y_low_end_blue <- df$low[df$t == x_blue]
  
  # Ensure no NA values for blue ribbon
  # if (is.na(y_high_end_blue) | is.na(y_low_end_blue)) {
  #   y_high_end_blue <- max(df$high, na.rm = TRUE)
  #   y_low_end_blue <- min(df$low, na.rm = TRUE)
  # }
  
  p <- ggplot(df, aes(x=t, y=x)) +  
    
    # Reverse the order of ribbons for the legend
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.3, fill="#56b4e9") +
    geom_ribbon(aes(ymin=low2, ymax=high2), alpha=0.3, fill="grey70") +
    geom_segment(aes(x = x_end, xend = x_end, 
                     y = y_low_end+0.01, yend = y_high_end), 
                 arrow = arrow(length = unit(0.2, "cm"), ends="both"), 
                 colour = "grey50", linewidth = 0.7)+
    # geom_segment(aes(x = x_blue-5, xend = x_blue-5,  # Slight offset to avoid overlap
    #                  y = y_low_end_blue, yend = y_high_end_blue), 
    #              arrow = arrow(length = unit(0.2, "cm"), ends = "both"), 
    #              colour = "#5E96B2", linewidth = 0.7) +
    geom_text(data = df[2, ], aes(x = x_end-10, y = y_high_end-0.5, label = paste0("Accepted variability of the measurement relative to the mean\nWidth determined by ", greeks("alpha"), "-outlier")), vjust = 0, colour = "grey50", hjust = "inward", size=5)+
    # add the cut-off line
    geom_line(aes(y=cut, x=t), linewidth=0.70, linetype="dashed", color = "black") +
    geom_text(data = df[1, ], aes(x = t, y = cut-0.2, label = "Legal deferral threshold"), vjust = 0, colour = "black", hjust = "inward")+
    
    geom_line(aes(y=mean), linewidth=1, linetype="solid", color = "#0072b2") + 
    
    geom_line(aes(y=x, x=t), color="black", linewidth=0.5, linetype = 'solid') +
    
    # Points for deferred and deferred_new with legends
    geom_point(aes(shape = deferred, color = deferred_new), size=4, stroke=1) +
    scale_color_manual(name = "Alternative algorithm", values = c("#009e73", "#d55e00"), 
                       labels = c("Donation", "Deferral"), guide = guide_legend(order=4)) +
    scale_shape_manual(name = "Current policy", values = c(16, 1), labels = c("Donation", "Deferral"), guide=guide_legend(order=3)) +
    new_scale("shape")+
    theme_minimal() +
    labs(title=paste0("Simulated donor career, ", greeks("alpha"), "-outlier =", siglevel_outl, ", ", greeks("alpha"), "-mean = ", siglevel_mean), x="Time (days)", y="Hemoglobin (g/dL)") +
    theme(legend.position="right", legend.spacing.y = unit(5, "pt"),text = element_text(size = 14)) #+
  # guides(
  #   linetype=guide_legend(nrow=2, title=NULL),
  #   fill=guide_legend(nrow=2, title=NULL),
  #   color=guide_legend(nrow=2)
  # )
  
  # Adding exclamation mark for outliers in the legend
  if (outlier_lab){
    p <-  p + geom_text(data = df[df$deferred_outlier==T,], aes(x=t+60, y=x, label = "Outlier"), vjust = 0.5, colour = "#d55e00", hjust = 1)
  }
  
  if (siglevel_mean != 0){
    p <- p + geom_text(data = df[2, ], aes(x = x_blue-5, y = y_high_end_blue+0.1, label = paste0("Accepted variability of historical mean Hb\nWidth determined by ", greeks("alpha"), "-mean")), vjust = 0, colour = "#5E96B2", hjust = "inward", size = 5)
  } else {
    p<- p+geom_text(data = df[2, ], aes(x = t, y = mean+0.1, label = "  Historical mean Hb"), vjust = 0, colour = "#0072b2", hjust = "inward", size = 5)
  }
  
  
  
  print(seed)
  if (saveplot) {
    # Check if directory exists, if not create it
    if (!dir.exists("~/Amber/DonorDeferral/results/")) {
      dir.create("~/Amber/DonorDeferral/results/", recursive = TRUE)
    }
    # Save the plot with ggsave
    ggsave(filename = paste0("~/Amber/DonorDeferral/results/", seed, "_", siglevel_outl, "_",siglevel_mean, ".pdf"), plot = p, width = 12, height = 8, device = cairo_pdf)
  }
}


set.seed(sample(1:1e5, 1))

seed1 <- sample(1:1e5, 1)
print(seed1)

data <- c(145,135,138,130,146,138,124,135,130,120,133)/10
make_example(data=data, seed=seed1, n=length(data), saveplot=T, siglevel_outl=0.999, siglevel_mean = 0)

make_example <- function(data=F, siglevel_outl = 0.999, siglevel_mean = 0, sexM=T, n=20, seed=11, saveplot=F, outlier_lab =F) {
  
  set.seed(seed)
  
  x = data
  
  if (sexM){
    sig <- 10.03/10  #based on finnish
    cut <- 135/10
  } else {
    sig <- 9.7/10 #based on finnish
    cut <- 125/10
  }
  
  if (sexM) {
    dt <- round(runif((n-1), min=56, max = 85))
  } else {
    dt <- round(runif((n-1), min=122, max = 150))
  }
  
  t <- c(0,cumsum(dt))
  print(t)
  
  nsig_outl <- qnorm((1+siglevel_outl)/2)
  nsig_mean <- qnorm((1+siglevel_mean)/2)
  
  mean <- cumsum(x) / seq(1, n)
  mean[1]<- NA
  historicmean <- mean[1:n]
  std <- sig / sqrt(seq(1, n))
  low <- mean - nsig_mean * std
  high <- mean + nsig_mean * std
  high[2]<- NA
  low[2]<- NA
  
  low2 <- historicmean - nsig_outl * sig
  high2 <- historicmean + nsig_outl * sig
  
  deferred <- x < cut
  deferred_outlier <- x + nsig_outl * sig < historicmean
  deferred_new <- ifelse(deferred_outlier == TRUE & !is.na(deferred_outlier), deferred_outlier, high < cut)
  deferred_new[2]<- NA
  
  print(x)
  
  df <- data.frame(t=t, x=x, mean=mean, low=low, high=high, low2=low2, high2=high2, cut=cut,
                   deferred=deferred, deferred_new=deferred_new, deferred_outlier=deferred_outlier)
  
  x_end <- max(df$t, na.rm = TRUE)
  
  # Get corresponding high2 and low2 values at max(df$t)
  y_high_end <- df$high2[df$t == x_end]
  y_low_end <- df$low2[df$t == x_end]
  
  # Ensure there are no NA values
  if (is.na(y_high_end) | is.na(y_low_end)) {
    y_high_end <- max(df$high2, na.rm = TRUE)
    y_low_end <- min(df$low2, na.rm = TRUE)
  }
  
  x_blue <- df[3, ]$t
  
  # Get corresponding high and low values for the blue ribbon at max(df$t)
  y_high_end_blue <- df$high[df$t == x_blue]
  y_low_end_blue <- df$low[df$t == x_blue]
  
  # Ensure no NA values for blue ribbon
  if (is.na(y_high_end_blue) | is.na(y_low_end_blue)) {
    y_high_end_blue <- max(df$high, na.rm = TRUE)
    y_low_end_blue <- min(df$low, na.rm = TRUE)
  }
  
  p <- ggplot(df, aes(x=t, y=x)) +  
    
    # Reverse the order of ribbons for the legend
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.3, fill="#56b4e9") +
    geom_ribbon(aes(ymin=low2, ymax=high2), alpha=0.3, fill="grey70") +
    geom_segment(aes(x = x_end, xend = x_end, 
                     y = y_low_end+0.01, yend = y_high_end), 
                 arrow = arrow(length = unit(0.2, "cm"), ends="both"), 
                 colour = "grey50", linewidth = 0.7)+
    geom_segment(aes(x = x_blue-5, xend = x_blue-5,  # Slight offset to avoid overlap
                     y = y_low_end_blue, yend = y_high_end_blue), 
                 arrow = arrow(length = unit(0.2, "cm"), ends = "both"), 
                 colour = "#5E96B2", linewidth = 0.7) +
    geom_text(data = df[2, ], aes(x = x_end-10, y = y_high_end-0.5, label = paste0("Accepted variability of the measurement relative to the mean\nWidth determined by ", greeks("alpha"), "-outlier")), vjust = 0, colour = "grey50", hjust = "inward", size=5)+
    # add the cut-off line
    geom_line(aes(y=cut, x=t), linewidth=0.70, linetype="dashed", color = "black") +
    geom_text(data = df[1, ], aes(x = t, y = cut-0.2, label = "Legal deferral threshold"), vjust = 0, colour = "black", hjust = "inward")+
    
    geom_line(aes(y=mean), linewidth=1, linetype="solid", color = "#0072b2") + 
    
    geom_line(aes(y=x, x=t), color="black", linewidth=0.5, linetype = 'solid') +
    
    # Points for deferred and deferred_new with legends
    geom_point(aes(shape = deferred, color = deferred_new), size=4, stroke=1) +
    scale_color_manual(name = "Alternative algorithm", values = c("#009e73", "#d55e00"), 
                       labels = c("Donation", "Deferral"), guide = guide_legend(order=4)) +
    scale_shape_manual(name = "Current policy", values = c(16, 1), labels = c("Donation", "Deferral"), guide=guide_legend(order=3)) +
    new_scale("shape")+
    theme_minimal() +
    labs(title=paste0("Simulated donor career, ", greeks("alpha"), "-outlier =", siglevel_outl, ", ", greeks("alpha"), "-mean = ", siglevel_mean), x="Time (days)", y="Hemoglobin (g/dL)") +
    theme(legend.position="right", legend.spacing.y = unit(5, "pt"),text = element_text(size = 14)) #+
  # guides(
  #   linetype=guide_legend(nrow=2, title=NULL),
  #   fill=guide_legend(nrow=2, title=NULL),
  #   color=guide_legend(nrow=2)
  # )
  
  # Adding exclamation mark for outliers in the legend
  if (outlier_lab){
    p <-  p + geom_text(data = df[df$deferred_outlier==T,], aes(x=t+60, y=x, label = "Outlier"), vjust = 0.5, colour = "#d55e00", hjust = 1)
  }
  
  if (siglevel_mean != 0){
    p <- p + geom_text(data = df[2, ], aes(x = x_blue-5, y = y_high_end_blue+0.1, label = paste0("Accepted variability of historical mean Hb\nWidth determined by ", greeks("alpha"), "-mean")), vjust = 0, colour = "#5E96B2", hjust = "inward", size = 5)
  } else {
    p<- p+geom_text(data = df[2, ], aes(x = t, y = mean+0.1, label = "  Historical mean Hb"), vjust = 0, colour = "#0072b2", hjust = "inward", size = 5)
  }
  
  
  
  print(seed)
  if (saveplot) {
    # Check if directory exists, if not create it
    if (!dir.exists("~/Amber/DonorDeferral/results/")) {
      dir.create("~/Amber/DonorDeferral/results/", recursive = TRUE)
    }
    # Save the plot with ggsave
    ggsave(filename = paste0("~/Amber/DonorDeferral/results/", seed, "_", siglevel_outl, "_",siglevel_mean, ".pdf"), plot = p, width = 12, height = 8, device = cairo_pdf)
  }
}


set.seed(sample(1:1e5, 1))
seed3 <-sample(1:1e5, 1)
print(seed3)

data <- c(9.1, 8.9, 9.0, 8.7, 9.1, 6.4, 9.1, 9.3, 9.1, 9.0, 9.2, 9.1)/0.6206

make_example(data = data, seed=seed3, n=length(data), saveplot=T, siglevel_outl=0.999, siglevel_mean = 0.5, outlier_lab = T)
