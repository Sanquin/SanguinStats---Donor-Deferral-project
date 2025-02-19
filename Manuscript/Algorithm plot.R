library(ggplot2)
library(dplyr)
library(ggnewscale)
library(greekLetters)

#everything in mmol/L
mu <- 8.4 
sig <- 0.45
cut <- 8.4 

# make_example <- function(mu, siglevel = 0.999, sexM=T, n=20, seed=11) {
#   
#   set.seed(seed)
#   
#   if (sexM){
#     #mu <- 8.4 / 0.0626
#     sig <- 10.03  #based on finnish
#     cut <- 8.4 / 0.0626
#   } else {
#     #mu <- 7.8 / 0.0626
#     sig <- 9.7 #based on finnish
#     cut <- 7.8 / 0.0626
#   }
#   
#   if (sexM) {
#     dt <- round(runif((n-1), min=56, max = 85))
#   } else {
#     dt <- round(runif((n-1), min=122, max = 150))
#   }
#   
#   t <- c(0,cumsum(dt))
#   print(t)
#   
#   nsig <- qnorm(siglevel)
#   x <- rnorm(n, mean=mu, sd=sig)
#   x[sample(1:n, 1)] <- x[sample(1:n, 1)] - sig*(nsig+1)
#   
#   mean <- cumsum(x) / seq(1, n)
#   historicmean <- c(NA,mean[1:(n-1)])
#   std <- sig / sqrt(seq(1, n))
#   low <- mean - nsig * std
#   high <- mean + nsig * std
#   
#   low2 <- historicmean - nsig * sig
#   high2 <- historicmean + nsig * sig
#   
#   #low3 <- low - nsig * sig
#   #high3 <- high + nsig * sig
#   
#   deferred <- x < cut
#   deferred_outlier <- x + nsig * sig < historicmean
#   deferred_new <- ifelse(deferred_outlier==TRUE & ! is.na(deferred_outlier), deferred_outlier, high < cut)
#   
#   df <<- data.frame(t=t, x=x, mean=mean, low=low, high=high, low2=low2, high2=high2,
#                    deferred=deferred, deferred_new=deferred_new, deferred_outlier=deferred_outlier)
#   
#   ggplot(df, aes(x=t, y=x)) +  
#     
#     geom_ribbon(aes(ymin=low2, ymax=high2), fill="grey70", alpha=0.3) +
#     geom_ribbon(aes(ymin=low, ymax=high), fill="#56b4e9", alpha=0.3) +
#     
#     geom_hline(yintercept=cut, linetype="dashed", color="#e69f00") +
#     geom_line(aes(y=mean), color="#0072b2") +
#     geom_line(color="black") +
#     
#     geom_point(data=df, aes(shape = deferred, col = deferred_new), size=2, stroke=1) +
#     scale_shape_manual(deferred_new, values = c(16, 1), labels = c("Donation", "Deferral"), name = "Current algorithm") +
#     scale_color_manual(values = c("#009e73", "#d55e00"), labels = c("Donation", "Deferral"), name = "New algorithm") +
#     
#     geom_text(data=subset(df, deferred_outlier), aes(x=t+15, y=x, label="!"), color="#d55e00", size=6, fontface="bold")+
#     
#     theme_minimal() +
#     labs(title="Deferral algorithm on simulated donor career", x="Time (days)", y="Hemoglobin (g/L)") +
#     theme(legend.position="right") +
#     guides(color=guide_legend(nrow=2))
# }

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
  
  p <- ggplot(df, aes(x=t, y=x)) +  
    
    # Reverse the order of ribbons for the legend
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.3, fill="#56b4e9") +
    geom_ribbon(aes(ymin=low2, ymax=high2), alpha=0.3, fill="grey70") +
    geom_text(data = df[2, ], aes(x = t, y = high2+0.1, label = "Expected variability of the measurement relative to the mean"), vjust = 0, colour = "grey50", hjust = "inward")+
    # add the cut-off line
    geom_line(aes(y=cut, x=t), linewidth=0.70, linetype="dashed", color = "black") +
    geom_text(data = df[1, ], aes(x = t, y = cut-0.2, label = "Legal deferral cut off"), vjust = 0, colour = "black", hjust = "inward")+
    
    geom_line(aes(y=mean), linewidth=1, linetype="solid", color = "#0072b2") + 
    geom_text(data = df[2, ], aes(x = t, y = mean+0.1, label = "  Historical mean Hb"), vjust = 0, colour = "#0072b2", hjust = "inward")+
    geom_line(aes(y=x, x=t), color="black", linewidth=0.5, linetype = 'solid') +
    
    # Points for deferred and deferred_new with legends
    geom_point(aes(shape = deferred, color = deferred_new), size=4, stroke=1) +
    scale_color_manual(name = "New algorithm", values = c("#009e73", "#d55e00"), 
                       labels = c("Donation", "Deferral"), guide = guide_legend(order=4)) +
    scale_shape_manual(name = "Current algorithm", values = c(16, 1), labels = c("Donation", "Deferral"), guide=guide_legend(order=3)) +
    new_scale("shape")+
    theme_minimal() +
    labs(title=paste("Simulated donor career, ", greeks("alpha"), " outlier =", siglevel_outl, ", ", greeks("alpha"), " mean = ", siglevel_mean), x="Time (days)", y="Hemoglobin (g/dL)") +
    theme(legend.position="right", legend.spacing.y = unit(5, "pt")) #+
    # guides(
    #   linetype=guide_legend(nrow=2, title=NULL),
    #   fill=guide_legend(nrow=2, title=NULL),
    #   color=guide_legend(nrow=2)
    # )
  
  # Adding exclamation mark for outliers in the legend
  if (outlier_lab){
    p <-  p + geom_text(data = df[df$deferred_outlier==T,], aes(x=t+17, y=x, label = "Outlier"), vjust = 0.5, colour = "#d55e00", hjust = "outward")
  }
  
  if (siglevel_mean != 0){
    p <- p + geom_text(data = df[2, ], aes(x = t, y = mean+0.1, label = "  Historical mean Hb + CI"), vjust = 0, colour = "#0072b2", hjust = "inward")
  }
  
  print(seed)
  if (saveplot) {
    # Check if directory exists, if not create it
    if (!dir.exists("~/Amber/DonorDeferral/results/")) {
      dir.create("~/Amber/DonorDeferral/results/", recursive = TRUE)
    }
    # Save the plot with ggsave
    ggsave(filename = paste0("~/Amber/DonorDeferral/results/", seed, "_", siglevel_outl, "_",siglevel_mean, ".jpg"), plot = p, width = 12, height = 8)
  }
}

set.seed(sample(1:1e5, 1))

seed1 <- sample(1:1e5, 1)
print(seed1)

data <- c(145,135,138,130,146,138,124,135,130,120,133)/10
make_example(data=data, seed=seed1, n=length(data), saveplot=T, siglevel_outl=0.999, siglevel_mean = 0)



set.seed(sample(1:1e5, 1))
seed3 <-sample(1:1e5, 1)
print(seed3)

data <- c(9.5, 9, 8.5, 8.6, 7.8, 8.1, 8.3, 7.6, 8, 6.1, 7.8, 7.85, 8)/0.6206

make_example(data = data, seed=seed3, n=length(data), saveplot=T, siglevel_outl=0.999, siglevel_mean = 0.5, sexM = F, outlier_lab = T)
