library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Define paths
file_dir <- "~/Amber/alt_Hb_deferral/results"
results_path <- "~/Amber/alt_Hb_deferral/results"

# List all country files
cf_list <- list.files(file_dir, full.names = TRUE, include.dirs=F)
cf_list <- cf_list[grepl("confusion.*\\.csv$", cf_list)]

hist_list <- list.files(file_dir, full.names = TRUE)
hist_list <- hist_list[grepl("\\.rds$", hist_list)]

char_list <- list.files(file_dir, full.names = TRUE, include.dirs=F)
char_list <- char_list[grepl("characteristics.*\\.csv$", char_list, ignore.case = TRUE)]

# Define functions outside the loop
get_hist <- function(hists, sex, label, alpha_mean, alpha_outlier) {
  h <- hists[[paste0(sex, "_", alpha_mean, "_", alpha_outlier)]][[paste0("hist_output_", label)]]
  return(h)
}

get_mean_Hb <- function(hists, sex, label, alpha_mean, alpha_outlier) {
  return(get_hist(hists, sex, label, alpha_mean, alpha_outlier)$mean_Hb)
}

get_mean_mean_Hb <- function(hists, sex, label, alpha_mean, alpha_outlier) {
  return(get_hist(hists, sex, label, alpha_mean, alpha_outlier)$mean_mean_Hb)
}

v_get_mean_Hb <- Vectorize(get_mean_Hb, c("sex", "alpha_mean", "alpha_outlier"))
v_get_mean_mean_Hb <- Vectorize(get_mean_mean_Hb, c("sex", "alpha_mean", "alpha_outlier"))

# Function to create plots
template_plot <- function(data, x_var, y_var, y_label, title, y_limits = NULL) {

  if (is.null(y_limits)) {
    y_min <- min(data[[y_var]], na.rm = TRUE)  # get min value of y_var
    y_max <- max(data[[y_var]], na.rm = TRUE)  # get max value of y_var
    y_limits <- c(y_min, y_max)
  }
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var, color = "sex", shape = "alpha_outlier", linetype = "alpha_outlier")) +
    geom_hline(yintercept = 0, col = "black", size = 0.25) +
    geom_vline(xintercept = 0, col = "black", size = 0.25) +
    geom_line(size = 1) +
    theme_bw() +
    scale_color_manual(values = c(F="#d22020b9", M="#153fe8c6"), labels = c("Females", "Males")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    labs(
      x = expression(alpha[" mean"]),
      y = y_label,
      linetype = expression(alpha[" outlier"]),
      shape = expression(alpha[" outlier"]),
      title = title
    ) +
    scale_y_continuous(limits = y_limits) 
  
  if(grepl('Deferral rate', title)){
    current_deferral_male <- (mean(cf_processed$n_true_def[cf_processed$sex=="M"])/mean(cf_processed$total[cf_processed$sex=="M"]))*100
    current_deferral_female <- (mean(cf_processed$n_true_def[cf_processed$sex=="F"])/mean(cf_processed$total[cf_processed$sex=="F"]))*100
    p<- p +
      geom_hline(yintercept = current_deferral_male, color = "darkblue", linetype = "solid") +
      geom_hline(yintercept = current_deferral_female, color = "darkred", linetype = "solid") +
      annotate("text", x = max(cf_processed$alpha_mean), y = current_deferral_male+0.3, 
               label = "Current deferral rate males", color = "darkblue",hjust = 1) +
      annotate("text", x = max(cf_processed$alpha_mean), y = current_deferral_female+0.3, 
               label = "Current deferral rate females", color = "darkred", hjust = 1) 
  } else if (grepl('Increase in number of donations', title)){
    max_increase_males <- 100*(((mean(cf_processed$n_true_def[cf_processed$sex=="M"])+mean(cf_processed$total_donations_old[cf_processed$sex=="M"]))-cf_processed$total_donations_old[cf_processed$sex=="M"])/cf_processed$total_donations_old[cf_processed$sex=="M"])
    max_increase_females <- 100*(((mean(cf_processed$n_true_def[cf_processed$sex=="F"])+mean(cf_processed$total_donations_old[cf_processed$sex=="F"]))-cf_processed$total_donations_old[cf_processed$sex=="F"])/cf_processed$total_donations_old[cf_processed$sex=="F"])
    
    if(max_increase_females > max(y_limits)){
      y_limits <- c(min(y_limits), max_increase_females+5)
    } else if(max_increase_males > max(y_limits)& max_increase_males > max_increase_females){
      y_limits <- c(min(y_limits), max_increase_males+5)
    }
    
    p<- p +
      geom_hline(yintercept = max_increase_males, color = "darkblue", linetype = "solid") +
      geom_hline(yintercept = max_increase_females, color = "darkred", linetype = "solid") +
      annotate("text", x = max(cf_processed$alpha_mean), y = max_increase_males+0.7, 
               label = "Max. donation increase males", color = "darkblue",hjust = 1) +
      annotate("text", x = max(cf_processed$alpha_mean), y = max_increase_females+0.7, 
               label = "Max. donation increase females", color = "darkred", hjust = 1)+
      scale_y_continuous(limits = y_limits)
  }
  
  return(p)
}

# Loop through country files
for (file in cf_list) {
  cf <- read.csv(file, row.names = 1)
  country <- tail(strsplit(gsub(".csv", "", file), "_")[[1]], 1)
  country_path <- file.path(results_path, country)
  dir.create(country_path, showWarnings = FALSE)
  # Identify corresponding histogram file
  hist_file <- hist_list[grepl(paste0(country, "\\.rds$"), hist_list)]
  characteristics_file <- char_list[grepl(paste0(country), char_list)]
  characteristics <- read.csv(characteristics_file, row.names = 1)

  # Process data
  cf_processed <- cf %>% mutate(
    country=country,
    total = TN + FP + FN + TP,
    n_def_alt = n_mean_def + n_outlier_def - n_both_def,
    total_donations_old = total - n_true_def,
    total_donations_new = total - n_def_alt,
    deferral_reduction = n_true_def - n_def_alt,
    deferral_reduction_perc = deferral_reduction / n_true_def * 100,
    alpha_outlier = as.factor(alpha_outlier),
    wrong_donations_prevented_perc = 100 * FP / total_donations_old,
    extra_donations_perc = 100 * (total_donations_new - total_donations_old) / total_donations_old,
    perc_mean_def = (n_mean_def / n_def_alt) * 100,
    perc_outlier_def = (n_outlier_def / n_def_alt) * 100
  )
    if (length(hist_file) == 1) {
    hists <- readRDS(hist_file)
  for (label in c("FN", "FP", "TP", "TN", "old", "new")) {
    cf_processed[paste0("mean_Hb_", label)] <- v_get_mean_Hb(hists, cf_processed$sex, label, cf_processed$alpha_mean, cf_processed$alpha_outlier)
    cf_processed[paste0("mean_mean_Hb_", label)] <- v_get_mean_mean_Hb(hists, cf_processed$sex, label, cf_processed$alpha_mean, cf_processed$alpha_outlier)
  }
    cf_processed$mean_mean_Hb_donations <- ((cf_processed$TN*cf_processed$mean_mean_Hb_TN)+(cf_processed$FN*cf_processed$mean_mean_Hb_FN))/(cf_processed$TN+cf_processed$FN)
    
    cols_meanHb <- grepl("mean_Hb", colnames(cf_processed))
    if (characteristics[12,2]=="g/L"){
      cf_processed[cols_meanHb] <- cf_processed[cols_meanHb]/10
    } else if (characteristics[12,2]=='mmol/L'){
      cf_processed[cols_meanHb] <- cf_processed[cols_meanHb]/0.6206
    }
    }
  
  # Create and save plots
  plots <- list(
    "increase_donations" = template_plot(cf_processed, "alpha_mean", "extra_donations_perc","Percentage", paste0("Increase in number of donations - ", country)), #lijn voor maximale toename in donations
    "deferral_rate" = template_plot(cf_processed, "alpha_mean", "(n_def_alt / total) * 100", "Percentage", paste0("Deferral rate - ", country), c(0,10)), 
    "perc_mean_deferrals" = template_plot(cf_processed, "alpha_mean", "perc_mean_def", 
"Percentage", paste0("Percentage mean deferrals - ", country), c(0, 100)),
    "perc_outlier_deferrals" = template_plot(cf_processed, "alpha_mean", "perc_outlier_def","Percentage", paste0("Percentage outlier deferrals - ", country), c(0, 100)),
    "deferral_reduction" = template_plot(cf_processed, "alpha_mean", "deferral_reduction_perc", "Percentage", paste0("Deferral reduction - ", country), c(-50, 100)),
    "ineligible_donations" = template_plot(cf_processed, "alpha_mean", "wrong_donations_prevented_perc", "Percentage", paste0("Ineligible donations - ", country))
  )
  
  if (length(hist_file) == 1) {
    plots2 <- list(
      "mean mean hb all donors alt strat" = template_plot(cf_processed, "alpha_mean", "mean_mean_Hb_donations", "Hb (g/dL)" , paste0("Mean of mean Hb in all donors donating with alternative strategy - ", country)),
    "mean mean hb new donors alt strat" = template_plot(cf_processed, "alpha_mean", "mean_mean_Hb_FN", "Hb (g/dL)", paste0("Mean of mean Hb in new donors donating with alternative strategy - ", country))
    )
    plots <- append(plots, plots2)
  }
  
  
  for (plot_name in names(plots)) {
    ggsave(filename = file.path(country_path, paste0(plot_name, ".png")), plot = plots[[plot_name]], width = 7, height = 7)
  }
}

