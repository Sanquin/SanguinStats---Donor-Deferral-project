library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Define paths
file_dir <- "~/Amber/alt_Hb_deferral/results"
results_path <- "~/Amber/alt_Hb_deferral/results"

# List all country files
cf_list <- list.files(file_dir, full.names = TRUE)
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

# Initialize empty data frame for combined data
combined_data <- data.frame()

# Loop through country files and process data
for (file in cf_list) {
  cf <- read.csv(file, row.names = 1)
  country <- tail(strsplit(gsub(".csv", "", file), "_")[[1]], 1)
  cf$country <- country
  
  # Identify corresponding histogram file
  hist_file <- hist_list[grepl(country, hist_list)]
  characteristics_file <- char_list[grepl(paste0(country), char_list)]
  characteristics <- read.csv(characteristics_file, row.names = 1)
  
  
  if (length(hist_file) == 1) {
    hists <- readRDS(hist_file)
    
    cf_processed <- cf %>% mutate(
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
      perc_outlier_def = (n_outlier_def / n_def_alt) * 100,
      current_deferral_rate = (n_true_def/total)*100,
      max_increase_donations = 100*((((n_true_def)+(total_donations_old))-total_donations_old)/total_donations_old)
    )
    
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

    
    combined_data <- bind_rows(combined_data, cf_processed)
  } else {
    warning(paste("No unique histogram file found for country:", country))
  }
}

# Define color palettes
colors <- c("#6A3D9A", "#FFD700","#4A90E2", "#27408B", "#FF7F00",  "#33A02C","#E31A1C")

# Function to create facet plots with different color schemes
template_plot <- function(data, x_var, y_var, y_label, title, y_limits = NULL) {
  p <- ggplot(data, aes_string(x = x_var, y = y_var, color = "country", linetype = "'Alternative'")) +
    geom_hline(yintercept = 0, col = "black", size = 0.25) +
    geom_vline(xintercept = 0, col = "black", size = 0.25) +
    geom_line(size = 0.8) +
    theme_bw() +
    labs(
      x = expression(alpha["mean"]),
      y = y_label,
      title = title
    ) +
    facet_grid(sex ~ alpha_outlier, scales = "free", space = "free") +
    scale_color_manual(values = setNames(
      colors,
      unique(data$country)
    ), name = "Country") 
  
  # Apply ylim only if limits are provided
  if (!is.null(y_limits)) {
    p <- p + ylim(y_limits)
  }
  
  if(grepl('Deferral rate', title)){
    p <- p + geom_line(
      aes(y = current_deferral_rate, color = country, linetype = "Current"), size=0.5,alpha=0.5
    )+ scale_linetype_manual(
    name = "Algorithm",
    values = c("Alternative" = "solid", "Current" = "dotted")
  )
  } else if (grepl('Increase in number of donations', title)){
    p <- p + geom_line(
      aes(y = max_increase_donations, color = country, linetype = "Maximum"),
      size = 0.5, alpha=0.5
    ) +
      scale_linetype_manual(
        name = " ",
        values = c("Maximum" = "dashed", 'Alternative'="solid")
      )  
  } else{
    p <- p + guides(linetype = "none")
  }

  return(p)
}

# Create and save combined plots
plots <- list(
  "increase_donations" = template_plot(combined_data, "alpha_mean", "extra_donations_perc", "Percentage", "Increase in number of donations", c(-5, 10)),
  "deferral_rate" = template_plot(combined_data, "alpha_mean", "(n_def_alt / total) * 100", "Percentage", "Deferral rate", c(0, 10)),
  "perc_mean_deferrals" = template_plot(combined_data, "alpha_mean", "perc_mean_def", "Percentage", "Percentage mean deferrals", c(0, 100)),
  "perc_outlier_deferrals" = template_plot(combined_data, "alpha_mean", "perc_outlier_def", "Percentage", "Percentage outlier deferrals", c(0, 100)),
  "deferral_reduction" = template_plot(combined_data, "alpha_mean", "deferral_reduction_perc", "Percentage", "Deferral reduction", c(-50, 100)),
  "ineligible_donations" = template_plot(combined_data, "alpha_mean", "wrong_donations_prevented_perc", "Percentage", "Ineligible donations", c(0, 20)),
  "mean mean hb all donors alt strat" = template_plot(combined_data, "alpha_mean", "mean_mean_Hb_donations", "Hb (g/dL)", "Mean of mean Hb in all donors donating with alternative strategy", c(13, 17)),
  "mean mean hb new donors alt strat" = template_plot(combined_data, "alpha_mean", "mean_mean_Hb_FN", "Hb (g/dL)", "Mean of mean Hb in new donors donating with alternative strategy", c(12, 15))
)

for (plot_name in names(plots)) {
  ggsave(filename = file.path(results_path, paste0('combined/', plot_name, ".png")), plot = plots[[plot_name]], width = 10, height = 8)
}
