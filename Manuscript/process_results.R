library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(writexl)

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
combined_characteristics <-  data.frame()

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
    cf_processed$mean_Hb_all <- ((cf_processed$TN*cf_processed$mean_Hb_TN)+(cf_processed$FN*cf_processed$mean_Hb_FN)+(cf_processed$TP*cf_processed$mean_Hb_TP)+(cf_processed$FP*cf_processed$mean_Hb_FP))/(cf_processed$total)
    
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
  
  characteristics[nrow(characteristics) + 1,] = list("country",country)
  characteristics_t <- setNames(data.frame(t(characteristics[,-1])), characteristics[,1])
  combined_characteristics <- bind_rows(combined_characteristics, characteristics_t)
}

# process the characteristics file
conversion_factor_mmol <- 0.6206 
conversion_factor_gl <- 10

convert_cols <- c("Cutoff_m", "Cutoff_f", "mean_std_dev_meas_M", "mean_std_dev_meas_F")

combined_characteristics[convert_cols] <- lapply(combined_characteristics[convert_cols], function(col) {
  col <- as.numeric(col) 
  col <- ifelse(combined_characteristics$units == "mmol/L", col / conversion_factor_mmol, col)
  col <- ifelse(combined_characteristics$units == "g/L", col / conversion_factor_gl, col)
  return(col) 
})

combined_characteristics$units <- ifelse(combined_characteristics$units == "mmol/L", "g/dL", combined_characteristics$units)
combined_characteristics$units <- ifelse(combined_characteristics$units == "g/L", "g/dL", combined_characteristics$units)

#process the combined data file to include only the results from parameters we selected
combined_data <- combined_data[combined_data$alpha_mean==0.5 & combined_data$alpha_outlier==0.999,]

#save both files
write_xlsx(combined_data, paste0(results_path, "/combined/combined_data.xlsx"))
write_xlsx(combined_characteristics, paste0(results_path, "/combined/combined_characteristics.xlsx"))

#results table
# Get the list of unique countries
countries <- unique(combined_characteristics$country)

# Initialize an empty list to store results
country_data <- list()

# Loop over each country
for (country in countries) {
  daterange <- as.character(c(combined_characteristics$daterange_min[combined_characteristics$country == country],combined_characteristics$daterange_max[combined_characteristics$country == country]))
  
  donors <- round(sum(as.numeric(combined_characteristics$donors_M[combined_characteristics$country == country]),as.numeric(combined_characteristics$donors_F[combined_characteristics$country == country])),0)
  
  donations <- round(sum(as.numeric(combined_characteristics$donations_M[combined_characteristics$country == country]),as.numeric(combined_characteristics$donations_F[combined_characteristics$country == country])),0)
  
  current_deferral <- round((sum(combined_data$n_true_def[combined_data$country == country & combined_data$sex == "M"], combined_data$n_true_def[combined_data$country == country & combined_data$sex == "F"]) / sum(combined_data$total[combined_data$country == country & combined_data$sex == "M"], combined_data$total[combined_data$country == country & combined_data$sex == "F"])) * 100, 2)
  
  alt_deferral <- round((sum(combined_data$n_def_alt[combined_data$country == country & combined_data$sex == "M"], combined_data$n_def_alt[combined_data$country == country & combined_data$sex == "F"]) / sum(combined_data$total[combined_data$country == country & combined_data$sex == "M"], combined_data$total[combined_data$country == country & combined_data$sex == "F"])) * 100,2)
  
  change_def <- round(((alt_deferral - current_deferral) / current_deferral) * 100,2)
  
  ineligible_donations <- round((sum(combined_data$FP[combined_data$country == country & combined_data$sex == "M"], combined_data$FP[combined_data$country == country & combined_data$sex == "F"]) / sum(combined_data$total[combined_data$country == country & combined_data$sex == "M"], combined_data$total[combined_data$country == country & combined_data$sex == "F"])) * 100, 2)
  
  change_donations <- round((((sum(combined_data$total_donations_new[combined_data$country == country & combined_data$sex == "M"], combined_data$total_donations_new[combined_data$country == country & combined_data$sex == "F"]) - sum(combined_data$total_donations_old[combined_data$country == country & combined_data$sex == "M"], combined_data$total_donations_old[combined_data$country == country & combined_data$sex == "F"])) /sum(combined_data$total_donations_old[combined_data$country == country & combined_data$sex == "M"], combined_data$total_donations_old[combined_data$country == country & combined_data$sex == "F"]))) * 100,2)
  
  mean_males <- round(combined_data$mean_mean_Hb_new[combined_data$sex=="M"&combined_data$country==country],1)
  mean_females <- round(combined_data$mean_mean_Hb_new[combined_data$sex=="F"&combined_data$country==country],1)

# Store values in a named vector for the country
country_data[[country]] <- c(daterange, donors, donations, current_deferral, alt_deferral, change_def, ineligible_donations, change_donations, mean_males, mean_females)
}

# Convert the list into a dataframe
results_table <- as.data.frame(do.call(cbind, country_data))

# Add row names to specify the metrics
rownames(results_table) <- c("Date min", "Date max", "Donors", "Donations", "Current Deferral", "Alternative Deferral", "Change in Deferral", "Ineligible Donations", "Change in Donations", "Mean new males", "Mean new females")

# Print or return the dataframe
print(results_table)
write_xlsx(results_table, paste0(results_path, "/combined/results_table.xlsx"))
