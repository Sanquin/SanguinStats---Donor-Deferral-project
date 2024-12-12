library(dplyr) # for data processing
library(zoo) # for rolling mean
library(caret)
codeversion <- 20241210


# The code is structed at follows:
# 1. read in donation records from DATAFILE_NAME in FILE_DIR, this file should already be preprocessed to contain the correct column names 
#    this is done in the function read_data(), this takes a while because it needs to calculate the mean of all previous Hb measurements per donor
#    a further preprocessed file is stored in FILE_DIR so that this does not need to be run again in the future
# ... for each alpha_mean and alpha_outlier do:
#   2. The further preprocessed dataframe is going into the function deferral_algorithm for different alpha_means and alpha_outliers (specified below)
#      here based on the confidence interval around the mean it is decided if a donations is a deferral based on the new algorithm
#   3. The new deferrals and true (old) deferrals are compared using confusion_matrix to create TP,TN,FP,FN (positive is deferral) to assess
#      the number of deferrals in both new and old algorithm
# ...done
# 4. merge the confusion matrices for each alpha_mean and alpha_outlier and save them to csv
# 5. a few other population level characteristics are stored in a separate file (see word instructions for specifications)


# These values need to be changed for each country
FILE_DIR <- #datafile directory here (see word file with instructions)
DATAFILE_NAME <- "fullhistory.rds" #or change to own file name

datafile_name_withoutext <- gsub(".rds", "", DATAFILE_NAME)
DATAFILE_NAME_PREPROCESSED <- paste0(datafile_name_withoutext, "_preprocessed.rds")

# Change the deferral threshold here per country, in the same units as is in the data, i.e. mmol/L or g/dL or g/L or ...
CUTOFF_M <- 8.4 # mmol/L
CUTOFF_F <- 7.8 # mmol/L

TEST <- F # If test use only 10000 rows

# Do NOT change anything below

test_deferral <- function(df,
                          out_colname = "def",
                          Hb_test_col = "Hb",
                          cutoff_M = CUTOFF_M,
                          cutoff_F = CUTOFF_F) {
  # test if a donation is a deferral (but vectorized implementation)
  return(df %>% mutate(
    !!as.name(out_colname) := if_else(
      Sex == "M",
      if_else(!!as.name(Hb_test_col) < cutoff_M, 1, 0),
      if_else(!!as.name(Hb_test_col) < cutoff_F, 1, 0),
    )
  ))
}

read_data <- function(file_path,
                      test = T,
                      nmin_donations = 2) {
  # Read the data and do some preprocessing
  # Also does prev_min_4, but this is not used by default
  cat("Reading file", file_path, "\n")
  df <- readRDS(file_path)
  
  if (test) {
    df <- head(df, 10000)
  }
  print("Calculating true deferrals")
  df <- df %>%
    arrange(KeyID, DonDate) %>%
    filter(!is.na(Hb))
  df <- test_deferral(df, out_colname = "def_true", Hb_test_col = "Hb")
  
  print("Calculating prev means")
  # TODO: rolling mean optional?
  df <- df %>%
    group_by(KeyID) %>%
    mutate(
      nth_don = row_number(),
      prev_Hb = lag(Hb),
      mean_Hb = cumsum(Hb) / nth_don,
      prev_mean_Hb = lag(mean_Hb),
      mean_Hb_min_4 = rollapply(Hb, 4, mean,
                                align = "right",
                                fill = NA
      ),
      prev_mean_Hb_min_4 = rollapply(prev_Hb, 4, mean,
                                     align = "right",
                                     fill = NA
      ),
    )
  
  print("Calculating std dev of measurement (approx)")
  df <- df %>%
    group_by(Sex) %>%
    mutate(std_dev_meas = sd(Hb - prev_Hb,
                             na.rm =
                               T
    ) / sqrt(2)) %>%
    ungroup() %>%
    mutate(
      mean_Hb_std = std_dev_meas / sqrt(nth_don),
      mean_Hb_std_min_4 = std_dev_meas / sqrt(4),
      prev_mean_Hb_std = std_dev_meas / sqrt(nth_don - 1),
      prev_mean_Hb_min_4_std = std_dev_meas / sqrt(4)
    )
  
  print(df %>% group_by(Sex) %>% summarise(stddev_meas = mean(std_dev_meas)))
  return(df)
}

deferral_algorithm <- function(df,
                               alpha_mean,
                               alpha_outlier,
                               cutoff_M = CUTOFF_M,
                               cutoff_F = CUTOFF_F, verbose = T,
                               use_prev_4 = F) {
  # Function to do the new deferral threshold algorithm
  mean_ppf <- qnorm((1 + alpha_mean) / 2) # inverse of cdf
  outlier_ppf <- qnorm((1 + alpha_outlier) / 2)
  
  if (verbose) {
    print("Calculating deferrals for alt strategy")
    cat("Using alpha for mean = ", alpha_mean, "-> percent point frac = ", mean_ppf, "\n")
    cat("Using alpha for outlier = ", alpha_outlier, "-> percent point frac = ", outlier_ppf, "\n")
  }
  
  df <- df %>% mutate(prob_below_thres = if_else(Sex == "M", pnorm(CUTOFF_M, mean_Hb, mean_Hb_std), pnorm(CUTOFF_F, mean_Hb, mean_Hb_std)))
  # NOTE only downward for outliers
  if (use_prev_4) {
    df <- df %>% mutate(
      prev_mean_Hb_CI = mean_Hb_min_4 + mean_ppf * mean_Hb_min_4_std,
      outlier_Hb_CI = prev_mean_Hb_min_4 - outlier_ppf * std_dev_meas
    )
  } else {
    df <- df %>% mutate(
      prev_mean_Hb_CI = mean_Hb + mean_ppf * mean_Hb_std,
      outlier_Hb_CI = prev_mean_Hb - outlier_ppf * std_dev_meas
    )
  }
  
  df <- test_deferral(
    df,
    out_colname = "def_mean",
    Hb_test_col = "prev_mean_Hb_CI",
    cutoff_M = cutoff_M,
    cutoff_F = cutoff_F
  )
  df <- df %>% mutate(def_mean = if_else(nth_don == 1, def_true, def_mean))
  
  # outlier detect
  # only below
  df <- df %>% mutate(
    def_outlier = if_else(Hb < outlier_Hb_CI, 1, 0),
    def_outlier = if_else(nth_don == 1, 0, def_outlier)
  )
  
  df <- df %>% mutate(
    def_OR = if_else(def_mean | def_outlier, 1, 0),
    def_AND = if_else(def_mean & def_outlier, 1, 0)
  )
  
  df <- df %>% mutate(
    Hb_wrt_thres = ifelse(Sex == "M", Hb - CUTOFF_M, Hb - CUTOFF_F),
    FN = (def_true == 1) & (def_OR == 0),
    FP = (def_true == 0) & (def_OR == 1),
    TN = (def_true == 0) & (def_OR == 0),
    TP = (def_true == 1) & (def_OR == 1)
  )
  
  return(df)
}



create_confusion_matrix_data <- function(df, col, col_true = "def_true") {
  # create confusion matrix and store true/false positives/negatives counts in dataframe
  # Positive = deferral, negative = no deferral
  cf <- as.data.frame(confusionMatrix(as.factor(df[, col]), as.factor(df[, col_true]),
                                      dnn = c("alt", "true"), positive = "1"
  )$table)
  
  # Complicated way of transposing column to rows...
  # and convert to dataframe
  cf_df <- cbind(c("TN", "FP", "FN", "TP"), cf$Freq)
  cf_df <- as.data.frame(t(cf_df))
  colnames(cf_df) <- cf_df[1, ]
  return(cf_df[-1, ]) # drop column CF
}


create_and_merge_confusion_matrix_and_save_means <- function(df, alpha_mean, alpha_outlier, sex = "") {
  # Merge confusion matrix dataframes from create_confusion_matrix_data for different columns
  cfs <- list()
  def_col_names <- c("def_mean", "def_outlier", "def_OR", "def_AND")
  
  # loop over different types of deferral strageties based on mean or outlier and the combination OR/AND
  # calculate confusion matrix (FN, FP, TN, TP) for each strategy
  # Positive = deferral, negative = no deferral
  for (i in seq_along(def_col_names)) {
    cf <- create_confusion_matrix_data(df, col = def_col_names[i])
    cf[, "type"] <- def_col_names[i]
    cfs[[i]] <- cf
  }
  cf_df <- bind_rows(cfs)
  cf_df[, "alpha_mean"] <- alpha_mean
  cf_df[, "alpha_outlier"] <- alpha_outlier
  cf_df[, "sex"] <- sex
  # calculate mean wrt threshold for false negatives and store
  # And store a few more things why not
  #TODO: make functions of these. Sorry for all the repeated code now...
  cf_df$mean_mean_Hb_FN <- mean(df$mean_Hb[df$FN == T], na.rm = T)
  cf_df$std_mean_Hb_FN <- sd(df$mean_Hb[df$FN == T], na.rm = T)
  cf_df$mean_prob_below_thres_FN <- mean(df$prob_below_thres[df$FN == T], na.rm = T)
  cf_df$std_prob_below_thres_fn <- sd(df$prob_below_thres[df$FN == T], na.rm = T)
  
  cf_df$mean_mean_Hb_TN <- mean(df$mean_Hb[df$TN == T], na.rm = T)
  cf_df$std_mean_Hb_TN <- sd(df$mean_Hb[df$TN == T], na.rm = T)
  cf_df$mean_prob_below_thres_TN <- mean(df$prob_below_thres[df$TN == T], na.rm = T)
  cf_df$std_prob_below_thres_TN <- sd(df$prob_below_thres[df$TN == T], na.rm = T)
  
  cf_df$mean_mean_Hb_FP <- mean(df$mean_Hb[df$FP == T], na.rm = T)
  cf_df$std_mean_Hb_FP <- sd(df$mean_Hb[df$FP == T], na.rm = T)
  cf_df$mean_prob_below_thres_FP <- mean(df$prob_below_thres[df$FP == T], na.rm = T)
  cf_df$std_prob_below_thres_FP <- sd(df$prob_below_thres[df$FP == T], na.rm = T)
  
  cf_df$mean_mean_Hb_TP <- mean(df$mean_Hb[df$TP == T], na.rm = T)
  cf_df$std_mean_Hb_TP <- sd(df$mean_Hb[df$TP == T], na.rm = T)
  cf_df$mean_prob_below_thres_TP <- mean(df$prob_below_thres[df$TP == T], na.rm = T)
  cf_df$std_prob_below_thres_TP <- sd(df$prob_below_thres[df$TP == T], na.rm = T)
  
  if (sex == "M") {
    mask_below <- df$mean_Hb < CUTOFF_M
  } else {
    mask_below <- df$mean_Hb < CUTOFF_F
  }
  
  #below thres for FN,FP,TN,TP
  cf_df$n_mean_below_thres_FN <- sum(mask_below & (df$FN == T))
  cf_df$below_thres_mean_mean_Hb_FN <- mean(df$mean_Hb[(df$FN==T)&mask_below], na.rm=T)
  cf_df$below_thres_std_mean_Hb_FN <- sd(df$mean_Hb[(df$FN==T)&mask_below], na.rm=T)
  cf_df$n_mean_below_thres_FP <- sum(mask_below & (df$FP == T))
  cf_df$below_thres_mean_mean_Hb_FP <- mean(df$mean_Hb[(df$FP==T)&mask_below], na.rm=T)
  cf_df$below_thres_std_mean_Hb_FP <- sd(df$mean_Hb[(df$FP==T)&mask_below], na.rm=T)
  cf_df$n_mean_below_thres_TP <- sum(mask_below & (df$TP == T))
  cf_df$below_thres_mean_mean_Hb_TP <- mean(df$mean_Hb[(df$TP==T)&mask_below], na.rm=T)
  cf_df$below_thres_std_mean_Hb_TP <- sd(df$mean_Hb[(df$TP==T)&mask_below], na.rm=T)
  cf_df$n_mean_below_thres_TN <- sum(mask_below & (df$TN == T))
  cf_df$below_thres_mean_mean_Hb_TN <- mean(df$mean_Hb[(df$TN==T)&mask_below], na.rm=T)
  cf_df$below_thres_std_mean_Hb_TN <- sd(df$mean_Hb[(df$TN==T)&mask_below], na.rm=T)
  
  #below thres old
  cf_df$n_mean_below_thres_old <- sum(mask_below & ((df$TN == T) | (df$FP==T)))
  cf_df$below_thres_mean_mean_Hb_old <- mean(df$mean_Hb[((df$TN==T)|(df$FP==T))&mask_below], na.rm=T)
  cf_df$below_thres_std_mean_Hb_old <- sd(df$mean_Hb[((df$TN==T)|(df$FP==T))&mask_below], na.rm=T)
  #below thres for new
  cf_df$n_mean_below_thres_new <- sum(mask_below & ((df$TN == T) | (df$FN==T)))
  cf_df$below_thres_mean_mean_Hb_new <- mean(df$mean_Hb[((df$TN==T)|(df$FN==T))&mask_below], na.rm=T)
  cf_df$below_thres_std_mean_Hb_new <- sd(df$mean_Hb[((df$TN==T)|(df$FN==T))&mask_below], na.rm=T)
  
  return(cf_df)
}

run_and_save <- function(use_prev_4 = F, test = T, save_intermediate = F) {
  if (file.exists(file.path(FILE_DIR, DATAFILE_NAME_PREPROCESSED))) {
    df <- readRDS(file.path(FILE_DIR, DATAFILE_NAME_PREPROCESSED))
    if (test) {
      df <- head(df, 10000)
    }
  } else {
    data_file <- file.path(FILE_DIR, DATAFILE_NAME)
    df <- read_data(data_file, test = test)
    if (!test) {
      saveRDS(df, file.path(FILE_DIR, DATAFILE_NAME_PREPROCESSED))
    }
  }
  
  # filter out first donations?
  # df <- df %>% filter(nth_don > 1)
  
  if (test) {
    alpha_mean_test <- c(0.5, 0.99)
    alpha_outlier_test <- c(0.9, 0.999)
  } else {
    alpha_mean_test <- c(
      -0.999, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.,
      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999
    )
    alpha_outlier_test <- c(0.99, 0.999, 0.9999)
  }
  cfs <- list()
  i <- 1
  for (alpha_mean in alpha_mean_test) {
    for (alpha_outlier in alpha_outlier_test) {
      # TODO: summary stats for all combinations of deferrals. confusion matrices...
      
      # First copy dataframe and set deferral status for new algorithm using alpha mean and outlier
      df_ <- data.frame(df) # copy
      df_ <- deferral_algorithm(df_, alpha_mean, alpha_outlier, use_prev_4 = use_prev_4)
      if (!test) {
        if (use_prev_4) {
          savefile <- file.path(FILE_DIR, paste0("results_alt_def_a1=", alpha_mean, "_a2=", alpha_outlier, "_prev_4", ".rds"))
        } else {
          savefile <- file.path(FILE_DIR, paste0("results_alt_def_a1=", alpha_mean, "_a2=", alpha_outlier, ".rds"))
        }
        if (save_intermediate) {
          cat("Saving at", savefile, "\n")
          saveRDS(df_, savefile)
        }
      }
      
      # add to list and later rbind rows
      df_m <- df_ %>% filter(Sex == "M")
      df_f <- df_ %>% filter(Sex == "F")
      
      # Create confusion matrix (TP,TF,FN,FP) for each strategy for this alphas
      cat("Creating confusion matrix\n")
      cf_m <- create_and_merge_confusion_matrix_and_save_means(df_m, alpha_mean, alpha_outlier, sex = "M")
      cf_f <- create_and_merge_confusion_matrix_and_save_means(df_f, alpha_mean, alpha_outlier, sex = "F")
      cfs[[i]] <- cf_m
      cfs[[i + 1]] <- cf_f
      i <- i + 2
      cat("\n")
    }
  }
  
  result <- bind_rows(cfs)
  if (!test) {
    if (use_prev_4) {
      write.csv(result, file.path(FILE_DIR, "confusion_matrix_results_prev_4.csv"))
    } else {
      write.csv(result, file.path(FILE_DIR, "confusion_matrix_results.csv"))
    }
  }
  
  tosave <- data.frame(
    Metric = c(
      "daterange_min",
      "daterange_max",
      "donations_M",
      "donations_F",
      "donors_M",
      "donors_F",
      "mean_std_dev_meas_M",
      "mean_std_dev_meas_F",
      "InputFile",
      "codeversion"
    ),
    Value = c(
      as.character(min(df$DonDate)),
      as.character(max(df$DonDate)),
      sum(df$Sex == "M"),
      sum(df$Sex == "F"),
      length(unique(df$KeyID[df$Sex == "M"])),
      length(unique(df$KeyID[df$Sex == "F"])),
      mean(df$std_dev_meas[df$Sex == "M"]),
      mean(df$std_dev_meas[df$Sex == "F"]),
      as.character(file.path(FILE_DIR, DATAFILE_NAME)),
      as.character(codeversion)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(tosave, file.path(FILE_DIR, "results_characteristics.csv"))
  cat('... done \n')
  return(result)
}

# WARNING if use prev_4:
# How to handle the 1st, 2nd, 3rd donation is NOT implemented
# They are set to NA and so the total number of donations is NOT correct
# So either do not use this, or implement some strategy for this



## mean of mean Hb (previous) of false negatives
# prob of Hb below thres?
# totaal aantal donaties

result <- run_and_save(test = TEST, use_prev_4 = F)
