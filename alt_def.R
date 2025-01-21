library(dplyr) # for data processing
library(zoo) # for rolling mean
library(caret)
library(ggplot2)
library(ggnewscale)
library(patchwork)
codeversion <- 20250107


# The code is structed at follows:
# 1. read in donation records from DATAFILE_NAME in FILE_DIR, this file should already be preprocessed to contain the correct column names (instruction for this?)
#    this is done in the function read_data(), this takes a while because it needs to calculate the mean of all previous Hb measurements per donor
#    a further preprocessed file is stored in FILE_DIR so that this does not need to be run again in the future
# ... for each alpha_mean and alpha_outlier do:
#   2. The further preprocessed dataframe is going into the function deferral_algorithm for different alpha_means and alpha_outliers (specified below)
#      here based on the confidence interval around the mean it is decided if a donations is a deferral based on the new algorithm
#   3. The new deferrals and true (old) deferrals are compared using confusion_matrix to create TP,TN,FP,FN (positive is deferral) to assess
#      the number of deferrals in both new and old algorithm
# ...done
# 4. merge the confusion matrices for each alpha_mean and alpha_outlier and save them to csv
# 5. store a few data characterics in another file



# These values need to be changed for each country
FILE_DIR <- "./insert_path"
DATAFILE_NAME <- "fullhistory.rds"

datafile_name_withoutext <- gsub(".rds", "", DATAFILE_NAME)
DATAFILE_NAME_PREPROCESSED <- paste0(datafile_name_withoutext, "_preprocessed.rds")

# Change the deferral threshold here per country, in the same units as is in the data, i.e. mmol/L or g/dL or g/L or ...
#
CUTOFF_M <- 8.4 # mmol/L, change to cut off in your units of measurement
CUTOFF_F <- 7.8 # mmol/L, change to cut off in your units of measurement
UNITS <- "mmol/L"
# uncomment correct below
# UNITS <- 'g/dL'
# UNITS <- 'g/L'

# Do NOT change anything below, skip to the last line of the code to run the analysis.

TEST <- F # If test use only 10000 rows
NMIN_DONATIONS <- 2 # NMIN_DONATIONS indicates the number of donations disregarded from the analysis, set at 2 because at the first 2 donations there is no "mean of previous donations"
ROLLING_MEAN <- 0
SAVE_EXAMPLES <- T

# only used for outliers
if (UNITS == "mmol/L") {
  CONVERSION_FACTOR <- 1.
} else if (UNITS == "g/dL") {
  CONVERSION_FACTOR <- 1.610306
} else if (UNITS == "g/L") {
  CONVERSION_FACTOR <- 16.10306
} else {
  stop("UNITS should be one of [mmol/L, g/dL, g/L]")
}


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
                      nmin_donations = 0, rolling_mean = 0) {
  # Read the data and do some preprocessing
  cat("Reading file", file_path, "\n")
  df <- readRDS(file_path)
  if (test) {
    df <- head(df, 10000)
  }
  print("Preprocessing data, takes some time...")
  # remove donors with 2 or fewer donations, because they don't have a prev mean anyway
  cat("Removing donors with less than", nmin_donations, "donations\n")
  df <- df %>%
    group_by(KeyID) %>%
    mutate(ndons = n()) %>%
    ungroup() %>%
    filter(ndons > nmin_donations, !is.na(Hb), Hb < 16. * CONVERSION_FACTOR, Hb > 1 * CONVERSION_FACTOR, Sex %in% c("M", "F"))
  # remove really high Hb values (Hb > 16 mmol/L or 25 g/dL)

  print("Calculating true deferrals")
  df <- test_deferral(df, out_colname = "def_true", Hb_test_col = "Hb")
  print("Calculating means")
  df <- df %>%
    group_by(KeyID) %>%
    arrange(DonDate) %>%
    mutate(
      nth_don = row_number(),
      # prev_Hb = lag(Hb), #see below for much faster implement
      mean_Hb = cummean(Hb),
      # prev_mean_Hb = lag(mean_Hb),
    )
  print("Calculating lagged vars")
  df <- df %>% arrange(KeyID, DonDate)
  idx <- 1:nrow(df)
  precursor <- idx - 1
  precursor[1] <- nrow(df)
  df$prev_Hb <- NA
  df$prev_Hb <- df$Hb[precursor]
  df$prev_Hb[df$KeyID != df$KeyID[precursor]] <- NA
  df$prev_mean_Hb <- NA
  df$prev_mean_Hb <- df$mean_Hb[precursor]
  df$prev_mean_Hb[df$KeyID != df$KeyID[precursor]] <- NA

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
      prev_mean_Hb_std = std_dev_meas / sqrt(nth_don - 1),
      mean_Hb = if_else(nth_don == 1, NA, mean_Hb), # not defined, so set to NA
      prev_mean_Hb = if_else(nth_don == 2, NA, prev_mean_Hb),
    )

  if (rolling_mean > 1) {
    print("Calculating rolling means")
    df <- df %>%
      group_by(KeyID) %>%
      arrange(DonDate) %>%
      mutate(
        mean_Hb = rollapply(Hb, rolling_mean, mean,
          align = "right",
          fill = NA
        ),
        prev_mean_Hb = rollapply(prev_Hb, rolling_mean, mean,
          align = "right",
          fill = NA
        )
      )
    # need to divide by sqrt(rolling_mean) instead
    df <- df %>%
      mutate(
        mean_Hb_std = std_dev_meas / sqrt(rolling_mean),
        prev_mean_Hb_std = std_dev_meas / sqrt(rolling_mean - 1),
      )
  }

  print(df %>% group_by(Sex) %>% summarise(stddev_meas = mean(std_dev_meas)))
  print("... done with preprocessing")
  return(df)
}

sample_example_donor <- function(df, n = 1) {
  sample_donor <- df %>%
    ungroup() %>%
    filter(KeyID %in% sample(unique(KeyID), n))

  return(sample_donor)
}

plot_example <- function(sample_donor) {
  if (sample_donor$Sex[1] == "M") {
    sample_donor$cutoff <- CUTOFF_M
  } else {
    sample_donor$cutoff <- CUTOFF_F
  }
  sample_donor$low_mean <- sample_donor$mean_Hb_CI
  sample_donor$high_mean <- sample_donor$mean_Hb + (sample_donor$mean_Hb - sample_donor$mean_Hb_CI)
  sample_donor$low_outlier <- sample_donor$outlier_Hb_CI
  sample_donor$high_outlier <- sample_donor$mean_Hb + (sample_donor$prev_mean_Hb - sample_donor$outlier_Hb_CI)
  sample_donor$def_mean <- as.logical(sample_donor$def_mean)
  sample_donor$def_outlier <- as.logical(sample_donor$def_outlier)
  sample_donor$def_alt <- as.logical(sample_donor$def_alt)
  sample_donor$def_true <- as.logical(sample_donor$def_true)

  p <- ggplot(sample_donor, aes(x = nth_don, y = Hb)) +
    geom_ribbon(aes(ymin = low_mean, ymax = high_mean, fill = "Mean CI"), alpha = 0.3) +
    geom_ribbon(aes(ymin = low_outlier, ymax = high_outlier, fill = "Outlier CI"), alpha = 0.3) +
    geom_line(aes(y = cutoff, x = nth_don, linetype = "Legal Hb threshold", color = "Legal Hb threshold"), linewidth = 0.70) +
    geom_line(aes(y = mean_Hb, x = nth_don, linetype = "Mean", color = "Mean"), linewidth = 1) +
    geom_line(aes(y = Hb, x = nth_don), color = "black", linewidth = 0.5, linetype = "solid") +
    scale_color_manual(
      name = "", values = c("black", "#0072b2"),
      labels = c("Legal Hb threshold", "Mean"),
      guide = guide_legend(order = 1, title = NULL)
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Legal Hb threshold" = "solid", "Mean" = "solid"),
      guide = guide_legend(order = 1, title = NULL)
    ) +
    new_scale_color() +
    geom_point(aes(shape = def_true, color = def_alt), size = 4, stroke = 1) +
    scale_color_manual(
      name = "New algorithm", values = c("#009e73", "#d55e00"),
      labels = c("Donation", "Deferral"), guide = guide_legend(order = 4)
    ) +
    scale_shape_manual(name = "Current algorithm", values = c(16, 1), labels = c("Donation", "Deferral"), guide = guide_legend(order = 3)) +
    new_scale("shape") +
    geom_point(data = subset(sample_donor, def_outlier), aes(x = nth_don + 0.1, y = Hb, shape = "Outlier"), color = "#d55e00", size = 3, stroke = 4) +
    scale_shape_manual(values = c("Outlier" = "!"), labels = c("Outlier observation"), guide = guide_legend(order = 5, title = NULL)) +
    theme_minimal() +
    labs(x = "nth donation", y = paste("Hemoglobin", UNITS))

  return(p)
}

sample_and_plot <- function(df) {
  sample_donor <- sample_example_donor(df)
  p <- plot_example(sample_donor)
  return(p)
}

make_all_examples_and_plot <- function(df, savefile) {
  df_group <- df %>%
    filter(ndons > 3) %>%
    group_by(KeyID)

  p_FN_m <- sample_and_plot(df_group %>% filter(any(FN == T), Sex == "M")) + labs(title = "FN")
  p_FP_m <- sample_and_plot(df_group %>% filter(any(FP == T), Sex == "M")) + labs(title = "FP")
  p_TN_m <- sample_and_plot(df_group %>% filter(any(TN == T), Sex == "M")) + labs(title = "TN")
  p_TP_m <- sample_and_plot(df_group %>% filter(any(TP == T), Sex == "M")) + labs(title = "TP")

  p_FN_f <- sample_and_plot(df_group %>% filter(any(FN == T), Sex == "F")) + labs(title = "FN")
  p_FP_f <- sample_and_plot(df_group %>% filter(any(FP == T), Sex == "F")) + labs(title = "FP")
  p_TN_f <- sample_and_plot(df_group %>% filter(any(TN == T), Sex == "F")) + labs(title = "TN")
  p_TP_f <- sample_and_plot(df_group %>% filter(any(TP == T), Sex == "F")) + labs(title = "TP")

  plot_dir <- paste0(FILE_DIR, "/example_plots/")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  p_m <- p_FN_m + p_FP_m + p_TN_m + p_TP_m
  p_m <- p_m + plot_layout(guides = "collect") & theme(legend.position = "right", legend.spacing.y = unit(5, "pt"))
  ggsave(filename = paste0(plot_dir, paste0(savefile, "_m.png")), plot = p_m, height = 15, width = 20)
  p_f <- p_FN_f + p_FP_f + p_TN_f + p_TP_f
  p_f <- p_f + plot_layout(guides = "collect") & theme(legend.position = "right", legend.spacing.y = unit(5, "pt"))
  ggsave(filename = paste0(plot_dir, paste0(savefile, "_f.png")), plot = p_f, height = 15, width = 20)
}

deferral_algorithm <- function(df,
                               alpha_mean,
                               alpha_outlier,
                               cutoff_M = CUTOFF_M,
                               cutoff_F = CUTOFF_F,
                               rolling_mean = 0) {
  # Function to do the new deferral threshold algorithm
  # a fraction alpha is higher than this value
  mean_ppf <- qnorm((1 + alpha_mean) / 2) # inverse of cdf, aka ppf
  # TODO: alpha_outlier = 0, 1 does not work
  if (alpha_outlier == 0) {
    outlier_ppf <- 1000.
  } else if (alpha_outlier == 1) {
    outlier_ppf <- -1000.
  } else {
    outlier_ppf <- qnorm(alpha_outlier)
  }

  print("Calculating deferrals for alt strategy")
  cat("Using alpha for mean = ", alpha_mean, "-> percent point frac = ", mean_ppf, "\n")
  cat("Using alpha for outlier = ", alpha_outlier, "-> percent point frac = ", outlier_ppf, "\n")

  df <- df %>% mutate(prob_below_thres = if_else(Sex == "M", pnorm(CUTOFF_M, mean_Hb, mean_Hb_std), pnorm(CUTOFF_F, mean_Hb, mean_Hb_std)))
  # NOTE only downward for outliers
  df <- df %>% mutate(
    mean_Hb_CI = mean_Hb + mean_ppf * mean_Hb_std,
    outlier_Hb_CI = prev_mean_Hb - outlier_ppf * sqrt(std_dev_meas^2 + prev_mean_Hb_std^2)
  )

  df <- test_deferral(
    df,
    out_colname = "def_mean",
    Hb_test_col = "mean_Hb_CI",
    cutoff_M = cutoff_M,
    cutoff_F = cutoff_F
  )

  # outlier detect, only below
  df <- df %>% mutate(
    def_outlier = if_else(Hb < outlier_Hb_CI, 1, 0),
    def_outlier = if_else(is.na(def_outlier), 0, def_outlier),
    def_mean = if_else(is.na(def_mean), 0, def_mean)
    # def_outlier = if_else(nth_don == 1, 0, def_outlier) #no need for this but could be implemented like so
  )

  # for 1st and 2nd donation there is no prev mean, outlier possible for 1st
  # df <- df %>% mutate(def_mean = case_when(nth_don == 1 ~ def_true,
  #   nth_don == 2 ~ def_true,
  #   .default = def_mean
  # ))

  # if alpha_outliers == 0 this means only use outliers and no means
  # alpha_outliers == 1 means only use means and not outliers
  # TODO: this does not work as expected
  if (alpha_outlier == 0) {
    df <- df %>% mutate(def_alt = def_outlier)
  } else if (alpha_outlier == 1) {
    df <- df %>% mutate(def_alt = def_mean)
  } else {
    df <- df %>% mutate(
      def_alt = if_else(def_mean | def_outlier, 1, 0),
    )
  }

  df <- df %>% mutate(
    FN = (def_true == 1) & (def_alt == 0),
    FP = (def_true == 0) & (def_alt == 1),
    TN = (def_true == 0) & (def_alt == 0),
    TP = (def_true == 1) & (def_alt == 1)
  )

  # save some examples
  if (SAVE_EXAMPLES){
    make_all_examples_and_plot(df, savefile = paste0("examples_amean=", alpha_mean, "_alphaout=", alpha_outlier, "_rm=", rolling_mean))
  }

  if (rolling_mean > 1) {
    df %>% filter(nth_don > 4)
  } else {
    df <- df %>% filter(nth_don > 2)
  }

  return(df)
}


create_hist_Hb <- function(data, mask) {
  # Create histograms of Hb, mean Hb and prev mean Hb for masked data
  # also means
  data_masked <- data[mask, ]
  if (nrow(data_masked) == 0) {
    breaks <- seq(6, 12, 0.1) # Hb bins
  } else {
    breaks <- seq(
      round(min(data_masked$Hb, data_masked$mean_Hb, data_masked$prev_mean_Hb) - 0.1, 1),
      round(max(data_masked$Hb, data_masked$mean_Hb, data_masked$prev_mean_Hb) + 0.1, 1), 0.1
    )
  }
  hist_Hb <- hist(data_masked$Hb, breaks = breaks, plot = F)$counts
  mean_Hb <- mean(data_masked$Hb)
  hist_mean_Hb <- hist(data_masked$mean_Hb, breaks = breaks, plot = F)$counts
  mean_mean_Hb <- mean(data_masked$mean_Hb)
  hist_prev_mean_Hb <- hist(data_masked$prev_mean_Hb, breaks = breaks, plot = F)$counts
  mean_prev_mean_Hb <- mean(data_masked$prev_mean_Hb)
  output <- list(
    breaks = breaks, hist_Hb = hist_Hb, hist_mean_Hb = hist_mean_Hb, hist_prev_mean_Hb = hist_prev_mean_Hb,
    mean_Hb = mean_Hb, mean_mean_Hb = mean_mean_Hb, mean_prev_mean_Hb = mean_prev_mean_Hb
  )
  return(output)
}


create_aggregated_results_and_hists <- function(df) {
  mask_old <- as.logical(df$def_true)
  mask_new <- as.logical(df$def_alt)
  mask_TP <- df$TP
  mask_TN <- df$TN
  mask_FP <- df$FP
  mask_FN <- df$FN
  mask_def_mean <- as.logical(df$def_mean)
  mask_def_outlier <- as.logical(df$def_outlier)

  TP <- sum(mask_TP)
  TN <- sum(mask_TN)
  FP <- sum(mask_FP)
  FN <- sum(mask_FN)
  # checks
  # total <- TP+TN+FP+FN
  # print(total)
  n_true_def <- sum(mask_old)
  n_new_def <- sum(mask_new)
  n_mean_def <- sum(mask_def_mean)
  n_outlier_def <- sum(mask_def_outlier)
  n_both_def <- sum(mask_def_mean & mask_def_outlier)

  hist_output_TP <- create_hist_Hb(df, mask_TP)
  hist_output_TN <- create_hist_Hb(df, mask_TN)
  hist_output_FN <- create_hist_Hb(df, mask_FN)
  hist_output_FP <- create_hist_Hb(df, mask_FP)
  hist_output_new <- create_hist_Hb(df, mask_new)
  hist_output_old <- create_hist_Hb(df, mask_old)

  output <- list(
    TP = TP, TN = TN, FP = FP, FN = FN, n_true_def = n_true_def, n_new_def = n_new_def, n_mean_def = n_mean_def, n_outlier_def = n_outlier_def, n_both_def = n_both_def,
    hist_output_TP = hist_output_TP, hist_output_TN = hist_output_TN, hist_output_FN = hist_output_FN, hist_output_FP = hist_output_FP,
    hist_output_new = hist_output_new, hist_output_old = hist_output_old
  )
  return(output)
}


run_and_save <- function(test = T, save_intermediate = F, use_preprocessed = F, nmin_donations = 2, rolling_mean = 0) {
  if (file.exists(file.path(FILE_DIR, DATAFILE_NAME_PREPROCESSED)) & (!test) & use_preprocessed) {
    df <- readRDS(file.path(FILE_DIR, DATAFILE_NAME_PREPROCESSED))
    cat("Using preprocessed file ", DATAFILE_NAME_PREPROCESSED, "BEWARE!\n")
  } else {
    data_file <- file.path(FILE_DIR, DATAFILE_NAME)
    df <- read_data(data_file, test = test, nmin_donations = nmin_donations, rolling_mean = rolling_mean)
    if (!test) {
      cat("Saving preprocessed file at", DATAFILE_NAME_PREPROCESSED, "\n")
      saveRDS(df, file.path(FILE_DIR, DATAFILE_NAME_PREPROCESSED))
    }
  }

  if (test) {
    alpha_mean_test <- c(0.68)
    alpha_outlier_test <- c(0.999)
  } else {
    # the confidence (interval) that the mean is below the cutoff
    alpha_mean_test <- c(
      -0.999, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.,
      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999
    )
    # the prob of a measurement deviating from the mean
    alpha_outlier_test <- c(0.99, 0.999, 0.9999)
  }

  outputs <- list()
  outputs_hist <- list()
  i <- 1
  for (alpha_mean in alpha_mean_test) {
    for (alpha_outlier in alpha_outlier_test) {
      # First copy dataframe and set deferral status for new algorithm using alpha mean and outlier
      df_ <- data.frame(df) # copy
      df_ <- deferral_algorithm(df_, alpha_mean, alpha_outlier, rolling_mean = rolling_mean)
      if (!test) {
        if (rolling_mean > 1) {
          savefile <- file.path(FILE_DIR, paste0("results_alt_def_a1=", alpha_mean, "_a2=", alpha_outlier, "_rolling_mean", rolling_mean, ".rds"))
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
      cat("Creating results for alphas ", alpha_mean, ", ", alpha_outlier, "\n")
      output_m <- create_aggregated_results_and_hists(df_m)
      output_f <- create_aggregated_results_and_hists(df_f)

      agg_cols <- c("TP", "TN", "FP", "FN", "n_true_def", "n_mean_def", "n_outlier_def", "n_both_def")
      output_m_agg <- c(sex = "M", alpha_mean = alpha_mean, alpha_outlier = alpha_outlier, output_m[agg_cols])
      output_f_agg <- c(sex = "F", alpha_mean = alpha_mean, alpha_outlier = alpha_outlier, output_f[agg_cols])
      output_hist_m <- c(sex = "M", alpha_mean = alpha_mean, alpha_outlier = alpha_outlier, output_m[setdiff(names(output_m), agg_cols)])
      output_hist_f <- c(sex = "F", alpha_mean = alpha_mean, alpha_outlier = alpha_outlier, output_f[setdiff(names(output_f), agg_cols)])

      outputs[[i]] <- output_m_agg
      outputs[[i + 1]] <- output_f_agg
      i <- i + 2
      outputs_hist[[paste0("M_", alpha_mean, "_", alpha_outlier)]] <- output_hist_m
      outputs_hist[[paste0("F_", alpha_mean, "_", alpha_outlier)]] <- output_hist_f
      cat("\n")
    }
  }
  results <- bind_rows(outputs) # combine agg outputs to dataframe

  cat("Saving files at", FILE_DIR)
  if (!test) {
    if (rolling_mean > 1) {
      write.csv(results, file.path(FILE_DIR, paste0("confusion_matrix_results_rolling_mean_", rolling_mean, ".csv")))
      saveRDS(outputs_hist, file.path(FILE_DIR, paste0("hist_results_rolling_mean_", rolling_mean, ".rds")))
    } else {
      write.csv(results, file.path(FILE_DIR, "confusion_matrix_results.csv"))
      saveRDS(outputs_hist, file.path(FILE_DIR, "hist_results.rds"))
    }
  } else {
    write.csv(results, file.path(FILE_DIR, "test_results.csv"))
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
      "Cutoff_m",
      "Cutoff_f",
      "units",
      "nmin_donations",
      "rolling_mean",
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
      CUTOFF_M,
      CUTOFF_F,
      UNITS,
      NMIN_DONATIONS,
      ROLLING_MEAN,
      as.character(codeversion)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(tosave, file.path(FILE_DIR, "results_characteristics.csv"))
  cat("... done \n")
  return(results)
}

# WARNING if use rolling_mean:
# How to handle the 1st, 2nd, 3rd... etc. donation is NOT implemented
# They are set to NA and so the total number of donations is NOT correct
# So either do not use this, or implement some strategy for this

#Run the analysis using the line below
results <- run_and_save(test = TEST, use_preprocessed = F, save_intermediate = F, nmin_donations = NMIN_DONATIONS, rolling_mean = ROLLING_MEAN)
