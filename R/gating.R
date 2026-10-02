# Statistical gating and gate-table calculations.
# These functions intentionally retain the app's existing cutoff behavior,
# including the two-component GMM fallback for negatively skewed data.

remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
  low_threshold <- quantile(target, low_percentile / 100)
  high_threshold <- quantile(target, high_percentile / 100)
  target[target >= low_threshold & target <= high_threshold]
}

keep_marker_inlier_rows <- function(data, markers) {
  keep_rows <- rep(TRUE, nrow(data))
  markers <- unique(markers[markers %in% names(data)])

  for (marker in markers) {
    inlier_values <- remove_outliers2(data[[marker]])
    keep_rows <- keep_rows & data[[marker]] %in% inlier_values
  }

  keep_rows[is.na(keep_rows)] <- FALSE
  keep_rows
}

skew_gate <- function(x, alpha = 0.01) {
  sk <- moments::skewness(x)

  n <- length(x)
  a <- locmodes(x)$locations[1]
  b <- max(x)

  if (sk < 0) {
    message("The skewness is negative!")
    outputGMM <- Mclust(x, G = 2)
    gmm_gate <- mean(outputGMM$parameters$mean)

    n_removed <- sum(x > gmm_gate)
    perc_removed <- round(n_removed / n, 3)

    return(list(
      skewness = sk,
      cutoff = gmm_gate,
      N = n,
      N_removed = n_removed,
      percentage_removed = perc_removed,
      returnvalplot = x[which(x > a + (b - a) / 2)]
    ))
  }

  iteration <- 0
  while (abs(sk) > alpha & iteration <= 100) {
    if (sk >= 0) {
      b <- a + (b - a) / 2
    } else {
      a <- a + (b - a) / 2
    }
    a <- min(a, b)
    b <- max(a, b)

    sk <- skewness(x[which(x < b)])
    if (is.nan(sk)) {
      message("Warning: skewness is NaN")
      break
    }
    iteration <- iteration + 1
  }

  n_removed <- sum(x > b)
  perc_removed <- round(n_removed / n, 3)

  print(n_removed)
  print(perc_removed)
  list(
    skewness = sk,
    cutoff = b,
    N = n,
    N_removed = n_removed,
    percentage_removed = perc_removed,
    returnvalplot = x[which(x > a + (b - a) / 2)]
  )
}

# The multi-marker path removes percentile outliers; the single-marker path
# historically does not. Keep those distinct until the scientific policy changes.
get_gates_csv <- function(dataframe, session) {
  data <- dataframe
  unique_imageids <- unique(data$imageid)
  alpha <- 0.01
  results_df <- data.frame(
    Patient = character(), Marker = character(), Gate = numeric(),
    stringsAsFactors = FALSE
  )

  total_iterations <- length(unique_imageids) * ncol(data)
  progress <- Progress$new(session, min = 1, max = total_iterations)
  progress$set(message = "Randkluft in action...", value = 0)

  for (imageid in unique_imageids) {
    sub_data <- data[data$imageid == imageid, -1]
    for (col_idx in 1:ncol(sub_data)) {
      col_name <- colnames(sub_data)[col_idx]
      values <- sub_data[, col_idx]
      values <- remove_outliers2(values)
      values <- values[!is.infinite(values)]
      result <- skew_gate(values, alpha)$cutoff
      results_df <- rbind(results_df, data.frame(Patient = imageid, Marker = col_name, Gate = result))
      cat("Image ID:", imageid, "Marker:", col_name, "Result:", result, "\n")
      current_iteration <- (match(imageid, unique_imageids) - 1) * ncol(sub_data) + col_idx
      progress$set(message = "Randkluft in Action...", value = current_iteration)
    }
  }

  progress$close()
  results_df
}

get_gates_csv_single <- function(dataframe) {
  data <- dataframe
  markers <- colnames(data)[-1]
  results_df <- data.frame(
    Patient = character(), Marker = character(), Gate = numeric(),
    stringsAsFactors = FALSE
  )

  unique_imageids <- unique(data$imageid)
  for (imageid in unique_imageids) {
    for (marker in markers) {
      subset_data <- data[data$imageid == imageid, c("imageid", marker)]
      result_to_plot <- skew_gate(subset_data[[marker]])
      cutoff_value <- result_to_plot$cutoff
      results_df <- rbind(results_df, data.frame(Patient = imageid, Marker = marker, Gate = cutoff_value))
      cat("Image ID:", imageid, "Marker:", marker, "Result:", cutoff_value, "\n")
    }
  }

  results_df
}

determine_positivity <- function(marker_intensity, gate_value) {
  if (marker_intensity > gate_value) "+" else "-"
}
