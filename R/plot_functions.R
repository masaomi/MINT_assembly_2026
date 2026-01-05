# plot_functions.R
# Visualization functions for assembly results

library(ggplot2)
library(dplyr)

#' Create contig length histogram
#' @param lengths Vector of contig lengths
#' @param title Plot title
#' @param bins Number of bins for histogram
#' @return ggplot object
plot_length_histogram <- function(lengths, title = "Contig Length Distribution", bins = 30) {
  if (length(lengths) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available") +
             theme_void())
  }
  
  df <- data.frame(length = lengths)
  
  p <- ggplot(df, aes(x = length)) +
    geom_histogram(bins = bins, fill = "#3498db", color = "#2980b9", alpha = 0.7) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = title,
      x = "Contig Length (bp)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Create contig length histogram with log scale
#' @param lengths Vector of contig lengths
#' @param title Plot title
#' @param bins Number of bins
#' @return ggplot object
plot_length_histogram_log <- function(lengths, title = "Contig Length Distribution (Log Scale)", bins = 30) {
  if (length(lengths) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available") +
             theme_void())
  }
  
  df <- data.frame(length = lengths)
  
  p <- ggplot(df, aes(x = length)) +
    geom_histogram(bins = bins, fill = "#9b59b6", color = "#8e44ad", alpha = 0.7) +
    scale_x_log10(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = title,
      x = "Contig Length (bp, log scale)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Create cumulative length plot (Nx plot)
#' @param lengths Vector of contig lengths
#' @param title Plot title
#' @return ggplot object
plot_cumulative_length <- function(lengths, title = "Cumulative Assembly Length (Nx Plot)") {
  if (length(lengths) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available") +
             theme_void())
  }
  
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  cumsum_lengths <- cumsum(sorted_lengths)
  total_length <- sum(sorted_lengths)
  
  # Calculate percentage
  cumsum_percent <- cumsum_lengths / total_length * 100
  
  df <- data.frame(
    contig_index = seq_along(sorted_lengths),
    cumsum_length = cumsum_lengths,
    cumsum_percent = cumsum_percent,
    length = sorted_lengths
  )
  
  # Find N50 position
  n50_idx <- which(cumsum_percent >= 50)[1]
  n50_value <- sorted_lengths[n50_idx]
  
  p <- ggplot(df, aes(x = cumsum_percent, y = length)) +
    geom_line(color = "#e74c3c", linewidth = 1.2) +
    geom_hline(yintercept = n50_value, 
               linetype = "dashed", color = "#2ecc71", linewidth = 0.8) +
    annotate("text", x = 5, y = n50_value, 
             label = "N50", vjust = -0.5, color = "#2ecc71", fontface = "bold") +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = title,
      x = "Cumulative Percentage (%)",
      y = "Contig Length (bp)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Create contig count by length category
#' @param lengths Vector of contig lengths
#' @param title Plot title
#' @return ggplot object
plot_length_categories <- function(lengths, title = "Contigs by Length Category") {
  if (length(lengths) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available") +
             theme_void())
  }
  
  # Define categories
  categories <- c("< 1 kb", "1-5 kb", "5-10 kb", "10-50 kb", "50-100 kb", "> 100 kb")
  
  df <- data.frame(
    category = factor(
      dplyr::case_when(
        lengths < 1000 ~ "< 1 kb",
        lengths < 5000 ~ "1-5 kb",
        lengths < 10000 ~ "5-10 kb",
        lengths < 50000 ~ "10-50 kb",
        lengths < 100000 ~ "50-100 kb",
        TRUE ~ "> 100 kb"
      ),
      levels = categories
    )
  )
  
  category_counts <- as.data.frame(table(df$category))
  names(category_counts) <- c("Category", "Count")
  
  p <- ggplot(category_counts, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = title,
      x = "Length Category",
      y = "Number of Contigs"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Create comparison bar plot for multiple assemblies
#' @param comparison_df Data frame with comparison data
#' @param metric Which metric to plot
#' @param title Plot title
#' @return ggplot object
plot_comparison <- function(comparison_df, metric = "n50", title = NULL) {
  if (nrow(comparison_df) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No comparison data available") +
             theme_void())
  }
  
  metric_labels <- list(
    n50 = "N50 (bp)",
    l50 = "L50 (contigs)",
    total_contigs = "Total Contigs",
    total_length = "Total Length (bp)",
    filtered_length = "Filtered Length (bp)",
    filtered_contigs = "Filtered Contigs",
    largest_contig = "Largest Contig (bp)"
  )
  
  if (is.null(title)) {
    title <- paste("Comparison:", metric_labels[[metric]])
  }
  
  if (!metric %in% names(comparison_df)) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = paste("Metric not found:", metric)) +
             theme_void())
  }
  
  # Create plot data
  plot_data <- comparison_df
  plot_data$y_value <- plot_data[[metric]]
  
  p <- ggplot(plot_data, aes(x = run_id, y = y_value, fill = run_id)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = scales::comma(y_value)), vjust = -0.5, size = 3.5) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = title,
      x = "Run",
      y = metric_labels[[metric]]
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Create multi-metric comparison plot
#' @param comparison_df Data frame with comparison data
#' @return ggplot object (faceted)
plot_multi_comparison <- function(comparison_df) {
  if (nrow(comparison_df) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No comparison data available") +
             theme_void())
  }
  
  # Reshape data for faceting
  metrics_to_plot <- c("n50", "total_contigs", "filtered_length")
  metric_labels <- c("N50 (bp)", "Total Contigs", "Total Length (bp)")
  
  plot_data <- comparison_df %>%
    tidyr::pivot_longer(
      cols = all_of(metrics_to_plot),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = factor(metric, levels = metrics_to_plot, labels = metric_labels))
  
  p <- ggplot(plot_data, aes(x = run_id, y = value, fill = run_id)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Assembly Comparison",
      x = "Run",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}
