# parse_results.R
# Assembly result parsing and statistics calculation

#' Parse GFA file and extract sequences
parse_gfa <- function(gfa_file) {
  if (!file.exists(gfa_file)) {
    stop(paste("GFA file not found:", gfa_file))
  }
  
  lines <- readLines(gfa_file, warn = FALSE)
  segment_lines <- lines[startsWith(lines, "S\t")]
  
  if (length(segment_lines) == 0) {
    return(data.frame(name = character(0), sequence = character(0), length = integer(0)))
  }
  
  segments <- lapply(segment_lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      return(data.frame(name = parts[2], sequence = parts[3], 
                        length = nchar(parts[3]), stringsAsFactors = FALSE))
    }
    return(NULL)
  })
  
  segments <- do.call(rbind, segments[!sapply(segments, is.null)])
  
  if (nrow(segments) > 0) {
    segments <- segments[order(-segments$length), ]
    rownames(segments) <- NULL
  }
  
  return(segments)
}

#' Calculate Nx and Lx
calculate_nx_lx <- function(lengths, x = 50) {
  if (length(lengths) == 0) return(list(Nx = 0, Lx = 0))
  
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  total_length <- sum(sorted_lengths)
  target <- total_length * (x / 100)
  cumsum_lengths <- cumsum(sorted_lengths)
  idx <- which(cumsum_lengths >= target)[1]
  
  return(list(Nx = sorted_lengths[idx], Lx = idx))
}

#' Calculate comprehensive assembly statistics
calculate_assembly_stats <- function(gfa_file, min_length = 500) {
  segments <- parse_gfa(gfa_file)
  
  if (nrow(segments) == 0) {
    return(list(total_contigs = 0, total_length = 0, filtered_contigs = 0,
                filtered_length = 0, n50 = 0, l50 = 0, n90 = 0, l90 = 0,
                largest_contig = 0, smallest_contig = 0, mean_length = 0,
                median_length = 0, gc_content = NA))
  }
  
  all_lengths <- segments$length
  filtered <- segments[segments$length >= min_length, ]
  filtered_lengths <- filtered$length
  
  n50_l50 <- calculate_nx_lx(filtered_lengths, 50)
  n90_l90 <- calculate_nx_lx(filtered_lengths, 90)
  
  gc_content <- NA
  if (nrow(filtered) > 0) {
    tryCatch({
      all_seq <- paste(filtered$sequence, collapse = "")
      gc_count <- nchar(gsub("[^GCgc]", "", all_seq))
      total_bases <- nchar(gsub("[^ATGCatgc]", "", all_seq))
      if (total_bases > 0) gc_content <- round(gc_count / total_bases * 100, 2)
    }, error = function(e) { gc_content <- NA })
  }
  
  list(
    total_contigs = nrow(segments),
    total_length = sum(all_lengths),
    filtered_contigs = nrow(filtered),
    filtered_length = sum(filtered_lengths),
    n50 = n50_l50$Nx,
    l50 = n50_l50$Lx,
    n90 = n90_l90$Nx,
    l90 = n90_l90$Lx,
    largest_contig = ifelse(length(filtered_lengths) > 0, max(filtered_lengths), 0),
    smallest_contig = ifelse(length(filtered_lengths) > 0, min(filtered_lengths), 0),
    mean_length = ifelse(length(filtered_lengths) > 0, round(mean(filtered_lengths), 0), 0),
    median_length = ifelse(length(filtered_lengths) > 0, round(median(filtered_lengths), 0), 0),
    gc_content = gc_content
  )
}

#' Format statistics as a data frame for display
format_stats_table <- function(stats) {
  data.frame(
    Metric = c("Total Contigs", "Total Length (bp)", "Filtered Contigs (>= min length)",
               "Filtered Length (bp)", "N50 (bp)", "L50", "N90 (bp)", "L90",
               "Largest Contig (bp)", "Smallest Contig (bp)", "Mean Length (bp)",
               "Median Length (bp)", "GC Content (%)"),
    Value = c(format(stats$total_contigs, big.mark = ","),
              format(stats$total_length, big.mark = ","),
              format(stats$filtered_contigs, big.mark = ","),
              format(stats$filtered_length, big.mark = ","),
              format(stats$n50, big.mark = ","),
              stats$l50,
              format(stats$n90, big.mark = ","),
              stats$l90,
              format(stats$largest_contig, big.mark = ","),
              format(stats$smallest_contig, big.mark = ","),
              format(stats$mean_length, big.mark = ","),
              format(stats$median_length, big.mark = ","),
              ifelse(is.na(stats$gc_content), "N/A", paste0(stats$gc_content, "%"))),
    stringsAsFactors = FALSE
  )
}

#' Convert GFA to FASTA format
gfa_to_fasta <- function(gfa_file, output_fasta = NULL, min_length = 0) {
  segments <- parse_gfa(gfa_file)
  if (nrow(segments) == 0) stop("No sequences found in GFA file")
  
  if (min_length > 0) segments <- segments[segments$length >= min_length, ]
  if (nrow(segments) == 0) stop("No sequences remaining after length filter")
  
  if (is.null(output_fasta)) output_fasta <- sub("\\.gfa$", ".fasta", gfa_file)
  
  fasta_lines <- character(nrow(segments) * 2)
  for (i in seq_len(nrow(segments))) {
    fasta_lines[(i - 1) * 2 + 1] <- paste0(">", segments$name[i], " length=", segments$length[i])
    fasta_lines[(i - 1) * 2 + 2] <- segments$sequence[i]
  }
  
  writeLines(fasta_lines, output_fasta)
  return(output_fasta)
}

#' Get top contigs for preview
get_top_contigs <- function(gfa_file, n = 10, preview_length = 500) {
  segments <- parse_gfa(gfa_file)
  
  if (nrow(segments) == 0) {
    return(data.frame(Rank = integer(0), Name = character(0), 
                      Length = integer(0), Preview = character(0)))
  }
  
  top_n <- head(segments, n)
  
  data.frame(
    Rank = seq_len(nrow(top_n)),
    Name = top_n$name,
    Length = top_n$length,
    Preview = sapply(top_n$sequence, function(seq) {
      if (nchar(seq) > preview_length) paste0(substr(seq, 1, preview_length), "...")
      else seq
    }),
    stringsAsFactors = FALSE
  )
}

#' Get contig lengths for plotting
get_contig_lengths <- function(gfa_file, min_length = 0) {
  segments <- parse_gfa(gfa_file)
  if (nrow(segments) == 0) return(integer(0))
  
  lengths <- segments$length
  if (min_length > 0) lengths <- lengths[lengths >= min_length]
  return(lengths)
}

#' Create a summary string for the assembly
create_summary_string <- function(stats) {
  sprintf("Assembly: %s contigs, %s bp total, N50: %s bp",
          format(stats$filtered_contigs, big.mark = ","),
          format(stats$filtered_length, big.mark = ","),
          format(stats$n50, big.mark = ","))
}

#' Parse FASTA file and extract sequences
parse_fasta <- function(fasta_file) {
  if (!file.exists(fasta_file)) {
    stop(paste("FASTA file not found:", fasta_file))
  }
  
  lines <- readLines(fasta_file, warn = FALSE)
  
  if (length(lines) == 0) {
    return(data.frame(name = character(0), sequence = character(0), length = integer(0)))
  }
  
  # Find header lines
  header_idx <- which(startsWith(lines, ">"))
  
  if (length(header_idx) == 0) {
    return(data.frame(name = character(0), sequence = character(0), length = integer(0)))
  }
  
  # Parse each sequence
  sequences <- lapply(seq_along(header_idx), function(i) {
    start <- header_idx[i]
    end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    
    name <- sub("^>", "", lines[start])
    name <- sub(" .*", "", name)  # Take only the first part before space
    
    seq_lines <- lines[(start + 1):end]
    seq_lines <- seq_lines[!startsWith(seq_lines, ">")]
    sequence <- paste(seq_lines, collapse = "")
    
    data.frame(
      name = name,
      sequence = sequence,
      length = nchar(sequence),
      stringsAsFactors = FALSE
    )
  })
  
  result <- do.call(rbind, sequences)
  return(result)
}

#' Get top reads from FASTA for preview
get_top_reads <- function(fasta_file, n = 10, preview_length = 100) {
  tryCatch({
    reads <- parse_fasta(fasta_file)
    
    if (nrow(reads) == 0) {
      return(data.frame(Rank = integer(0), Name = character(0), 
                        Length = integer(0), Preview = character(0)))
    }
    
    # Sort by length (descending) and take top n
    reads <- reads[order(-reads$length), ]
    top_n <- head(reads, n)
    
    data.frame(
      Rank = seq_len(nrow(top_n)),
      Name = top_n$name,
      Length = top_n$length,
      Preview = sapply(top_n$sequence, function(seq) {
        if (nchar(seq) > preview_length) paste0(substr(seq, 1, preview_length), "...")
        else seq
      }),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    return(data.frame(Rank = integer(0), Name = character(0), 
                      Length = integer(0), Preview = character(0)))
  })
}

#' Calculate input reads statistics
calculate_input_stats <- function(fasta_file) {
  tryCatch({
    reads <- parse_fasta(fasta_file)
    
    if (nrow(reads) == 0) {
      return(list(total_reads = 0, total_bases = 0, mean_length = 0, 
                  median_length = 0, longest = 0, shortest = 0))
    }
    
    list(
      total_reads = nrow(reads),
      total_bases = sum(reads$length),
      mean_length = round(mean(reads$length), 0),
      median_length = round(median(reads$length), 0),
      longest = max(reads$length),
      shortest = min(reads$length)
    )
  }, error = function(e) {
    return(list(total_reads = 0, total_bases = 0, mean_length = 0, 
                median_length = 0, longest = 0, shortest = 0))
  })
}
