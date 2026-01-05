# run_hifiasm.R
# hifiasm execution logic and progress monitoring

#' Find hifiasm executable
#' @return Path to hifiasm or NULL if not found
find_hifiasm <- function() {
  possible_paths <- c(
    Sys.which("hifiasm"),
    "/usr/local/bin/hifiasm",
    "/usr/bin/hifiasm",
    "/opt/homebrew/bin/hifiasm",
    file.path(Sys.getenv("HOME"), "bin", "hifiasm"),
    file.path(Sys.getenv("HOME"), "software", "hifiasm", "hifiasm")
  )
  
  for (path in possible_paths) {
    if (nchar(path) > 0 && file.exists(path)) {
      return(path)
    }
  }
  return(NULL)
}

#' Run hifiasm assembly
run_hifiasm <- function(input_file, output_prefix, kmer = 51, threads = 4, 
                        error_rounds = 3, log_file = NULL) {
  hifiasm_path <- find_hifiasm()
  if (is.null(hifiasm_path)) {
    stop("hifiasm not found. Please install hifiasm and ensure it's in PATH.")
  }
  if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
  }
  if (is.null(log_file)) {
    log_file <- paste0(output_prefix, ".log")
  }
  
  args <- c("-o", output_prefix, "-t", as.character(threads),
            "-k", as.character(kmer), "-r", as.character(error_rounds), input_file)
  
  process <- processx::process$new(
    command = hifiasm_path, args = args,
    stdout = log_file, stderr = "2>&1",
    cleanup = TRUE, cleanup_tree = TRUE
  )
  
  return(list(process = process, log_file = log_file, output_prefix = output_prefix))
}

#' Parse hifiasm log file to estimate progress
parse_hifiasm_progress <- function(log_file) {
  if (!file.exists(log_file)) {
    return(list(progress = 0, stage = "Starting...", message = ""))
  }
  
  tryCatch({
    log_content <- readLines(log_file, warn = FALSE)
    if (length(log_content) == 0) {
      return(list(progress = 0, stage = "Starting...", message = ""))
    }
    
    stages <- list(
      "Reading" = list(pattern = "reading", progress = 10),
      "Counting k-mers" = list(pattern = "counting", progress = 20),
      "Building graph" = list(pattern = "graph|assembly graph", progress = 40),
      "Cleaning graph" = list(pattern = "clean", progress = 60),
      "Resolving haplotypes" = list(pattern = "haplotype|phasing", progress = 75),
      "Writing output" = list(pattern = "writing|output", progress = 90),
      "Complete" = list(pattern = "Real time:|CPU time:", progress = 100)
    )
    
    current_progress <- 5
    current_stage <- "Initializing..."
    full_log <- paste(tolower(log_content), collapse = " ")
    
    for (stage_name in names(stages)) {
      stage_info <- stages[[stage_name]]
      if (grepl(stage_info$pattern, full_log, ignore.case = TRUE)) {
        if (stage_info$progress > current_progress) {
          current_progress <- stage_info$progress
          current_stage <- stage_name
        }
      }
    }
    
    last_lines <- tail(log_content, 3)
    message <- paste(last_lines, collapse = "\n")
    
    return(list(progress = current_progress, stage = current_stage, message = message))
  }, error = function(e) {
    return(list(progress = 0, stage = "Error reading log", message = as.character(e)))
  })
}

#' Check if hifiasm process is still running
is_hifiasm_running <- function(process_info) {
  if (is.null(process_info$process)) return(FALSE)
  return(process_info$process$is_alive())
}

#' Get hifiasm output files
get_hifiasm_outputs <- function(output_prefix) {
  possible_outputs <- list(
    primary = paste0(output_prefix, ".bp.p_ctg.gfa"),
    alternate = paste0(output_prefix, ".bp.a_ctg.gfa"),
    hap1 = paste0(output_prefix, ".bp.hap1.p_ctg.gfa"),
    hap2 = paste0(output_prefix, ".bp.hap2.p_ctg.gfa"),
    raw_unitig = paste0(output_prefix, ".bp.r_utg.gfa")
  )
  
  existing_outputs <- list()
  for (name in names(possible_outputs)) {
    if (file.exists(possible_outputs[[name]])) {
      existing_outputs[[name]] <- possible_outputs[[name]]
    }
  }
  return(existing_outputs)
}

#' Clean up temporary files
cleanup_hifiasm_files <- function(output_prefix, keep_gfa = TRUE) {
  dir_path <- dirname(output_prefix)
  base_name <- basename(output_prefix)
  all_files <- list.files(dir_path, pattern = paste0("^", base_name), full.names = TRUE)
  
  if (keep_gfa) {
    files_to_remove <- all_files[!grepl("\\.gfa$", all_files)]
  } else {
    files_to_remove <- all_files
  }
  
  for (f in files_to_remove) {
    tryCatch(file.remove(f), error = function(e) NULL)
  }
}
