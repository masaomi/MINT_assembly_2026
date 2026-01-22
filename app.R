# app.R - MINT Assembly Demo Shiny Application
# hifiasm-based genome assembly demonstration tool

# Load required packages
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(plotly)
library(processx)
library(dplyr)
library(tidyr)
library(scales)

# Source helper functions
source("R/run_hifiasm.R")
source("R/parse_results.R")
source("R/plot_functions.R")

# Configuration
DATA_DIR <- "MINT_Tage"

# Get available sample files
get_sample_files <- function() {
  patterns <- c("\\.fasta$", "\\.fa$", "\\.fna$")
  files <- character(0)
  
  # Check MINT_Tage directory
  if (dir.exists(DATA_DIR)) {
    for (pattern in patterns) {
      found <- list.files(DATA_DIR, pattern = pattern, full.names = FALSE, recursive = FALSE)
      if (length(found) > 0) {
        files <- c(files, file.path(DATA_DIR, found))
      }
    }
  }
  
  # Also check data directory
  if (dir.exists("data")) {
    for (pattern in patterns) {
      found <- list.files("data", pattern = pattern, full.names = FALSE, recursive = TRUE)
      if (length(found) > 0) {
        files <- c(files, file.path("data", found))
      }
    }
  }
  
  if (length(files) == 0) {
    files <- c("No FASTA files found")
  }
  
  return(files)
}

# Get available pre-computed GFA files for demo mode
get_demo_gfa_files <- function() {
  gfa_files <- list()
  
  # Check known directories with pre-computed results
  demo_dirs <- c(
    file.path(DATA_DIR, "kokai_100K"),
    file.path(DATA_DIR, "kokai_1M_11x"),
    file.path(DATA_DIR, "kokai_1M_24x")
  )
  
  for (dir in demo_dirs) {
    if (dir.exists(dir)) {
      primary_gfa <- file.path(dir, "kokai.asm.bp.p_ctg.gfa")
      if (file.exists(primary_gfa)) {
        # Extract friendly name from directory
        name <- basename(dir)
        gfa_files[[name]] <- primary_gfa
      }
    }
  }
  
  return(gfa_files)
}

# Get corresponding input FASTA file for demo datasets
get_demo_fasta_file <- function(demo_name) {
  fasta_mapping <- list(
    "kokai_100K" = file.path(DATA_DIR, "Ckokai_chr1_a_and_b_100K.fasta"),
    "kokai_1M_11x" = file.path(DATA_DIR, "kokai_1M_11x", "Ckokai_chr1_a_and_b_1M_0.10_percent.fasta"),
    "kokai_1M_24x" = file.path(DATA_DIR, "Ckokai_chr1_a_and_b_1M.fasta")
  )
  
  fasta_file <- fasta_mapping[[demo_name]]
  if (!is.null(fasta_file) && file.exists(fasta_file)) {
    return(fasta_file)
  }
  return(NULL)
}

# Check if hifiasm is available
hifiasm_available <- function() {
  !is.null(find_hifiasm())
}

# UI Definition
ui <- dashboardPage(
  skin = "blue",
  
  # Header
  dashboardHeader(
    title = span(icon("dna"), "MINT Assembly Demo"),
    titleWidth = 280
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 280,
    
    sidebarMenu(
      id = "tabs",
      menuItem("Assembly", tabName = "assembly", icon = icon("play")),
      menuItem("Results", tabName = "results", icon = icon("chart-bar")),
      menuItem("Comparison", tabName = "comparison", icon = icon("balance-scale")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    ),
    
    hr(),
    
    # Parameter controls
    div(style = "padding: 10px;",
        h4(icon("sliders-h"), "Parameters", style = "color: #fff;"),
        
        # Mode selection
        radioGroupButtons(
          "mode",
          label = "Mode:",
          choices = c("Demo" = "demo", "Run" = "run"),
          selected = "demo",
          status = "primary",
          size = "normal",
          justified = TRUE
        ),
        
        # Demo mode: select pre-computed results
        conditionalPanel(
          condition = "input.mode == 'demo'",
          selectInput(
            "demo_gfa",
            label = "Select Dataset:",
            choices = names(get_demo_gfa_files()),
            selected = NULL,
            width = "100%"
          ),
          helpText("Load pre-computed assembly results")
        ),
        
        # Run mode: assembly parameters
        conditionalPanel(
          condition = "input.mode == 'run'",
          
          # Input file selection
          selectInput(
            "input_file",
            label = "Input FASTA File:",
            choices = get_sample_files(),
            selected = NULL,
            width = "100%"
          ),
          
          # K-mer size
          sliderInput(
            "kmer",
            label = "K-mer Size:",
            min = 15,
            max = 63,
            value = 51,
            step = 2
          ),
          helpText("Larger k-mer = more specific but may miss overlaps"),
          
          # Error correction rounds
          sliderInput(
            "error_rounds",
            label = "Error Correction Rounds:",
            min = 0,
            max = 10,
            value = 3,
            step = 1
          ),
          helpText("More rounds = better accuracy but slower"),
          
          # Thread count (fixed)
          div(
            style = "pointer-events: none; opacity: 0.6;",
            numericInput(
              "threads",
              label = "CPU Threads:",
              value = 2,
              min = 2,
              max = 2,
              step = 1,
              width = "100%"
            )
          )
        ),
        
        # Minimum contig length (applies to both modes)
        numericInput(
          "min_length",
          label = "Min Contig Length (bp):",
          value = 500,
          min = 0,
          max = 10000,
          step = 100
        ),
        
        hr(),
        
        # Action button (changes based on mode)
        conditionalPanel(
          condition = "input.mode == 'demo'",
          actionButton(
            "load_demo",
            label = "Load Results",
            icon = icon("folder-open"),
            class = "btn-success action-button"
          )
        ),
        
        conditionalPanel(
          condition = "input.mode == 'run'",
          actionButton(
            "run_assembly",
            label = "Run Assembly",
            icon = icon("rocket"),
            class = "btn-primary action-button"
          )
        ),
        
        br(), br(),
        
        # Status indicator
        div(
          id = "status_container",
          style = "text-align: center;",
          uiOutput("status_indicator")
        ),
        
        # hifiasm availability warning
        conditionalPanel(
          condition = "input.mode == 'run'",
          uiOutput("hifiasm_status")
        )
    )
  ),
  
  # Body
  dashboardBody(
    # Include custom CSS
    tags$head(
      # Viewport meta tag for mobile responsiveness
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0, maximum-scale=5.0, user-scalable=yes"),
      tags$meta(name = "mobile-web-app-capable", content = "yes"),
      tags$meta(name = "apple-mobile-web-app-capable", content = "yes"),
      tags$meta(name = "apple-mobile-web-app-status-bar-style", content = "black-translucent"),
      
      # CSS
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$style(HTML("
        .content-wrapper { background-color: #0a1628; }
        .box { background-color: rgba(17, 34, 64, 0.8); color: #e6f1ff; }
      ")),
      
      # JavaScript for mobile detection and sidebar auto-close
      tags$script(HTML("
        $(document).ready(function() {
          // Auto-close sidebar on mobile after menu click
          if (window.innerWidth <= 768) {
            $('.sidebar-menu a').on('click', function() {
              setTimeout(function() {
                $('body').removeClass('sidebar-open');
                $('body').addClass('sidebar-collapse');
              }, 100);
            });
          }
          
          // Handle orientation change
          window.addEventListener('orientationchange', function() {
            setTimeout(function() {
              window.dispatchEvent(new Event('resize'));
            }, 100);
          });
          
          // Close sidebar when clicking outside on mobile
          $(document).on('click', '.content-wrapper', function() {
            if (window.innerWidth <= 768 && $('body').hasClass('sidebar-open')) {
              $('body').removeClass('sidebar-open');
              $('body').addClass('sidebar-collapse');
            }
          });
        });
      "))
    ),
    
    tabItems(
      # Assembly Tab
      tabItem(
        tabName = "assembly",
        
        fluidRow(
          box(
            title = span(icon("terminal"), "Assembly Progress"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            # Progress bar
            div(
              style = "margin-bottom: 20px;",
              progressBar(
                id = "assembly_progress",
                value = 0,
                total = 100,
                status = "primary",
                display_pct = TRUE,
                striped = TRUE
              )
            ),
            
            # Current stage
            div(
              style = "text-align: center; margin-bottom: 15px;",
              h4(textOutput("current_stage"))
            ),
            
            # Log output
            div(
              class = "log-output",
              style = "height: 250px; overflow-y: auto;",
              verbatimTextOutput("log_output")
            )
          )
        ),
        
        fluidRow(
          # Quick stats boxes
          valueBoxOutput("stat_contigs", width = 3),
          valueBoxOutput("stat_total_length", width = 3),
          valueBoxOutput("stat_n50", width = 3),
          valueBoxOutput("stat_largest", width = 3)
        )
      ),
      
      # Results Tab
      tabItem(
        tabName = "results",
        
        fluidRow(
          # Comparison Statistics Table
          box(
            title = span(icon("scale-balanced"), "Input / Assembly Comparison"),
            status = "primary",
            solidHeader = TRUE,
            width = 5,
            
            tableOutput("comparison_stats_table")
          ),
          
          # Plots
          box(
            title = span(icon("chart-area"), "Visualizations"),
            status = "primary",
            solidHeader = TRUE,
            width = 7,
            
            tabsetPanel(
              tabPanel(
                "Length Dist. (Input)",
                br(),
                plotlyOutput("plot_histogram_input", height = "350px")
              ),
              tabPanel(
                "Length Dist. (Assembly)",
                br(),
                plotlyOutput("plot_histogram_assembly", height = "350px")
              )
            )
          )
        ),
        
        fluidRow(
          # Input Reads Preview
          box(
            title = span(icon("file-lines"), "Input Reads Preview"),
            status = "info",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            
            fluidRow(
              column(3, valueBoxOutput("input_reads_count", width = 12)),
              column(3, valueBoxOutput("input_total_bases", width = 12)),
              column(3, valueBoxOutput("input_mean_length", width = 12)),
              column(3, valueBoxOutput("input_longest", width = 12))
            ),
            DTOutput("reads_preview_table")
          )
        ),
        
        fluidRow(
          # Assembled Contigs Preview
          box(
            title = span(icon("dna"), "Assembled Contigs Preview"),
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            
            DTOutput("contig_preview_table")
          )
        ),
        
        fluidRow(
          # Download section
          box(
            title = span(icon("download"), "Download Results"),
            status = "success",
            solidHeader = TRUE,
            width = 12,
            
            div(
              class = "download-section",
              downloadButton("download_fasta", "Download FASTA", class = "btn-success"),
              downloadButton("download_stats", "Download Statistics (CSV)", class = "btn-info"),
              downloadButton("download_gfa", "Download GFA", class = "btn-warning")
            )
          )
        )
      ),
      
      # Comparison Tab
      tabItem(
        tabName = "comparison",
        
        fluidRow(
          box(
            title = span(icon("balance-scale"), "Assembly Run Comparison"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            div(
              style = "margin-bottom: 20px;",
              actionButton("clear_history", "Clear History", 
                           icon = icon("trash"), class = "btn-danger btn-sm")
            ),
            
            DTOutput("comparison_table")
          )
        ),
        
        fluidRow(
          box(
            title = span(icon("chart-bar"), "Comparison Charts"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            tabsetPanel(
              tabPanel(
                "N50 Comparison",
                br(),
                plotlyOutput("comparison_n50", height = "300px")
              ),
              tabPanel(
                "Total Length",
                br(),
                plotlyOutput("comparison_length", height = "300px")
              ),
              tabPanel(
                "Contig Count",
                br(),
                plotlyOutput("comparison_contigs", height = "300px")
              )
            )
          )
        )
      ),
      
      # Help Tab
      tabItem(
        tabName = "help",
        
        fluidRow(
          box(
            title = span(icon("info-circle"), "About This Tool"),
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            h4("MINT Assembly Demo"),
            p("This tool demonstrates genome assembly using hifiasm, a fast haplotype-resolved 
               de novo assembler for PacBio HiFi reads."),
            
            hr(),
            
            h4("Modes"),
            tags$ul(
              tags$li(tags$strong("Demo Mode:"), " Load pre-computed assembly results from the 
                      MINT_Tage sample data. No hifiasm installation required."),
              tags$li(tags$strong("Run Mode:"), " Execute hifiasm assembly with custom parameters. 
                      Requires hifiasm to be installed on the server.")
            ),
            
            hr(),
            
            h4("Parameters"),
            tags$ul(
              tags$li(tags$strong("K-mer Size:"), " The length of k-mers used for overlap detection. 
                      Larger values increase specificity but may miss some overlaps. Default: 51"),
              tags$li(tags$strong("Error Correction Rounds:"), " Number of rounds for correcting 
                      sequencing errors. More rounds improve accuracy but increase runtime. Default: 3"),
              tags$li(tags$strong("Min Contig Length:"), " Minimum length threshold for filtering 
                      output contigs in statistics. Default: 500 bp"),
              tags$li(tags$strong("CPU Threads:"), " Fixed to 1 for demo consistency.")
            ),
            
            hr(),
            
            h4("Assembly Statistics"),
            tags$ul(
              tags$li(tags$strong("N50:"), " The sequence length of the shortest contig at 50% 
                      of the total genome length. Higher is better."),
              tags$li(tags$strong("L50:"), " The smallest number of contigs whose length sum 
                      produces N50. Lower is better."),
              tags$li(tags$strong("GC Content:"), " Percentage of guanine-cytosine base pairs 
                      in the assembly.")
            ),
            
            hr(),
            
            h4("Sample Datasets"),
            tags$ul(
              tags$li(tags$strong("kokai_100K:"), " 100K coverage assembly"),
              tags$li(tags$strong("kokai_1M_11x:"), " 1M dataset with 11x coverage"),
              tags$li(tags$strong("kokai_1M_24x:"), " 1M dataset with 24x coverage")
            ),
            
            hr(),
            
            h4("References"),
            p("Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., & Li, H. (2021). 
               Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. 
               Nature Methods, 18(2), 170-175."),
            
            tags$a(href = "https://github.com/chhylp123/hifiasm", 
                   target = "_blank", 
                   icon("github"), " hifiasm GitHub Repository")
          )
        )
      )
    )
  )
)

# Server Definition
server <- function(input, output, session) {
  
  # Store GFA file paths
  demo_gfa_files <- get_demo_gfa_files()
  
  # Reactive values for state management
  rv <- reactiveValues(
    assembly_running = FALSE,
    assembly_complete = FALSE,
    process_info = NULL,
    current_gfa = NULL,
    current_stats = NULL,
    current_source = NULL,  # Track data source (demo name or run params)
    input_fasta = NULL,     # Path to input FASTA file
    input_stats = NULL,     # Input reads statistics
    run_history = data.frame(
      run_id = character(0),
      timestamp = character(0),
      source = character(0),
      kmer = integer(0),
      error_rounds = integer(0),
      min_length = integer(0),
      n50 = integer(0),
      l50 = integer(0),
      total_contigs = integer(0),
      filtered_contigs = integer(0),
      filtered_length = integer(0),
      largest_contig = integer(0),
      stringsAsFactors = FALSE
    ),
    log_content = ""
  )
  
  # Timer for progress updates
  autoInvalidate <- reactiveTimer(1000)
  
  # hifiasm availability status
  output$hifiasm_status <- renderUI({
    if (!hifiasm_available()) {
      div(
        class = "alert alert-warning",
        style = "margin-top: 10px; padding: 10px; font-size: 12px;",
        icon("exclamation-triangle"),
        " hifiasm not found. Use Demo mode instead."
      )
    }
  })
  
  # Status indicator
  output$status_indicator <- renderUI({
    if (rv$assembly_running) {
      div(
        class = "status-running",
        icon("spinner", class = "fa-spin"),
        " Assembly in progress..."
      )
    } else if (rv$assembly_complete) {
      div(
        class = "status-ready",
        icon("check-circle"),
        " Data loaded!"
      )
    } else {
      div(
        class = "status-ready",
        icon("circle"),
        " Ready"
      )
    }
  })
  
  # Load demo data
  observeEvent(input$load_demo, {
    req(input$demo_gfa)
    
    gfa_path <- demo_gfa_files[[input$demo_gfa]]
    
    if (is.null(gfa_path) || !file.exists(gfa_path)) {
      showNotification("GFA file not found", type = "error")
      return()
    }
    
    # Update state
    rv$current_gfa <- gfa_path
    rv$current_source <- input$demo_gfa
    rv$log_content <- paste0("Loaded pre-computed results from: ", gfa_path, "\n",
                             "Dataset: ", input$demo_gfa, "\n",
                             "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    
    # Get corresponding input FASTA file
    fasta_path <- get_demo_fasta_file(input$demo_gfa)
    rv$input_fasta <- fasta_path
    if (!is.null(fasta_path) && file.exists(fasta_path)) {
      rv$input_stats <- calculate_input_stats(fasta_path)
    } else {
      rv$input_stats <- NULL
    }
    
    # Calculate statistics
    tryCatch({
      rv$current_stats <- calculate_assembly_stats(rv$current_gfa, input$min_length)
      rv$assembly_complete <- TRUE
      updateProgressBar(session, "assembly_progress", value = 100)
      
      # Add to history
      new_row <- data.frame(
        run_id = paste0("demo_", format(Sys.time(), "%H%M%S")),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        source = input$demo_gfa,
        kmer = NA,
        error_rounds = NA,
        min_length = input$min_length,
        n50 = rv$current_stats$n50,
        l50 = rv$current_stats$l50,
        total_contigs = rv$current_stats$total_contigs,
        filtered_contigs = rv$current_stats$filtered_contigs,
        filtered_length = rv$current_stats$filtered_length,
        largest_contig = rv$current_stats$largest_contig,
        stringsAsFactors = FALSE
      )
      rv$run_history <- rbind(rv$run_history, new_row)
      
      showNotification(paste("Loaded:", input$demo_gfa), type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error")
    })
  })
  
  # Run assembly when button clicked
  observeEvent(input$run_assembly, {
    req(input$input_file)
    
    # Validate input file
    if (grepl("No FASTA files found", input$input_file)) {
      showNotification("Please add FASTA files to the data/ directory", type = "error")
      return()
    }
    
    # Check if file exists
    if (!file.exists(input$input_file)) {
      showNotification(paste("File not found:", input$input_file), type = "error")
      return()
    }
    
    # Check if hifiasm is available
    if (!hifiasm_available()) {
      showNotification("hifiasm not found. Please use Demo mode instead.", type = "error")
      return()
    }
    
    # Create session-specific output directory for user isolation
    session_id <- session$token
    output_dir <- file.path(tempdir(), paste0("session_", session_id))
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_prefix <- file.path(output_dir, paste0("assembly_", run_id))
    
    # Reset state
    rv$assembly_running <- TRUE
    rv$assembly_complete <- FALSE
    rv$log_content <- ""
    rv$current_source <- basename(input$input_file)
    
    # Store input FASTA info
    rv$input_fasta <- input$input_file
    rv$input_stats <- calculate_input_stats(input$input_file)
    
    # Update progress bar
    updateProgressBar(session, "assembly_progress", value = 0)
    
    # Run hifiasm
    tryCatch({
      rv$process_info <- run_hifiasm(
        input_file = input$input_file,
        output_prefix = output_prefix,
        kmer = input$kmer,
        threads = 2,
        error_rounds = input$error_rounds
      )
      
      showNotification("Assembly started!", type = "message")
      
    }, error = function(e) {
      rv$assembly_running <- FALSE
      showNotification(paste("Error starting assembly:", e$message), type = "error")
    })
  })
  
  # Progress monitoring
  observe({
    autoInvalidate()
    
    if (rv$assembly_running && !is.null(rv$process_info)) {
      # Check if process is still running
      if (is_hifiasm_running(rv$process_info)) {
        # Update progress
        progress <- parse_hifiasm_progress(rv$process_info$log_file)
        updateProgressBar(session, "assembly_progress", value = progress$progress)
        
        # Update log
        if (file.exists(rv$process_info$log_file)) {
          rv$log_content <- paste(readLines(rv$process_info$log_file, warn = FALSE), 
                                  collapse = "\n")
        }
      } else {
        # Assembly finished
        rv$assembly_running <- FALSE
        rv$assembly_complete <- TRUE
        updateProgressBar(session, "assembly_progress", value = 100)
        
        # Get output files
        outputs <- get_hifiasm_outputs(rv$process_info$output_prefix)
        
        if (length(outputs) > 0 && !is.null(outputs$primary)) {
          rv$current_gfa <- outputs$primary
          
          # Calculate statistics
          rv$current_stats <- calculate_assembly_stats(rv$current_gfa, input$min_length)
          
          # Add to history
          new_row <- data.frame(
            run_id = basename(rv$process_info$output_prefix),
            timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            source = basename(input$input_file),
            kmer = input$kmer,
            error_rounds = input$error_rounds,
            min_length = input$min_length,
            n50 = rv$current_stats$n50,
            l50 = rv$current_stats$l50,
            total_contigs = rv$current_stats$total_contigs,
            filtered_contigs = rv$current_stats$filtered_contigs,
            filtered_length = rv$current_stats$filtered_length,
            largest_contig = rv$current_stats$largest_contig,
            stringsAsFactors = FALSE
          )
          rv$run_history <- rbind(rv$run_history, new_row)
          
          showNotification("Assembly completed successfully!", type = "message")
        } else {
          showNotification("Assembly completed but no output files found", type = "warning")
        }
      }
    }
  })
  
  # Log output
  output$log_output <- renderText({
    rv$log_content
  })
  
  # Current stage
  output$current_stage <- renderText({
    if (rv$assembly_running && !is.null(rv$process_info)) {
      progress <- parse_hifiasm_progress(rv$process_info$log_file)
      progress$stage
    } else if (rv$assembly_complete) {
      paste("Complete -", rv$current_source)
    } else {
      "Ready to start"
    }
  })
  
  # Value boxes
  output$stat_contigs <- renderValueBox({
    val <- if (!is.null(rv$current_stats)) rv$current_stats$filtered_contigs else 0
    valueBox(
      value = format(val, big.mark = ","),
      subtitle = "Contigs",
      icon = icon("layer-group"),
      color = "blue"
    )
  })
  
  output$stat_total_length <- renderValueBox({
    val <- if (!is.null(rv$current_stats)) rv$current_stats$filtered_length else 0
    valueBox(
      value = format(val, big.mark = ","),
      subtitle = "Total Length (bp)",
      icon = icon("ruler"),
      color = "green"
    )
  })
  
  output$stat_n50 <- renderValueBox({
    val <- if (!is.null(rv$current_stats)) rv$current_stats$n50 else 0
    valueBox(
      value = format(val, big.mark = ","),
      subtitle = "N50 (bp)",
      icon = icon("chart-line"),
      color = "purple"
    )
  })
  
  output$stat_largest <- renderValueBox({
    val <- if (!is.null(rv$current_stats)) rv$current_stats$largest_contig else 0
    valueBox(
      value = format(val, big.mark = ","),
      subtitle = "Largest Contig (bp)",
      icon = icon("trophy"),
      color = "yellow"
    )
  })
  
  # Comparison statistics table (Input vs Assembly)
  output$comparison_stats_table <- renderTable({
    req(rv$current_stats)
    
    # Get input stats (may be NULL)
    input_stats <- rv$input_stats
    assembly_stats <- rv$current_stats
    
    # Build comparison table
    comparison_df <- data.frame(
      Metric = c(
        "Sequences",
        "Total Length (bp)",
        "Mean Length (bp)",
        "Median Length (bp)",
        "Longest (bp)",
        "Shortest (bp)",
        "N50 (bp)",
        "L50",
        "GC Content (%)"
      ),
      Input = c(
        if (!is.null(input_stats)) format(input_stats$total_reads, big.mark = ",") else "-",
        if (!is.null(input_stats)) format(input_stats$total_bases, big.mark = ",") else "-",
        if (!is.null(input_stats)) format(input_stats$mean_length, big.mark = ",") else "-",
        if (!is.null(input_stats)) format(input_stats$median_length, big.mark = ",") else "-",
        if (!is.null(input_stats)) format(input_stats$longest, big.mark = ",") else "-",
        if (!is.null(input_stats)) format(input_stats$shortest, big.mark = ",") else "-",
        "-",  # N50 not applicable for input reads
        "-",  # L50 not applicable for input reads
        "-"   # GC not calculated for input
      ),
      Assembly = c(
        format(assembly_stats$filtered_contigs, big.mark = ","),
        format(assembly_stats$filtered_length, big.mark = ","),
        format(assembly_stats$mean_length, big.mark = ","),
        format(assembly_stats$median_length, big.mark = ","),
        format(assembly_stats$largest_contig, big.mark = ","),
        format(assembly_stats$smallest_contig, big.mark = ","),
        format(assembly_stats$n50, big.mark = ","),
        assembly_stats$l50,
        ifelse(is.na(assembly_stats$gc_content), "-", paste0(assembly_stats$gc_content, "%"))
      ),
      stringsAsFactors = FALSE
    )
    
    comparison_df
  }, striped = TRUE, hover = TRUE, bordered = TRUE, align = "lrr")
  
  # Plots
  output$plot_histogram_input <- renderPlotly({
    req(rv$input_fasta)
    
    input_lengths <- c()
    if (file.exists(rv$input_fasta)) {
      tryCatch({
        input_reads <- parse_fasta(rv$input_fasta)
        input_lengths <- input_reads$length
      }, error = function(e) {
        input_lengths <- c()
      })
    }
    
    p <- plot_length_histogram(input_lengths, title = "Input Reads Length Distribution")
    ggplotly(p)
  })
  
  output$plot_histogram_assembly <- renderPlotly({
    req(rv$current_gfa)
    lengths <- get_contig_lengths(rv$current_gfa, input$min_length)
    p <- plot_length_histogram(lengths, title = "Assembled Contigs Length Distribution")
    ggplotly(p)
  })
  
  # Input reads value boxes
  output$input_reads_count <- renderValueBox({
    count <- if (!is.null(rv$input_stats)) rv$input_stats$total_reads else 0
    valueBox(
      format(count, big.mark = ","),
      "Reads",
      icon = icon("list"),
      color = "light-blue"
    )
  })
  
  output$input_total_bases <- renderValueBox({
    bases <- if (!is.null(rv$input_stats)) rv$input_stats$total_bases else 0
    valueBox(
      format(bases, big.mark = ","),
      "Total Bases",
      icon = icon("ruler"),
      color = "light-blue"
    )
  })
  
  output$input_mean_length <- renderValueBox({
    mean_len <- if (!is.null(rv$input_stats)) rv$input_stats$mean_length else 0
    valueBox(
      format(mean_len, big.mark = ","),
      "Mean Length",
      icon = icon("chart-line"),
      color = "light-blue"
    )
  })
  
  output$input_longest <- renderValueBox({
    longest <- if (!is.null(rv$input_stats)) rv$input_stats$longest else 0
    valueBox(
      format(longest, big.mark = ","),
      "Longest Read",
      icon = icon("arrow-up"),
      color = "light-blue"
    )
  })
  
  # Input reads preview table
  output$reads_preview_table <- renderDT({
    req(rv$input_fasta)
    
    if (!file.exists(rv$input_fasta)) {
      return(NULL)
    }
    
    preview <- get_top_reads(rv$input_fasta, n = 10, preview_length = 100)
    
    if (nrow(preview) == 0) {
      return(NULL)
    }
    
    # Truncate long read names for display
    preview$Name <- sapply(preview$Name, function(name) {
      if (nchar(name) > 20) paste0(substr(name, 1, 17), "...") else name
    })
    
    datatable(
      preview,
      options = list(
        pageLength = 5,
        dom = 'tip',
        autoWidth = FALSE,
        scrollX = FALSE,
        columnDefs = list(
          list(width = '40px', targets = 0),   # Rank - narrow
          list(width = '100px', targets = 1),  # Name - slightly wider for truncated names
          list(width = '70px', targets = 2),   # Length - narrow
          list(width = '60%', targets = 3, className = 'sequence-preview')  # Preview - wide
        )
      ),
      rownames = FALSE
    ) %>%
      formatStyle('Preview', fontFamily = 'JetBrains Mono, monospace', fontSize = '11px')
  })
  
  # Contig preview table
  output$contig_preview_table <- renderDT({
    req(rv$current_gfa)
    preview <- get_top_contigs(rv$current_gfa, n = 10, preview_length = 100)
    
    datatable(
      preview,
      options = list(
        pageLength = 5,
        dom = 'tip',
        autoWidth = FALSE,
        columnDefs = list(
          list(width = '40px', targets = 0),   # Rank - narrow
          list(width = '80px', targets = 1),   # Name - narrow
          list(width = '70px', targets = 2),   # Length - narrow
          list(width = '60%', targets = 3, className = 'sequence-preview')  # Preview - wide
        )
      ),
      rownames = FALSE
    ) %>%
      formatStyle('Preview', fontFamily = 'JetBrains Mono, monospace', fontSize = '11px')
  })
  
  # Comparison table
  output$comparison_table <- renderDT({
    req(nrow(rv$run_history) > 0)
    
    display_df <- rv$run_history %>%
      select(run_id, timestamp, source, kmer, error_rounds, 
             n50, l50, filtered_contigs, filtered_length, largest_contig) %>%
      rename(
        "Run ID" = run_id,
        "Time" = timestamp,
        "Source" = source,
        "K-mer" = kmer,
        "Error Rounds" = error_rounds,
        "N50" = n50,
        "L50" = l50,
        "Contigs" = filtered_contigs,
        "Length" = filtered_length,
        "Largest" = largest_contig
      )
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        dom = 'tip',
        order = list(list(1, 'desc'))
      ),
      rownames = FALSE
    ) %>%
      formatStyle('N50', backgroundColor = styleInterval(
        quantile(rv$run_history$n50, probs = c(0.33, 0.66), na.rm = TRUE),
        c('#553333', '#555533', '#335533')
      ))
  })
  
  # Comparison plots
  output$comparison_n50 <- renderPlotly({
    req(nrow(rv$run_history) > 0)
    p <- plot_comparison(rv$run_history, "n50", "N50 Comparison")
    ggplotly(p)
  })
  
  output$comparison_length <- renderPlotly({
    req(nrow(rv$run_history) > 0)
    p <- plot_comparison(rv$run_history, "filtered_length", "Total Length Comparison")
    ggplotly(p)
  })
  
  output$comparison_contigs <- renderPlotly({
    req(nrow(rv$run_history) > 0)
    p <- plot_comparison(rv$run_history, "filtered_contigs", "Contig Count Comparison")
    ggplotly(p)
  })
  
  # Clear history
  observeEvent(input$clear_history, {
    rv$run_history <- data.frame(
      run_id = character(0),
      timestamp = character(0),
      source = character(0),
      kmer = integer(0),
      error_rounds = integer(0),
      min_length = integer(0),
      n50 = integer(0),
      l50 = integer(0),
      total_contigs = integer(0),
      filtered_contigs = integer(0),
      filtered_length = integer(0),
      largest_contig = integer(0),
      stringsAsFactors = FALSE
    )
    showNotification("History cleared", type = "message")
  })
  
  # Download handlers
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0("assembly_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta")
    },
    content = function(file) {
      req(rv$current_gfa)
      gfa_to_fasta(rv$current_gfa, file, input$min_length)
    }
  )
  
  output$download_stats <- downloadHandler(
    filename = function() {
      paste0("assembly_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(rv$current_stats)
      stats_df <- format_stats_table(rv$current_stats)
      write.csv(stats_df, file, row.names = FALSE)
    }
  )
  
  output$download_gfa <- downloadHandler(
    filename = function() {
      paste0("assembly_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".gfa")
    },
    content = function(file) {
      req(rv$current_gfa)
      file.copy(rv$current_gfa, file)
    }
  )
  
  # Cleanup on session end
  session$onSessionEnded(function() {
    # Kill running process
    if (!is.null(rv$process_info) && !is.null(rv$process_info$process)) {
      tryCatch({
        rv$process_info$process$kill()
      }, error = function(e) NULL)
    }
    
    # Remove session-specific temp directory
    session_dir <- file.path(tempdir(), paste0("session_", session$token))
    if (dir.exists(session_dir)) {
      tryCatch({
        unlink(session_dir, recursive = TRUE)
      }, error = function(e) NULL)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
