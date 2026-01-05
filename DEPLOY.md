# Deployment Guide for FGCZ Shiny Server

## Overview

This document describes how to deploy the MINT Assembly Demo to the FGCZ Shiny proxy server (https://fgcz-shiny.uzh.ch/).

## Prerequisites

### Server Requirements

1. **R version**: 4.0 or higher
2. **Required R packages** (install via `install_packages.R`):
   - shiny, shinydashboard, shinyWidgets
   - DT, ggplot2, plotly
   - processx, dplyr, tidyr, scales

3. **Optional**: hifiasm binary (for Run mode)

### Files to Deploy

```
MINT_assembly_2026/
├── app.R                    # Main application (required)
├── R/                       # Helper functions (required)
│   ├── run_hifiasm.R
│   ├── parse_results.R
│   └── plot_functions.R
├── www/                     # Static assets (required)
│   └── custom.css
├── MINT_Tage/              # Sample data (required for demo)
│   ├── *.fasta
│   ├── kokai_100K/
│   ├── kokai_1M_11x/
│   └── kokai_1M_24x/
├── install_packages.R       # Package installer (optional)
└── README.md               # Documentation (optional)
```

## Deployment Steps

### 1. Prepare the Package

```bash
# Create deployment archive
cd /path/to/MINT_assembly_2026
tar -czvf mint_assembly_demo.tar.gz \
  app.R R/ www/ MINT_Tage/ README.md install_packages.R
```

### 2. Upload to FGCZ Server

Contact the FGCZ bioinformatics team or use your allocated upload method:
- Email: bioinformatics@fgcz.uzh.ch
- SFTP to designated server path

### 3. Install Dependencies

On the server, run:
```r
source("install_packages.R")
```

### 4. Configure App Settings

If needed, create a `shiny-server.conf` entry or app configuration:

```
# Example Shiny Server configuration
location /mint-assembly {
  app_dir /srv/shiny-server/mint-assembly;
  log_dir /var/log/shiny-server;
  
  # Increase timeout for assembly operations
  app_idle_timeout 600;
  
  # Resource limits
  reconnect true;
  sanitize_errors true;
}
```

## Configuration Options

### Timeout Settings

For long-running assemblies, configure appropriate timeouts:

```r
# In app.R or global.R
options(shiny.maxRequestSize = 100*1024^2)  # 100MB max upload
```

### hifiasm Path

If hifiasm is installed in a non-standard location, update `R/run_hifiasm.R`:

```r
# Add custom path to find_hifiasm()
possible_paths <- c(
  "/path/to/custom/hifiasm",
  Sys.which("hifiasm"),
  ...
)
```

## Testing the Deployment

### 1. Verify App Loads

Navigate to: https://fgcz-shiny.uzh.ch/your-app-path

### 2. Test Demo Mode

1. Select "Demo (Pre-computed)" mode
2. Choose a dataset
3. Click "Load Results"
4. Verify statistics and plots display correctly

### 3. Test Download Functions

1. Load any dataset
2. Try downloading FASTA, GFA, and CSV files

## Troubleshooting

### Common Issues

1. **Package not found**
   - Run `install_packages.R` again
   - Check R library path

2. **GFA files not found**
   - Verify MINT_Tage directory is uploaded
   - Check file permissions

3. **hifiasm not found** (Run mode)
   - Normal for servers without hifiasm
   - Demo mode will still work

4. **Timeout errors**
   - Increase server timeout settings
   - Use smaller datasets for testing

### Logs

Check Shiny server logs:
```bash
tail -f /var/log/shiny-server/*.log
```

## Contact

For FGCZ-specific deployment issues:
- Email: bioinformatics@fgcz.uzh.ch
- Web: https://fgcz.ch/

## Version History

- v1.0.0: Initial release with demo mode and run mode support

