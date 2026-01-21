# MINT Assembly Demo

A Shiny application for demonstrating genome assembly using hifiasm with interactive visualization and parameter comparison.

## Features

- **Demo Mode**: Load pre-computed assembly results instantly
- **Run Mode**: Execute hifiasm assembly with custom parameters
- **Interactive Visualization**: 
  - Contig length distribution histograms
  - Cumulative length (Nx) plots
  - Length category breakdown
- **Assembly Statistics**: N50, L50, GC content, and more
- **Parameter Comparison**: Compare results across multiple runs
- **Download Results**: Export FASTA, GFA, and statistics CSV

## Requirements

### R Packages

```r
# CRAN packages
install.packages(c(
  "shiny",
  "shinydashboard", 
  "shinyWidgets",
  "DT",
  "ggplot2",
  "plotly",
  "processx",
  "dplyr",
  "tidyr",
  "scales"
))
```

### Optional (for Run Mode)

- **hifiasm**: Required only for running new assemblies
  - Installation: https://github.com/chhylp123/hifiasm

## Quick Start

### Local Development

```bash
# Clone the repository
cd MINT_assembly_2026

# Start the Shiny app
R -e "shiny::runApp('app.R', port = 3838)"
```

Then open http://localhost:3838 in your browser.

### Remote Access (Server Deployment)

To allow access from other machines on the network:

```bash
R -e "shiny::runApp('app.R', port = 3838, host = '0.0.0.0')"
```

Then access via `http://<server-hostname>:3838` from a remote machine.

#### Important Notes

- **`host = '0.0.0.0'`**: This binds the server to all network interfaces, allowing external access. Without this option, Shiny defaults to `127.0.0.1` (localhost only).
- **Firewall**: Ensure port 3838 is open on the server's firewall. Contact your system administrator if needed.
- **Security Warning**: Exposing Shiny apps directly to the internet is not recommended for production use. Consider:
  - Using a reverse proxy (nginx/Apache) with HTTPS
  - Deploying to a managed Shiny Server or ShinyProxy
  - For FGCZ: Use the official proxy at `https://fgcz-shiny.uzh.ch/`
- **Background Execution**: For long-running sessions, consider using `nohup` or `screen`:
  ```bash
  nohup R -e "shiny::runApp('app.R', port = 3838, host = '0.0.0.0')" &
  ```

### Using Demo Mode

1. Select "Demo (Pre-computed)" mode
2. Choose a dataset from the dropdown:
   - `kokai_100K` - 100K coverage assembly
   - `kokai_1M_11x` - 1M dataset with 11x coverage
   - `kokai_1M_24x` - 1M dataset with 24x coverage
3. Click "Load Results"
4. Explore the Assembly and Results tabs

### Running New Assemblies

1. Select "Run Assembly" mode
2. Choose an input FASTA file
3. Adjust parameters:
   - **K-mer Size**: 15-63 (default: 51)
   - **Error Correction Rounds**: 0-10 (default: 3)
   - **CPU Threads**: 1-8 (default: 4)
4. Click "Run Assembly"
5. Monitor progress in real-time

## Project Structure

```
MINT_assembly_2026/
├── app.R                    # Main Shiny application
├── R/
│   ├── run_hifiasm.R       # hifiasm execution logic
│   ├── parse_results.R     # Result parsing and statistics
│   └── plot_functions.R    # Visualization functions
├── www/
│   └── custom.css          # Custom styling
├── MINT_Tage/              # Sample data
│   ├── *.fasta             # Input FASTA files
│   ├── kokai_100K/         # Pre-computed results
│   ├── kokai_1M_11x/       # Pre-computed results
│   └── kokai_1M_24x/       # Pre-computed results
└── README.md
```

## Deployment to FGCZ Shiny Server

### Prerequisites

1. Ensure all R packages are installed on the server
2. Verify hifiasm is available (optional, for Run mode)
3. Upload the project files

### Deployment Steps

1. Upload the entire project directory to the FGCZ Shiny server
2. Configure the app in the server's apps directory
3. Access via https://fgcz-shiny.uzh.ch/ after authentication

### Server Configuration

For FGCZ proxy server deployment, you may need to:

1. Set appropriate timeout values for long-running assemblies
2. Configure resource limits (memory, CPU)
3. Ensure the MINT_Tage data directory is accessible

## Assembly Parameters

| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| K-mer Size | 15-63 | 51 | Larger = more specific, may miss overlaps |
| Error Rounds | 0-10 | 3 | More rounds = better accuracy, slower |
| Min Contig | 0-10000 | 500 | Filter threshold for statistics |
| Threads | 1-8 | 4 | Parallel processing threads |

## Output Statistics

- **N50**: Sequence length at 50% of total assembly
- **L50**: Number of contigs to reach N50
- **N90/L90**: Similar metrics at 90%
- **GC Content**: Percentage of G+C bases
- **Total/Filtered Contigs**: With length filter applied

## References

Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., & Li, H. (2021). 
Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. 
*Nature Methods*, 18(2), 170-175.

## License

MIT License
