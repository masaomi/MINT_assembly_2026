# install_packages.R
# Script to install required R packages for MINT Assembly Demo

# CRAN packages
cran_packages <- c(
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
)

# Install missing CRAN packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    cat("Installing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, repos = "https://cloud.r-project.org")
  } else {
    cat("All packages are already installed.\n")
  }
}

cat("Checking and installing required packages...\n")
install_if_missing(cran_packages)

cat("\nLoading packages to verify installation...\n")
for (pkg in cran_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  [OK] %s\n", pkg))
  } else {
    cat(sprintf("  [ERROR] %s failed to load\n", pkg))
  }
}

cat("\nPackage installation complete!\n")
cat("You can now run the app with: shiny::runApp('.')\n")

