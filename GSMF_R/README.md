#Galaxy Stellar Mass Function (GSMF) Analysis
This repository contains R scripts to analyze and fit Schechter functions (single and double) to galaxy stellar mass function data. It includes galaxy selection, volume calculations, maximum likelihood fitting, and plotting results.

##Features
- Fits single and double Schechter functions to star-forming galaxy stellar mass data
- Supports star-forming and passive galaxy populations
- Uses maximum likelihood estimation
- Generates diagnostic plots for fits and likelihood surfaces

##Installation
###R
Ensure you have R (version 4.0 or higher recommended) installed.

Install required R packages:
```r
install.packages(c("data.table", "mvtnorm", "magicaxis", "plotrix", "xtable"))
#For dftools and celestial packages, install from GitHub or relevant sources if not available on CRAN
# Example:
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#devtools::install_github("author/dftools")
#devtools::install_github("author/celestial")

For `dftools` and `celestial` packages, install from GitHub or other sources if not available on CRAN:

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("author/dftools")
devtools::install_github("author/celestial")

##Usage
Run the main R script to perform fits and generate plots. Modify parameters inside the script to change fitting ranges, population selections, and output options.
```r
source("scripts/main_fit.R")

## Requirements
Install the dependencies with:

```bash
pip install -r requirements.txt

##License
This project is intended for academic and educational use. For full terms and conditions, please refer to the LICENSE file.
