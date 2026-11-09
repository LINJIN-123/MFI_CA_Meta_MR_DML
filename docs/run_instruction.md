# Reproduction Instructions

This document explains how to reproduce every result in the manuscript from a clean computer.
It is written for researchers with no prior command-line experience. Two software stacks are used:
**R** (meta-analysis and two-sample Mendelian randomization) and **Python** (double machine
learning). Run the R and Python parts independently; they do not depend on each other.

---

## Part A — R analyses (meta-analysis and two-sample MR)

### A1. Install the software (once)

1. Download and install **R** for Windows from <https://cran.r-project.org/bin/windows/base/>.
2. Download and install **RStudio Desktop (Free)** from <https://posit.co/download/rstudio-desktop/>.
3. Open RStudio.

### A2. Restore the locked R environment

The project uses **`renv`** so that every package is installed at the exact version used in the
manuscript. The lock file is `R_analysis/renv.lock`.

In the RStudio console, run the following once. (If you prefer a graphical workflow, use
*File → Open Project* to open the project, then accept the prompt to install the locked packages.)

```r
# Install the renv package if it is missing
install.packages("renv")

# Point R at the folder that contains the project (edit the path to your own location)
setwd("D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML")

# Activate renv and install every package recorded in renv.lock
renv::activate()
renv::restore()
```

#### How the lock file was created (for reference / future updates)

```r
# Run inside the project after all required packages are installed
renv::init()          # initialises renv and snapshots the loaded packages
renv::snapshot()      # updates renv.lock after adding/updating a package
renv::status()        # checks that the library matches the lock file
```

Key R packages recorded: `meta`, `metafor`, `robumeta`, `clubSandwich`, `readxl`,
`TwoSampleMR`, `ieugwasr`, `data.table`, `remotes`.

### A3. Provide the OpenGWAS token (do not hard-code it)

The 2SMR scripts call the OpenGWAS API. Provide your personal token through an environment
variable rather than pasting it into the code.

```r
# Set for the current R session (replace with your own token; do not save this line in a committed file)
Sys.setenv(OPENGWAS_TOKEN = "paste-your-token-here")

# The scripts authenticate with:
ieugwasr::get_access_token()
```

To make the variable available permanently on Windows, search *"Edit environment variables for your
account"* in the Start menu and add a variable named `OPENGWAS_TOKEN`.

### A4. Run the scripts in numbered order

Open each file in `R_analysis/` and click **Source** (top-right of the editor), or run:

```r
setwd("D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/R_analysis")
source("01_meta_analysis.R")          # pooled aHR, heterogeneity, bias, subgroup RVE
source("02_forward_MR.R")             # MFI/AS -> cancers
source("03_reverse_MR.R")             # cancers -> MFI
source("04_sensitivity_FDR_MDOR.R")   # sensitivity, BH-FDR, MDOR
source("05_forest_plot.R")            # forest and funnel plots
```

Tabular outputs are written to `output_results/`; the observational meta-analysis input table
(an Excel/CSV of extracted aHR values) is read by `01_meta_analysis.R` — update its file path in
the script header if you store it elsewhere.

### A5. Notes on data access

- All GWAS data are retrieved online; no large summary-statistic file is stored in the repository.
- If a FinnGen outcome cannot be read online, place the corresponding `.gz` file locally and use
  the local-fallback branch in the MR script. The raw `.gz` file stays on your computer and is
  excluded from Git by `.gitignore`.

---

## Part B — Python analyses (double machine learning)

### B1. Install the software (once)

1. Download and install **Miniconda** for Windows from
   <https://docs.conda.io/en/latest/miniconda.html> (accept the default options).
2. Open **Anaconda Prompt** from the Start menu.

### B2. Recreate the environment from the lock file

```bat
:: Move to the project folder (edit the path as needed)
cd /d "D:\研究生科研小组\小论文\MFI-cancer-meta-MR-DML\Python_DML"

:: Create the environment exactly as recorded in environment.yml
conda env create -f environment.yml

:: Activate it
conda activate dml-repro
```

Alternative using pip only (if you do not use conda):

```bat
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

#### How the environment files were created (for reference / future updates)

```bat
:: conda export
conda env export --no-builds > environment.yml

:: pip freeze
pip freeze > requirements.txt
```

### B3. Run the scripts in order

```bat
cd /d "D:\研究生科研小组\小论文\MFI-cancer-meta-MR-DML\Python_DML"
conda activate dml-repro
python 01_dml_main.py        # 5-fold cross-fitting DML (LassoCV) + random forest sensitivity
python 02_figure3_plot.py    # Figure 3 visualisation panel (300 dpi)
```

The DML scripts read the harmonised instrument table (a small CSV) produced within the workflow.
Raw FinnGen `.gz` inputs are read only if present locally and are never committed.

---

## Expected run time

- Meta-analysis: seconds.
- Two-sample MR (all datasets, online): depends on network speed, typically tens of minutes.
- DML: under one minute on a standard laptop.

## Troubleshooting

- **API authentication error**: confirm `OPENGWAS_TOKEN` is set and not expired.
- **Package installation fails**: keep R and RStudio updated; on a corporate network, install
  packages through the RStudio *Tools → Install Packages* menu.
- **`conda` not recognised**: open the Anaconda Prompt rather than the ordinary command prompt.
- **File-not-found errors**: check that the working directory matches the path in each script.
