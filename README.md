# dbn: Dynamic Bilinear Network Models

An R package for analyzing temporal ordinal relational data using Dynamic Bilinear Network (DBN) models.

Here's a revised version without the broken link:

## Installation

### Prerequisites

This package uses OpenMP for parallel computing. Setup by platform:

**macOS:**

1. Install Xcode Command Line Tools: 

```bash
xcode-select --install
```

2. Install gfortran (includes OpenMP support):

- For R 4.5.0+: Download [gfortran-14.2-universal.pkg](https://github.com/R-macos/gfortran-for-macOS/releases)
- For R 4.4.x and earlier: Download [gfortran-12.2-universal.pkg](https://mac.r-project.org/tools/gfortran-12.2-universal.pkg)
- Double-click the downloaded .pkg file to install

**Linux:**

- OpenMP is typically pre-installed with GCC
- If compilation fails, install: 

```bash
# Ubuntu/Debian
sudo apt-get install libomp-dev

# RHEL/CentOS/Fedora
sudo yum install libomp-devel
```

**Windows:**

- Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version
- OpenMP is included automatically

### Install Package

```r
# Install from GitHub
devtools::install_github("netify-dev/dbn")
```

### Troubleshooting

If installation fails, verify your setup:

```r
# Test if compilers work
install.packages("minpack.lm", type = "source")

# Check your R version
R.version.string
```

**Common macOS Issues:**

- If you see "gfortran not found", add to your `~/.zshrc`:

```bash
export PATH="/opt/gfortran/bin:$PATH"
```

Then restart your terminal and R.

- For persistent OpenMP errors on macOS, you may need to create `~/.R/Makevars` with appropriate compiler flags. Here an example:

```bash
# gfortran 14.2 configuration for R 4.5.1 on ARM64 Mac
FC = /opt/gfortran/bin/gfortran
F77 = /opt/gfortran/bin/gfortran

# Use the actual path where libgfortran.dylib is located
FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0 -L/opt/gfortran/lib -lgfortran -lquadmath -lm

# Compiler flags
FCFLAGS = -mmacosx-version-min=11.0
FFLAGS = -mmacosx-version-min=11.0

# brew
CPPFLAGS = -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib
```

## Quick Start

```r
library(dbn)

# Set number of threads for parallel computation (optional)
# By default, dbn uses 1
# Enable parallelization for better performance:
set_dbn_threads(4)  # Use 4 threads

# Check current thread setting
get_dbn_threads()

# Load example data
data(example_data)

# Run static model
results_static <- dbn(
  data = example_data,
  model = "static",
  nscan = 1000,  # iterations after burn-in
  burn = 500,    # burn-in period
  odens = 1     # keep every iteration
)

# View results
summary(results_static)
plot(results_static)

# Run dynamic model (requires more data)
results_dynamic <- dbn(
  data = example_data,
  model = "dynamic", 
  nscan = 2000,  # iterations after burn-in
  burn = 1000,   # burn-in period
  odens = 5      # keep every 5th iteration
)
```

## Models

### Static DBN

- Fixed sender/receiver effects over time
- Bilinear structure: `Θ_t = A Θ_{t-1} B' + σE_t`
- More stable, good for smaller datasets

### Dynamic DBN

- Time-varying A_t and B_t matrices
- Evolution: `Θ_t = A_t Θ_{t-1} B_t' + σE_t`
- Requires larger datasets (20+ actors, 50+ time points)

## Main Functions

- `dbn()`: Main wrapper for analysis
- `plot()`: Diagnostic plots
- `summary()`: Model summary statistics
- `check_convergence()`: MCMC diagnostics
- `compare_dbn()`: Compare multiple models

## Data Format

Input data should be a 4-dimensional array:

- Dimension 1: Sender actors
- Dimension 2: Receiver actors  
- Dimension 3: Relation types
- Dimension 4: Time points

## Author

- [Tosin Salau](https://polisci.msu.edu/people/directory/salau-tosin.html) <salaubol@msu.edu>  
- [Shahryar Minhas](https://s7minhas.com) <sminhas@msu.edu>


## License

MIT License - see [LICENSE](LICENSE) file for details.
