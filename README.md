# Pareto

Pareto is an R package providing methods and tools for the Pareto and the piecewise Pareto distribution which are 
useful for pricing of reinsurance treaties:

- Distribution functions, densities and quantile functions
- Layer mean and variance
- Simulation of Pareto and piecewise Pareto distributions
- Pareto extrapolation
- Finding the Pareto alpha between the expected losses of two layers
- Finding the Pareto alpha between an excess frequency and the expected loss of a layer
- Fitting a Piecewise Pareto distribution to a tower of layer losses

## Installation

To install the current development version from github you need the [devtools package](https://cran.r-project.org/package=devtools) and the other packages on which Pareto depends and links to:

```s
install.packages(c("knitr", "rmarkdown"))
```

To install Pareto run:
```s
library(devtools)
install_github("ulrichriegel/Pareto", build_vignettes = TRUE)
```

## Usage

```s
library(Pareto)
```

## License

This package is free and open source software, licensed under [GPL](https://www.gnu.org/copyleft/gpl.html).
