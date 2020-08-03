# Pareto

Pareto is an R package providing methods and tools for the Pareto, the piecewise Pareto and the generalized Pareto
distribution which are useful for pricing of reinsurance treaties:

- Distribution functions, densities and quantile functions
- Layer mean and variance
- Simulation of Pareto and piecewise Pareto distributions
- Pareto extrapolation
- Finding the Pareto alpha between the expected losses of two layers
- Finding the Pareto alpha between an excess frequency and the expected loss of a layer
- Maximum likelihood estimation of the alpha(s) of a (piecewise) Pareto distribution
- Calculation of local Pareto alphas for normal, lognormal and gamma distributions
- Fitting a Piecewise Pareto distribution to the expected losses of an arbitrary number of reference layers
  and the excess frequencies at given thresholds

Moreover, the package provides some functions for collective models with a claim count distribution from the Panjer class 
(i.e. Binomial, Poisson and Negative Binomial) and a piecewise Pareto distributed severity:

- Layer mean, variance and standard deviation
- Simulation of losses with the collective model

All methods are also available for truncated versions of the (piecewise) Pareto distribution.

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

## Thanks

To Stefan Foerster for his NPAddins which I used a lot in my daily work and which inspired me to create this package.
