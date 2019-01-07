library(Pareto)
library(ggplot2)

ggplot(data.frame(x = c(0,6000)), aes(x)) + stat_function(fun = function(x) PiecewisePareto_PDF(x, c(1000,3000), c(1,2), truncation = 5000, truncation_type = "wd"), n=1000, color = "red") +
  stat_function(fun = function(x) PiecewisePareto_PDF(x, c(1000,3000), c(1,2), truncation = 5000, truncation_type = "lp"), n=1000, color = "green") +
  stat_function(fun = function(x) PiecewisePareto_PDF(x, c(1000,3000), c(1,2)), n=1000)

plot(cumsum(PiecewisePareto_PDF(1:6000, c(1000,3000), c(1,2), truncation = 5000, truncation_type = "lp")) - PiecewisePareto_CDF(1:6000, c(1000,3000), c(1,2), truncation = 5000, truncation_type = "lp"))



