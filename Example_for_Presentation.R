





















library(Pareto)
library(ggplot2)

thresholds <- c(1000,1500,2000,2500,3000,4000,5000,10000)
expexted_losses <- c(100,90,50,40,60,40,50,50)
lowest_frequency <- 0.25

PP_Fit <- PiecewisePareto_Match_Layer_Losses(thresholds, expexted_losses, FQ_at_lowest_AttPt = lowest_frequency)

ggplot(data.frame(x = c(0, 15000)), aes(x = x)) +
  stat_function(fun = function(x) PiecewisePareto_CDF(x, PP_Fit$t, PP_Fit$alpha))

PiecewisePareto_Layer_Mean(500, 1000, PP_Fit$t, PP_Fit$alpha)
