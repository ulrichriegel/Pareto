library(Pareto)
library(ggplot2)

ap <- (1:20) * 100
el <- c(rep(100,10), rep(50, 10))
el[1] <- 200
el[4] <- 110


res <- PiecewisePareto_Match_Layer_Losses(ap, el)
res

ggplot(data.frame(x = c(0,3000)), aes(x = x)) + stat_function(fun = function(x) PiecewisePareto_CDF(x, res$t, res$alpha), n=1000)

PiecewisePareto_Layer_Mean(1800, 100, res$t, res$alpha) * res$FQ
sum(el[1:18])
PiecewisePareto_Layer_Mean(1900, 100, res$t, res$alpha) * res$FQ
sum(el[1:19])
PiecewisePareto_Layer_Mean(1800, 200, res$t, res$alpha) * res$FQ
sum(el[2:19])
PiecewisePareto_Layer_Mean(Inf, 2000, res$t, res$alpha) * res$FQ
sum(el[20])


res <- PiecewisePareto_Match_Layer_Losses(ap, el, truncation = 2300)
ggplot(data.frame(x = c(0,3000)), aes(x = x)) + stat_function(fun = function(x) PiecewisePareto_CDF(x, res$t, res$alpha, truncation = 2300), n=1000)

PiecewisePareto_Layer_Mean(1800, 100, res$t, res$alpha, truncation = 2300) * res$FQ
sum(el[1:18])
PiecewisePareto_Layer_Mean(1900, 100, res$t, res$alpha, truncation = 2300) * res$FQ
sum(el[1:19])
PiecewisePareto_Layer_Mean(1800, 200, res$t, res$alpha, truncation = 2300) * res$FQ
sum(el[2:19])
PiecewisePareto_Layer_Mean(Inf, 2000, res$t, res$alpha, truncation = 2300) * res$FQ
sum(el[20])

