library(Pareto)
t <- 1000
alpha <- 2
truncation <- 10000
#truncation <- Inf
x <- 1000*(1:12)-1000

y <- Pareto_CDF(x,t,alpha,truncation)
y

Pareto_Inverse_CDF(y,t,alpha,truncation)



t <- c(1000,2000,3000)
alpha <- c(2,3,4)
truncation <- 10000
#truncation <- Inf
x <- 1000*(1:12)-1000

y <- PiecewisePareto_CDF(x,t,alpha,truncation, truncation_type = "wd")
y

PiecewisePareto_Inverse_CDF(y,t,alpha,truncation,truncation_type = "wd")





t <- c(1000,2000,3000)
alpha <- c(2,3,4)
truncation <- 10000
#truncation <- Inf
x <- 1000*(1:12)-1000

y <- PiecewisePareto_CDF(x,t,alpha,truncation, truncation_type = "lp")
y

PiecewisePareto_Inverse_CDF(y,t,alpha,truncation,truncation_type = "lp")

