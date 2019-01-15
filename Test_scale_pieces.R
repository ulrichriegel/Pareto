library(Pareto)
library(ggplot2)

number_of_simulations <- 10000

ap <- (1:5) * 1000
alpha <- c(1,3,2,4,2)

scale_pieces <- c(0.5,1,1,4,1)

losses_wo_scaling <- rPiecewisePareto(number_of_simulations, ap, alpha)
losses_with_scaling <- rPiecewisePareto(number_of_simulations, ap, alpha, scale_pieces = scale_pieces)

df_plot <- data.frame(losses_wo_scaling = losses_wo_scaling, losses_with_scaling = losses_with_scaling)

ggplot(df_plot) + stat_ecdf(aes(x = losses_with_scaling), color = "red") +
  stat_ecdf(aes(x = losses_wo_scaling)) +
  coord_cartesian(xlim = c(0,6000))
