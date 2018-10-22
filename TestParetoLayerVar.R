t <- 5000
Cover <- Inf
AttPoint <- 1000
truncation <- 20000
alpha <- 0.7

NumberOfSimulations <- 1000000

losses <- rPareto(NumberOfSimulations, t, alpha, truncation = truncation)
xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
var(xs_losses)

Pareto_Layer_Var(Cover, AttPoint, alpha, t, truncation)

var(xs_losses) / Pareto_Layer_Var(Cover, AttPoint, alpha, t, truncation)




