t <- c(1000,3000,5000)
t <- 1000
Cover <- Inf
AttPoint <- 9000
truncation <- NULL
alpha <- c(1.5,1.5,1)
alpha <- 1
truncation_type <- "lp"

NumberOfSimulations <- 1000000

losses <- rPiecewisePareto(NumberOfSimulations, t, alpha, truncation = truncation, truncation_type = truncation_type)
xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
mean(xs_losses^2)

PiecewisePareto_Layer_SM(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)

mean(xs_losses^2) / PiecewisePareto_Layer_SM(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)

var(xs_losses) / PiecewisePareto_Layer_Var(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)



