t <- c(1000,3000,5000)
#t <- 1000
Cover <- 1000
AttPoint <- 9000
truncation <- 10000
alpha <- c(0.7,1.5,2)
#alpha <- 1.5
truncation_type <- "lp"

NumberOfSimulations <- 1000000

losses <- rPiecewisePareto(NumberOfSimulations, t, alpha, truncation = truncation, truncation_type = truncation_type)
xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
mean(xs_losses^2)

PiecewisePareto_Layer_SM(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)

mean(xs_losses^2) / PiecewisePareto_Layer_SM(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)

var(xs_losses) / PiecewisePareto_Layer_Var(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)



