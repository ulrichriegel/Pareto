t <- c(1000,3000,5000)
Cover <- 8000
AttPoint <- 0
truncation <- 6000
alpha <- c(0.7,1.5,2)
truncation_type ="wd"

NumberOfSimulations <- 1000000

losses <- rPiecewisePareto(NumberOfSimulations, t, alpha, truncation = truncation, truncation_type = truncation_type)
xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
mean(xs_losses)

PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)

mean(xs_losses) / PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type)




