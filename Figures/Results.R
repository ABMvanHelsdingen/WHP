# baseline integral
integral_vectorized <- function(parameters, times){
  
  # times must be ascending order, 1st entry can be zero
  integrals <- numeric(length(times))
  slots <- ceiling(rincr*max(times))
  if (slots == 0){ # Only 1 time, which is zero
    return(0)
  }
  lambdas <- numeric(slots)
  
  
  # y = exp(b0 + b%*%x)
  b0 <- parameters[[1]] # state 0
  b1 <- parameters[[2]] # state 1/2
  b2 <- parameters[[3]] # depth coefficient
  b3 <- parameters[[4]] # rate coefficient
  
  for(i in 1:slots){
    depth <- depths[i]
    rate <- rates[i]
    state <- states[i]
    if (state == 0){
      lambdas[i] <- exp(b0 + b2*depth + b3*rate)
    } else{
      lambdas[i] <- exp(b0 + b1 + b2*depth + b3*rate)
    }
  }
  
  for(j in 1:length(times)){
    slots <- ceiling(rincr * times[j])
    if (slots > 1){
      integrals[j] <- (1/rincr) * sum(lambdas[1:(slots-1)]) + (times[j] - (1/rincr)*(slots-1)) * lambdas[slots]
    } else if (times[j] == 0){
      integrals[j] <- 0
    } else { # 0 < times[i] < (1/rincr)
      integrals[j] <- times[j] * lambda[1]
    }
    
  }
  
  return(integrals)
}

show_GOF <-  function(coefs, times, mu = NULL, plot = TRUE, return_values = FALSE) {
  # For constructing GOF plot with custom coefficients
  # use show_Hawkes_GOF2 if TMB output is available
  # mu is the integral of the CBF
    
  k <- coefs[length(coefs)]
  beta <- coefs[(length(coefs)-1)]
  alpha <- coefs[(length(coefs)-2)]
  background_parameters <- coefs[1:(length(coefs)-3)]
  marks <- rep(1, length(times))
  
  A <- numeric(length = length(times))
  for(i in 2:length(times)) {
    A[i] <- exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1])
  }
  
  compensator <- numeric(length = length(times))
  baselines <- mu(background_parameters, times)
  for(i in 1:length(times)){
    compensator[i] <- baselines[i] + ((alpha / beta) * (sum(marks[1:i]) - marks[i] - A[i]))
  }
    
  compensator <- compensator - mu(background_parameters,0) ## Subtract integral at zero

  interarrivals <- compensator[2:length(compensator)] - compensator[1:(length(compensator)-1)]
  ## Kolmogorov-Smirnov Test

  weibull_lambda = 1 / gamma(1+(1/k))
  print(stats::ks.test(interarrivals, "pweibull", shape = k, scale = weibull_lambda))

  ## Ljung-Box Test
  print(stats::Box.test(interarrivals, type = "Ljung"))
  if (plot) {
    data <- data.frame(xs = times, observed = 1:length(times), compensator = compensator)
    ## Plot of compensator and actual events
    data <- reshape(data, direction = "long", idvar = "xs",
                    varying = c("observed", "compensator"), v.names = "val",
                    times = c("observed", "compensator"),
                    new.row.names = 1:(2*length(times)))
    lineplot <- ggplot2::ggplot(data = data,
                                ggplot2::aes(x = .data$xs, y = .data$val, colour = .data$time)) +
      ggplot2::xlab("Time") +
      ggplot2::ylab("Events") +
      ggplot2::geom_line() +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position=c(0.8,0.2),  text = element_text(size = 25)) +
      ggplot2::ggtitle("Actual Events and Compensator")
    ## Histogram of transformed interarrival times with exponential distribution overlaid
    data <- data.frame(data = interarrivals)
    hist <-  ggplot2::ggplot(data = data,  ggplot2::aes(x = .data$data)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 0.02) + ggplot2::theme_minimal() +
      ggplot2::xlab("Compensator Differences") +  ggplot2::ylab("Density") +
      ggplot2::stat_function(fun = dweibull, args = list(shape = k, scale = weibull_lambda), color = "red") +
      ggplot2::xlim(0,max(1.25, quantile(interarrivals, probs=0.99))) +
      ggplot2::theme(text = element_text(size = 25)) +
      ggplot2::ggtitle("Compensator Differences")
    
    p <- ppoints(100)    ## 100 equally spaced points on (0,1), excluding endpoints
    q <- quantile(interarrivals,p = p) ## percentiles of the sample distribution
    ## Q-Q plot of transformed interarrival times
    data <- data.frame(x = qweibull(p, shape = k, scale = weibull_lambda), y = q)

    qqplot <- ggplot2::ggplot(data =  data,
                              ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::xlab("Theoretical Quantiles") +
      ggplot2::ylab("Observed Quantiles") + 
      ggplot2::geom_point() +  ggplot2::theme_minimal() + 
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
      ggplot2::theme(text = element_text(size = 25)) +
      ggplot2::ggtitle("Compensator Differences")
    
    U <- pweibull(interarrivals, shape = k, scale = weibull_lambda)

    ## Scatterplot of the CDF of consecutive interarrival times
    data <-  data.frame(x = U[1:(length(U)-1)], y = U[2:length(U)])
    scatter <- ggplot2::ggplot(data = data,
                               ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::xlab("F(Interarrival time k)") +
      ggplot2::ylab("F(Interarrival time k+1)") + 
      ggplot2::geom_point() +  ggplot2::theme_minimal() +
      ggplot2::theme(text = element_text(size = 25)) +
      ggplot2::ggtitle("Consecutive Interarrival Times")
    gridExtra::grid.arrange(lineplot, qqplot, hist, scatter, ncol = 2)
  }
  if(return_values) {
    return(list(interarrivals = interarrivals))
  }
}

# Show compensator differences
show_cd <-  function(coefs, times, mu = NULL, plot = TRUE, return_values = FALSE) {
  # Same as the Bottom Left plot in show_GOF()
  
  k <- coefs[length(coefs)]
  beta <- coefs[(length(coefs)-1)]
  alpha <- coefs[(length(coefs)-2)]
  background_parameters <- coefs[1:(length(coefs)-3)]
  marks <- rep(1, length(times))
  
  A <- numeric(length = length(times))
  for(i in 2:length(times)) {
    A[i] <- exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1])
  }
  
  compensator <- numeric(length = length(times))
  baselines <- mu(background_parameters, times)
  for(i in 1:length(times)){
    compensator[i] <- baselines[i] + ((alpha / beta) * (sum(marks[1:i]) - marks[i] - A[i]))
  }
  
  compensator <- compensator - mu(background_parameters,0) ## Subtract integral at zero
  
  interarrivals <- compensator[2:length(compensator)] - compensator[1:(length(compensator)-1)]
  ## Kolmogorov-Smirnov Test
  
  weibull_lambda = 1 / gamma(1+(1/k))
  print(stats::ks.test(interarrivals, "pweibull", shape = k, scale = weibull_lambda))
  
  ## Ljung-Box Test
  print(stats::Box.test(interarrivals, type = "Ljung"))
  if (plot) {
    ## Histogram of transformed interarrival times with Weibull distribution overlaid
    data <- data.frame(data = interarrivals)
    hist <-  ggplot2::ggplot(data = data,  ggplot2::aes(x = .data$data)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 0.02) + ggplot2::theme_minimal() +
      ggplot2::xlab("Compensator Differences") +  ggplot2::ylab("Density") +
      ggplot2::stat_function(fun = dweibull, args = list(shape = k, scale = weibull_lambda), color = "red") +
      ggplot2::xlim(0,max(1.25, quantile(interarrivals, probs=0.99))) +
      ggplot2::theme(text = element_text(size = 25)) +
      ggplot2::ggtitle("Compensator Differences")
    
    gridExtra::grid.arrange(hist, ncol = 1)
  }
  if(return_values) {
    return(list(interarrivals = interarrivals))
  }
}

# Presentation of modeled intensities, depths and clicks
show_lambda <-  function(alpha, beta, coefs, X, depths, times, incr, n = length(depths), log = TRUE) {
  p <- seq(0, incr * (nrow(X) - 1), length.out = n)
  idx <- ceiling(p / incr); idx[1] = 1
  baseline <- X %*% coefs
  baseline <- baseline[idx]
  lambda <- numeric(n)
  for(i in 1:length(p)){
    A <- exp(-beta * (p[i] - times))[times < p[i]]
    if(length(A) == 0){
      A <- 0
    }
    if (log){
      lambda[i] = log(exp(baseline[i]) + alpha * sum(A))
    } else {
      lambda[i] = exp(baseline[i]) + alpha * sum(A)
      baseline[i] = exp(baseline[i])
    }
  }
  
  if (log){
    lmax <- max(lambda)
    ylab <- expression("log("*lambda*")")
  } else {
    lmax <- quantile(lambda, 0.95)
    ylab <- expression(lambda)
  }
  lmin <- min(lambda)
  
  data <-  data.frame(xp = p, lambda = lambda, mu = baseline)
  data$xp <- data$xp/3600 # seconds to hours
  data <- reshape(data, direction = "long", idvar = "p",
                  varying = c("lambda", "mu"), v.names = "val",
                  times = c("lambda", "mu"),
                  new.row.names = 1:(2*length(lambda)))
  line <- ggplot2::ggplot(data = data,
                          ggplot2::aes(x = .data$xp, y = .data$val, color = .data$time)) +
    ggplot2::xlab("Time (hr)") +
    ggplot2::ylab(ylab) + 
    ggplot2::geom_line(lwd = 2.5) +  ggplot2::theme_minimal() +
    ggplot2::ylim(c(lmin, lmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25))
  gridExtra::grid.arrange(line, ncol = 1)
}

# Presentation of modeled intensities, depths and clicks
show_mu <-  function(alpha, beta, coefs, X, depths, times, incr, n = length(depths), log = TRUE) {
  p <- seq(0, incr * (nrow(X) - 1), length.out = n)
  idx <- ceiling(p / incr); idx[1] = 1
  baseline <- X %*% coefs
  baseline <- baseline[idx]
  
  if (log){
    lmax <- max(baseline)
    ylab <- expression("log("*lambda*")")
  } else {
    baseline <- exp(baseline)
    lmax <- max(baseline)
    ylab <- expression(lambda)
  }
  lmin <- min(baseline)
  
  data <-  data.frame(xp = p, mu = baseline, depths = 100*depths[idx], rates = 100*X[idx,4])
  data$xp <- data$xp/3600 # seconds to hours
  line <- ggplot2::ggplot(data = data,
                          ggplot2::aes(x = .data$xp, y = .data$mu)) +
    ggplot2::xlab("Time (hr)") +
    ggplot2::ylab(ylab) + 
    ggplot2::geom_line(lwd = 2.5) +  ggplot2::theme_minimal() +
    ggplot2::ylim(c(lmin, lmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25))
  
  line2 <- ggplot2::ggplot(data = data,
                          ggplot2::aes(x = .data$xp, y = .data$depths)) +
    ggplot2::xlab("Time (hr)") +
    ggplot2::ylab("Depth (m)") + 
    ggplot2::geom_line(lwd = 2.5) +  ggplot2::theme_minimal() +
    #ggplot2::ylim(c(lmin, lmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25))
  
  line3 <- ggplot2::ggplot(data = data,
                          ggplot2::aes(x = .data$xp, y = .data$rates)) +
    ggplot2::xlab("Time (hr)") +
    ggplot2::ylab("Rate of Descent (m/s)") + 
    ggplot2::geom_line(lwd = 2.5) +  ggplot2::theme_minimal() +
    #ggplot2::ylim(c(lmin, lmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25))
  
  gridExtra::grid.arrange(line, line2, line3, nrow = 3)
}

show_fig4 <-  function(alpha, beta, coefs, X, depths, times, incr, n = length(depths)) {
  p <- seq(0, incr * (nrow(X) - 1), length.out = n)
  idx <- ceiling(p / incr); idx[1] = 1
  baseline <- X %*% coefs
  baseline <- baseline[idx]
  lambda <- numeric(n)
  for(i in 1:length(p)){
    A <- exp(-beta * (p[i] - times))[times < p[i]]
    if(length(A) == 0){
      A <- 0
    }
      lambda[i] = exp(baseline[i]) + alpha * sum(A)
      baseline[i] = exp(baseline[i])
    }
  

  mmax <- max(baseline); mmin <- min(baseline)
  lmax <- max(lambda); lmin <- min(baseline)
  ylab0 <- expression(lambda(t))
  ylab1 <- expression(mu(t))
  
  data <-  data.frame(xp = p, lambda = lambda, mu = baseline, depths = -100*depths[idx], rates = -100*X[idx,4])
  #data$xp <- data$xp/3600 # seconds to hours
  
  line0 <- ggplot2::ggplot(data = data,
                          ggplot2::aes(x = .data$xp, y = .data$lambda)) +
    #ggplot2::xlab("Time (hr)") +
    ggplot2::ylab(ylab0) + 
    ggplot2::geom_line(lwd = 2.25) +  ggplot2::theme_minimal() +
    ggplot2::ylim(c(lmin, lmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25),
                   axis.title.x = element_blank())
  
  line1 <- ggplot2::ggplot(data = data,
                          ggplot2::aes(x = .data$xp, y = .data$mu)) +
    #ggplot2::xlab("Time (hr)") +
    ggplot2::ylab(ylab1) + 
    ggplot2::geom_line(lwd = 2.25) +  ggplot2::theme_minimal() +
    ggplot2::ylim(c(mmin, mmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25),
                   axis.title.x = element_blank())
  
  line2 <- ggplot2::ggplot(data = data,
                           ggplot2::aes(x = .data$xp, y = .data$depths)) +
    ggplot2::xlab("Time (s)") +
    ggplot2::ylab("Depth (m)") + 
    ggplot2::geom_line(lwd = 2.25) +  ggplot2::theme_minimal() +
    #ggplot2::ylim(c(lmin, lmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25))
  
  line3 <- ggplot2::ggplot(data = data,
                           ggplot2::aes(x = .data$xp, y = .data$rates)) +
    #ggplot2::xlab("Time (hr)") +
    ggplot2::ylab("Rate of Descent (m/s)") + 
    ggplot2::geom_line(lwd = 2.25) +  ggplot2::theme_minimal() +
    #ggplot2::ylim(c(lmin, lmax)) +
    ggplot2::theme(legend.position = "none", text = element_text(size = 25))
  
  image <- ggpubr::ggarrange(line0, line1, line2, nrow = 3)
  return(image)
}

# START OF SCRIPT
load("PhD-oversize/cleaned_whale_data.Rdata")
i = 53


times <- cleaned_data[[i]]$times
depths <- cleaned_data[[i]]$depths
rates <- cleaned_data[[i]]$rates
rincr <- 1/cleaned_data[[i]]$incr # reciprocal of the increment
  
# assign states
states <- numeric(length(depths))
for(j in 1:length(depths)){
  if (depths[j] > 0.2 && rates[j] >= 0){
    states[j] = 1 # descent
  } else if (depths[j] > 0.2 && rates[j] < 0){
    states[j] = 2 # ascent
  } # if depth < 20m, state = 0 (surface)
}

# create matrix to calculate mu(t)
X = matrix(0, nrow = length(depths), ncol = 4)
X[,3] = depths; X[,4] = rates;
X[,1] = rep(1, length(depths))
for(i in 1:length(depths)){
  if (depths[i] > 0.2){
    X[i, 2] = 1
  }
}

# Import coefs
Bayesian <- T
NIMBLEcoefs <- read.csv("IWP/NIMBLE_Coefs53.csv")[,2]

if (Bayesian){
  coefs <- NIMBLEcoefs[2:5]
  beta <- NIMBLEcoefs[6]; alpha <- NIMBLEcoefs[1] * NIMBLEcoefs[6]
}
  
  


# Create Figure 4

image <- show_fig4(alpha,beta,coefs,X,depths,times,1/rincr)
ggsave("fig4.jpeg", height = 15, width = 12)
  
  
out <- tryCatch(
  {
    jpeg(paste("PhD/Whales/GOF_",names(cleaned_data)[i],"_",i, ".jpeg", sep=""), width = 2000, height = 2000, quality=100)
    zzz = show_GOF(coefs,times, mu = integral_vectorized, return_values = TRUE)
    dev.off()
  }, 
  error=function(cond){
    print(paste("Whale ",i," Error creating cbf GOF/TIAT plots",sep=""))
    print(cond)
    }
  )


write.csv(cbfCoefs,paste("PhD/Whales/cbfCoefs_",i, ".csv", sep=""))
write.csv(cbfNLLs,paste("PhD/Whales/cbfNLLs_",i, ".csv", sep=""))


coefs <- mcmc.out2$summary[c(2,3,4,5,1,6,7),1]
coefs[5] <- coefs[6] * coefs[5]


# Plot Intensities

for(j in c(i)){
  
  
  times <- cleaned_data[[j]]$times
  depths <- cleaned_data[[j]]$depths
  rates <- cleaned_data[[j]]$rates
  incr <- cleaned_data[[j]]$incr
  
  
  
  X = matrix(0, nrow = length(depths), ncol = 4)
  X[,3] = depths; X[,4] = rates;
  X[,1] = rep(1, length(depths))
  for(i in 1:length(depths)){
    if (depths[i] > 0.2){
      X[i, 2] = 1
    }
  }
  
  alpha = coefs[5]; beta = coefs[6]
  coefs2 = coefs[1:4]

  # Only for model F: exp coef column 2
  #coefs[2] = exp(coefs[2])
  
  jpeg(paste("PhD/Whales/",names(cleaned_data)[j],"_",j, ".jpeg", sep=""), width = 2000, height = 2000, quality=100)
  zzz <- show_lambda(alpha, beta, coefs2, X, depths, times, cleaned_data[[j]]$incr, n = 4000, log = FALSE)
  dev.off()
  
  jpeg(paste("PhD/Whales/",names(cleaned_data)[j],"_",j, "_logged",".jpeg", sep=""), width = 2000, height = 2000, quality=100)
  zzz <- show_lambda(alpha, beta, coefs2, X, depths, times, cleaned_data[[j]]$incr, n = 4000, log = TRUE)
  dev.off()
  
}

jpeg(paste("PhD/Whales/",names(cleaned_data)[j],"_",j, "_mu.jpeg", sep=""), width = 2000, height = 2000, quality=100)
zzz <- show_mu(alpha, beta, coefs2, X, depths, times, cleaned_data[[j]]$incr, n = 4000, log = FALSE)
dev.off()

jpeg(paste("PhD/Whales/",names(cleaned_data)[j],"_",j, "_lambda.jpeg", sep=""), width = 2000, height = 2000, quality=100)
zzz <- show_fig4(alpha, beta, coefs2, X, depths, times, cleaned_data[[j]]$incr, n = 4000, log = FALSE)
dev.off()

# plot showing effect of underwater, depth and rate