# Functions to simulate Weibull-Hawkes processes

WHPintegral <- function(t,mu,A,alpha,beta){
  return((t*mu) + (alpha*A/beta)*(1 - exp(-beta*t)))
}

WHPsolve <- function(min, max, x, mu,A,alpha,beta, iter = 10){
  
  bounds <- matrix(0, nrow = iter, ncol = 2)
  bounds[1, 1] <- min; bounds[1, 2] <- max
  for(j in 2:iter){
    t_min <- bounds[(j-1),1]
    t_max <- bounds[(j-1),2]
    x_min <- WHPintegral(t_min, mu,A,alpha,beta)
    x_max <- WHPintegral(t_max, mu,A,alpha,beta)
    bounds[j, 1] <- ((x - x_min)/(mu + alpha*A*exp(-beta*t_min))) + t_min
    bounds[j, 2] <- ((x - x_max)/(mu + alpha*A*exp(-beta*t_max))) + t_max
  }
  
  return(list(min=bounds[iter,1],max=bounds[iter,2]))
}

simWHP <- function(mu, alpha, beta, k, t, tol = 1.e-8){
  errors <- c(0)
  Ai <- 0
  
  # First point can be calculated explicitly
  gk <- 1/gamma(1+(1/k))
  u <- runif(1)
  w <- gk * (-1 * log(u))^(1/k)
  times <- c(w/mu)
  n <- 1
  
  while (times[n] < t){
    u <- runif(1)
    w <- gk * (-1 * log(u))^(1/k)
    y <- WHPsolve((w/(mu+(alpha*(Ai+1)))), (w/mu), w, mu, Ai+1, alpha, beta)
    
    min <- y$min; max <- y$max
    while ((max - min) > tol){
      y <- WHPsolve(min, max, w, mu, Ai+1, alpha, beta)
      min <- y$min; max <- y$max
    }
    
    times <- append(times, times[n] + mean(min,max))
    errors <- append(errors, (max-min))
    n <- n + 1
    Ai <- (Ai + 1) * exp(-beta*(times[n] - times[(n-1)]))
    
    
  }
  
  if (n == 0 || (n == 1 && max(times) > t)){
    return(list(times=numeric(0), errors=numeric(0)))
  } else if (max(times)<t){
    return(list(times=times, errors=errors))
  } else {
    return(list(times=times[1:(n-1)], errors=errors[1:(n-1)]))
  }
  
  
}

WIHPsolve <- function(min, max, x, t0, Mu,A,alpha,beta, pars, iter = 10){
  
  minbound <- min
  maxbound <- max
  newt <- mean(c(minbound,maxbound))
  Mu_t0 <- Mu(t0, pars) # integral at the previous event
  
  for(i in 1:iter){
    integral <- Mu(newt + t0, pars) + (alpha*A/beta)*(1 - exp(-beta*newt))
    if (integral - Mu_t0 > x){
      maxbound <- newt
    } else {
      minbound <- newt
    }
    newt <- mean(c(minbound,maxbound))
  }
  
  return(list(min=minbound, max = maxbound))
}

simWIHP <- function(Mu, mu_max, mu_min, alpha, beta, pars, k, t, tol = 1.e-8){
  errors <- c(0)
  Ai <- 0
  
  # First event
  gk <- 1/gamma(1+(1/k))
  u <- runif(1)
  w <- gk * (-1 * log(u))^(1/k)
  minbound <- w/(mu_max +(alpha*(Ai+1)))
  maxbound <- w/(mu_min)
  iter <- ceiling((log(maxbound-minbound) - log(tol))/log(2))
  y <- WIHPsolve(minbound, maxbound, w, 0, Mu, 0, alpha, beta, pars, iter = iter)
  
  times <- c(mean(y$min,y$max))
  n <- 1
  
  while (times[n] < t){
    u <- runif(1)
    w <- gk * (-1 * log(u))^(1/k)
    minbound <- w/(mu_max +(alpha*(Ai+1)))
    maxbound <- w/(mu_min)
    iter <- ceiling((log(maxbound-minbound) - log(tol))/log(2))
    y <- WIHPsolve(minbound, maxbound, w, times[n], Mu, Ai+1, alpha, beta, pars, iter)
    
    times <- append(times, times[n] + mean(c(y$min,y$max)))
    errors <- append(errors, (y$max-y$min))
    n <- n + 1
    Ai <- (Ai + 1) * exp(-1*beta*(times[n] - times[(n-1)]))
    
    
  }
  
  if (n == 0 || (n == 1 && max(times) > t)){
    return(list(times=numeric(0), errors=numeric(0)))
  } else if (max(times)<t){
    return(list(times=times, errors=errors))
  } else {
    return(list(times=times[1:(n-1)], errors=errors[1:(n-1)]))
  }
  
  
}

simMWIHP <- function(Mu, mu_max, mu_min, alpha, beta, pars, 
                     k1, k2, p, mu1, t, tol = 1.e-8){
  errors <- c(0)
  Ai <- 0
  
  # Weibull means
  mean1 <- mu1/p
  mean2 <- (1 - mu1)/(1-p)
  
  # Weibull g(k)
  gk1 <- mean1/gamma(1+(1/k1))
  gk2 <- mean2/gamma(1+(1/k2))
  
  
  # First event
  m <- runif(1)
  u <- runif(1)
  if (m < p){ # Weibull 1
    w <- gk1 * (-1 * log(u))^(1/k1)
  } else { # Weibull 2
    w <- gk2 * (-1 * log(u))^(1/k2)
  }
    
  minbound <- w/(mu_max +(alpha*(Ai+1)))
  maxbound <- w/(mu_min)
  iter <- ceiling((log(maxbound-minbound) - log(tol))/log(2))
  y <- WIHPsolve(minbound, maxbound, w, 0, Mu, 0, alpha, beta, pars, iter = iter)
  
  times <- c(mean(y$min,y$max))
  n <- 1
  
  while (times[n] < t){
    m <- runif(1)
    u <- runif(1)
    if (m < p){ # Weibull 1
      w <- gk1 * (-1 * log(u))^(1/k1)
    } else { # Weibull 2
      w <- gk2 * (-1 * log(u))^(1/k2)
    }

    minbound <- w/(mu_max +(alpha*(Ai+1)))
    maxbound <- w/(mu_min)
    iter <- ceiling((log(maxbound-minbound) - log(tol))/log(2))
    y <- WIHPsolve(minbound, maxbound, w, times[n], Mu, Ai+1, alpha, beta, pars, iter)
    
    times <- append(times, times[n] + mean(c(y$min,y$max)))
    errors <- append(errors, (y$max-y$min))
    n <- n + 1
    Ai <- (Ai + 1) * exp(-1*beta*(times[n] - times[(n-1)]))
    
    
  }
  
  if (n == 0 || (n == 1 && max(times) > t)){
    return(list(times=numeric(0), errors=numeric(0)))
  } else if (max(times)<t){
    return(list(times=times, errors=errors))
  } else {
    return(list(times=times[1:(n-1)], errors=errors[1:(n-1)]))
  }
  
  
}


