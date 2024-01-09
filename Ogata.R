# Functions to simulate univariate Hawkes Process and univariate Weibull-Hawkes process
# Chen (2016a): https://www.math.fsu.edu/~ychen/research/Thinning%20algorithm.pdf
# Chen (2016b): https://www.math.fsu.edu/~ychen/research/multiHawkes.pdf

# Ogata's Modified Thinning Algorithm for univariate HP

UOgata <- function(mu, alpha, beta, t){
  # following Chen (2016a)
  # using the definition of A in Laub (2014) and stelfi
  times = numeric(0)
  s = 0
  n = 0
  A = 0
  lambda_bar = mu
  lambda_s = mu
  while (s < t){
    if (n > 0){
      # lambda_bar is intensity immediately after last accepted event
      lambda_bar = mu + alpha*(A + 1)
    }
    u <- runif(1)
    w <- -1*log(u)/lambda_bar
    s <- s + w
    D <- runif(1)
    if (n > 0){
      lambda_s <- mu + alpha*(A + 1)*exp(-beta*(s-max(times)))
    }
    if (D*lambda_bar <= lambda_s){
      if (n > 0){
        A <- (A + 1) * (exp(-beta*(s-max(times))))
      }
      n <- n + 1
      times <- append(times, s)
    }
  }
  # Check if last event > t, and return output
  if(max(times)<t && n > 0){
    return(times)
  } else if (n < 2){
    return(numeric(0))
  } else {
    return(times[1:(n-1)])
  }
}


UWOgata <- function(mu, alpha, beta, k, t){
  # Simulate a Hawkes process with Weibull IAT
  # following Chen (2016a)
  # using the definition of A in Laub (2014) and R package stelfi
  times = numeric(0)
  s = 0
  n = 0
  A = 0
  lambda_bar = mu
  lambda_s = mu
  weibull_factor = gamma(1+(1/k))
  while (s < t){
    if (n > 0){
      # lambda_bar is intensity immediately after last accepted event
      lambda_bar = mu + alpha*(A + 1)
    }
    u <- runif(1)
    #w <- -1*log(u)/lambda_bar
    weibull_lambda <- 1/(lambda_bar * weibull_factor)
    w <- weibull_lambda * (-1 * log(u))^(1/k)
    s <- s + w
    D <- runif(1)
    if (n > 0){
      lambda_s <- mu + alpha*(A + 1)*exp(-beta*(s-max(times)))
    }
    if (D*lambda_bar <= lambda_s){
      if (n > 0){
        A <- (A + 1) * (exp(-beta*(s-max(times))))
      }
      n <- n + 1
      times <- append(times, s)
    }
  }
  # Check if last event > t, and return output
  if(max(times)<t && n > 0){
    return(times)
  } else if (n < 2){
    return(numeric(0))
  } else {
    return(times[1:(n-1)])
  }
}

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

WIHPsolve <- function(min, max, x, Mu,A,alpha,beta, pars, iter = 10){
  
  minbound <- min
  maxbound <- max
  newt <- mean(c(minbound,maxbound))
  #print(minbound); print(maxbound); print(newt)
  #print(pars)
  
  for(i in 1:iter){
    integral <- Mu(newt, pars) + (alpha*A/beta)*(1 - exp(-beta*newt))
    #print(integral); print(x)
    if (integral > x){
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
  #print(mu_max); print(mu_min); print(alpha); print(beta)
  
  # First event
  gk <- 1/gamma(1+(1/k))
  u <- runif(1)
  w <- gk * (-1 * log(u))^(1/k)
  minbound <- w/(mu_max +(alpha*(Ai+1)))
  maxbound <- w/(mu_min)
  iter <- ceiling((log(maxbound-minbound) - log(tol))/log(2))
  y <- WIHPsolve(minbound, maxbound, w, Mu, Ai+1, alpha, beta, pars, iter = iter)
  
  times <- c(mean(y$min,y$max))
  n <- 1
  
  while (times[n] < t){
    u <- runif(1)
    w <- gk * (-1 * log(u))^(1/k)
    minbound <- w/(mu_max +(alpha*(Ai+1)))
    maxbound <- w/(mu_min)
    iter <- ceiling((log(maxbound-minbound) - log(tol))/log(2))
    y <- WIHPsolve(minbound, maxbound, w, Mu, Ai+1, alpha, beta, pars, iter)
    
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


