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


