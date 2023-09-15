# Functions to simulate univariate, multivariate and SESCR HPs
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


# Ogata's Modified Thinning Algorithm for multivariate HP, fixed beta
MOgata <- function(mu, alpha, beta, d, t){
  # check arguments are correct size
  if(length(mu) != d){
    stop("mu must be length d")
  }
  if((nrow(alpha) != d) || (ncol(alpha) != d)){
    stop("alpha must be a dxd matrix")
  }
  # following Chen (2016b)
  # using the definition of A in Laub (2014) and stelfi
  times = numeric(0)
  streams = numeric(0)
  s = 0
  n = 0
  A = numeric(d)
  lambda_bar = sum(mu)
  lambda_s = mu
  last_event = 0 # stream of last event
  while (s < t){
    if (n > 0){
      # lambda_bar is sum of intensities immediately after last accepted event
      lambda_bar = sum(mu)
      for(i in 1:d){
        for(k in 1:d){
          lambda_bar = lambda_bar + alpha[i,k]*A[k]
        }
      }
      # Add 1 to A for the last event
      lambda_bar = lambda_bar + sum(alpha[i,])
    }
    u <- runif(1)
    w <- -1*log(u)/lambda_bar
    s <- s + w
    D <- runif(1)
    if (n > 0){
      # lambda_s is intensity at candidate point s
      lambda_s = mu
      decay = exp(-beta*(s-max(times)))
      for(i in 1:d){
        for(k in 1:d){
          lambda_s[i] = lambda_s[i] + alpha[i,k]*A[k]*decay
        }
        # Add 1 to A for the last event
        lambda_s[i] = lambda_s[i] + alpha[i,last_event]*decay
      }
    }
    if (D*lambda_bar <= sum(lambda_s)){
      k = 1
      while(D*lambda_bar > sum(lambda_s[1:k])){
        k = k + 1
      }
      # re-calculate A
      if (n > 0){
        for(i in 1:d){
          if(i == last_event){
            A[i] <- (A[i] + 1) * decay
          } else {
            A[i] <- A[i] * decay
          }
        }
      }
      n <- n + 1
      last_event <- k
      times <- append(times, s)
      streams <- append(streams, k)
    }
  }
  # Check if last event > t, and return output
  if(max(times)<t && n > 0){
    return(list(times=times, streams=streams))
  } else if (n < 2){
    return(list(times=numeric(0), streams=numeric(0)))
  } else {
    return(list(times=times[1:(n-1)], streams=streams[1:(n-1)]))
  }
}



UWOgata <- function(mu, alpha, beta, k, t){
  # Simulate a Hawkes process with Weibull IAT
  # following Chen (2016a)
  # using the definition of A in Laub (2014) and stelfi
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


