# include <TMB.hpp>
# include <math.h>
# include <vector>
# include <iostream>
/* Experimental model, intensity is as Hawkes process but IATs are Weibull (with shape parameter k) */
template<class Type>
Type objective_function<Type>::operator() () {
  using namespace Eigen;
  // vector of time
  DATA_VECTOR(times);
  DATA_VECTOR(marks);
  Type marks_mean = marks.sum()/marks.size(); // Average mark
  // parameters of the Hawkes process
  PARAMETER(log_mu);
  PARAMETER(logit_abratio);
  PARAMETER(log_beta);
  PARAMETER(log_k);
  Type mu = exp(log_mu);
  Type beta = exp(log_beta);
  Type alpha = exp(logit_abratio) / (Type(1.) + exp(logit_abratio)) * (beta/marks_mean); // enforcing 0<=alpha<=beta
  Type k = exp(log_k); // k > 0

  vector<Type> A = vector<Type>::Zero(times.size());
  
  int N = times.size();
  for(int i = 1; i < N; ++i){
    // Page 28 of https://pat-laub.github.io/pdfs/honours_thesis.pdf
    A[i] = exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1]);
  }
  vector<Type> lambdas = mu + alpha * A;
  vector<Type> compensators(N);
  Type marks_sum = 0;
  for(int i = 0; i < N; ++i){
    marks_sum += marks[i];
    compensators[i] = (mu * times[i]) + (alpha/beta) * (marks_sum - marks[i] - A[i]);
  }
  
  vector<Type> diff_compensators(N);
  diff_compensators[0] = compensators[0];
  for(int i = 1; i < N; ++i){
    diff_compensators[i] = compensators[i] - compensators[(i-1)];
  }


  Type nll = 0;
  nll -= sum(log(lambdas));
  nll -= N * log(k);
  nll -= N * k * lgamma(1+(1/k));
  nll -= (k - 1) * sum(log(diff_compensators));
  for(int i = 0; i < N; ++i){
    nll += pow((diff_compensators[i] * exp(lgamma(1+(1/k)))),k);
  }
  

  SIMULATE {
    Type eps = 1e-10, t = 0, M = mu, U;
    int index = 0;
    while (index < times.size()){
      M = mu + alpha * (-beta * (t + eps - times.array().head(index))).exp().sum();
      t += rexp(Type(1.) / M); U = runif(Type(0.), M); // There is currently a bug as at TMB-1.7.20, 14/05/2021.
      if (U <= mu + alpha * (-beta * (t - times.array().head(index))).exp().sum()){
	times[index] = t;
	index++;
      }
    }
    REPORT(times);
  }

  ADREPORT(mu);
  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(k);

  return nll;
}
