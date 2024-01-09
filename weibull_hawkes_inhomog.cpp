# include <TMB.hpp>
# include <math.h>
# include <vector>
# include <iostream>
/* baseline intensity is mu + B*s(t) + C*sin(t) */
/* where B, C >= 0, B < mu */
/* s(t) alternates between 0 and 1 every pi time units, so that s(t) = when 0<t<pi */
/* and s(t) = 1 when pi<t<(2*pi) etc */
template<class Type>
Type objective_function<Type>::operator() () {
  using namespace Eigen;
  DATA_VECTOR(times); // vector of event times
  DATA_IVECTOR(slots); // slot number (0-indexing) for each event
  DATA_VECTOR(marks);
  DATA_IVECTOR(states);
  DATA_VECTOR(covs); // sin(t)
  DATA_SCALAR(interval); // first x value is at interval, then 2*interval etc
  Type marks_mean = marks.sum()/marks.size(); // Average mark
  // parameters of the Hawkes process
  //PARAMETER(log_mu);
  PARAMETER_VECTOR(baseline); // log mu, log B and logit(C/mu)
  PARAMETER(logit_abratio);
  PARAMETER(log_beta);
  PARAMETER(log_k);

  Type beta = exp(log_beta);
  Type alpha = exp(logit_abratio) / (Type(1.) + exp(logit_abratio)) * (beta/marks_mean); // enforcing 0<=alpha<=beta
  Type k = exp(log_k); // k > 0
  
  Type mu = exp(baseline[0]);
  Type B = exp(baseline[1]);
  Type C = exp(baseline[2]) / (Type(1.) + exp(baseline[2])) * mu;
  
  int N = times.size();
  int n_slots = covs.size();
  
  vector<Type> baselines = vector<Type>::Zero(n_slots);
  for(int j = 0; j< n_slots; ++j){
    if (states[j] > 0){
      baselines[j] = mu + B + (C *covs[j]);
    } else {
      baselines[j] = mu + (C * covs[j]);
    }
  }
  
  vector<Type> baselines_sums = vector<Type>::Zero(n_slots);
  baselines_sums[0] = baselines[0];
  for(int j = 1; j < n_slots; ++j){
    baselines_sums[j] = baselines_sums[(j-1)] + baselines[j];
  }
  
  vector<Type> A = vector<Type>::Zero(times.size());
  for(int i = 1; i < N; ++i){
    // Page 28 of https://pat-laub.github.io/pdfs/honours_thesis.pdf
    A[i] = exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1]);
  }
  
  vector<Type> compensators(N);
  Type marks_sum = 0;
  vector<Type> lambdas(N);
  Type baseline_integral = 0;
  
  for(int i = 0; i < N; ++i){
    marks_sum += marks[i]; // sum of marks up to time i
    int slot = slots[i];
    lambdas[i] = baselines[slot] + alpha * A[i]; // lambda at time i
    if (slot > 0){ // integral of baseline from 0 to time i
      baseline_integral = baselines_sums[(slot-1)] * interval;
      baseline_integral += (times[i] - interval*slot) * baselines[slot];
    } else {
      baseline_integral = times[i] * baselines[0];
    }
    compensators[i] = baseline_integral + (alpha/beta) * (marks_sum - marks[i] - A[i]);
  }
  
  vector<Type> diff_compensators(N);
  diff_compensators[0] = compensators[0];
  for(int i = 1; i < N; ++i){
    diff_compensators[i] = compensators[i] - compensators[(i-1)];
  }
  // likelihood is density of standardized IATs (from a Weibull distribution)
  // multiplied by rates at each event

  Type nll = 0;
  Type lgk = lgamma(1+(1/k)); // Avoid having to recompute
  Type gk = exp(lgk);
  
  nll -= sum(log(lambdas));
  nll -= N * log(k);
  nll -= N * k * lgk;
  nll -= (k - 1) * sum(log(diff_compensators));
  for(int i = 0; i < N; ++i){
    nll += pow((diff_compensators[i] * gk),k);
  }
  


  ADREPORT(mu);
  ADREPORT(B);
  ADREPORT(C);
  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(k);

  return nll;
}
