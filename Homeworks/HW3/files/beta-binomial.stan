data{
  int N;
  int y;
  real<lower=0> alpha;
  real<lower=0> beta;
}
parameters{
  real theta;
}
model{
  target+=binomial_lpmf(y|N, theta); 
  target+=beta_lpdf(theta|alpha, beta);
}
