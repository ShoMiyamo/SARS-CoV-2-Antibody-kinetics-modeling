data {
  int<lower=0> N;
  int K;
   vector[N] Y;
    int<lower=1, upper=K>  VAX[N];
    vector[N] X;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real <lower=0, upper=6> a0;
  real <lower=-0.005, upper=0> b0;
  vector<lower=0, upper=7>[K] a;
  vector<lower=-0.01, upper=0>[K]b;
  real<lower=0, upper=5> s_a;
  real<lower=0, upper=0.03> s_b;
  real<lower=0> s_y;
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
   target += normal_lpdf(a0|3,3); // a0 ~ normal (3,3);
   target += normal_lpdf(b0|-0.005, 0.005); //b0 ~ normal (-0.005, 0.005);
   target += student_t_lpdf(s_b|4, 0, 0.03) ; //s_b ~ student_t (4, 0, 0.03);
   target += student_t_lpdf(s_a|4, 1.5, 5) ; //s_a ~student_t(4, 1.5, 5);

  target += normal_lpdf(a|a0,s_a); //a~normal(a0,s_a);
  target += normal_lpdf(b|b0,s_b);//b~normal(b0,s_b);
 
  target += normal_lpdf(Y|a[VAX]+b[VAX].* X,s_y); //Y ~ normal(a[VAX]+b[VAX].* X,s_y);
 
}

//generated quantities{
//vector[N] log_likelihood;
//for (n in 1:N)
//  log_likelihood[n] = normal_lpdf(Y[n]|a[VAX[n]]+b[VAX[n]]* X[n],s_y);

//}
