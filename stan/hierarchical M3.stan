data {
  int<lower=0> N;
  int K;
  int S;
  vector[N] Y;
  vector[N] X;
  array[N] int<lower=1, upper=K> VAX; // Updated array syntax
  array[K] int<lower=1, upper=S> AGE; // Updated array syntax
  array[N] int<lower=1, upper=S> AGE0; // Updated array syntax
}

parameters {
  real<lower=0.0, upper=5.0> a0;
  real<lower=0.0, upper=0.02> b0;
  real<lower=0.0, upper=7.0> c0;
  vector<lower=0.0, upper=5.0>[K] a;
  vector<lower=0.0, upper=0.02>[K] b;
  vector<lower=0.0, upper=7.0>[K] c;
  vector<lower=0.0, upper=5.0>[S] a1;
  vector<lower=0.0, upper=0.02>[S] b1;
  vector<lower=0.0, upper=7.0>[S] c1;
  real<lower=0, upper=5> s_ag;
  real<lower=0, upper=0.025> s_bg;
  real<lower=0, upper=4> s_cg;
  real<lower=0, upper=5> s_a;
  real<lower=0, upper=0.025> s_b;
  real<lower=0, upper=4> s_c;
  vector<lower=0, upper=1.0>[S] s_y;
}

model {
  target += normal_lpdf(a0 | 2.5, 2.5);
  target += normal_lpdf(b0 | 0.005, 0.01);
  target += normal_lpdf(c0 | 2.25, 2.0);
  
  target += student_t_lpdf(s_a | 4, 1.5, 5);
  target += student_t_lpdf(s_b | 4, 0.005, 0.01);
  target += student_t_lpdf(s_c | 4, 2, 2);
  target += student_t_lpdf(s_ag | 4, 1.5, 5);
  target += student_t_lpdf(s_bg | 4, 0.005, 0.01);
  target += student_t_lpdf(s_cg | 4, 2, 2);
  target += student_t_lpdf(s_y | 4, 0.4, 0.05);

  target += normal_lpdf(a1 | a0, s_ag);
  target += normal_lpdf(b1 | b0, s_bg);
  target += normal_lpdf(c1 | c0, s_cg);

  target += normal_lpdf(a | a1[AGE], s_a);
  target += normal_lpdf(b | b1[AGE], s_b);
  target += normal_lpdf(c | c1[AGE], s_c);

  target += normal_lpdf(Y | c[VAX] + a[VAX] .* exp(-b[VAX] .* X), s_y[AGE0]);
}

generated quantities {
  vector[K] intercept = a + c;
}
