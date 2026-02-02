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
  real <lower=1.0, upper=6.0> a0;
  real <lower=0.0,upper=0.015> b0;
  vector<lower=1.0, upper=6.0>[K]a;
  vector<lower=0, upper=0.01666667>[K]b;
  real<lower=0, upper=5> s_a;  //real<lower=0, upper=5> s_b;
  real<lower=0, upper=0.025> s_b;
  real<lower=0,upper=1.0> s_y;
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.


model {
    target += normal_lpdf(a0|3,3); //a0 ~normal (2.25,1.2);
    target += normal_lpdf(b0|0.005,0.005); // b0 ~ normal (3.25,2.25);
   
    target += student_t_lpdf(s_a|4, 1.5, 5);//s_a ~student_t(4, 1.5, 1.5);
    target += student_t_lpdf(s_b|4,0.005,0.01);//s_b ~ student_t (4, 2.5, 3.5);
    
    target += student_t_lpdf(s_y|4,0.4,0.05);//s_y ~student_t(4, 0.4, 0.05);

  target += normal_lpdf(a|a0,s_a);//a~normal(a0,s_a);
  target += normal_lpdf(b|b0,s_b);//b~normal(b0,s_b);

  
  target +=normal_lpdf(Y|a[VAX].*exp(-1*b[VAX].*X), s_y);
  //Y ~ normal(a[VAX] + b[VAX] ./ (1 + exp(c[VAX].*X)), s_y);
 
 }

//generated quantities{
//vector[K] intercept;
//vector[K] st;
 //intercept= a+b/2;
 //st=log(b*10^(4)-1) ./c;
 
//}