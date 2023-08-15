data {
  int N; //Number of galaxies
  vector[N] x; //x values (absolute magnitude)
  vector[N] y; //y values (colour)
}

parameters {

  vector[2] slopes;
  ordered[2] intercepts;
  vector<lower=0>[2] scatter; 

  real<lower=0, upper=1> mixture_prob;

  real<lower=1> nu;
}

model{

  slopes[1] ~ normal(0.07, 0.1);
  slopes[2] ~ normal(0.00, 0.05);
  intercepts[1] ~ normal(0.6, 0.05);
  intercepts[2] ~ normal(0.8, 0.05);
  
  scatter[1] ~ normal(0.1, 0.2);
  scatter[2] ~ normal(0, 0.1);

  mixture_prob ~ beta(2, 2);

  nu ~ gamma(2,0.1);

  for (i in 1:N){
  target += log_mix(mixture_prob, 
                    student_t_lpdf(y[i] | nu, slopes[1] * x[i] + intercepts[1], scatter[1]),
                    student_t_lpdf(y[i] | nu, slopes[2] * x[i] + intercepts[2], scatter[2])
                    );

  }

}

generated quantities{

  vector[N] y_red;
  vector[N] y_blue;
  vector[N] log_Pr;

  for (i in 1:N){

    real lp1;
    real lp2;

    y_red[i] = normal_rng(slopes[1] * x[i] + intercepts[1], scatter[1]);
    y_blue[i] = normal_rng(slopes[2] * x[i] + intercepts[2], scatter[2]);
  
    // From https://discourse.mc-stan.org/t/mixture-models-how-to-get-p-data-cluster/662
    lp1 = log(mixture_prob) + student_t_lpdf(y[i] | nu, slopes[1] * x[i] + intercepts[1], scatter[1]);
    lp2 = log1m(mixture_prob) + student_t_lpdf(y[i] | nu, slopes[2] * x[i] + intercepts[2], scatter[2]);
    log_Pr[i] = lp1 - log_sum_exp(lp1, lp2);
  }
}