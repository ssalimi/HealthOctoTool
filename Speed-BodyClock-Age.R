
#####Speed-Body Clock using Gamma distribution in Bayesian model in Stan and brms


stan_code<-list("
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_sigma;  // number of population-level effects
  matrix[N, K_sigma] X_sigma;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int Kc_sigma = K_sigma - 1;
  matrix[N, Kc_sigma] Xc_sigma;  // centered version of X_sigma without an intercept
  vector[Kc_sigma] means_X_sigma;  // column means of X_sigma before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_sigma) {
    means_X_sigma[i - 1] = mean(X_sigma[, i]);
    Xc_sigma[, i - 1] = X_sigma[, i] - means_X_sigma[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Kc_sigma] b_sigma;  // population-level effects
  real Intercept_sigma;  // temporary intercept for centered predictors
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b | 0, 10);
  lprior += normal_lpdf(Intercept | 0, 10);
  lprior += student_t_lpdf(Intercept_sigma | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    // initialize linear predictor term
    vector[N] sigma = Intercept_sigma + Xc_sigma * b_sigma;
    for (n in 1:N) {
      // apply the inverse link function
      sigma[n] = exp(sigma[n]);
    }
    target += normal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_sigma_Intercept = Intercept_sigma - dot_product(means_X_sigma, b_sigma);
}
"
)

fit_Speed_BodyClock <- stan(
  model_code = stan_code,  # Stan code as a string
  data = BLSA_data)


posterior_predict_fit <- stan("fit_Speed_BodyClock", data=blsa_data,
                              fit = fit_Speed_BodyClock, sample_file = "post_pred_SpeedBodyClock.csv")

#######Speed-Body Clock 



prior_a = c(prior(normal(0, 10), class = Intercept),
            prior(normal(0, 10), class = b),
            prior(cauchy(0, 1), class = sds))

fit_Speed-BodyClock<-brm(bf(walkingspeed~bodyclock,sigma~bodyclock),family = gaussian(),data =blsa_data,
                         prior =prior_a,seed=123,init=1,chains=5,iter=20000,
                         control=list(adapt_delta=0.99,max_treedepth=12),
                         cores=future::availableCores())




Speed-BodyClock<-posterior_predict(fit_Speed-BodyClock, newdata=data_blsa, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)


Speed-BodyClock<-as.data.frame(posterior_summary(Speed-BodyClock))

write.csv(Speed-BodyClock,"Speed-BodyClock.csv",row.names=FALSE)


#####Speed-Body Age Stan code 


stan_code<-list("
                
      functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b | 0, 10);
  lprior += normal_lpdf(Intercept | 0, 10);
  lprior += gamma_lpdf(shape | 0.01, 0.01);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = shape * exp(-(mu[n]));
    }
    target += gamma_lpdf(Y | shape, mu);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}          
"
)



fit_Speed_BodyAge <- stan(
  model_code = stan_code,  # Stan code as a string
  data = BLSA_data)


posterior_predict <- stan("fit_Speed_BodyAge", data=blsa_data,
                          fit = fit_Speed_BodyAge, sample_file = "post_pred_SpeedBodyAge.csv")



######Speed-Body Age 
prior_a = c(prior(student_t(3, 1.1, 10), class = Intercept),
            prior(normal(0, 10), class = b))


fit_Speed-BodyAge<-brm(bf(age~(Speed-BodyClock),sigma~Speed-BodyClock),family = gaussian(),data =blsa_data,
                       seed=123,init=1,chains=5,iter=8000,
                       control=list(adapt_delta=0.999,max_treedepth=12))



Speed-BodyAge<-posterior_predict(fit_Speed-BodyAge, newdata=data_blsa, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)


Speed-BodyAge<-as.data.frame(posterior_summary( Speed-BodyAge))

write.csv( Speed-BodyAge," Speed-BodyAge.csv",row.names=FALSE)
