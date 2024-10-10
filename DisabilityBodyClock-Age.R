
###### Development of Disability Index in Stan 


stan_code<-list("functions {
  /* zero-inflated beta-binomial log-PDF of a single response
  * Args:
    *   y: the response value
  *   trials: number of trials of the binomial part
  *   mu: mean parameter of the beta distribution
  *   phi: precision parameter of the beta distribution
  *   zi: zero-inflation probability
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real zero_inflated_beta_binomial_lpmf(int y, int trials,
                                          real mu, real phi, real zi) {
      if (y == 0) {
        return log_sum_exp(bernoulli_lpmf(1 | zi),
                           bernoulli_lpmf(0 | zi) +
                             beta_binomial_lpmf(0 | trials,
                                                mu * phi,
                                                (1 - mu) * phi));
      } else {
        return bernoulli_lpmf(0 | zi) +
          beta_binomial_lpmf(y | trials, mu * phi, (1 - mu) * phi);
      }
    }
  /* zero-inflated beta-binomial log-PDF of a single response
  * logit parameterization of the zero-inflation part
  * Args:
    *   y: the response value
  *   trials: number of trials of the binomial part
  *   mu: mean parameter of the beta distribution
  *   phi: precision parameter of the beta distribution
  *   zi: linear predictor for zero-inflation part
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real zero_inflated_beta_binomial_logit_lpmf(int y, int trials,
                                                real mu, real phi, real zi) {
      if (y == 0) {
        return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                           bernoulli_logit_lpmf(0 | zi) +
                             beta_binomial_lpmf(0 | trials,
                                                mu * phi,
                                                (1 - mu) * phi));
      } else {
        return bernoulli_logit_lpmf(0 | zi) +
          beta_binomial_lpmf(y | trials, mu * phi, (1 - mu) * phi);
      }
    }
  // zero-inflated beta-binomial log-CCDF and log-CDF functions
  real zero_inflated_beta_binomial_lccdf(int y, int trials, real mu, real phi,
                                         real zi) {
    return bernoulli_lpmf(0 | zi) + beta_binomial_lccdf(y | trials, 
                                                        mu * phi,
                                                        (1 - mu) * phi);
  }
  real zero_inflated_beta_binomial_lcdf(int y, int trials, real mu, real phi,
                                        real zi) {
    return log1m_exp(zero_inflated_beta_binomial_lccdf(y | trials, mu, phi, zi));
  }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int trials[N];  // number of trials
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> phi;  // precision parameter
  real<lower=0,upper=1> zi;  // zero-inflation probability
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_1 = (sd_1[1] * (z_1[1]));
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += gamma_lpdf(phi | 0.04, 0.04);
  lprior += beta_lpdf(zi | 2,2);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = inv_logit(mu[n]);
    }
    for (n in 1:N) {
      target += zero_inflated_beta_binomial_lpmf(Y[n] | trials[n], mu[n], phi, zi);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
"
) 


### Disability Index 

prior_DI<-c(prior("student_t(3, 0, 2.5)",class=Intercept),prior("gamma(0.04, 0.04)",class="phi"),prior("student_t(3, 0, 2.5)",class="sd"),prior("beta(2,2)",class="zi"))

fit_DI<- brm(sumevent | trials(DItrial) ~ 1 + (1|id),
             data = blsa_data, family = zero_inflated_beta_binomial( link = "logit", link_phi = "log", link_zi = "logit" ),
             prior= prior_DI,seed=1,init=1,chains=6,iter=8000,
             control=list(adapt_delta=0.99,max_treedepth=11))


save(fit_DI,file="fit_DI.RDATA")



DisabilityIndex<-posterior_predict(fit_DI, newdata=data_blsa, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)


DisabilityIndex<-as.data.frame(posterior_summary(Disability))

write.csv(DisabilityIndex,"DisabilityIndex.csv",row.names=FALSE)


######Stan code Disability-Body Clock 


stan_code<-list("functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int trials[N];  // number of trials
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
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
  real<lower=0> phi;  // precision parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  lprior += normal_lpdf(b | 0, 5);
  lprior += normal_lpdf(Intercept | 0,1);
  lprior += gamma_lpdf(phi | 0.01,0.01);
  lprior += cauchy_lpdf(sd_1 | 0, 5)
    - 2 * cauchy_lccdf(0 | 0, 5);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
    }
    mu = inv_logit(mu);
    for (n in 1:N) {
      target += beta_binomial_lpmf(Y[n] | trials[n], mu[n] * phi, (1 - mu[n]) * phi);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
                
               
 "
)



fit_Disability_BodyClock <- stan(
  model_code = stan_code,  # Stan code as a string
  data = BLSA_data)


posterior_predict_fit <- stan("fit_Disability_BodyClock", data=blsa_data,
                              fit = fit_Disability_BodyClock, sample_file = "post_pred_Disability_BodyClock.csv")



###Disability-Body Clock 
## sumevent=total number of having deficits
###DItrial: Total number of disability measured

fit_Disability-BodyClock <-brm(sumevent | trials(DItrial) ~ bodyclock+(1+yrs|id),
                               prior = c(prior("cauchy(0, 5)", class = "sd"),prior("normal(0, 5)", class = "b"),
                                         prior("gamma(0.01,0.01)", class = "phi"),
                                         prior("normal(0,1)", class = "Intercept"),prior("lkj(1)", class = "cor")),     
                               save_pars=save_pars(all= TRUE),data = blsa_data,family =beta_binomial(),seed=1,init=1,chains=8,iter=10000,
                               control=list(adapt_delta=0.99,max_treedepth=12))




Disability-BodyClock<-posterior_predict(fit_Disability-BodyClock, newdata=data_blsa, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)


Disability-BodyClock<-as.data.frame(posterior_summary(Disability-BodyClock))

write.csv(Disability-BodyClock,"Disability-BodyClock.csv",row.names=FALSE)



#############################Disability-Body Age Stan code 


stan_code<-list("functions {
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
  lprior += normal_lpdf(b |  0, 5);
  lprior += normal_lpdf(Intercept |  0, 10);
  lprior += gamma_lpdf(shape | 0.04, 0.04);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    mu = exp(mu);
    target += gamma_lpdf(Y | shape, shape ./ mu);
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



fit_Disability_BodyAge <- stan(
  model_code = stan_code,  # Stan code as a string
  data = BLSA_data)


posterior_predict_fit <- stan("fit_Disability_BodyAge", data=blsa_data,
                              fit = fit_Disability_BodyAge, sample_file = "post_pred_Disability_BodyAge.csv")



######Disability-Body Age 

prior_DI_BodyAge<-c(prior("normal( 0, 5)",class="b"),prior("normal( 0, 10)",class="Intercept"),prior("gamma(0.04, 0.04)",class="shape"))


fit_Disability-BodyAge <- brm(age ~ (Disability-BodyClock),
                              data = blsa_data,family = Gamma(link = "log"),
                              save_pars=save_pars(all= TRUE)
                              ,seed=1,init=1,chains=8,iter=80000,prior=prior_DI_BodyAge,
                              control=list(adapt_delta=0.99,max_treedepth=12))




Disability-BodyAge<-posterior_predict(fit_Disability-BodyAge, newdata=data_blsa, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)


Disability-BodyAge<-as.data.frame(posterior_summary(Disability-BodyAge))

write.csv(Disability-BodyAge,"Disability-BodyAge.csv",row.names=FALSE)
