######Shabnam Salimi. Summer 2021


library(tidyverse)
library(Rcpp)
library(brms)
library(ggplot2)
library(bayesplot)
library(loo)
library(tidybayes)
library(rstan)



data_blsa=read.csv("data.csv",header=TRUE,sep=",")


##mo function that is detaled in STAN syntax in this page, which considers the maximum effect of levels of each disease, known as monotonic effect, 
#for each person (each person is modeled as random effect of the model(1+1|id))

p_all <-get_prior(bodn~mo(Hypertension)+mo(congestiveHeartFailure)+mo(IschemicHeartDisease)+mo(Arrhythmia)+mo(Kidney)+mo(Diabetes)+mo(Hyperlipidemia)+
                    PrepheralArteryDisease+mo(Stroke)+mo(Anemia)+Thrombocytopnia+mo(GastrointestinalDisease)+mo(Liver)+mo(COPD)+Asthma+mo(OralHealth)+mo(Hypothyroisism)+Hyperthyroisism+mo(OsteoArtheristis)+
                    mo(Osteoporosis)+mo(Hearing)+mo(Eye)+mo(Depression)+mo(sParkinsons)+mo(Cognition)+Cancer+yrs+(1+1|id),
                  family=cumulative(),data =data_blsa)


##mo function consideres the levels of each disease. The diseases with levesls have Dirichlet distribution 
##Cauchy distribution is for variance

p_all$prior[1]<-"normal(0,10)" #
p_all$prior[29]<-"student_t(3, 0, 10)" #
p_all$prior[41]<-"cauchy(0, 10)" #
p_all$prior[44]<-"dirichlet(2,2,2,2)" #
p_all$prior[45]<-"dirichlet(2,2)" #
p_all$prior[46]<-"dirichlet(2,2,2,2,2)" #
p_all$prior[47]<-"dirichlet(2,2)" #
p_all$prior[48]<-"dirichlet(2,2)" #
p_all$prior[49]<-"dirichlet(2,2,2,2)" #
p_all$prior[50]<-"dirichlet(2,2)" #
p_all$prior[50]<-"dirichlet(2,2)" #
p_all$prior[51]<-"dirichlet(2,2,2,2,2)" #
p_all$prior[52]<-"dirichlet(2,2)" #
p_all$prior[53]<-"dirichlet(2,2,2)" #
p_all$prior[54]<-"dirichlet(2,2)" #
p_all$prior[55]<-"dirichlet(2,2)" #
p_all$prior[56]<-"dirichlet(2,2)" #
p_all$prior[57]<-"dirichlet(2,2,2)" #
p_all$prior[58]<-"dirichlet(2,2)" #
p_all$prior[59]<-"dirichlet(2,2,2)" #
p_all$prior[60]<-"dirichlet(2,2)" #
p_all$prior[61]<-"dirichlet(2,2)" #
p_all$prior[62]<-"dirichlet(2,2)" #
p_all$prior[63]<-"dirichlet(2,2,2)" #

fit_blsa<-brm(bodn~mo(Hypertension)+mo(congestiveHeartFailure)+mo(IschemicHeartDisease)+mo(Arrhythmia)+mo(Kidney)+mo(Diabetes)+mo(Hyperlipidemia)+
                PrepheralArteryDisease+mo(Stroke)+mo(Anemia)+Thrombocytopnia+mo(GastrointestinalDisease)+mo(Liver)+mo(COPD)+Asthma+mo(OralHealth)+mo(Hypothyroisism)+Hyperthyroisism+mo(OsteoArtheristis)+
                mo(Osteoporosis)+mo(Hearing)+mo(Eye)+mo(Depression)+mo(sParkinsons)+mo(Cognition)+Cancer+yrs+(1+1|id),
              family=cumulative(),data =data_blsa,prior= p_all,seed=1,inits=0,chains=5,iter=5000,control=list(adapt_delta=0.99,max_treedepth=12),
              cores=future::availableCores())

save(fit_blsa,file="fit_blsa.RDATA")


#####Graphically checking the model performance 


pp_check(fit_blsa,newdata=data_blsa,type="bars")


##In data validation using the BLSA model to apply to the same data and re-perform the model 


BodyClock_blsa<-posterior_predict(fit_blsa, newdata=data_blsa, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)


BodyClock_blsa<-as.data.frame(posterior_summary(BodyClock))

write.csv(BodyClock_blsa,"BodyClock_blsa.csv",row.names=FALSE)




#Read Inchianti Data
data_inchianti<-read.csv(data_inch.csv,header=TRUE)



######Using BLSA model to created  Body Clock in InCHIANTI Data (new data validation )
BodyClock_Inchianti<-posterior_predict(fit_blsa,newdata=data_inchianti, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)



BodyClock_Inchianti<-as.data.frame(posterior_summary(BodyClock_Inchianti))

write.csv(BodyClock_Inchianti,"BodyClock_Inchianti.csv",row.names=FALSE)




#####Replication using the same algorithm using InCHIANTI data 



p_all <-get_prior(bodn~mo(Hypertension)+mo(congestiveHeartFailure)+mo(IschemicHeartDisease)+mo(Arrhythmia)+mo(Kidney)+mo(Diabetes)+mo(Hyperlipidemia)+
                    PrepheralArteryDisease+mo(Stroke)+mo(Anemia)+Thrombocytopnia+mo(GastrointestinalDisease)+mo(Liver)+mo(COPD)+Asthma+mo(OralHealth)+mo(Hypothyroisism)+Hyperthyroisism+mo(OsteoArtheristis)+
                    mo(Osteoporosis)+mo(Hearing)+mo(Eye)+mo(Depression)+mo(sParkinsons)+mo(Cognition)+Cancer+yrs+(1+1|id),
                  family=cumulative(),data =data_blsa)


######Setting up the prior for the model
p_all$prior[1]<-"normal(0,10)" #
p_all$prior[29]<-"student_t(3, 0, 10)" #
p_all$prior[41]<-"cauchy(0, 10)" #
p_all$prior[44]<-"dirichlet(2,2,2,2)" #
p_all$prior[45]<-"dirichlet(2,2)" #
p_all$prior[46]<-"dirichlet(2,2,2,2,2)" #
p_all$prior[47]<-"dirichlet(2,2)" #
p_all$prior[48]<-"dirichlet(2,2)" #
p_all$prior[49]<-"dirichlet(2,2,2,2)" #
p_all$prior[50]<-"dirichlet(2,2)" #
p_all$prior[50]<-"dirichlet(2,2)" #
p_all$prior[51]<-"dirichlet(2,2,2,2,2)" #
p_all$prior[52]<-"dirichlet(2,2)" #
p_all$prior[53]<-"dirichlet(2,2,2)" #
p_all$prior[54]<-"dirichlet(2,2)" #
p_all$prior[55]<-"dirichlet(2,2)" #
p_all$prior[56]<-"dirichlet(2,2)" #
p_all$prior[57]<-"dirichlet(2,2,2)" #
p_all$prior[58]<-"dirichlet(2,2)" #
p_all$prior[59]<-"dirichlet(2,2,2)" #
p_all$prior[60]<-"dirichlet(2,2)" #
p_all$prior[61]<-"dirichlet(2,2)" #
p_all$prior[62]<-"dirichlet(2,2)" #
p_all$prior[63]<-"dirichlet(2,2,2)" #

fit_InChianti<-brm(bodn~mo(Hypertension)+mo(congestiveHeartFailure)+mo(IschemicHeartDisease)+mo(Arrhythmia)+mo(Kidney)+mo(Diabetes)+mo(Hyperlipidemia)+
                     PrepheralArteryDisease+mo(Stroke)+mo(Anemia)+Thrombocytopnia+mo(GastrointestinalDisease)+mo(Liver)+mo(COPD)+Asthma+mo(OralHealth)+mo(Hypothyroisism)+Hyperthyroisism+mo(OsteoArtheristis)+
                     mo(Osteoporosis)+mo(Hearing)+mo(Eye)+mo(Depression)+mo(sParkinsons)+mo(Cognition)+Cancer+yrs+(1+1|id),
                   family=cumulative(),data =data_inchianti,
                   prior= p_all,seed=1,inits=0,chains=5,iter=5000,control=list(adapt_delta=0.99,max_treedepth=12),cores=future::availableCores())



BodyClock_Inch<-posterior_predict(fit_InChianti, newdata=data_inchianti, re_formula = NA, allow_new_levels = TRUE, summary=TRUE, scale = "response",draws=3000)



BodyClock_Inch<-as.data.frame(posterior_summary(BodyClock_Inch))

write.csv(BodyClock_Inch,"BodyClock_Inch.csv",row.names=FALSE)




#### For each system, diseases of that specific systems inlcuded in the model and the process is the same as Body Clock. 






##################################Code called in R 
install.packages("rstan")

library(rstan)

####The model to use BODN and add diseases levels of all organs into the model and the predicted values of the model is called Body Clock detailed in STAN 

stan_code <- list("
functions {
  real cumulative_logit_lpmf(int y, real mu, real disc, vector thres) {
    int nthres = num_elements(thres);
    real p;
    if (y == 1) {
      p = inv_logit(disc * (thres[1] - mu));
    } else if (y == nthres + 1) {
      p = 1 - inv_logit(disc * (thres[nthres] - mu));
    } else {
      p = inv_logit(disc * (thres[y] - mu)) -
          inv_logit(disc * (thres[y - 1] - mu));
    }
    return log(p);
  }
  real cumulative_logit_merged_lpmf(int y, real mu, real disc, vector thres, int[] j) {
    return cumulative_logit_lpmf(y | mu, disc, thres[j[1]:j[2]]);
  }
  real mo(vector scale, int i) {
    if (i == 0) {
      return 0;
    } else {
      return rows(scale) * sum(scale[1:i]);
    }
  }
}
data {
  int<lower=1> N;
  int Y[N];
  int<lower=2> nthres;
  int<lower=1> K;
  matrix[N, K] X;
  int<lower=1> Ksp;
  int<lower=1> Imo;
  int<lower=1> Jmo[Imo];
  int Xmo_1[N];
  int Xmo_2[N];
  int Xmo_3[N];
  int Xmo_4[N];
  int Xmo_5[N];
  int Xmo_6[N];
  int Xmo_7[N];
  int Xmo_8[N];
  int Xmo_9[N];
  int Xmo_10[N];
  int Xmo_11[N];
  int Xmo_12[N];
  int Xmo_13[N];
  int Xmo_14[N];
  int Xmo_15[N];
  int Xmo_16[N];
  int Xmo_17[N];
  int Xmo_18[N];
  int Xmo_19[N];
  int Xmo_20[N];
  int Xmo_21[N];
  vector[Jmo[1]] con_simo_1;
  vector[Jmo[2]] con_simo_2;
  vector[Jmo[3]] con_simo_3;
  vector[Jmo[4]] con_simo_4;
  vector[Jmo[5]] con_simo_5;
  vector[Jmo[6]] con_simo_6;
  vector[Jmo[7]] con_simo_7;
  vector[Jmo[8]] con_simo_8;
  vector[Jmo[9]] con_simo_9;
  vector[Jmo[10]] con_simo_10;
  vector[Jmo[11]] con_simo_11;
  vector[Jmo[12]] con_simo_12;
  vector[Jmo[13]] con_simo_13;
  vector[Jmo[14]] con_simo_14;
  vector[Jmo[15]] con_simo_15;
  vector[Jmo[16]] con_simo_16;
  vector[Jmo[17]] con_simo_17;
  vector[Jmo[18]] con_simo_18;
  vector[Jmo[19]] con_simo_19;
  vector[Jmo[20]] con_simo_20;
  vector[Jmo[21]] con_simo_21;
  real<lower=0> disc;
  int<lower=1> N_1;
  int<lower=1> M_1;
  int<lower=1> J_1[N];
  vector[N] Z_1_1;
  int prior_only;
}
transformed data {
  int Kc = K;
  matrix[N, Kc] Xc;
  vector[Kc] means_X;
  for (i in 1:K) {
    means_X[i] = mean(X[, i]);
    Xc[, i] = X[, i] - means_X[i];
  }
}
parameters {
  vector[Kc] b;
  ordered[nthres] Intercept;
  vector[Ksp] bsp;
  simplex[Jmo[1]] simo_1;
  simplex[Jmo[2]] simo_2;
  simplex[Jmo[3]] simo_3;
  simplex[Jmo[4]] simo_4;
  simplex[Jmo[5]] simo_5;
  simplex[Jmo[6]] simo_6;
  simplex[Jmo[7]] simo_7;
  simplex[Jmo[8]] simo_8;
  simplex[Jmo[9]] simo_9;
  simplex[Jmo[10]] simo_10;
  simplex[Jmo[11]] simo_11;
  simplex[Jmo[12]] simo_12;
  simplex[Jmo[13]] simo_13;
  simplex[Jmo[14]] simo_14;
  simplex[Jmo[15]] simo_15;
  simplex[Jmo[16]] simo_16;
  simplex[Jmo[17]] simo_17;
  simplex[Jmo[18]] simo_18;
  simplex[Jmo[19]] simo_19;
  simplex[Jmo[20]] simo_20;
  simplex[Jmo[21]] simo_21;
  vector<lower=0>[M_1] sd_1;
  vector[N_1] z_1[M_1];
}
transformed parameters {
  vector[N_1] r_1_1;
  r_1_1 = (sd_1[1] * (z_1[1]));
}
model {
  vector[N] mu = Xc * b;
  for (n in 1:N) {
    mu[n] += (bsp[1]) * mo(simo_1, Xmo_1[n]) + (bsp[2]) * mo(simo_2, Xmo_2[n]) +
             (bsp[3]) * mo(simo_3, Xmo_3[n]) + r_1_1[J_1[n]];
  }
  target += student_t_lpdf(Intercept | 3, 0, 2.5);
  target += normal_lpdf(b | 0, 1);
  target += normal_lpdf(bsp | 0, 1);
  target += normal_lpdf(sd_1 | 0, 1);
  if (!prior_only) {
    for (n in 1:N) {
      target += cumulative_logit_lpmf(Y[n] | mu[n], disc, Intercept);
    }
  }
}
generated quantities {
  vector[nthres] b_Intercept = Intercept + dot_product(means_X, b);
  int y_pred[N];
  for (n in 1:N) {
    y_pred[n] = ordered_logistic_rng(mu[n], b_Intercept);
  }
}
"
)


#######example of the data 


# Run the Stan model
fit <- stan(
  model_code = stan_code,  # Stan code as a string
  data = BLSA_data)


########################################Validity experiment 


stan_code_validity<-"

data {
  int<lower=1> N_new;          // number of new observations
  matrix[N_new, K] X_new;      // new population-level design matrix
  int Xmo_1_new[N_new];        // new monotonic variables
  int Xmo_2_new[N_new];
  // ... (more monotonic variables as needed)
  int<lower=1> J_1_new[N_new]; // new grouping indicator
  vector[N_new] Z_1_1_new;     // new group-level predictor values
}

parameters {
  // These parameters are already estimated from the previous model,
  // but we include them here to be able to sample posterior values for prediction
  vector[Kc] b;               // population-level effects (posterior samples)
  ordered[nthres] Intercept;   // thresholds (posterior samples)
  vector[Ksp] bsp;            // special effects coefficients (posterior samples)
  // simplexes of monotonic effects
  simplex[Jmo[1]] simo_1;     // posterior simplexes
  // ... (more simplexes for monotonic effects as needed)
  vector[N_1] r_1_1;          // random effects (posterior samples)
}

generated quantities {
  vector[N_new] mu_new;        // new linear predictor
  int<lower=1, upper=nthres + 1> y_new[N_new];  // predicted response for new data
  
  // Compute linear predictor for new data
  for (n in 1:N_new) {
    mu_new[n] = dot_product(X_new[n], b);
    mu_new[n] += bsp[1] * mo(simo_1, Xmo_1_new[n]);
    mu_new[n] += bsp[2] * mo(simo_2, Xmo_2_new[n]);
    // ... (more terms for monotonic effects as needed)
    mu_new[n] += r_1_1[J_1_new[n]] * Z_1_1_new[n];
  }
  
  // Generate predictions for y using the ordered logistic distribution
  for (n in 1:N_new) {
    y_new[n] = ordered_logistic_rng(mu_new[n], Intercept);
  }
}
"


# After fitting your original model, extract the posterior samples

fit_blsa <- sampling(fit, data=BLSA_data)



# Generate posterior predictions using the posterior samples. 
# Posterior_predict_fit is Body Clcok in InCHIANTI using BLSA model parameters and stored 
#in post_pred_InChianti.csv 
posterior_predict_fit <- stan("fit_blsa", data=InCHIANTI_data,
                              fit = fit_blsa, sample_file = "post_pred_InChianti.csv")