# Health Octo Tool (HOT)
## Global Statistical Frame Work
### The BODN, Disability Index, and Health Octo Tool definitions and its components, algorithms, and codes are on the provisional patent. Copy or claiming similar algorithms are prone to copy right and legal responses. 

**Body Organ Disease Number (BODN):** BODN is determined as the number of organs with at least one impairment. Pathology at the specific organ level is established based on predefined impairment criteria for clinical practice disease diagnoses. An organ system contributes to BODN  if it has at least one positive criterion. 
- BODN behaves as a progressive irregular arithmetic sequence and, hence, behaves ordinal in the statistical models. 
- The BODN includes 11 organs, stroke, and cancer. it ranges  from 1 indicating no organ disease to 14 representing 13 organs.

 ### The organs and their diseases and disease levels include:
   
 **1)Cardiovascular System (CV):** Level_1) No CV, Disease_2) Hypertension (3 l3v3ls), Disease_2) Ischemic Heart Disease (3 level), Disease_3) Congestive Heart Failure (3 levels), Disease_4) Arrhythmia (6 levels), Disease_5) 
  Peripheral Artery Disease (2 levels)
  
 **2)Cerebrovascular Accident (CA):** 4 levels
 
 **3)Renal (Ren):** 5 levels 
 
 **4)Metabolic (Met):** Disease_1) Diabetes (6 levels), Diseae_2) Hyperlipidemia (3 levels)

 **5)Gastrointestinal and Liver (GL):** Disease_1) Liver (3 levels) Disease_2)Gastrointestinal disease (3 levels)
 
 **6)Respiratory (Res):** Disease_1) Chronic Obstructive Pulmonary Disease (3 levels), Disease_2) Asthma (2 levels)

 **7)Dsythyroidism (Th):** Disease_1) Hypothyroidism (3 levels), Disease_2)Hyperthyroidism (2 levels)

 **8)Hemaotopoietic System(He):** Disease_1) Anemia (4 levels), Disease_2) Thrombocytopenia (2 levels), Disease_3) White blood cells (3 levels) 

 **9)Oral Health(Pe):** 3 levels 

 **10) Musckeloskeletal System(MS):** Disease_1) Osteoartheritis (OA) (3 levels), Disease_2) Osteoporosis (3 levels), Disease_3) Gout (3 levels)

 **11) MSensory System(Se):** Disease_1) Hearing (4 levels), Disease_2)Eye diseases (5 levels)

 **12) Central Nervous System(CNS):** Disease_1) Depression (3 levels), Disease_2)Parkinson's Disease (3 levels), Disease_3) Cognitive function (3 levels)
 
 **13) Malignant Cancer(Ca):** 2 levels

### 1- Body Clock using the Bayesian model 

BODN~Hypertension+congestiveHeartFailure+IschemicHeartDisease+Arrhythmia+Kidney+Diabetes+Hyperlipidemia+PrepheralArteryDisease+Stroke+Anemia+Thrombocytopnia+GastrointestinalDisease+Liver+COPD+Asthma+OralHealth+Hypothyroisism+Hyperthyroisism+OsteoArtheristis+Osteoporosis+Hearing+Eye+Depression+Parkinsons+Cognition+Cancer+yrs+(1+1|id)

### 2- Body Age using the Bayesian model 

Age ~ Body Clock 

### 3- Bodily Organ Specific Clock (BSC)
The N diseases of specific organs are included in the model so that BODN ~ Diseases1+...+DiseaseN
The codes are similar to the Body Clock code, limited to specific organs and the model's predicted value is an organ-specific Clock. 




### 4- Bodily Organ Specific Age (BSA) the Bayesian model 
Age ~ BSC 


### 5- Speed-Body Clock using the Bayesian model 

Walking Speed (m/S) ~ Body Clock 

### 6- Speed-Body Age 

Age ~ Speed-Body Clock


### Disability Index using the Bayesian model 

Each item of Daily Living Disability, Instrumental Daily Living Disability, and items of cognitive functions measured by MMS  are coded as yes or no. The probability of deficits over the total number of deficits using the Zero Inflated Beta Binomial in the Bayesian framework has been developed. The predicted values of deficits are defined as Disability Index. 


### 7- Disability-Body Clock using the Bayesian model 

 DI ~ Body Clock 

### 8- Disability-Body Age using the Bayesian model 

Age ~ Disability-Body Clock
 
### Description of the codes:
 
**1. Functions Block**
The functions block defines custom functions that are reused in the model. There are two main custom functions here:
  
**cumulative_logit_lpmf(int y, real mu, real disc, vector thres):**
  This function computes the log probability mass function (log-PMF) for a cumulative logit model for a given observation y. It uses the latent mean mu, discrimination parameter disc, and a vector of thresholds thres.
  
**cumulative_logit_merged_lpmf(int y, real mu, real disc, vector thres, int[] j):**
  A wrapper function that computes the cumulative logit log-PMF for ordinal variables by using a subset of the thresholds thres[j[1]:j[2]].

**mo(vector scale, int i):**
  A helper function to compute the monotonic effects of predictors. This function sums over the simplex and returns a value between 0 and 1.

**2. Data Block**
The data block defines the inputs that the model needs to fit the data. It includes:
  
  N: The number of observations (sample size).
Y: The ordinal response variable, with integer values representing categories.
nthres: The number of thresholds for the ordinal logistic model.
K: The number of population-level effects (fixed effects).
X: The population-level design matrix for fixed effects (N by K matrix).
Monotonic Variables (Xmo_1 to Xmo_21): These variables represent the monotonic effects predictors.
disc: The discrimination parameter for the ordinal logistic model.
J_1: A grouping indicator for group-level effects.
prior_only: A flag indicating if only prior information should be used (no likelihood).

**3. Transformed Data Block**
In this block:
  
  The design matrix X is centered to create Xc by subtracting the mean of each column. Centering the data helps with the convergence of the model.
  
**4. Parameters Block**
This block defines the model parameters:
  
 b: A vector of regression coefficients for the population-level (fixed) effects.
Intercept: Ordered thresholds for the ordinal logistic regression model.
bsp: Coefficients for the monotonic predictors.
simo_*: Simplex parameters for the monotonic effects (for each monotonic variable).
sd_1: Group-level standard deviations for the random effects.
z_1: Standardized group-level effects (random effects) for each group.

**5. Transformed Parameters Block**
Here, the actual group-level effects (r_1_1) are computed by scaling the standardized random effects z_1 by the group-level standard deviations sd_1.

**6. Model Block**
This is the core of the model, where the likelihood is defined:
  
  Linear Predictor (mu):
  The model starts by computing the linear predictor mu for each observation using the population-level effects (Xc * b), monotonic effects (mo(simo_*, Xmo_*)), and group-level random effects (r_1_1).

**Priors:**
  
Intercept (thresholds) follows a Student-t prior.
Population-level effects (b) and monotonic effects (bsp) follow a normal prior.
The simplex parameters (simo_*) follow Dirichlet priors.
The group-level effects (z_1) follow a standard normal prior.
The group-level standard deviations (sd_1) also follow a normal prior (with truncation to enforce positivity).
Likelihood:
  The cumulative_logit_lpmf function is used to define the likelihood of the data based on the cumulative logit model.

**7. Generated Quantities Block**
This block computes additional quantities after the model has been fit:
  
  b_Intercept: The thresholds are adjusted by adding the dot product of the means of the predictors and the estimated coefficients (b). This gives the "actual" thresholds for the ordinal categories.

Predicted Response (y_pred):
  The predicted category for each observation is generated using the ordered_logistic_rng function, which samples from the ordinal logistic distribution using the posterior linear predictor (mu) and the computed thresholds (b_Intercept).

**Key Sections:**
  Custom Functions: Define custom cumulative logit functions for ordinal regression and compute monotonic effects.
Data & Parameters: The data block specifies the observed data, while the parameters block defines the unknown quantities (such as regression coefficients) to estimate.
Model Block: The core model that defines how the predictors relate to the response variable through a cumulative logit model. It also specifies the priors for parameters.
Generated Quantities: This section generates posterior predictions of the response variable, allowing you to see what the model would predict for each observation.
**Model Summary:**
  This is a hierarchical ordinal logistic regression model with monotonic effects and 
group-level random effects. It allows for predictors that are treated as monotonic 
variables and estimates both population-level (fixed) effects and group-level (random)
effects. The final generated quantities block generates posterior predictions for 
the ordinal response variable based on the estimated parameters.



**Validity model using fit model parameters in the BLSA to develop Body Clock in InChianti data**


**Data Block:**
  
  This contains the new data you want to predict on, including the design matrix (X_new), monotonic variables (Xmo_*_new), and the group-level indicator (J_1_new).
N_new is the number of new observations.

**Parameters Block:**
  
  Posterior Parameters from the previously fitted model are declared here (b, Intercept, bsp, simo_*, and r_1_1). You are not estimating these parameters again but using them from a previous fit. They will come from the posterior samples loaded into the model.
  
**Generated Quantities Block:**
  
  mu_new: The linear predictor for the new data is computed using the same formula as in the original model.
  
y_new: The new predicted response categories are generated using the ordered_logistic_rng function, which simulates responses from the posterior predictive distribution based on mu_new and the previously estimated thresholds (Intercept).

**Using This Code in Practice:**

  Step 1: You will need to extract the posterior samples from the original model for the parameters b, Intercept, bsp, simo_*, and r_1_1. This can be done using a tool like RStan or PyStan, where you save the posterior samples after the first fit.

Step 2: For the new data, pass in the new design matrix (X_new), new monotonic variables, and any group-level information.

Step 3: Run the Stan model in "generated quantities" mode to generate new predictions using the y_new vector.
