
# Load
library(tidyverse)
library(lme4)      # lmer/glmer
library(lmerTest)  # p-values for lmer
library(sjPlot)    # pretty tables/plots for mixed models
library(patchwork) # combine plots

#_____________________________Simulation of data
set.seed(2025)
# design
n_sites <- 10
patients_per_site <- 100
timepoints <- 5   # e.g., baseline + 4 follow-ups

# Total counts
n_patients <- n_sites * patients_per_site;n_patients
n_rows <- n_patients * timepoints;n_rows

# Create site and patient IDs
sites <- tibble(site_id = 1:n_sites,
                site_random_intercept = rnorm(n_sites, 0, 0.4));sites  # site-level SD ~0.4

patients <- expand.grid(site_id = sites$site_id,
                        patient_within_site = 1:patients_per_site) %>%
  as_tibble() %>%
  mutate(patient_id = row_number(),
         age = round(rnorm(n_patients, mean = 45, sd = 12)),
         sex = rbinom(n_patients, 1, 0.55),         # 1 = female
         drug = sample(c("A","B"), n_patients, replace = TRUE, prob = c(0.5,0.5)),
         dose = if_else(drug == "A",
                        round(rnorm(n_patients, 20, 2)),   # mg
                        round(rnorm(n_patients, 50, 5))));patients
#View(patients)

#For instance, the patient in row 1 is the first patient (patient_within_site = 1) from site_id 1, while the patient in row 11 is the second patient (patient_within_site = 2) from that same site_id 1.

#______________________________________________________

# Expand to repeated measures
data_long <- patients %>%
  uncount(timepoints, .id = "timepoint") %>%
  mutate(time = timepoint - 1) ;data_long # 0,1,2,3,4

#uncount() is designed to create a sequence starting from 1. It will always produce 1, 2, 3, 4, 5. It doesn't have an option to start counting from 0.
#mutate() -Use that intermediate column to create the final, analysis-ready time column that starts at 0, following statistical conventions.

#View(data_long)

# Create random effects: site & patient
site_re <- rnorm(n_sites, 0, 0.4)
patient_re <- rnorm(n_patients, 0, 0.6)

#rnorm:The whole mathematical framework for linear mixed-effects models—the type of model you would use to analyze this simulated data—is built on the assumption that the random effects follow a normal distribution.

# Map random intercepts
data_long <- data_long %>%
  mutate(site_re = site_re[site_id],
         patient_re = patient_re[patient_id])
#  View(data_long)

# Fixed effects (true values used for simulation)
beta_0 <- 1.5              # intercept (baseline side-effect severity)
beta_time <- 0.25         # increase per timepoint (maybe side-effects accumulate)
beta_drugB <- 0.6         # drug B increases severity relative to A
beta_dose <- 0.01         # small effect per mg
beta_age <- 0.005         #For every one-year increase in a patient's age, their predicted side-effect score will increase by 0.005 points, assuming all other factors (like drug, dose, etc.) remain constant.
beta_sex <- 0.15          # female slightly higher reported severity [Male = 0 (the baseline group),Female=1]


# Residual SD
sigma_eps <- 1.2

# Simulate continuous side-effect severity (0-10 scale, will clip)
data_long <- data_long %>%
  mutate(
    linear_pred = beta_0 +
      beta_time * time +
      (drug == "B") * beta_drugB +
      dose * beta_dose +
      age * beta_age +
      sex * beta_sex +
      site_re +
      patient_re,
    side_effect_score = linear_pred + rnorm(n(), 0, sigma_eps),
    side_effect_score = pmin(pmax(side_effect_score, 0), 10) # clip 0-10
  );data_long 

#side_effect_score: It then adds one final bit of purely random noise (rnorm(...)) to the "perfect" score. This represents the unexplainable, measurement-to-measurement randomness that always exists.
#pmin(pmax(...)): The side-effect score is supposed to be on a 0-10 scale. # Replaces any score below 0 with 0 and any score above 10 with 10.

#View(data_long)

# Make a binary present/absent variable (for alternative)
data_long <- data_long %>%
  mutate(
    prob_present = plogis((side_effect_score - 3)/2), # map score to probability
    side_effect_present = rbinom(n(), 1, prob_present)
  )

#plogis() Function: This is the logistic function, which is a standard way to convert any number into a probability,is designed so that it always returns 0.5 (or 50%) when its input is exactly 0.
#If a patient's side_effect_score is exactly 3, their probability of the side effect being present is 50%.
#If their score is higher than 3, the probability of a "yes" becomes greater than 50%.
#If their score is lower than 3, the probability of a "yes" is less than 50%.
#The choice of / 2 creates a moderately steep curve. It means the probability changes at a reasonable, smooth rate as the side effect score increases or decreases around 3.
#The simulation is set up to treat a score of 3 as the threshold, and any score higher than that is increasingly likely to be flagged as a "yes" for a side effect being present.

# Quick look
data_long %>% 
  dplyr::select(site_id, patient_id, time, drug, dose, age, sex, side_effect_score, side_effect_present) %>% 
  slice_sample(n = 10)


#_________________________________________________EDA

# overall distribution
p1 <- ggplot(data_long, aes(side_effect_score)) +
  geom_histogram(bins = 30) +
  labs(title = "Side-effect score distribution (0-10)")
p1

# mean by time and drug
p2 <- data_long %>%
  group_by(time, drug) %>%
  summarise(mean_score = mean(side_effect_score), .groups = "drop") %>%
  ggplot(aes(time, mean_score, color = drug)) +
  geom_line() + geom_point() +
  labs(title = "Mean side-effect score by time and drug")
p2

p1 + p2 #both plots in side by side view

#Interpretation:

#Graph 1: Histogram (on the left)
#This plot shows the overall distribution of the side_effect_score across all patients and all time points.
#Most Common Scores: The most frequent side-effect scores are concentrated between approximately 2.0 and 4.0, with the peak of the distribution around a score of 2.5.
#Distribution Shape: The distribution is unimodal (has one peak) and is skewed to the right. This indicates that while most scores are in the lower range, there are a smaller number of instances with higher scores, and very high scores (above 7.5) are rare.

#Graph 2: Line Plot (on the right)
#This plot shows how the average (mean_score) side-effect score changes over time for two different drugs, A and B. There are two main takeaways:
#Effect of Drug: At every time point, the blue line (Drug B) is consistently higher than the red line (Drug A). This means that, on average, patients taking Drug B reported more severe side effects than patients taking Drug A.
#Effect of Time: Both lines trend upwards from left to right. This indicates that for both drugs, the average side-effect score increases over time, from a baseline at time 0 to a peak at time 4.


#___________________________________________Three-level mixed model (continuous outcome)___________________________________________________________________________________________________________________________________________________________________________________

#The main purpose of a longitudinal study (tracking the same individuals over time) is to understand how things change over time. A central and very reasonable scientific question is: "Does everyone change at the same rate?
#time: This is a "within-patient" variable
#age, sex, and drug: These are "between-patient" variables.

# _________________Model:Random intercepts: site and patient

# Basic 3-level random intercepts model
model1 <- lmer(side_effect_score ~ time + drug + dose + age + sex + (1 | site_id) + (1 | patient_id),
               data = data_long, REML = TRUE)
summary(model1)


#model1: Random Intercept Model (1 | patient_id)-This code fits a random intercept for each patient.
#Assumption: It assumes that each patient has their own unique baseline level of side effects (their own starting point, or "intercept"). However, it assumes that the effect of time is the same for everyone.

#model2: Random Intercept and Random Slope Model (1 + time | patient_id)-This code fits both a random intercept and a random slope for each patient.
#Assumption: This is a more flexible and complex model. It assumes that each patient has their own unique baseline AND that the effect of time is also unique to each patient.



# _________________Model:Random slope for time at patient level
model2 <- lmer(side_effect_score ~ time + drug + dose + age + sex + 
                 (1 | site_id) + (1 + time | patient_id),
               data = data_long, REML = TRUE)
#🔸 “boundary (singular) fit”-This warning means the random slope variance for time is estimated as almost zero, so the model thinks that allowing random slopes doesn’t add meaningful improvement.
#This is not an error, just a signal that the model simplified itself — effectively behaving like the random-intercepts model.
summary(model2)

#Both models show that the side-effect score significantly increases with time (Estimate ≈ 0.21, p < .001) and for patients on drugB (Estimate ≈ 0.82, p < .001), while dose has no significant effect (p ≈ 0.77). 
#The more complex model2 produces a "singular fit" warning, indicating it is overfitted because it is trying to estimate a random slope for time that doesn't exist in the data.
#This is confirmed by the near-zero variance for the time slope (7.3e-05) and its perfect correlation (-1.00) with the random intercept in the model2 output. 
#Therefore, the simpler model1 is the correct choice, showing significant random variation between patients (Variance = 0.34) and sites (Variance = 0.17) but a fixed, consistent effect of time for everyone.
#----->   [A perfect -1 or +1 correlation is a mathematical red flag. It is the primary symptom of the boundary (singular) fit warning you received. It means the model is too complex and had to force an artificial, perfect relationship between the starting point (intercept) and the trend (slope) to find a solution.]
#✅ So, the final conclusion is: The data does not actually support the complex idea that every patient has a unique time slope. The simpler model1, which assumes the effect of time is the same for everyone, is the more reliable and appropriate model to use.



#__________________________________________________________________________________________________________________________________________________________________
#dose: This is usually a "between-patient" variable as well, unless the study was designed to have each patient's dose change at different time points.
#     A model that allows each patient to have a unique slope for TIME and DOSE
#-->              model3 <- lmer(side_effect_score ~ time + drug + dose + age + sex + 
#                       (1 | site_id) + (1 + time + dose | patient_id),
#                        data = data_long, REML = TRUE)
#Site-level random slope for drug requires enough sites to estimate the slope variance reliably. With only ~10 sites this estimate will often be unstable or singular. If you have 30+ sites, it’s much safer.
#__________________________________________________________________________________________________________________________________________________________________



#============================================================================
# --- Test for site_id ---

# 1. Refit your original model1 using Maximum Likelihood (ML)
model1_ml <- lmer(side_effect_score ~ time + drug + dose + age + sex + 
                    (1 | site_id) + (1 | patient_id),
                  data = data_long, REML = FALSE) # Note: REML = FALSE

# 2. Fit a simpler model without the site_id random effect
model_no_site <- lmer(side_effect_score ~ time + drug + dose + age + sex + 
                        (1 | patient_id),
                      data = data_long, REML = FALSE)

# 3. Compare the two models
#The ANOVA table is performing a Likelihood Ratio Test. Its purpose is to compare two "nested" models to see if the more complex one is significantly better than the simpler one.
# A small p-value (< 0.05) means model1_ml is significantly better, 
# so the site_id effect is significant.
anova(model_no_site, model1_ml)

# --- Test for patient_id ---

# 1. We already have our full model fit with ML from the step above (model1_ml)

# 2. Fit a simpler model without the patient_id random effect
model_no_patient <- lmer(side_effect_score ~ time + drug + dose + age + sex + 
                           (1 | site_id),
                         data = data_long, REML = FALSE)

# 3. Compare the two models
anova(model_no_patient, model1_ml)
#==============================================================================

#CONCLUSION:Since this p-value is much smaller than 0.05, it means that model1_ml (the model that includes the site_id random effect) is a significantly better fit for the data.


#__________________Extract variance components(For -ICC Calculations)
vc <- as.data.frame(VarCorr(model1))
vc
sigma_site <- vc$vcov[vc$grp == "site_id"]
sigma_patient <- vc$vcov[vc$grp == "patient_id"]
sigma_resid <- attr(VarCorr(model1), "sc")^2


#How ICC Helps: The sample size formula for these types of studies requires you to plug in an estimated ICC value.
#A high estimated ICC means patients within a site are very similar. This tells you that adding more patients to the same site is not very efficient. To get enough statistical power, you will need to recruit from more sites.
#A low estimated ICC means patients within a site are very different. This tells you that the sites themselves don't matter as much, and you might not need as many sites.



icc_site <- sigma_site / (sigma_site + sigma_patient + sigma_resid)
icc_patient <- sigma_patient / (sigma_site + sigma_patient + sigma_resid)
data.frame(icc_site, icc_patient)

#------------------Interpretation

#1. icc_site ≈ 0.092 (or 9.2%)
# This means that approximately 9.2% of the total variation in the side_effect_score is due to systematic differences between the hospital sites.
#In simple terms: Knowing which site a patient attended helps explain about 9.2% of the variability in their side-effect scores.
#This suggests that the sites are somewhat different from each other, but it's not the largest source of variation.

#2. icc_patient ≈ 0.182 (or 18.2%)
# This means that approximately 18.2% of the total variation is due to stable, consistent differences between the patients themselves.
#In simple terms: This effect is larger than the site effect.It indicates that individual patient characteristics (e.g., genetics, general health)
#that make one person consistently different from another are a more significant source of variation in the data than which site they belong to.

#Overall Conclusion
#Combined, about 27.4% (9.2% + 18.2%) of the total variability in the side-effect scores is accounted for by the clustering structure (patients being grouped within sites). 
#The remaining ~73% is the within-patient (or residual) variation, which is the random fluctuation in a single patient's score from one measurement to the next.
#This confirms that accounting for both the site and patient random effects was an important step, as they explain a meaningful portion of the data's structure.


#___________________________________________Three-level mixed model (Binary outcome)___________________________________________________________________________________________________________________________________________________________________________________

model_bin <- glmer(side_effect_present ~ time + drug + dose + age + sex + 
                     (1 | site_id) + (1 | patient_id),
                   data = data_long, family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
#control = glmerControl(...): This is a technical argument to help the model's complex calculations converge to a stable solution, which is sometimes necessary for glmer models.

## This is the model you would fit in that hypothetical scenario

summary(model_bin)
# Convert log-odds to odds ratios
exp(fixef(model_bin))

#Most important thing to notice are the warning messages at the end:
#Model failed to converge... and Model is nearly unidentifiable....
#Interpretation: This means the model had trouble finding a stable, reliable solution. You should be very cautious with these results, as they may not be trustworthy. The model is likely too complex for this binary data.

#Significant Predictors: Based on the p-values (the Pr(>|z|) column), the only factor that is a significant predictor is time (p-value is 3.77e-09, which is very small).
#Non-Significant Predictors: The effects of drugB, dose, age, and sex are not statistically significant. Their p-values are all much larger than the standard 0.05 threshold.

#ODDS RATIO:
#time (OR ≈ 1.13): For each one-unit increase in time (e.g., one week or month), the odds of a patient having a side effect increase by about 13%.
#Other Variables (ORs ≈ 1.0): The odds ratios for drugB, dose, age, and sex are all very close to 1. 
#This is consistent with what we saw above: they have little to no significant effect on the odds of a side effect being present. 
#For example, the OR for drugB is 1.09, meaning it increases the odds by only 9%, which the model found was not a statistically significant difference from Drug A.


# 1. Extract the variance components from your fitted model_bin
vc_bin <- as.data.frame(VarCorr(model_bin))

# Get the variance for site and patient
var_site <- vc_bin$vcov[vc_bin$grp == "site_id"]
var_patient <- vc_bin$vcov[vc_bin$grp == "patient_id"]

# The fixed residual variance for the logit link
var_residual_logit <- (pi^2) / 3


# 2. Calculate the Median Odds Ratio (MOR) for site and patient
mor_site <- exp(sqrt(2 * var_site) * qnorm(0.75))
mor_patient <- exp(sqrt(2 * var_patient) * qnorm(0.75))

# 3. Calculate the Conditional ICC for site and patient
icc_site_bin <- var_site / (var_site + var_patient + var_residual_logit)
icc_patient_bin <- var_patient / (var_site + var_patient + var_residual_logit)


# 4. Display the results
print(paste("MOR for Site:", round(mor_site, 2)))
print(paste("MOR for Patient:", round(mor_patient, 2)))

print(paste("Conditional ICC for Site:", round(icc_site_bin, 3)))
print(paste("Conditional ICC for Patient:", round(icc_patient_bin, 3)))


#Of course. Here is the interpretation of those results in simple terms.


#================================================
#Interpretation of the Median Odds Ratios (MOR):
#An MOR of 1 means there is no variation between clusters.
#MOR for Site: 1.16: This means if you randomly picked two patients from two different sites, the odds of one having a side effect are only 1.16 times (or 16%) higher than the other, just because they are in different hospitals "MOR for Site: 1.16"]. This is a very small effect.
#MOR for Patient: 1.15: Similarly, the random variation between individual patients is also very small "MOR for Patient: 1.15"].

#Interpretation of the Conditional ICC
#The ICC tells you what percentage of the total variation is due to each level.

#Conditional ICC for Site: 0.007: This means that less than 1% (specifically 0.7%) of the variation in whether a side effect is present or absent is due to differences between the sites "Conditional ICC for Site: 0.007"].This is a negligible amount.
#Conditional ICC for Patient: 0.007: Likewise, less than 1% of the variation is due to stable differences between patients "Conditional ICC for Patient: 0.007"].

#Final Conclusion
#While your first model (model1) on the continuous score showed that the magnitude of side effects varied significantly between patients and sites, this analysis shows that the probability of crossing the threshold to have a side effect present is not very dependent on which site or which patient you are. For this binary outcome, almost all the variation is random (residual) noise.


