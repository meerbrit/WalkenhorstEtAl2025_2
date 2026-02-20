#### Analysis script for DTamide study: DIG (DIGGING) calls ###################
###### Bayesian Multilevel models ##################################
############## BWalkenhorst 2024 ###############################################

#### SETUP ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(ggplot2) 
library(ggpubr) 
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(rstan)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # emtrends
library(extrafont)# use font_import() on first use

set.seed(23)

DIG_data <- readRDS('../data/DIG_data.rds')

#Generic weakly informative prior: normal(0, 1);
priors <- c(set_prior("normal(0,1)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal(base_family='Calibri') +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(1.75)),
      axis.text = element_text(face = "bold", size= rel(1.25)),
      strip.text = element_text(face = "bold", size = rel(1.5), color='white'),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold", size = rel(1.25)),
      legend.text = element_text(face = 'italic', size = rel(1)))
}

#### FUNCTIONS ####
get_age_vars <- function(){
  age_vars <- c((30- mean(DIG_data$REC_AGE))/sd(DIG_data$REC_AGE), 
                (75- mean(DIG_data$REC_AGE))/sd(DIG_data$REC_AGE), 
                (120- mean(DIG_data$REC_AGE))/sd(DIG_data$REC_AGE))
  return(age_vars)
}

################################################################################
######################## Call duration ###########################################
################################################################################
ggplot(DIG_data, aes(x = DIG_avg_Len)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call length", x = "Call length (/s)")

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_DIG_len_TA <- brms::brm(formula = DIG_avg_Len  ~ TREATMENT * AGE_z * SEX + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_len_TA")
# summary(B_DIG_len_A)
# plot(B_DIG_len_A)
# pp_check(B_DIG_len_A, ndraws=100)

B_DIG_len_TA2 <- brms::brm(formula = DIG_avg_Len  ~  TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + (1|LITTER_CODE/ID),
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_len_TA2")
# summary(B_DIG_len_A2)
# plot(B_DIG_len_A2)
# pp_check(B_DIG_len_A2, ndraws=100)

loo(B_DIG_len_TA, B_DIG_len_TA2)

# 2) Define model for DIG DURATION ####
B_DIG_len <- brms::brm(formula = DIG_avg_Len  ~ TREATMENT * AGE_z * SEX + 
                         WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z +  (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(2),
                             file="B_DIG_len")

#### RESULTS: DIG DURATION ####
B_DIG_len <- readRDS("../models/B_DIG_len.rds")

summary(B_DIG_len)


plot(B_DIG_len)
pp_check(B_DIG_len, ndraws = 100)

describe_posterior(
  B_DIG_len,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_len),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_DIG_len, moment_match=T) 

bayes_R2(B_DIG_len)

performance::variance_decomposition(B_DIG_len)


#### EMMs: TAS ####
(emm_DIG_LEN <- emtrends(B_DIG_len, pairwise ~ TREATMENT:SEX, var='AGE_z'))

pd(emm_DIG_LEN) 


p_significance(emm_DIG_LEN, threshold = rope_range(B_DIG_len))


### EMMs at different ages TAS ####
(emm_DIG_LEN <- emmeans(B_DIG_len,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_DIG_LEN)

p_significance(emm_DIG_LEN, threshold = rope_range(B_DIG_len))


(pairs_within_sex <- # Contrasts between treatments within sex at each age
      +     contrast(emm_DIG_LEN, method = "pairwise", by = c("SEX", 'AGE_z')))

pd(pairs_within_sex)

p_significance(pairs_within_sex, threshold = rope_range(B_DIG_len))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
      +     contrast(emm_DIG_LEN, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

p_direction(pairs_within_treatment)


p_significance(pairs_within_treatment, threshold = rope_range(B_DIG_len))


rm(emm_DIG_LEN, pairs_within_sex, pairs_within_treatment)

### PLOTS: DIG LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_len,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_len),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size^2',
                  'Group size',
                  'Competition', 
                  'Body mass offset', 
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_LEN 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(DIG_data$REC_AGE))#25.36677
(mean_age <- mean(DIG_data$REC_AGE))#86.41091

range(DIG_data$REC_AGE)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE))/sd(DIG_data$REC_AGE)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_LEN_pred <- B_DIG_len %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
            re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_LEN_pred$REC_AGE <- DIG_LEN_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_LEN_pred$DIG_len <- DIG_LEN_pred$.epred
DIG_LEN_pred$TREATMENT <- factor(DIG_LEN_pred$TREATMENT, levels = c("DC", "SC", "DT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("DC", "SC", "DT"))

#800*500: DIG_LEN_TAS
ggplot(DIG_LEN_pred, aes(x = REC_AGE, y = DIG_len, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_avg_Len))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 7), alpha = 0.2, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call duration (s) \n", n.breaks = 10) +
  theme_clean()+ 
  facet_wrap(~SEX)

rm(DIG_LEN_pred, B_DIG_len)

################################################################################
######################## Call interval duration ################################
################################################################################
ggplot(DIG_data, aes(x = DIG_avg_Int)) + 
     geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
     labs(title = "Histogram of Call interval duration", x = "Call interval (s)")

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_DIG_int_TA <- brms::brm(formula = DIG_avg_Int  ~ TREATMENT * AGE_z * SEX  + (1|LITTER_CODE/ID),
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_int_TA")
# summary(B_DIG_int_TA)
# plot(B_DIG_int_TA)
# pp_check(B_DIG_int_TA, ndraws=100)

B_DIG_int_TA2 <- brms::brm(formula = DIG_avg_Int  ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX  + (1|LITTER_CODE/ID), 
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_int_TA2")
# summary(B_DIG_int_TA2)
# plot(B_DIG_int_TA2)
# pp_check(B_DIG_int_TA2, ndraws=100)

loo(B_DIG_int_TA, B_DIG_int_TA2)

# 2) Define models for DIG interval duration ####
B_DIG_int <- brms::brm(formula = DIG_avg_Int  ~ TREATMENT * AGE_z * SEX + 
                         WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(2),
                             file="B_DIG_int")

#### RESULTS: DIG interval duration ####
B_DIG_int <- readRDS("../models/B_DIG_int.rds")

summary(B_DIG_int)


plot(B_DIG_int)
pp_check(B_DIG_int, ndraws = 100)

describe_posterior(
  B_DIG_int,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_int),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_DIG_int, moment_match=T) 

bayes_R2(B_DIG_int)


performance::variance_decomposition(B_DIG_int)

### EMMS: DIG interval duration: TAS ####
(emm_DIG_INT <- emtrends(B_DIG_int, pairwise ~ TREATMENT:SEX, var = 'AGE_z'))

pd(emm_DIG_INT) 

p_significance(emm_DIG_INT, threshold = rope_range(B_DIG_int)) 

rm(emm_DIG_INT)

### EMMs at different ages TAS ####
(emm_DIG_INT <- emmeans(B_DIG_int,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_DIG_INT)


p_significance(emm_DIG_INT, threshold = rope_range(B_DIG_int))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
  +    contrast(emm_DIG_INT, method = "pairwise", by = c("SEX", 'AGE_z')))
pd(pairs_within_sex)

p_significance(pairs_within_sex, threshold = rope_range(B_DIG_int))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
    +   contrast(emm_DIG_INT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

p_direction(pairs_within_treatment)

p_significance(pairs_within_treatment, threshold = rope_range(B_DIG_int))

rm(emm_DIG_INT, pairs_within_sex, pairs_within_treatment)

### PLOTS: DIG INTERVAL DURATION ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_int,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_int),
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18,1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size^2',
                  'Group size',
                  'Competition', 
                  'Body mass offset', 
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_INT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('DIG call interval duration')+
  theme_clean()
  #theme(legend.position="none")

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(DIG_data$REC_AGE))
(mean_age <- mean(DIG_data$REC_AGE))
range(DIG_data$REC_AGE)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE))/sd(DIG_data$REC_AGE)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_INT_pred <- B_DIG_int %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
              re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_INT_pred$REC_AGE <- DIG_INT_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_INT_pred$DIG_int <- DIG_INT_pred$.epred
DIG_INT_pred$TREATMENT <- factor(DIG_INT_pred$TREATMENT, levels = c("DC", "SC", "DT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("DC", "SC", "DT"))

#800*500: DIG_INT_TAS
ggplot(DIG_INT_pred, aes(x = REC_AGE, y = DIG_int, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_avg_Int))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 7), alpha = 0.2, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call interval duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(DIG_INT_pred, B_DIG_int)

################################################################################
######################## Call rate #############################################
################################################################################
ggplot(DIG_data, aes(x = DIG_rate)) + 
       geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
       labs(title = "Histogram of Call Rates", x = "Call Rate (calls/s)")


# define prior for DIG rate
priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity
B_DIG_rat_TA <- brms::brm(formula = DIG_rate  ~ TREATMENT * AGE_z * SEX  + (1|LITTER_CODE/ID),
                          data = DIG_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_rat_TA")
# summary(B_DIG_rat_A)
# plot(B_DIG_rat_A)
# pp_check(B_DIG_rat_A, ndraws=100) # slight left shift

B_DIG_rat_TA2 <- brms::brm(formula = DIG_rate  ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + (1|LITTER_CODE/ID), 
                          data = DIG_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_rat_TA2")
# summary(B_DIG_rat_A2)
# plot(B_DIG_rat_A2)
# pp_check(B_DIG_rat_A2, ndraws=100)

loo(B_DIG_rat_TA, B_DIG_rat_TA2)
# elpd_diff se_diff
# B_DIG_rat_TA   0.0       0.0   
# B_DIG_rat_TA2 -4.6       2.8   

# 2) Define models for DIG rate ####
priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))
B_DIG_rat <- brms::brm(formula = DIG_rate  ~ TREATMENT * AGE_z * SEX + 
                         WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z +(1|LITTER_CODE/ID),
                       data = DIG_data, family = gaussian(link='identity'),
                       chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                       save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                       prior = priors, threads = threading(2),
                       file="B_DIG_rat")


#### RESULTS: DIG rate ####
B_DIG_rat <- readRDS("../models/B_DIG_rat.rds")
summary(B_DIG_rat)


plot(B_DIG_rat)
pp_check(B_DIG_rat, ndraws = 100) # shifted to left!

describe_posterior(
  B_DIG_rat,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_rat),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_DIG_rat, moment_match=T) 

bayes_R2(B_DIG_rat)

performance::variance_decomposition(B_DIG_rat)

#### EMMs DIG RAT TAS ####
(emm_DIG_RAT <- emtrends(B_DIG_rat, pairwise ~ TREATMENT:SEX, var = 'AGE_z'))

pd(emm_DIG_RAT) 

p_significance(emm_DIG_RAT, threshold = rope_range(B_DIG_rat))

rm(emm_DIG_RAT)

### EMMs at different ages TAS ####
(emm_DIG_RAT <- emmeans(B_DIG_rat,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_DIG_RAT)

p_significance(emm_DIG_RAT, threshold = rope_range(B_DIG_rat))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
   +   contrast(emm_DIG_RAT, method = "pairwise", by = c("SEX", 'AGE_z')))

pd(pairs_within_sex)

p_significance(pairs_within_sex, threshold = rope_range(B_DIG_rat))


(pairs_within_treatment <- # Contrasts between sex within treatment at each age
   +  contrast(emm_DIG_RAT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

p_direction(pairs_within_treatment)


p_significance(pairs_within_treatment, threshold = rope_range(B_DIG_rat))

rm(emm_DIG_RAT, pairs_within_sex, pairs_within_treatment)

### PLOTS: DIG RATE ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_rat,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_rat),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size^2',
                  'Group size',
                  'Competition', 
                  'Body mass offset', 
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC')  

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_RAT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

rm(posterior_desc)

# Ontogeny 30 - 130 ####
# get all needed values
sd_age <- sd(DIG_data$REC_AGE)#
mean_age <- mean(DIG_data$REC_AGE)#

range(DIG_data$REC_AGE)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE))/sd(DIG_data$REC_AGE)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_RAT_pred <- B_DIG_rat %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
             re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_RAT_pred$REC_AGE <- DIG_RAT_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_RAT_pred$DIG_rate <- DIG_RAT_pred$.epred
DIG_RAT_pred$TREATMENT <- factor(DIG_RAT_pred$TREATMENT, levels = c("DC", "SC", "DT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("DC", "SC", "DT"))

# TAS
ggplot(DIG_RAT_pred, aes(x = REC_AGE, y = DIG_rate, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_rate))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 7), alpha = 0.2, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call rate (calls/s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(DIG_RAT_pred, B_DIG_rat)

# Cleanup ####
rm(DIG_data, priors)
