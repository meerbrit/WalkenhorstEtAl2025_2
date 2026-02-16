#### Analysis script for Flutamide study: REP (REPEAT) calls ###################
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
library(extrafont)
# use font_import() on first use

set.seed(23)

REP_data <- readRDS('../data/REP_data.rds')

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
  age_vars <- c((30- mean(REP_data$REC_AGE))/sd(REP_data$REC_AGE), 
                (75- mean(REP_data$REC_AGE))/sd(REP_data$REC_AGE), 
                (120- mean(REP_data$REC_AGE))/sd(REP_data$REC_AGE))
  return(age_vars)
}

################################################################################
######################## Call Duration #########################################
################################################################################

ggplot(REP_data, aes(x = BEG_avg_Len)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call duration", x = "Call duration (/s)")

# 1) Check for non-linearity ####
B_REP_len_TA <- brms::brm(formula = BEG_avg_Len ~ TREATMENT * AGE_z * SEX + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4), 
                          file="B_REP_len_TA")
# summary(B_REP_len_TA)
# plot(B_REP_len_TA)
# pp_check(B_REP_len_TA, ndraws=100)
B_REP_len_TA2 <- brms::brm(formula = BEG_avg_Len ~  TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_len_TA2")
# plot(B_REP_len_TA2)
# pp_check(B_REP_len_TA2, ndraws=100)
loo(B_REP_len_TA, B_REP_len_TA2)


# 2) Define model for REP DURATION ####
B_REP_len <- brms::brm(formula = BEG_avg_Len  ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + 
                         WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),
                       data = REP_data, family = lognormal(link='identity'),
                       chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                       save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                       prior = priors, threads = threading(2),
                       file="B_REP_len")

#### RESULTS: REP DURATION ####
B_REP_len <- readRDS("../models/B_REP_len.rds")

summary(B_REP_len)

plot(B_REP_len)
pp_check(B_REP_len, ndraws = 100)

describe_posterior(
  B_REP_len,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_len),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_REP_len, moment_match=T) 

bayes_R2(B_REP_len)

performance::variance_decomposition(B_REP_len)

### EMMs at different ages TAS ####
(emm_REP_LEN <- emmeans(B_REP_len,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_REP_LEN)

p_significance(emm_REP_LEN, threshold = rope_range(B_REP_len))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
    contrast(emm_REP_LEN, method = "pairwise", by = c("SEX", 'AGE_z')))

pd(pairs_within_sex)


p_significance(pairs_within_sex, threshold = rope_range(B_REP_len))


(pairs_within_treatment <- # Contrasts between sex within treatment at each age
contrast(emm_REP_LEN, method = "pairwise", by = c("TREATMENT", 'AGE_z')))


p_direction(pairs_within_treatment)


p_significance(pairs_within_treatment, threshold = rope_range(B_REP_len))


rm(emm_REP_LEN, pairs_within_sex, pairs_within_treatment)

### PLOTS: REP DURATION ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_len,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_len),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(1,24),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
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
                  'DT:Age^2:Male','SC:Age^2:Male','Age^2:Male', 
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))
# Coeff_REP_LEN 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# Ontogeny 30 - 130 plot ####
# get all needed values
(sd_age <- sd(REP_data$REC_AGE))
(mean_age <- mean(REP_data$REC_AGE))

range(REP_data$REC_AGE)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE))/sd(REP_data$REC_AGE)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_LEN_pred <- B_REP_len %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_LEN_pred$REC_AGE <- REP_LEN_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_LEN_pred$REP_len <- REP_LEN_pred$.epred
REP_LEN_pred$TREATMENT <- factor(REP_LEN_pred$TREATMENT, levels = c("DC", "SC", "DT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("DC", "SC", "DT"))

#800*500: REP_LEN_TAS
ggplot(REP_LEN_pred, aes(x = REC_AGE, y = REP_len, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Len))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(5, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call length (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

# REP LEN at specific ages
REP_LEN_pred <- B_REP_len %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = get_all_ages(),
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
              re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_LEN_pred$REC_AGE <- REP_LEN_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_LEN_pred$REP_len <- REP_LEN_pred$.epred
REP_LEN_pred$TREATMENT <- factor(REP_LEN_pred$TREATMENT, levels = c("DC", "SC", "DT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("DC", "SC", "DT"))

ggplot(REP_LEN_pred, aes(x = REC_AGE, y = REP_len, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Len))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(5, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call length (s)\n") +
  scale_x_continuous(breaks = c(30, 45, 60, 75, 90, 105, 120), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(REP_LEN_pred,  B_REP_len)

################################################################################
######################## Call interval duration ##################################
################################################################################

ggplot(REP_data, aes(x = BEG_avg_Int)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call interval duration", x = "Call interval")

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_REP_int_TA <- brms::brm(formula = BEG_avg_Int  ~ TREATMENT * AGE_z * SEX + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_int_TA")
# summary(B_REP_int_A)
# plot(B_REP_int_A)
# pp_check(B_REP_int_A, ndraws=100)
B_REP_int_TA2 <- brms::brm(formula = BEG_avg_Int  ~  TREATMENT * (AGE_z + I(AGE_z^2))* SEX  + (1|LITTER_CODE/ID), 
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_int_TA2")
# summary(B_REP_int_A2)
# plot(B_REP_int_A2)
# pp_check(B_REP_int_A2, ndraws=100)

loo(B_REP_int_TA, B_REP_int_TA2)

# 2) Define models for REP interval length ####
B_REP_int<- brms::brm(formula = BEG_avg_Int  ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + 
                        WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(2),
                             file="B_REP_int")

#### RESULTS: REP interval duration ####
B_REP_int <- readRDS("../models/B_REP_int.rds")
summary(B_REP_int)

plot(B_REP_int)
pp_check(B_REP_int, ndraws = 100)

describe_posterior(
  B_REP_int,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_int),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_REP_int, moment_match = T) 

bayes_R2(B_REP_int)

performance::variance_decomposition(B_REP_int)

### EMMs at different ages TAS ####
(emm_REP_INT <- emmeans(B_REP_int,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))


pd(emm_REP_INT)

p_significance(emm_REP_INT, threshold = rope_range(B_REP_int))

# Contrasts between treatments within sex at each age
(pairs_within_sex <- contrast(emm_REP_INT, method = "pairwise", by = c("SEX", 'AGE_z')))

pd(pairs_within_sex)

p_significance(pairs_within_sex, threshold = rope_range(B_REP_int))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
    +         contrast(emm_REP_INT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

pd(pairs_within_treatment)


p_significance(pairs_within_treatment, threshold = rope_range(B_REP_int))

rm(emm_REP_INT, pairs_within_sex, pairs_within_treatment)

### PLOTS: REP INTERVAL LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_int,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_int),
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(24, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall','SC:Monthly rainfall','Monthly rainfall',
                  'DT:Group size^2','SC:Group size^2','Group size^2',
                  'DT:Group size','SC:Group size','Group size',
                  'DT:Competition','SC:Competition','Competition', 
                  'DT:Body mass offset','SC:Body mass offset','Body mass offset', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_REP_INT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()#+

# Ontogeny 30 - 130 ####
# get all needed values
sd_age <- sd(REP_data$REC_AGE)#
mean_age <- mean(REP_data$REC_AGE)#

range(REP_data$REC_AGE)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE))/sd(REP_data$REC_AGE)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_INT_pred <- B_REP_int %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 


#unscale AGE_z values:
REP_INT_pred$REC_AGE <- REP_INT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_INT_pred$REP_int <- REP_INT_pred$.epred
REP_INT_pred$TREATMENT <- factor(REP_INT_pred$TREATMENT, levels = c("DC", "SC", "DT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("DC", "SC", "DT"))

#800*500: REP_INT_TAS
ggplot(REP_INT_pred, aes(x = REC_AGE, y = REP_int, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Int))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(5, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call interval duration (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call interval duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~ SEX)

# REP INT at specific ages TAS
REP_INT_pred <- B_REP_int %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = get_all_ages(),
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
              re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_INT_pred$REC_AGE <- REP_INT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_INT_pred$REP_int <- REP_INT_pred$.epred
REP_INT_pred$TREATMENT <- factor(REP_INT_pred$TREATMENT, levels = c("DC", "SC", "DT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("DC", "SC", "DT"))

ggplot(REP_INT_pred, aes(x = REC_AGE, y = REP_int, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Int))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(5, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call interval duration (s)\n") +
  scale_x_continuous(breaks = c(30, 45, 60, 75, 90, 105, 120), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call interval duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~ SEX)

rm(REP_INT_pred, B_REP_int)

################################################################################
######################## Call rate #############################################
################################################################################
ggplot(REP_data, aes(x = BEG_rate)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call rate", x = "Call rate (calls/s)")

priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_REP_rat_TA <- brms::brm(formula = BEG_rate  ~ TREATMENT * AGE_z * SEX  + (1|LITTER_CODE/ID),
                          data = REP_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_rat_TA")
# summary(B_REP_rat_TA)
# plot(B_REP_rat_TA)
# pp_check(B_REP_rat_TA, ndraws=100) # slight left shift

B_REP_rat_TA2 <- brms::brm(formula = BEG_rate  ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX  + (1|LITTER_CODE/ID), 
                          data = REP_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(2),
                          file="B_REP_rat_TA2")
# summary(B_REP_rat_A2)
# plot(B_REP_rat_A2)
# pp_check(B_REP_rat_A2, ndraws=100)

loo(B_REP_rat_TA, B_REP_rat_TA2)

# 2) Define models for REP rate ####
priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))
B_REP_rat <- brms::brm(formula = BEG_rate  ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + 
                         WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),
                       data = REP_data, family = gaussian(link='identity'),
                       chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                       save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                       prior = priors, threads = threading(4),
                       file="B_REP_rat")

#### RESULTS: REP rate ####
B_REP_rat <- readRDS("../models/B_REP_rat.rds")
summary(B_REP_rat)

plot(B_REP_rat)
pp_check(B_REP_rat, ndraws = 100)

describe_posterior(
  B_REP_rat,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_rat),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_REP_rat, moment_match=T) 

bayes_R2(B_REP_rat)

performance::variance_decomposition(B_REP_rat)

### EMMs at different ages TAS ####
(emm_REP_RAT <- emmeans(B_REP_rat,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))


pd(emm_REP_RAT)

p_significance(emm_REP_RAT, threshold = rope_range(B_REP_rat))


(pairs_within_sex <- # Contrasts between treatments within sex at each age
    +         contrasts_treatments_within_sex <- contrast(emm_REP_RAT, method = "pairwise", by = c("SEX", 'AGE_z')))

pd(pairs_within_sex)


p_significance(pairs_within_sex, threshold = rope_range(B_REP_rat))


(pairs_within_treatment <- # Contrasts between sex within treatment at each age
    +          contrast(emm_REP_RAT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))


pd(pairs_within_treatment)

p_significance(pairs_within_treatment, threshold = rope_range(B_REP_rat))

rm(emm_REP_RAT, pairs_within_sex, pairs_within_treatment)

### PLOTS: REP RATE ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_rat,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_rat),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(24, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall','SC:Monthly rainfall','Monthly rainfall',
                  'DT:Group size^2','SC:Group size^2','Group size^2',
                  'DT:Group size','SC:Group size','Group size',
                  'DT:Competition','SC:Competition','Competition', 
                  'DT:Body mass offset','SC:Body mass offset','Body mass offset', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_REP_RAT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()#+

rm(posterior_desc)

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(REP_data$REC_AGE))#
(mean_age <- mean(REP_data$REC_AGE))#

range(REP_data$REC_AGE)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE))/sd(REP_data$REC_AGE)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_RAT_pred <- B_REP_rat %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_RAT_pred$REC_AGE <- REP_RAT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_RAT_pred$REP_rate <- REP_RAT_pred$.epred
REP_RAT_pred$TREATMENT <- factor(REP_RAT_pred$TREATMENT, levels = c("DC", "SC", "DT"))

#800*500: REP_RAT_TAS
ggplot(REP_RAT_pred, aes(x = REC_AGE, y = REP_rate, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_rate))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(5, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call rate (calls/s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call rate (calls/s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

# REP RAT at specific ages TAS
REP_RAT_pred <- B_REP_rat %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = get_all_ages(),
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
              re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_RAT_pred$REC_AGE <- REP_RAT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_RAT_pred$BEG_rate <- REP_RAT_pred$.epred
REP_RAT_pred$TREATMENT <- factor(REP_RAT_pred$TREATMENT, levels = c("DC", "SC", "DT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("DC", "SC", "DT"))

ggplot(REP_RAT_pred, aes(x = REC_AGE, y = BEG_rate, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_rate))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(5, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(5, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call rate (calls/s)\n") +
  scale_x_continuous(breaks = c(30, 45, 60, 75, 90, 105, 120), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call rate (calls/s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~ SEX)


rm(REP_RAT_pred, B_REP_rat)

# Cleanup ####
rm(priors, REP_data)
