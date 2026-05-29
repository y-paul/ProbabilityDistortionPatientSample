library(tidyverse)
library(bayestestR)
library(BayesFactor)
rm(list=ls())


dat <- read.csv("data/fit/pop_level/tk_pop_level.csv")

# calculate differences in posterior distributions
diff_description_lambda <- dat$mu_lambda_de_con - dat$mu_lambda_de_pat
diff_experience_lambda <- dat$mu_lambda_ex_con - dat$mu_lambda_ex_pat
diff_description_gamma <- dat$mu_gamma_de_con - dat$mu_gamma_de_pat
diff_experience_gamma <- dat$mu_gamma_ex_con - dat$mu_gamma_ex_pat
diff_description_theta <- dat$mu_theta_de_con - dat$mu_theta_de_pat
diff_experience_theta <- dat$mu_theta_ex_con - dat$mu_theta_ex_pat

de_gap_pat <- dat$mu_gamma_de_pat - dat$mu_gamma_ex_pat
de_gap_con <- dat$mu_gamma_de_con - dat$mu_gamma_ex_con
diff_de_gap <- de_gap_pat - de_gap_con


dat <- cbind(dat,
      data.frame(
        #diff_description_lambda,
        #diff_experience_lambda,
        diff_description_gamma,
        diff_experience_gamma,
        diff_description_theta,
        diff_experience_theta,
        de_gap_pat,
        de_gap_con,
        diff_de_gap))





# HDI and Probability of Direction

# Simple Parameters
hdi(dat$mu_gamma_ex_con - dat$mu_gamma_ex_pat)
median(dat$mu_gamma_ex_con - dat$mu_gamma_ex_pat)
p_direction(dat$mu_gamma_ex_con - dat$mu_gamma_ex_pat)



hdi(dat$mu_gamma_de_pat - dat$mu_gamma_de_con)
median(dat$mu_gamma_de_pat - dat$mu_gamma_de_con)
p_direction(dat$mu_gamma_de_pat - dat$mu_gamma_de_con)



hdi(dat$mu_lambda_de_pat - dat$mu_lambda_de_con)
median(dat$mu_lambda_de_pat - dat$mu_lambda_de_con)
p_direction(dat$mu_lambda_de_pat - dat$mu_lambda_de_con)



hdi(dat$mu_lambda_ex_pat - dat$mu_lambda_ex_con)
median(dat$mu_lambda_ex_pat - dat$mu_lambda_ex_con)
p_direction(dat$mu_lambda_ex_pat - dat$mu_lambda_ex_con)



hdi(dat$de_gap_pat)
median(dat$de_gap_pat)
p_direction(dat$de_gap_pat)



hdi(dat$de_gap_con)
median(dat$de_gap_con)
p_direction(dat$de_gap_con)






# DE GAP
hdi(de_gap_pat)
p_direction(dat$de_gap_pat)

hdi(de_gap_con)
p_direction(dat$de_gap_con)


hdi(diff_de_gap)
p_direction(dat$diff_de_gap)
median(diff_de_gap)







# behavioral results
b1 <- read.csv("data/processed/patient/behave_dat.csv")
b1$group <- "Patient"
b2 <- read.csv("data/processed/control/overall/behave_dat.csv")[,c(-1)]
b2$group <- "Control"

behave_dat <- rbind(b1, b2)


ttestBF(data=behave_dat, formula=relative_switching~group)
tapply(behave_dat$relative_switching, behave_dat$group, FUN=median)
tapply(behave_dat$relative_switching, behave_dat$group, FUN=quantile)


ttestBF(data=behave_dat, formula=sample_time~group,)
tapply(behave_dat$sample_time, behave_dat$group, FUN=median)
tapply(behave_dat$sample_time, behave_dat$group, FUN=quantile)

ttestBF(data=behave_dat, formula=desc_time~group)
tapply(behave_dat$desc_time, behave_dat$group, FUN=median)
tapply(behave_dat$desc_time, behave_dat$group, FUN=quantile)



ttestBF(data=behave_dat, formula=sample_mean~group)
t.test(formula=sample_mean~group, data=behave_dat )

tapply(behave_dat$sample_mean, behave_dat$group, FUN=median)


# choice strategies


# Optimize EV
ttestBF(data=behave_dat, formula=optimal_ev_choice_desc~group)
tapply(behave_dat$optimal_ev_choice_desc, behave_dat$group, FUN=mean)
tapply(behave_dat$optimal_ev_choice_desc, behave_dat$group, FUN=sd)


ttestBF(data=behave_dat, formula=optimal_ev_choice_samp~group)
tapply(behave_dat$optimal_ev_choice_samp, behave_dat$group, FUN=mean)
tapply(behave_dat$optimal_ev_choice_samp, behave_dat$group, FUN=sd)



# optimize SE

ttestBF(data=behave_dat, formula=minimax_choice_desc~group)
tapply(behave_dat$minimax_choice_desc, behave_dat$group, FUN=median)

ttestBF(data=behave_dat, formula=minimax_choice_samp~group)
tapply(behave_dat$minimax_choice_samp, behave_dat$group, FUN=median)



# optimize Treatment

ttestBF(data=behave_dat, formula=minimax_treatment_desc~group)
tapply(behave_dat$minimax_treatment_desc, behave_dat$group, FUN=median)

ttestBF(data=behave_dat, formula=minimax_treatment_samp~group)
tapply(behave_dat$minimax_treatment_samp, behave_dat$group, FUN=median)





# de gap rule

d1 <- read.csv("data/processed/patient/exp_dat.csv")
d1$group <- "patient"
d2 <- read.csv("data/processed/control/overall/exp_dat.csv")
d2$group <- "control"
d <- rbind(d1, d2[,-c(1)])
d <- d %>% slice(which(row_number() %% 120 == 1)) 

ttestBF(data=d, formula=de_gap_rule~group)
tapply(d$de_gap_rule, d$group, FUN=mean)




