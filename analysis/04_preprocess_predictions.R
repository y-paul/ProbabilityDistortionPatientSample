rm(list=ls())
library(tidyverse)




# this file is for preprocessing trial by trial data on:
# a) choice rules,
# b) probabilities of acting within these choice rules by individual parameter posteriors


# we are dropping the LIN model here
ALL_FILES <- c("Prelec_parameters.csv",
               "WG_parameters.csv",
               "TK_parameters.csv")

# stays the same, trial by trial data
exp_data <- rbind(read.csv("data/processed/control/overall/exp_dat.csv")[,-c(1)] %>% mutate(group="control"),
                  read.csv("data/processed/patient/exp_dat.csv") %>% mutate(group="patient"))


# value functions, all the same
value_functions <- list(
  POW_gam_LA = function(v, lambda) {
    return(ifelse(v<0, lambda * v, v))
  },
  Prelec_parameters.csv = function(v, lambda) {
    return(ifelse(v<0, lambda * v, v))
  },
  WG_parameters.csv = function(v, lambda) {
    return(ifelse(v<0, lambda * v, v))
  },
  TK_parameters.csv = function(v, lambda) {
    return(ifelse(v<0, lambda * v, v))
  }
)



# probability weighting functions taken from the paper/supplement
probability_weighting_functions <- list(
  POW_gam_LA = function(p, gamma) {
    return(p^gamma)
  },
  Prelec_parameters.csv = function(p, gamma) {
    return(    exp(1) ^ (-  (-  (log(p))    )^gamma   )    )
  },
  WG_parameters.csv = function(p, gamma) {
    return( (p^gamma) / (   (p^gamma) + ((1-p)^gamma)   )   )
  },
  TK_parameters.csv = function(p, gamma) {
    
    return( (p^gamma) / (  (  p^gamma  + ((1-p)^gamma)  )^(gamma^(-1))   )   )
  }
)



# softmax necessary for calculating P(A)
softmax <- function(a,b, theta) {
  return(1 / (1 + exp(-theta*(a-b))))
}


a <- runif(1000, 1, 5)
b <- rep(2, times=1000)
theta = 3


plot(a, softmax(a,b,theta))

# no loop through all models
all_mods <- list()

for (mod in ALL_FILES) {
  
  
  current_model <- mod # remove _parameters.csv
  
  # read in data and join with experiment trial by trial data
  dat <- read.csv(paste0("data/fit/", mod))
  combined_dat <- left_join(exp_data, dat, by=c("pro_id" = "id"))
  
  # we loop through all trials, calculating via PT functions the probability
  # of choosing A
  prob_a_vec <- c()
  for (i in 1:length(combined_dat$participant_id)) {
    
    row <- combined_dat[i,]
    
    # collapsing across conditions to make it more elegant
    if (is.na(row$a.b.f)) {
      lambda = row$lambda_de
      gamma = row$gamma_de
      theta = row$theta_de
      pa_be <- row$pa_be
      pa_se <- row$pa_se
      pb_be <- row$pb_be
      pb_se <- row$pb_se
    } else {
      lambda = row$lambda_ex
      gamma=row$gamma_ex
      theta = row$theta_ex
      pa_be <- row$a.b.f 
      pa_se <- row$a.se.f
      pb_be <- row$b.b.f
      pb_se <- row$b.se.f
    }
    
    
    a_be <- row$a_be
    b_be <- row$b_be
    a_se <- -row$a_se
    b_se <- -row$b_se
    

    
    # calc terms of subjective evaluations
    vab <- value_functions[[current_model]](a_be, lambda)
    vbb <- value_functions[[current_model]](b_be, lambda)
    vase <-value_functions[[current_model]](a_se, lambda)
    vbse <-value_functions[[current_model]](b_se, lambda)
    pab <- probability_weighting_functions[[current_model]](pa_be, gamma)
    pase <- probability_weighting_functions[[current_model]](pa_se, gamma)
    pbb <- probability_weighting_functions[[current_model]](pb_be, gamma)
    pbse <- probability_weighting_functions[[current_model]](pb_se, gamma)
    subj_value_a <- vab * pab + vase * pase
    subj_value_b <- vbb * pbb + vbse * pbse
    
    # calc softmax
    prob_a <- softmax(subj_value_a, subj_value_b, theta)
    
    
    prob_a_vec[length(prob_a_vec)+1] <- prob_a # append to vector
  }
  
  combined_dat$prob_a <- prob_a_vec
  
  # now split files
  desc_data <- combined_dat %>% filter(cond=="desc")
  exp_dat <- combined_dat %>% filter(cond=="exp")
  
  # calculate expected values to various choice rules
  desc_data <- desc_data %>% 
    mutate(ev_a = a_be * pa_be - a_se * pa_se,
           ev_b = b_be * pb_be - b_se * pb_se,
           se_ev_a = a_se * pa_se,
           se_ev_b = b_se * pb_se,
           be_ev_a = a_be * pa_be,
           be_ev_b = b_be * pb_be,
           pase = pa_se,
           pbse = pb_se,
           option_chosen = recode(option_chosen, `0` ="B", `1`="A"))
  exp_dat <- exp_dat %>% 
    mutate(ev_a = a_be * a.b.f - a_se * a.se.f,
           ev_b = b_be * b.b.f - b_se * b.se.f,
           se_ev_a = a_se * a.se.f,
           se_ev_b = b_se * b.se.f,
           be_ev_a = a_be * a.b.f,
           be_ev_b = b_be * b.b.f,
           pase = a.se.f,
           pbse = b.se.f,
           option_chosen = recode(option_chosen, `0` ="B", `1`="A"))
  
  
  # combine back
  combined_dat <- rbind(desc_data, exp_dat)
  
  
  # Calculate wether A or B would be within the choice rule
  combined_dat$optim_ev <- ifelse(combined_dat$ev_a == combined_dat$ev_b,
                                  NA, ifelse(combined_dat$ev_a > combined_dat$ev_b, "A", "B"))
  combined_dat$optim_se_ev <- ifelse(combined_dat$se_ev_a == combined_dat$se_ev_b,
                                     NA, ifelse(combined_dat$se_ev_a < combined_dat$se_ev_b, "A", "B"))
  combined_dat$optim_be_ev <- ifelse(combined_dat$be_ev_a == combined_dat$be_ev_b,
                                     NA, ifelse(combined_dat$be_ev_a > combined_dat$be_ev_b, "A", "B"))
  combined_dat$optim_se_raw <- ifelse(combined_dat$a_se == combined_dat$b_se,
                                      NA, ifelse(combined_dat$a_se < combined_dat$b_se, "A", "B"))
  combined_dat$optim_se_freq <- ifelse(combined_dat$pase == combined_dat$pbse,
                                       NA, ifelse(combined_dat$pase < combined_dat$pbse, "A", "B"))
  
  # check if participant chose within this choice rule or not
  combined_dat$chosen_optim_ev <- ifelse(combined_dat$option_chosen == combined_dat$optim_ev,T,F)
  combined_dat$chosen_optim_se_ev <- ifelse(combined_dat$option_chosen == combined_dat$optim_se_ev,T,F)
  combined_dat$chosen_optim_be_ev <- ifelse(combined_dat$option_chosen == combined_dat$optim_be_ev,T,F)
  combined_dat$chosen_optim_se_raw <- ifelse(combined_dat$option_chosen == combined_dat$optim_se_raw,T,F)
  combined_dat$chosen_option_a <- recode(combined_dat$option_chosen, A=1, B=0)
  combined_dat$chosen_se_freq <- ifelse(combined_dat$option_chosen == combined_dat$optim_se_freq,T,F)
  # now calculate the model predicted probabiliy of acting within this choice rule
  combined_dat$modprob_optim_ev <- ifelse(combined_dat$ev_a == combined_dat$ev_b,
                                          NA,
                                          ifelse(combined_dat$ev_a > combined_dat$ev_b, combined_dat$prob_a, 1-combined_dat$prob_a))
  combined_dat$modprob_optim_se_ev <- ifelse(combined_dat$se_ev_a == combined_dat$se_ev_b,
                                             NA,
                                             ifelse(combined_dat$se_ev_a < combined_dat$se_ev_b, combined_dat$prob_a, 1-combined_dat$prob_a))
  combined_dat$modprob_optim_be_ev <- ifelse(combined_dat$be_ev_a == combined_dat$be_ev_b,
                                             NA,
                                             ifelse(combined_dat$be_ev_a > combined_dat$be_ev_b, combined_dat$prob_a, 1-combined_dat$prob_a))
  combined_dat$modprob_optim_se_raw <- ifelse(combined_dat$a_se == combined_dat$b_se,
                                              NA,
                                              ifelse(combined_dat$a_se < combined_dat$b_se, combined_dat$prob_a, 1-combined_dat$prob_a))
  combined_dat$modprob_chosen_option_a <- combined_dat$prob_a
  combined_dat$modpob_se_freq <- ifelse(combined_dat$pase == combined_dat$pbse,
                                        NA,
                                        ifelse(combined_dat$pase < combined_dat$pbse, combined_dat$prob_a, 1-combined_dat$prob_a))
  # add to list for saving
  all_mods[[mod]] <- combined_dat
}

# dumb everything into the folder "data/prediction_check"
for (mod_name in names(all_mods)) {
  current_mod <- mod_name
  
  current_file <- all_mods[[mod_name]]
  

  write.csv(current_file,
            paste0("data/prediction_check/probs_", current_mod))
}




