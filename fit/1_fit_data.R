library(rstan)
library(loo)
library(bayesplot)
library(tidyverse)

rm(list=ls())


# helper functions --------------------------------------------------------

# convergence diagnostics
source('fit/stan_estimation/functions/00_stan_diags.R')

# loo-accuracy
source('fit/stan_estimation/functions/binary_accuracy_loo.R')


binary_accuracy_from_ll <- function(log_lik_mat, y, cutoff = 0.5) {
  S <- nrow(log_lik_mat)
  n <- ncol(log_lik_mat)
  

  y_mat <- matrix(rep(y, each = S), nrow = S)  # S x n
  lik   <- exp(log_lik_mat)
  p_mat <- ifelse(y_mat == 1, lik, 1 - lik)   # S x n
  p_hat <- colMeans(p_mat)
  y_pred <- as.integer(p_hat >= cutoff)
  
  TP <- sum(y == 1 & y_pred == 1)
  TN <- sum(y == 0 & y_pred == 0)
  FP <- sum(y == 0 & y_pred == 1)
  FN <- sum(y == 1 & y_pred == 0)
  
  acc <- (TP + TN) / (TP + TN + FP + FN)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  
  ba <- mean(c(sens, spec), na.rm = TRUE)
  
  list(
    acc = acc,
    ba  = ba,
    sens = sens,
    spec = spec,
    p_hat = p_hat,
    y_pred = y_pred
  )
}


# data set with  ID, can be used to reconstruct stan data
data_patient <- read.csv("data/processed/patient/exp_dat.csv")
data_control <- read.csv("data/processed/control/exp_dat.csv")


# for convenience - can be reconstructed from dat
d_stan_patient <- readRDS("data/processed/patient/d_stan.rds")
d_stan_control <- readRDS("data/processed/control/d_stan.rds")


names(d_stan_patient) <- paste0(names(d_stan_patient), "_pat")
names(d_stan_control) <- paste0(names(d_stan_control), "_con")


d_stan <- c(d_stan_patient, d_stan_control)


# data --------------------------------------------------------------------



d_stan[c("xa_pat",
         "xb_pat",
         "pa_pat",
         "pb_pat",
         "fa_pat",
         "fb_pat",
         "choices_pat",
         "xa_con",
         "xb_con",
         "pa_con",
         "pb_con",
         "fa_con",
         "fb_con",
         "choices_con")] = lapply(d_stan[c("xa_pat",
                                           "xb_pat",
                                           "pa_pat",
                                           "pb_pat",
                                           "fa_pat",
                                           "fb_pat",
                                           "choices_pat",
                                           "xa_con",
                                           "xb_con",
                                           "pa_con",
                                           "pb_con",
                                           "fa_con",
                                           "fb_con",
                                           "choices_con")], function(x) { 
                                             
                                             # r = as.matrix(x[sub,]); return(r) 
                                             r = as.matrix(x); return(r) 
                                             
                                           } )


# stan set up -------------------------------------------------------------

# sampler parameters
n_chains = 8
n_iter = 5e3
n_burn = 1e3
n_thin = 2
n_cores = n_chains

## ----- group-level: theta -----
gp_t <- c(
  # patient means (0..5 scale, in generated quantities)
  "mu_theta_de_pat", "mu_theta_ex_pat",
  # control means
  "mu_theta_de_con", "mu_theta_ex_con",
  # SDs (per group & condition)
  "sig_theta_ex_pat", "sig_theta_de_pat",
  "sig_theta_ex_con", "sig_theta_de_con",
  # correlations between conditions (per group)
  "theta_r_pat", "theta_r_con"
)

## ----- group-level: gamma -----
gp_g1 <- c(
  "mu_gamma_de_pat", "mu_gamma_ex_pat",
  "mu_gamma_de_con", "mu_gamma_ex_con",
  "sig_gamma_ex_pat", "sig_gamma_de_pat",
  "sig_gamma_ex_con", "sig_gamma_de_con",
  "gam_r_pat", "gam_r_con"
)

## ----- group-level: lambda -----
gp_l1 <- c(
  "mu_lambda_de_pat", "mu_lambda_ex_pat",
  "mu_lambda_de_con", "mu_lambda_ex_con",
  "sig_lambda_ex_pat", "sig_lambda_de_pat",
  "sig_lambda_ex_con", "sig_lambda_de_con",
  "lam_r_pat", "lam_r_con"
)


## ----- individual-level parameters -----
ip_t <- c(
  "theta_de_pat", "theta_ex_pat",
  "theta_de_con", "theta_ex_con"
)

ip_g1 <- c(
  "gamma_de_pat", "gamma_ex_pat",
  "gamma_de_con", "gamma_ex_con"
)

ip_l1 <- c(
  "lambda_de_pat", "lambda_ex_pat",
  "lambda_de_con", "lambda_ex_con"
)

## ----- log-likelihoods (for LOO etc.) -----
l_pars <- c(
  "log_lik_de_pat", "log_lik_ex_pat",
  "log_lik_de_con", "log_lik_ex_con"
)

gap_pars <- c(
  "de_gap_gamma_pat", "de_gap_gamma_con", "group_diff_gamma_gap",
  "de_gap_lambda_pat", "de_gap_lambda_con", "group_diff_lambda_gap",
  "de_gap_theta_pat", "de_gap_theta_con", "group_diff_theta_gap"
)

gp_s  <- c(gp_t, gp_g1, gp_l1)      
ip_s  <- c(ip_t, ip_g1, ip_l1)
all_pars <- c(gp_s, ip_s, l_pars)

# stan models -------------------------------------------------------------

stan_mods = list.files('fit/stan_estimation/stan_mods_groups/no_delta',
                       full.names = T)[1:4]
names(stan_mods) = gsub('^0[^_]*_|.stan', '', list.files('fit/stan_estimation/stan_mods_groups/no_delta')[1:4])
# estimations -------------------------------------------------------------

# list for storing the results
est_mods = list()


# loop over the models
for(i in names(stan_mods)) {
  # set up parameters for monitoring
  if (i == "lin") {
    gp_s <- c(gp_t, gp_l1)
    ip_s <- c(ip_t, ip_l1)
    gap_pars_i <- gap_pars[!grepl("gamma", gap_pars)]  # Exclude gamma for lin
  } else {
    gp_s <- c(gp_t, gp_g1, gp_l1)
    ip_s <- c(ip_t, ip_g1, ip_l1)

    gap_pars_i <- gap_pars
  }
  
  trans = stanc(file = stan_mods[[i]])
  compiled = stan_model(stanc_ret = trans, verbose = F)
  
  stanfit = sampling(object = compiled,
                      data = d_stan,
                      pars = c(gp_s, ip_s, l_pars, gap_pars_i),  
                      init = '0',
                      chains = n_chains,
                      iter = n_iter,
                      warmup = n_burn,
                      thin = n_thin,
                      cores = n_chains,
                      control = NULL)
  
  
  
  
  
  # model performance -------------------------------------------------------
  
  perf_patient = lapply(list(de = 1, ex = 2), function(ii) {
    
    ll = ifelse(ii == 1, 'log_lik_de_pat', 'log_lik_ex_pat')
    
    # approximate loo 
    looE = rstan::loo(stanfit,
                      pars = ll,
                      moment_match = F)
    
    # loo balanced accuracy
    ba = binary_accuracy_loo(stanfit,
                             parameter_name = ll,
                             d_stan$choices_pat[,ii],
                             binary_cutoff = .5)
    
    # loo ind accuracy
    ind_ba_loo = ID_binary_accuracy_loo(stanfit,
                                        parameter_name = ll,
                                        y = d_stan$choices_pat[,ii],
                                        N = d_stan$N_pat, # number of participants 
                                        ncp = 60 # number of choices per participant & condition
    )
    
    return(list(looE = looE, ba = ba, ind_perf = ind_ba_loo, ll=ll))
    
  })
  
  perf_control = lapply(list(de = 1, ex = 2), function(ii) {
    
    ll = ifelse(ii == 1, 'log_lik_de_con', 'log_lik_ex_con')
    
    # approximate loo 
    looE = rstan::loo(stanfit,
                      pars = ll,
                      moment_match = F)
    
    # loo balanced accuracy
    ba = binary_accuracy_loo(stanfit,
                             parameter_name = ll,
                             d_stan$choices_con[,ii],
                             binary_cutoff = .5)
    
    # loo ind accuracy
    ind_ba_loo = ID_binary_accuracy_loo(stanfit,
                                        parameter_name = ll,
                                        y = d_stan$choices_con[,ii],
                                        N = d_stan$N_con, # number of participants 
                                        ncp = 60 # number of choices per participant & condition
    )
    
    return(list(looE = looE, ba = ba, ind_perf = ind_ba_loo, ll=ll))
    
  })
  
  
  
  
  # ----- overall performance (patients + controls together) -----
  
  # extract log-likelihood matrices
  ll_de_pat <- rstan::extract(stanfit, "log_lik_de_pat")[[1]]  # S x n_pat
  ll_ex_pat <- rstan::extract(stanfit, "log_lik_ex_pat")[[1]]  # S x n_pat
  ll_de_con <- rstan::extract(stanfit, "log_lik_de_con")[[1]]  # S x n_con
  ll_ex_con <- rstan::extract(stanfit, "log_lik_ex_con")[[1]]  # S x n_con
  
  # stack trials: columns = all trials from both groups
  ll_de_all <- cbind(ll_de_pat, ll_de_con)  
  ll_ex_all <- cbind(ll_ex_pat, ll_ex_con)
  
  # LOO on combined data
  reff_de_all <- loo::relative_eff(exp(ll_de_all),
                                   chain_id = rep(1:n_chains,
                                                  each = nrow(ll_de_all) / n_chains))
  reff_ex_all <- loo::relative_eff(exp(ll_ex_all),
                                   chain_id = rep(1:n_chains,
                                                  each = nrow(ll_ex_all) / n_chains))
  
  loo_de_all <- loo::loo(ll_de_all, r_eff = reff_de_all)
  loo_ex_all <- loo::loo(ll_ex_all, r_eff = reff_ex_all)
  
  
  
  
  y_de_all <- c(d_stan$choices_pat[,1], d_stan$choices_con[,1])
  y_ex_all <- c(d_stan$choices_pat[,2], d_stan$choices_con[,2])
  
  
  
  # This is a chatgpt function and I am just not sure its corrrect
  ba_de_all <- binary_accuracy_from_ll(ll_de_all, y_de_all, .5)
  ba_ex_all <- binary_accuracy_from_ll(ll_ex_all, y_ex_all, .5)
  
  
  perf_overall <- list(
    de = list(looE = loo_de_all, ba = ba_de_all, ll = ll_de_all),
    ex = list(looE = loo_ex_all, ba = ba_ex_all, ll = ll_ex_all)
  )
  
  
  
  
  
  
  
  # parameters --------------------------------------------------------------
  
  # print diagnostic plots
  
  try(stan_diag(stanFit = stanfit,
                n = d_stan$n_pat,
                ind_p = c("theta_de_pat", "theta_ex_pat", "lambda_de_pat", "lambda_ex_pat", "gamma_de_pat", "gamma_ex_pat"),
                group_p = c("mu_theta_de_pat", "mu_theta_ex_pat", "mu_lambda_de_pat", "mu_lambda_ex_pat", "mu_gamma_de_pat", "mu_gamma_ex_pat",
                            "theta_r_pat", "gam_r_pat", "lam_r_pat"), 
                write_path = paste0('fit/stan_estimation/stan_diags/patient/', i)
  ))
  
  # print diagnostic plots
  try(stan_diag(stanFit = stanfit,
                n = d_stan$n_con,
                ind_p = c("theta_de_con", "theta_ex_con", "lambda_de_con", "lambda_ex_con", "gamma_de_con", "gamma_ex_con"),
                group_p = c("mu_theta_de_con", "mu_theta_ex_con", "mu_lambda_de_con", "mu_lambda_ex_con", "mu_gamma_de_con", "mu_gamma_ex_con",
                            "theta_r_con", "gam_r_con", "lam_r_con"), 
                write_path = paste0('fit/stan_estimation/stan_diags/control/', i)
  ))
  
  # posterior samples
  p_pars = rstan::extract(stanfit)[c(gp_s)]
  gap_samples <- rstan::extract(stanfit)[gap_pars_i]
  

  
  i_pars = rstan::extract(stanfit)[c(ip_s)]
  
  # summary table
  fit_summary = summary(stanfit,
                        pars = c(gp_s, 'lp__'))[[1]]
  fit_summary = round(fit_summary, 3)
  
  # output ------------------------------------------------------------------
  
  # sampling info
  sampling_info = list(
    data = d_stan,
    pars = c(gp_s, ip_s),
    chains = n_chains,
    iter = n_iter,
    warmup = n_burn,
    thin = n_thin,
    control = NULL,
    model = compiled
  )
  
  
  # list with results
  est_mods[[i]] = list(pars = list(p_pars = p_pars, # posterior samples of pop-lvl pars
                                   i_pars = i_pars), # posterior samples of ind-lvl pars
                       fit_summary = fit_summary, # quick summary of pop-lvl pars
                       performance_pat = perf_patient,
                       performance_con = perf_control,# predictive performance metrics
                       sampling_info = sampling_info, # estimation info, including data,
                       performance_de_pat = perf_patient$de$ind_perf$ba,
                       performance_ex_pat = perf_patient$ex$ind_perf$ba,
                       performance_de_con = perf_control$de$ind_perf$ba,
                       performance_ex_con = perf_control$ex$ind_perf$ba,
                       perf_overall,
                       gap_samples = gap_samples
                       
  )
  

  
  print(i)
  
}

saveRDS(est_mods, file = 'data/fit/est_mods_fin.rds')
est_mods <- readRDS('data/fit/est_mods_fin.rds')

for (i in names(est_mods)) {
  stan_fit <- est_mods[[i]]
  
  
  
  
  performance_de_pat <- stan_fit$performance_pat$de$ba$ba
  performance_ex_pat <- stan_fit$performance_pat$ex$ba$ba
  performance_de_con <- stan_fit$performance_con$de$ba$ba
  performance_ex_con <- stan_fit$performance_con$ex$ba$ba
  stan_dat <- stan_fit$pars$i_pars

  tryCatch({
    pop_level <- stan_fit$pars$p_pars
    pop_level <- as.data.frame(pop_level)
    
    write.csv(pop_level, str_c("data/fit/pop_level/", i, "_pop_level.csv"))
  },
  error = function(cond) {
    print(i)
    print("Error saving pop-level")
    print(cond)
  },
  warning=function(cond){print(cond)})
  
  
  if (i!="lin") {
    stan_dat_pat <- stan_dat[c("theta_de_pat", "theta_ex_pat",
                               "gamma_de_pat", "gamma_ex_pat",
                               "lambda_de_pat", "lambda_ex_pat")]
    
    stan_dat_con <- stan_dat[c("theta_de_con", "theta_ex_con",
                               "gamma_de_con", "gamma_ex_con",
                               "lambda_de_con", "lambda_ex_con")]
    stan_dat_pat <- as.data.frame( sapply(stan_dat_pat, function(x) apply(x, 2, median) ) ) %>% 
      rename(theta_de = theta_de_pat,
             theta_ex = theta_ex_pat,
             lambda_de = lambda_de_pat,
             lambda_ex = lambda_ex_pat,
             gamma_de = gamma_de_pat,
             gamma_ex = gamma_ex_pat)
    
    stan_dat_con <-  as.data.frame( sapply(stan_dat_con, function(x) apply(x, 2, median) ) ) %>% 
      rename(theta_de = theta_de_con,
             theta_ex = theta_ex_con,
             lambda_de = lambda_de_con,
             lambda_ex = lambda_ex_con,
             gamma_de = gamma_de_con,
             gamma_ex = gamma_ex_con)
    
  } else {
    
    stan_dat_pat <- stan_dat[c("theta_de_pat", "theta_ex_pat",
                               "lambda_de_pat", "lambda_ex_pat")]
    
    stan_dat_con <- stan_dat[c("theta_de_con", "theta_ex_con",
                               "lambda_de_con", "lambda_ex_con")]
    stan_dat_pat <- as.data.frame( sapply(stan_dat_pat, function(x) apply(x, 2, median) ) ) %>% 
      rename(theta_de = theta_de_pat,
             theta_ex = theta_ex_pat,
             lambda_de = lambda_de_pat,
             lambda_ex = lambda_ex_pat)
    
    stan_dat_con <-  as.data.frame( sapply(stan_dat_con, function(x) apply(x, 2, median) ) ) %>% 
      rename(theta_de = theta_de_con,
             theta_ex = theta_ex_con,
             lambda_de = lambda_de_con,
             lambda_ex = lambda_ex_con)
    
  }
  
  
  
  
  
  stan_dat_con$id <- d_stan$pro_id_con
  stan_dat_pat$id <- d_stan$pro_id_pat
  stan_dat_con$group <- "Control"
  stan_dat_pat$group <- "patient"
  stan_dat <- rbind(stan_dat_con, stan_dat_pat)
  
  
  stan_dat$performance_de_g_level_pat <- rep(performance_de_pat, times=length(stan_dat$id))
  stan_dat$performance_ex_g_level_pat <- rep(performance_ex_pat, times=length(stan_dat$id))
  stan_dat$performance_de_g_level_con <- rep(performance_de_con, times=length(stan_dat$id))
  stan_dat$performance_ex_g_level_con <- rep(performance_ex_con, times=length(stan_dat$id))
  
  stan_dat$iperf_de <- c(stan_fit$performance_de_con, stan_fit$performance_de_pat ) 
  stan_dat$iperf_ex <- c(stan_fit$performance_ex_con, stan_fit$performance_ex_pat ) 
  write.csv(stan_dat, str_c("data/fit/", i, "_parameters.csv"))
  
  fit_stats <- stan_fit[c("performance_pat", "performance_con", "performance_de_pat",
                          "performance_ex_pat", "performance_de_con", "performance_ex_con")]
  fit_stats[["overall_perf"]] <- stan_fit[[10]]
  
  fit_stats$performance_pat$de$looE <- NULL
  fit_stats$performance_pat$ex$looE <- NULL
  fit_stats$performance_con$de$looE <- NULL
  fit_stats$performance_con$ex$looE <- NULL
  
  saveRDS(fit_stats, paste0("data/fit/fit_stats_", i, ".rds"))
  
}



