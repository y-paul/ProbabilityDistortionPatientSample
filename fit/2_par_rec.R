library(rstan)

rm(list=ls())


# CPT functions -----------------------------------------------------------

source('fit/functions/CPT.R')

# data --------------------------------------------------------------------



# for convenience - can be reconstructed from dat
d_stan_patient <- readRDS("data/processed/patient/d_stan.rds")
d_stan_control <- readRDS("data/processed/control/overall/d_stan.rds")


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



XA_pat = d_stan$xa_pat; XA_pat[,'a_se'] = -1 * XA_pat[,'a_se']
XB_pat = d_stan$xb_pat; XB_pat[,'b_se'] = -1 * XB_pat[,'b_se']
PA_pat = d_stan$pa_pat; PB_pat = d_stan$pb_pat
FA_pat = d_stan$fa_pat; FB_pat = d_stan$fb_pat

sub_pat = d_stan$sub_pat; N_pat = d_stan$N_pat; n_pat = d_stan$n_pat


XA_con = d_stan$xa_con; XA_con[,'a_se'] = -1 * XA_con[,'a_se']
XB_con = d_stan$xb_con; XB_con[,'b_se'] = -1 * XB_con[,'b_se']
PA_con = d_stan$pa_con; PB_con = d_stan$pb_con
FA_con = d_stan$fa_con; FB_con = d_stan$fb_con

sub_con = d_stan$sub_con; N_con = d_stan$N_con; n_con = d_stan$n_con

# modeling results --------------------------------------------------------

est_mods = readRDS('data/fit/est_mods.rds')

# main model's id pars
ip = sapply(est_mods$TK$pars$i_pars, function(x) apply(x, 2, median))
ip <- as.data.frame(ip)

# choice a probs ----------------------------------------------------------

pa_desc_pat = c(); pa_exp_pat = c()
pa_desc_con = c(); pa_exp_con = c()

for(i in 1:n_pat) {
  
  s = sub_pat[i]
  
  # simulate choices of pa
  pa_desc_pat[i] = cpt_mod_pr(XA = XA_pat[i,], XB = XB_pat[i,],
                          PA = PA_pat[i,], PB = PB_pat[i,],
                          g = ip$gamma_de_pat[s],
                          l = ip$lambda_de_pat[s],
                          t = ip$theta_de_pat[s],
                          out = 'pa')
  
  # simulate choices of pa
  pa_exp_pat[i] = cpt_mod_pr(XA = XA_pat[i,], XB = XB_pat[i,],
                          PA = FA_pat[i,], PB = FB_pat[i,],
                          g = ip$gamma_ex_pat[s],
                          l = ip$lambda_ex_pat[s],
                          t = ip$theta_ex_pat[s],
                          out = 'pa')
}


for(i in 1:n_con) {
  
  s = sub_con[i]
  
  # simulate choices of pa
  pa_desc_con[i] = cpt_mod_pr(XA = XA_con[i,], XB = XB_con[i,],
                              PA = PA_con[i,], PB = PB_con[i,],
                              g = ip$gamma_de_con[s],
                              l = ip$lambda_de_con[s],
                              t = ip$theta_de_con[s],
                              out = 'pa')
  
  # simulate choices of pa
  pa_exp_con[i] = cpt_mod_pr(XA = XA_con[i,], XB = XB_con[i,],
                             PA = FA_con[i,], PB = FB_con[i,],
                             g = ip$gamma_ex_con[s],
                             l = ip$lambda_ex_con[s],
                             t = ip$theta_ex_con[s],
                             out = 'pa')
}


# simulate choices --------------------------------------------------------

# no of replications
L = 20
co_sim_pat = replicate(L, apply(cbind(pa_desc_pat, pa_exp_pat), 1:2,
                          rbinom, n = 1, size = 1),
                   simplify = F)
co_sim_con = replicate(L, apply(cbind(pa_desc_con, pa_exp_con), 1:2,
                                rbinom, n = 1, size = 1),
                       simplify = F)

# fit replicated data -----------------------------------------------------

# transalate and compile the model
trans = stanc(file = 'fit/stan_estimation/stan_mods_groups/no_delta/03_TK.stan')
compiled = stan_model(stanc_ret = trans, verbose = F)

#
i_pars = c('theta_de_pat', "theta_ex_pat", 'gamma_de_pat', 'gamma_ex_pat', 'lambda_de_pat', 'lambda_ex_pat',
           'theta_de_con', "theta_ex_con", 'gamma_de_con', 'gamma_ex_con', 'lambda_de_con', 'lambda_ex_con')
p_pars = c('mu_theta_de_pat', "mu_theta_ex_pat", 'mu_gamma_de_pat', 'mu_gamma_ex_pat', 
           'mu_lambda_de_pat', 'mu_lambda_ex_pat',
           'mu_theta_de_con', "mu_theta_ex_con", 'mu_gamma_de_con', 'mu_gamma_ex_con', 
           'mu_lambda_de_con', 'mu_lambda_ex_con')

# sampler parameters
n_chains = 8
n_iter = 2e3
n_burn = 1e3
n_thin = 4
n_cores = n_chains

#
REC_PARS = list()

# loop over the models
for(l in 1:L) {
  
  # set data
  d_stan$choices_pat = co_sim_pat[[l]]
  d_stan$choices_con = co_sim_con[[l]]
  
  # get the posterior samples (estimation!!)
  stanfit = sampling(object = compiled,
                     data = d_stan,
                     pars = c(i_pars, p_pars),
                     init = '0',
                     chains = n_chains,
                     iter = n_iter,
                     warmup = n_burn,
                     thin = n_thin,
                     cores = n_chains,
                     control = NULL)
  
  
  # posterior samples
  p_pars_post = rstan::extract(stanfit)[p_pars]
  i_pars_post = rstan::extract(stanfit)[i_pars]
  
  # list with results
  REC_PARS[[l]] = list(p_pars = p_pars_post,
                       i_pars = i_pars_post)
  
  rm(p_pars_post, i_pars_post)
  print(l)
  
}

saveRDS(REC_PARS, file = 'data/fit/parameter_recovery.rds')

# results -----------------------------------------------------------------

REC_PARS = readRDS('data/fit/parameter_recovery.rds')




summarized <- list()
for (el in 1:length(REC_PARS)) {
  
  current_replication <- REC_PARS[[el]]$i_pars
  lambda_de_pat <- apply(current_replication$lambda_de_pat, 2, median)
  lambda_ex_pat <- apply(current_replication$lambda_ex_pat, 2, median)
  gamma_de_pat <- apply(current_replication$gamma_de_pat, 2, median)
  gamma_ex_pat <- apply(current_replication$gamma_ex_pat, 2, median)
  theta_de_pat <- apply(current_replication$theta_de_pat, 2, median)
  theta_ex_pat <- apply(current_replication$theta_ex_pat, 2, median)
  
  lambda_de_con <- apply(current_replication$lambda_de_con, 2, median)
  lambda_ex_con <- apply(current_replication$lambda_ex_con, 2, median)
  gamma_de_con <- apply(current_replication$gamma_de_con, 2, median)
  gamma_ex_con <- apply(current_replication$gamma_ex_con, 2, median)
  theta_de_con <- apply(current_replication$theta_de_con, 2, median)
  theta_ex_con <- apply(current_replication$theta_ex_con, 2, median)
  n_con <- length(lambda_de_con)
  n_pat <- length(lambda_de_pat)
  
  lambda_de <- c(lambda_de_pat, lambda_de_con)
  lambda_ex <- c(lambda_ex_pat, lambda_ex_con)
  gamma_de <- c(gamma_de_pat, gamma_de_con)
  gamma_ex <- c(gamma_ex_pat, gamma_ex_con)
  theta_de <- c(theta_de_pat, theta_de_con)
  theta_ex <- c(theta_ex_pat, theta_ex_con)
  group <- rep(c("pat", "con"), times=c(n_pat, n_con))
  
  summarized[[el]] <- 
    list(lambda_de_rec=lambda_de, lambda_ex_rec=lambda_ex,
         gamma_de_rec=gamma_de, gamma_ex_rec=gamma_ex,
         theta_de_rec = theta_de, theta_ex_rec = theta_ex,
         group = group)
}
num_fields <- c("lambda_de_rec","lambda_ex_rec","gamma_de_rec","gamma_ex_rec","theta_de_rec","theta_ex_rec")

means_between_repeats <- lapply(num_fields, function(nm) {
  rowMeans(do.call(cbind, lapply(summarized, `[[`, nm)))
})
names(means_between_repeats) <- num_fields

# keep group as-is (from first replication)
means_between_repeats$group <- summarized[[1]]$group

means_between_repeats <- as.data.frame(means_between_repeats)
gen_list <- as.data.frame(list(
lambda_de_gen = c(ip$lambda_de_pat, ip$lambda_de_con),
lambda_ex_gen = c(ip$lambda_ex_pat, ip$lambda_ex_con),
gamma_de_gen = c(ip$gamma_de_pat, ip$gamma_de_con),
gamma_ex_gen = c(ip$gamma_ex_pat, ip$gamma_ex_con),
theta_de_gen = c(ip$theta_de_pat, ip$theta_de_con),
theta_ex_gen = c(ip$theta_ex_pat, ip$theta_ex_con),
group = rep(c("pat", "con"), times=c(length(ip$lambda_de_pat), length(ip$lambda_de_con)))
))


rec_data <- cbind(means_between_repeats, gen_list)

all_cors <- c(
cor(rec_data$lambda_de_rec, rec_data$lambda_de_gen),
cor(rec_data$lambda_ex_rec, rec_data$lambda_ex_gen),
cor(rec_data$gamma_de_rec, rec_data$gamma_de_gen),
cor(rec_data$gamma_ex_rec, rec_data$gamma_ex_gen),
cor(rec_data$theta_de_rec, rec_data$theta_de_gen),
cor(rec_data$theta_ex_rec, rec_data$theta_ex_gen)
)
mean(all_cors)
median(all_cors)


cor(rec_data$lambda_de_rec, rec_data$lambda_de_gen)
cor(rec_data$lambda_ex_rec, rec_data$lambda_ex_gen)
cor(rec_data$gamma_de_rec, rec_data$gamma_de_gen)
cor(rec_data$gamma_ex_rec, rec_data$gamma_ex_gen)
cor(rec_data$theta_de_rec, rec_data$theta_de_gen)
cor(rec_data$theta_ex_rec, rec_data$theta_ex_gen)



write.csv2(rec_data, 'fit/parameter_recovery_medians.csv',
           row.names = F)


n <- length(rec_data$lambda_de_rec)
rec_data$id <- 1:n
rec_plot <- rec_data  %>%  pivot_longer(cols=c(
  lambda_de_rec, lambda_de_gen,
  lambda_ex_rec, lambda_ex_gen,
  gamma_de_rec, gamma_de_gen,
  gamma_ex_rec, gamma_ex_gen,
  theta_de_rec, theta_de_gen,
  theta_ex_rec, theta_ex_gen,
)) %>% separate(col=name, into=c("parameter", "condition", "type"), sep="_") %>% 
  pivot_wider(id_cols =  c(id, parameter, condition, group), names_from=c(type), values_from=value) %>% 
  ggplot(aes(x=gen, y=rec, color=condition, shape=group)) + 
  geom_line(stat="smooth", method="lm", alpha=.5) +
  geom_point(alpha=.4) + 
  facet_wrap(~parameter, scales="free") +
  theme_classic() +
  theme(text=element_text(size=12)) +
  scale_color_manual(values=c("darkolivegreen", "violet"))
rec_plot  

ggsave("rec_plot.png", plot=rec_plot, dpi=300, width=16, height=6)



