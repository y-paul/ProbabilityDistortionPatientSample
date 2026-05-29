library(tidyverse)
library(bayestestR)
library(BayesFactor)
library(latex2exp)
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
               diff_description_lambda,
               diff_experience_lambda,
               diff_description_gamma,
               diff_experience_gamma,
               diff_description_theta,
               diff_experience_theta,
               de_gap_pat,
               de_gap_con,
               diff_de_gap))

dat_ind <- read.csv("data/fit/tk_parameters.csv")



probability_weighting_function <- function(p, gamma) {
  return(  (p^gamma) / (  (  p^gamma  + ((1-p)^gamma)  )^(gamma^(-1))   )  )
}



probs <- seq(from=0, to=1, by=0.001)

gamma_data <- data.frame(
  pw_pat_de = probability_weighting_function(probs, median(dat$mu_gamma_de_pat)),
  pw_pat_ex = probability_weighting_function(probs, median(dat$mu_gamma_ex_pat)),
  pw_con_de = probability_weighting_function(probs, median(dat$mu_gamma_de_con)),
  pw_con_ex = probability_weighting_function(probs, median(dat$mu_gamma_ex_con)),
  p = probs
)



#====================================================#
#=========== DESCRIPTION GAMMA ======================#
#====================================================#

# make plots gamma




dat_pat <- dat_ind %>% filter(group == "patient") 

median_gamma_de_patient <- median(dat$mu_gamma_de_pat)
median_gamma_ex_patient <- median(dat$mu_gamma_ex_pat)
ps_de_pat <- probability_weighting_function(probs, median_gamma_de_patient)
ps_ex_pat <- probability_weighting_function(probs, median_gamma_ex_patient)
overall_pat_gamma <- data.frame(prob = c(probs, probs),
                                ps = c(ps_de_pat, ps_ex_pat),
                                condition = rep(c("description", "experience"), each=length(c(probs))))

ind_pat_gamma <- data.frame(id=c(), prob = c(), ps = c())
for (i in 1:length(dat_pat$gamma_de)) {
  id <- i
  prob <- probs
  ps_de <-  probability_weighting_function(probs, dat_pat$gamma_de[i])
  ps_ex <-  probability_weighting_function(probs, dat_pat$gamma_ex[i])
  
  ind_pat_gamma <- rbind(ind_pat_gamma,
                         data.frame(
                           id = id,
                           prob = c(prob, prob),
                           ps = c(ps_de, ps_ex),
                           condition = rep(c("description", "experience"), each=length(c(prob)))
                         ))
}





plot_gamma_patient <- ggplot() +
  geom_line(data=ind_pat_gamma, aes(x=prob, y= ps, group=interaction(id, condition), color=condition), alpha=.04, linewidth=1) +
  geom_line(linewidth=1.1, alpha=.7, data=overall_pat_gamma, aes(x=prob, y=ps, color=condition)) +
  theme_classic() + 
  theme(text=element_text(size=12),
        plot.title= element_text(size=12)) +
  scale_color_manual(values=c("slateblue1", "firebrick")) + 
  labs( x = "Objective Probability", y = "Decision Weight",
        title="Patients",
        color="") +
  guides(color="none")







dat_con <- dat_ind %>% filter(group == "Control") 

median_gamma_de_control <- median(dat$mu_gamma_de_con)
median_gamma_ex_control <- median(dat$mu_gamma_ex_con)
ps_de_con <- probability_weighting_function(probs, median_gamma_de_control)
ps_ex_con <- probability_weighting_function(probs, median_gamma_ex_control)
overall_con_gamma <- data.frame(prob = c(probs, probs),
                                ps = c(ps_de_con, ps_ex_con),
                                condition = rep(c("description", "experience"), each=length(c(probs))))

ind_con_gamma <- data.frame(id=c(), prob = c(), ps = c())
for (i in 1:length(dat_con$gamma_de)) {
  id <- i
  prob <- probs
  ps_de <-  probability_weighting_function(probs, dat_con$gamma_de[i])
  ps_ex <-  probability_weighting_function(probs, dat_con$gamma_ex[i])
  
  ind_con_gamma <- rbind(ind_con_gamma,
                         data.frame(
                           id = id,
                           prob = c(prob, prob),
                           ps = c(ps_de, ps_ex),
                           condition = rep(c("description", "experience"), each=length(c(prob)))
                         ))
}





plot_gamma_control <- ggplot() +
  geom_line(data=ind_con_gamma, aes(x=prob, y= ps, group=interaction(id, condition), color=condition), alpha=.04, linewidth=1) +
  geom_line(linewidth=1.1, alpha=.7, data=overall_con_gamma, aes(x=prob, y=ps, color=condition)) +
  theme_classic() + 
  theme(text=element_text(size=12),
        plot.title = element_text(size=12)) +
  scale_color_manual(values=c("slateblue1", "firebrick")) + 
  labs( x = "Objective Probability", y = "Decision Weight",
        title="Controls",
        color="") +
  guides(color="none")














#====================================================#
#=========== Description Lambda =====================#
#====================================================#


value_function <- function (v, lambda) {
  
  subj_v <- sapply(v, FUN= function(x) {
    if (x > 0) {
      subjective_value <- x
    } else if (x < 0) {
      subjective_value <- lambda * x
    } else {
      subjective_value <- 0
    }
    return(subjective_value)
  })
  
  return(subj_v)
  
}
v <- seq(from=-4, to=4, by=.005)

median_lambda_de_patient <- median(dat$mu_lambda_de_pat)
median_lambda_ex_patient <- median(dat$mu_lambda_ex_pat)
v_de_pat <- value_function(v, median_lambda_de_patient)
v_ex_pat <- value_function(v, median_lambda_ex_patient)
overall_pat_lambda <- data.frame(v = c(v, v),
                                 vs = c(v_de_pat, v_ex_pat),
                                 condition = rep(c("description", "experience"), each=length(c(v))))

ind_pat_lambda <- data.frame(id=c(), v = c(), vs = c())
for (i in 1:length(dat_pat$lambda_de)) {
  id <- i
  vv <- v
  vs_de <-  value_function(vv, dat_pat$lambda_de[i])
  vs_ex <-  value_function(vv, dat_pat$lambda_ex[i])
  
  ind_pat_lambda <- rbind(ind_pat_lambda,
                          data.frame(
                            id = id,
                            v = c(v, v),
                            vs = c(vs_de, vs_ex),
                            condition = rep(c("description", "experience"), each=length(c(v)))
                          ))
}



plot_lambda_patients <- ggplot() +
  geom_line(data=ind_pat_lambda, aes(x=v, y = vs, group=interaction(id, condition), color=condition), alpha=.05, linewidth=1) + 
  geom_line(linewidth=1.1, alpha=.7,data=overall_pat_lambda, aes(x = v, y = vs, color=condition)) + 
  theme_classic() + 
  theme(text=element_text(size=12),
        plot.title = element_text(size=12)) +
  labs(x= "Objective Value", y= "Subjective Value", title="Patients", color="") +
  scale_color_manual(values=c("slateblue1", "firebrick")) +
  xlim(c(-4, 4)) +
  ylim(c(-10, 4)) +
  scale_y_continuous(limits=c(-10, 4),
                     breaks=c(4,0,-4,-8)) +
  guides(color="none")







median_lambda_de_control <- median(dat$mu_lambda_de_con)
median_lambda_ex_control <- median(dat$mu_lambda_ex_con)
v_de_con <- value_function(v, median_lambda_de_control)
v_ex_con <- value_function(v, median_lambda_ex_control)
overall_con_lambda <- data.frame(v = c(v, v),
                                 vs = c(v_de_con, v_ex_con),
                                 condition = rep(c("description", "experience"), each=length(c(v))))

ind_con_lambda <- data.frame(id=c(), v = c(), vs = c())
for (i in 1:length(dat_con$lambda_de)) {
  id <- i
  vv <- v
  vs_de <-  value_function(vv, dat_con$lambda_de[i])
  vs_ex <-  value_function(vv, dat_con$lambda_ex[i])
  
  ind_con_lambda <- rbind(ind_con_lambda,
                          data.frame(
                            id = id,
                            v = c(v, v),
                            vs = c(vs_de, vs_ex),
                            condition = rep(c("description", "experience"), each=length(c(v)))
                          ))
}



plot_lambda_controls <- ggplot() +
  geom_line(data=ind_con_lambda, aes(x=v, y = vs, group=interaction(id, condition), color=condition), alpha=.05, linewidth=1) + 
  geom_line(linewidth=1.1, alpha=.7,data=overall_con_lambda, aes(x = v, y = vs, color=condition)) + 
  theme_classic() + 
  theme(text=element_text(size=12),
        plot.title = element_text(size=12)) +
  labs(x= "Objective Value", y= "Subjective Value", title="Controls", color="") +
  scale_color_manual(values=c("slateblue1", "firebrick")) +
  xlim(c(-4, 4)) +
  ylim(c(-10, 4)) +
  scale_y_continuous(limits=c(-10, 4),
                     breaks=c(4,0,-4,-8)) +
  guides(color="none")











#====================================================#
#====================== DE Gap ======================#
#====================================================#
de_gap_data <- data.frame(
  probs=probs,
  de_gap_de_pat = probability_weighting_function(probs, median(dat$mu_gamma_de_pat)),
  de_gap_ex_pat = probability_weighting_function(probs, median(dat$mu_gamma_ex_pat)),
  de_gap_de_con = probability_weighting_function(probs, median(dat$mu_gamma_de_con)),
  de_gap_ex_con = probability_weighting_function(probs, median(dat$mu_gamma_ex_con))
)

plot_de_gap <- ggplot(data=de_gap_data) +
  geom_ribbon(data=de_gap_data, aes(x = probs, ymin = de_gap_de_pat, ymax = de_gap_ex_pat), alpha=.1) +
  geom_ribbon(data=de_gap_data, aes(x = probs, ymin = de_gap_de_con, ymax = de_gap_ex_con), alpha=.1) +
  geom_line(aes(x = probs, y = de_gap_de_pat), linetype=1, color="slateblue1", linewidth=0.8, alpha=.5) +
  geom_line(aes(x = probs, y = de_gap_ex_pat), linetype=1, color="firebrick", linewidth=0.8, alpha=.5) +
  geom_line(aes(x = probs, y = de_gap_de_con), linetype=2, color="slateblue1", linewidth=0.8, alpha=.5) +
  geom_line(aes(x = probs, y = de_gap_ex_con), linetype=2, color="firebrick", linewidth=0.8, alpha=.5) +
  theme_classic() +
  theme(text = element_text(size=12), plot.title= element_text(size=12)) + 
  labs(x= "Objective Probabilities", y="Decision Weight", title="") +
  guides(color="none")





#====================================================#
#============== Individual Parameters ===============#
#====================================================#




# gamma
max_gammas <- c(
  max(c(dat_ind$gamma_de, dat_ind$gamma_ex)),
  min(c(dat_ind$gamma_de, dat_ind$gamma_ex))
)
max_gammas[1] <- round(max_gammas[1] + 0.2,1)
max_gammas[2] <- round(max_gammas[2] - 0.2,1)
max_gammas <- c(max_gammas[2], max_gammas[1])



individual_gammas <- dat_ind %>% select(id, group, gamma_de, gamma_ex) %>% 
  pivot_longer(cols=c(gamma_de, gamma_ex)) %>% filter(value<1.6) %>% 
  mutate(name = factor(as.character(name), levels=c("gamma_ex", "gamma_de"))) %>% 
  mutate(group = recode(group, Control="Controls", patient="Patients"))


ind_gamma_plot_patients <- ggplot(individual_gammas %>% filter(group=="Patients"), aes(y=group, x=value, fill=name)) +
  geom_boxplot(notch=T, alpha=.8, outlier.alpha=0) +
  geom_point(alpha=.3, position=position_jitterdodge(jitter.width=0.3, jitter.height=0)) +
  theme_classic() +
  theme(text = element_text(size=12)) +
  scale_fill_manual(values=c("firebrick", "slateblue1")) +
  guides(fill="none") +
  theme(axis.text.y = element_text(size=12, angle=90, hjust=0.5),
        plot.title= element_text(size=12)) +
  labs(x = TeX(r'($\gamma$)'), y="", title="") +
  scale_x_continuous(breaks=c(0.2, 0.6, 1, 1.4), limits=c(0.2, 1.6))
ind_gamma_plot_patients


ind_gamma_plot_controls <- ggplot(individual_gammas %>% filter(group=="Controls"), aes(y=group, x=value, fill=name)) +
  geom_boxplot(notch=T, alpha=.8, outlier.alpha=0) +
  geom_point(alpha=.3, position=position_jitterdodge(jitter.width=0.3, jitter.height=0)) +
  theme_classic() +
  theme(text = element_text(size=12)) +
  scale_fill_manual(values=c("firebrick", "slateblue1")) +
  guides(fill="none") +
  theme(axis.text.y = element_text(size=12, angle=90, hjust=0.5),
        plot.title = element_text(size=12)) +
  labs(x = TeX(r'($\gamma$)'), y="", title="") +
  scale_x_continuous(breaks=c(0.2, 0.6, 1,  1.4), limits=c(0.2, 1.6))
ind_gamma_plot_controls



individual_lambdas <- dat_ind %>% select(id, group, lambda_de, lambda_ex) %>% 
  pivot_longer(cols=c(lambda_de, lambda_ex)) %>% 
  mutate(group = recode(group, Control = "Controls", patient = "Patients"))


ind_lambda_plot_patients  <- ggplot(individual_lambdas %>% filter(group=="Patients"), aes(y=group, x=value, fill=name)) +
  geom_boxplot(notch=T, alpha=.8, outlier.alpha=0) +
  geom_point(alpha=.3, position=position_jitterdodge(jitter.width=0.3, jitter.height=0)) +
  theme_classic() +
  theme(text = element_text(size=12)) +
  scale_fill_manual(values=c("slateblue1", "firebrick")) +
  guides(fill="none") +
  theme(axis.text.y = element_text(size=12, angle=90, hjust=0.5),
        plot.title = element_text(size=12)) +
  labs(x = TeX(r'($\lambda$)'), y="",
       title="Patients") +xlim(0,4)+
  theme(plot.title = element_text(color = "transparent", size = 12, margin = margin(b = 2)))

ind_lambda_plot_controls  <- ggplot(individual_lambdas %>% filter(group=="Controls"), aes(y=group, x=value, fill=name)) +
  geom_boxplot(notch=T, alpha=.8, outlier.alpha=0) +
  geom_point(alpha=.3, position=position_jitterdodge(jitter.width=0.3, jitter.height=0)) +
  theme_classic() +
  theme(text = element_text(size=12)) +
  scale_fill_manual(values=c("slateblue1", "firebrick")) +
  guides(fill="none") +
  theme(axis.text.y = element_text(size=12, angle=90, hjust=0.5),
        plot.title = element_text(size=12)) +
  labs(x = TeX(r'($\lambda$)'), y="",
       title="")  +xlim(0,4)













# DE-Gap

dat_ind$de_gap = dat_ind$gamma_de - dat_ind$gamma_ex

plot_pat_control_de_gap <- ggplot(data = dat_ind %>% 
                                    mutate(group = recode(group, Control = "Controls", patient = "Patients")), aes(y = group, x = de_gap)) +
  geom_boxplot(notch=T, outlier.alpha=0, fill="lightgrey") +
  theme_classic() + 
  theme(text=element_text(size=12)) +
  geom_point(alpha=.3, position=position_jitterdodge(jitter.width=0.3, jitter.height=0)) + 
  labs(x = TeX(r'($\gamma_D - \gamma_E$)'), y = "") +
  theme(axis.text.y = element_text(size=12, angle=90, hjust=c(0.7, 0.3)),
        plot.title = element_text(size=12)) +
  scale_x_continuous(breaks=c(-1, -0.5, 0, 0.5), limits=c(-1.1, 0.55))

plot_pat_control_de_gap











#====================================================#
#============== Plot Differences ====================#
#====================================================#
scaler <- -0.08
hjuster = -.1
max_dens1 <- max(graphics::hist(dat$diff_description_gamma, plot = FALSE)$density)



hdi_plot1 <- hdi(dat$diff_description_gamma)
median_plot1 <- median(dat$diff_description_gamma)
plot_posterior_gamma_de <- ggplot(data=dat, aes(x=diff_description_gamma)) +
  geom_histogram(fill="slateblue1", color="black",  aes(y=stat(density))) +
  theme_classic() +
  theme(text=element_text(size=12),
        plot.title = element_text(size=12, hjust = hjuster)) +
  labs(x=TeX(r'($\mu_{\gamma,D;P}-\mu_{\gamma,D;C}$)'), y="%",
       title="Δ Description") +
  xlim(-0.2, 0.6) +
  annotate("point", y = max_dens1 * scaler, x = median_plot1, size=4)  +
  annotate("linerange", xmin = hdi_plot1[[2]], xmax = hdi_plot1[[3]], y=max_dens1 * scaler, size=1.2)

# geom_point(aes(y = -0.2, x=median_plot),size=4)




max_dens2 <- max(graphics::hist(dat$diff_experience_gamma, plot = FALSE)$density)
hdi_plot2 <- hdi(dat$diff_experience_gamma)
median_plot2 <- median(dat$diff_experience_gamma)
plot_posterior_gamma_ex <- ggplot(data=dat, aes(x=diff_experience_gamma)) +
  geom_histogram(fill="firebrick", color="black",  aes(y=stat(density))) +
  theme_classic() +
  theme(text=element_text(size=12),
        plot.title = element_text(size=12, hjust = hjuster)) +
  labs(x=TeX(r'($\mu_{\gamma,E;P}-\mu_{\gamma,E;C}$)'), y="%",
       title = "Δ Experience") +
  xlim(-0.2, 0.6) +
  annotate("point", y = max_dens2 * scaler, x = median_plot2, size=4)  +
  annotate("linerange", xmin = hdi_plot2[[2]], xmax = hdi_plot2[[3]], y=max_dens2 * scaler, size=1.2) 


max_dens3 <- max(graphics::hist(dat$diff_description_lambda, plot = FALSE)$density)
hdi_plot3 <- hdi(dat$diff_description_lambda)
median_plot3 <- median(dat$diff_description_lambda)
plot_posterior_lambda_de <- ggplot(data=dat, aes(x=diff_description_lambda)) +
  geom_histogram(fill="slateblue1", color="black",  aes(y=stat(density))) +
  theme_classic() +
  theme(text=element_text(size=12),
        plot.title = element_text(size=12, hjust = hjuster)) +
  labs(x=TeX(r'($\mu_{\lambda,D;P}-\mu_{\lambda,D;C}$)'), y="%",
       title="Δ Description") +
  annotate("point", y = max_dens3 * scaler, x = median_plot3, size=4)  +
  annotate("linerange", xmin = hdi_plot3[[2]], xmax = hdi_plot3[[3]], y=max_dens3 * scaler, size=1.2) +
  xlim(-1, 2)


max_dens4 <- max(graphics::hist(dat$diff_experience_lambda, plot = FALSE)$density)
hdi_plot4 <- hdi(dat$diff_experience_lambda)
median_plot4 <- median(dat$diff_experience_lambda)
plot_posterior_lambda_ex <- ggplot(data=dat, aes(x=diff_experience_lambda)) +
  geom_histogram(fill="firebrick", color="black",  aes(y=stat(density))) +
  theme_classic() +
  theme(text=element_text(size=12),
        plot.title = element_text(size=12, hjust = hjuster)) +
  labs(x=TeX(r'($\mu_{\lambda,E;P}-\mu_{\lambda,E;C}$)'), y="%",
       title="Δ Experience") +
  annotate("point", y = max_dens4 * scaler, x = median_plot4, size=4)  +
  annotate("linerange", xmin = hdi_plot4[[2]], xmax = hdi_plot4[[3]], y=max_dens4 * scaler, size=1.2)  +
  xlim(-1, 2)



max_dens5 <- max(graphics::hist(dat$diff_de_gap, plot = FALSE)$density)
hdi_plot5 <- hdi(dat$diff_de_gap)
median_plot5 <- median(dat$diff_de_gap)
plot_posterior_de_gap <- ggplot(data=dat, aes(x=diff_de_gap)) +
  geom_histogram(fill="lightgrey", color="black",  aes(y=stat(density))) +
  theme_classic() +
  theme(text=element_text(size=12),
        plot.title = element_text(size=12, hjust = hjuster)) +
  labs(x=TeX(r'($\mu_{DEG;P}-\mu_{DEG;C}$)'), y="%",
       title="Δ DE Gap") +
  annotate("point", y = max_dens5 * scaler, x = median_plot5, size=4)  +
  annotate("linerange", xmin = hdi_plot5[[2]], xmax = hdi_plot5[[3]], y=max_dens5 * scaler, size=1.2) 





























dummy_legend <- data.frame(label1 = factor(c("Description", "Experience", "DE gap"),
                                           levels=c("Description", "Experience", "DE gap")), label2 = factor(c("Controls", "Patients", "Patients"), levels=c("Patients", "Controls")),
                           x = c(0, 0, 0), y = c(0,0, 0)) %>% 
  ggplot(aes(color=label1, linetype = label2,x=x,y=y)) + 
  geom_line() + 
  theme_void() +
  scale_color_manual(values=c("slateblue1", "firebrick", "lightgrey")) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        text=element_text(size=12),
        legend.text=element_text(size=10),
        legend.key.size=unit(1, "cm")) +
  guides(color=guide_legend(nrow=1, override.aes=list(size=2, linewidth=2)),
         linetype=guide_legend(nrow=1, override.aes=list(size=2, linewidth=2))) +
  theme(plot.margin=margin(-5,0,0,0),
        legend.margin=margin(-5,0,0,0))



library(patchwork)
library(grid)


remove_plot_tags <- function(p) {
  p +
    labs(tag = NULL) +
    theme(
      plot.tag = element_blank(),
      plot.tag.position = "topleft"
    )
}

plot_lambda_patients      <- remove_plot_tags(plot_lambda_patients)
plot_lambda_controls      <- remove_plot_tags(plot_lambda_controls)
ind_lambda_plot_patients  <- remove_plot_tags(ind_lambda_plot_patients)
ind_lambda_plot_controls  <- remove_plot_tags(ind_lambda_plot_controls)
plot_posterior_lambda_de  <- remove_plot_tags(plot_posterior_lambda_de)
plot_posterior_lambda_ex  <- remove_plot_tags(plot_posterior_lambda_ex)

plot_gamma_patient        <- remove_plot_tags(plot_gamma_patient)
plot_gamma_control        <- remove_plot_tags(plot_gamma_control)
ind_gamma_plot_patients   <- remove_plot_tags(ind_gamma_plot_patients)
ind_gamma_plot_controls   <- remove_plot_tags(ind_gamma_plot_controls)
plot_posterior_gamma_de   <- remove_plot_tags(plot_posterior_gamma_de)
plot_posterior_gamma_ex   <- remove_plot_tags(plot_posterior_gamma_ex)

plot_de_gap               <- remove_plot_tags(plot_de_gap)
plot_pat_control_de_gap   <- remove_plot_tags(plot_pat_control_de_gap)
plot_posterior_de_gap     <- remove_plot_tags(plot_posterior_de_gap)


compact_theme <- theme(
  text = element_text(size = 10),
  
  plot.title = element_text(
    size = 10,
    margin = margin(t = 0, r = 0, b = 1, l = 0)
  ),
  
  axis.title.x = element_text(
    size = 10,
    margin = margin(t = 1, r = 0, b = 0, l = 0)
  ),
  
  axis.title.y = element_text(
    size = 10,
    margin = margin(t = 0, r = 1, b = 0, l = 0)
  ),
  
  axis.text = element_text(size = 10),
  
  axis.text.x = element_text(
    size = 10,
    margin = margin(t = 1)
  ),
  
  axis.text.y = element_text(
    size = 10,
    margin = margin(r = 1)
  ),
  
  strip.text = element_text(
    size = 10,
    margin = margin(t = 1, r = 1, b = 1, l = 1)
  ),
  
  legend.title = element_text(size = 10),
  legend.text  = element_text(size = 8),
  
  legend.key.size   = unit(0.35, "cm"),
  legend.key.width  = unit(0.45, "cm"),
  legend.key.height = unit(0.25, "cm"),
  
  legend.spacing.x = unit(0.05, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
  legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
  
  plot.margin = margin(t = 1, r = 2, b = 1, l = 2)
)

make_compact <- function(p) {
  p + compact_theme
}

plot_lambda_patients      <- make_compact(plot_lambda_patients)
plot_lambda_controls      <- make_compact(plot_lambda_controls)
ind_lambda_plot_patients  <- make_compact(ind_lambda_plot_patients)
ind_lambda_plot_controls  <- make_compact(ind_lambda_plot_controls)
plot_posterior_lambda_de  <- make_compact(plot_posterior_lambda_de)
plot_posterior_lambda_ex  <- make_compact(plot_posterior_lambda_ex)

plot_gamma_patient        <- make_compact(plot_gamma_patient)
plot_gamma_control        <- make_compact(plot_gamma_control)
ind_gamma_plot_patients   <- make_compact(ind_gamma_plot_patients)
ind_gamma_plot_controls   <- make_compact(ind_gamma_plot_controls)
plot_posterior_gamma_de   <- make_compact(plot_posterior_gamma_de)
plot_posterior_gamma_ex   <- make_compact(plot_posterior_gamma_ex)

plot_de_gap               <- make_compact(plot_de_gap)
plot_pat_control_de_gap   <- make_compact(plot_pat_control_de_gap)
plot_posterior_de_gap     <- make_compact(plot_posterior_de_gap)


add_invisible_title <- function(p) {
  p +
    labs(title = "placeholder") +
    theme(
      plot.title = element_text(
        size = 8,
        color = "transparent",
        margin = margin(t = 0, r = 0, b = 1, l = 0)
      )
    )
}

plot_de_gap             <- add_invisible_title(plot_de_gap)
plot_pat_control_de_gap <- add_invisible_title(plot_pat_control_de_gap)
plot_posterior_de_gap   <- add_invisible_title(plot_posterior_de_gap)


section_title <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = -0.1,
      y = 0.5,
      label = label,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 3
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(t = 1, r = 0, b = 1, l = 2)
    )
}



dummy_legend <- dummy_legend +
  theme(
    text = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.spacing.x = unit(0.05, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(t = -2, r = 0, b = 0, l = 0)
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(size = 1, linewidth = 1)
    ),
    linetype = guide_legend(
      nrow = 1,
      override.aes = list(size = 1, linewidth = 1)
    )
  )

sp1 <- (
  (plot_lambda_patients / plot_lambda_controls) |
    (ind_lambda_plot_patients / ind_lambda_plot_controls) |
    (plot_posterior_lambda_de / plot_posterior_lambda_ex)
) +
  plot_layout(widths = c(1, 1, 1))

sp2 <- (
  (plot_gamma_patient / plot_gamma_control) |
    (ind_gamma_plot_patients / ind_gamma_plot_controls) |
    (plot_posterior_gamma_de / plot_posterior_gamma_ex)
) +
  plot_layout(widths = c(1, 1, 1))

sp3 <- (
  plot_de_gap |
    plot_pat_control_de_gap |
    plot_posterior_de_gap
) +
  plot_layout(widths = c(1, 1, 1))


new_plot <-
  section_title("A) Loss Aversion") /
  sp1 /
  section_title("B) Nonlinear Probability Weighting") /
  sp2 /
  section_title("C) DE Gap") /
  sp3 /
  dummy_legend +
  plot_layout(
    heights = c(
      0.10,  # A title
      2.00,  # A plots
      0.10,  # B title
      2.00,  # B plots
      0.08,  # C title
      0.70,  # C plots
      0.16   # legend
    )
  )


# ------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------

ggsave(
  "main_res_plot.png",
  plot = new_plot,
  dpi = 300,
  units = "cm",
  width = 14,
  height = 18
)
