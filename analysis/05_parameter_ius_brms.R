library(data.table)
library(brms)
library(ggplot2)
library(patchwork)

RECALC_BAYES <- T
N_CORES <- 4
N_CHAINS <- 4
N_ITER <- 4000
N_WARMUP <- 1000


if (RECALC_BAYES) {
  
  default_prior(x~y,
                data = dd,
                family = student,
                cores = N_CORES,
                chains = N_CHAINS,
                iter = N_ITER,
                warmup = N_WARMUP,
                file = "brms_mods",
                file_refit = "on_change")
  

# data prep ---------------------------------------------------------------

dd = data.frame(fread("data_cor.csv"))
dd$V1 = NULL


x_vars <- c("gamma_de", "gamma_ex", "lambda_de", "lambda_ex", "de_gap", "sample_mean")
y_vars <- c("UIS_GESAMT", "UIS_A", "UIS_B", "UIS_C")



data <- data.frame(xvar=c(),
                   yvar=c(),
                   Estimate=c(),
                   Est.Error=c(),
                   Q2.5=c(),
                   Q97.5=c(),
                   R2.Estimate=c(),
                   R2.Est.Error=c(),
                   R2.Q2.5=c(),
                   R2.Q97.5=c())

for (x in x_vars) {
  for (y in y_vars) {
    
    if (x == "de_gap") {
    dd$x <- scale(dd[[x]])
    dd$y <- scale(log(dd[[y]]))
    } else {
      dd$x <- scale(log(dd[[x]]))
      dd$y <- scale(log(dd[[y]]))
    }

    
    m <- brm(x~y,
             data = dd,
             family = student,
             cores = N_CORES,
             chains = N_CHAINS,
             iter = N_ITER,
             warmup = N_WARMUP,
             file = "brms_mods",
             file_refit = "on_change")
    
    
    effects <- as.data.frame(t(as.data.frame(fixef(m)[2,])))
    r2 <- bayes_R2(m)
    colnames(r2) <- paste0("R2.", colnames(r2))
    
    temp <- cbind(
      data.frame(xvar = x, yvar = y),
      effects,
      r2
    )
    
   data <- rbind(data, temp)
    
    
  }
}




write.csv(data, "brms_data.csv")
}


data <- read.csv("brms_data.csv")


data$credible <- data$Q2.5 > 0 | data$Q97.5 < 0
library(tidyverse)
library(latex2exp)

brms_plot <- data %>% mutate(xvar = fct_relevel(xvar, "lambda_de", "lambda_ex", "gamma_de", "gamma_ex", "de_gap", "sample_mean")) %>% 

  mutate(yvar = fct_relevel(yvar, "UIS_GESAMT", "UIS_A", "UIS_B", "UIS_C")) %>% 
  mutate(yvar = recode(yvar, UIS_GESAMT="IUS-18\nTotal Score",
                       UIS_A = "IUS-18\nImpaired Ability",
                       UIS_B = "IUS-18\nDistress",
                       UIS_C = "IUS-18\nVigilance")) %>%
  mutate(r2_label = paste0("R^2 == '", sprintf("%.3f", R2.Estimate), "'")) %>% 
  mutate(credible = recode(as.factor(credible), `TRUE` = "1", `FALSE` = "2")) %>% 
  mutate(credible = ifelse((xvar == "sample_mean") & (credible == "1"), "3", credible)) %>% 
  mutate(credible = factor(as.character(credible), levels=c("1", "2", "3"))) %>% 
  mutate(credible = recode(credible, `2` = "credible", `1` = "not credible", `3` = "credible with 𝘙² < 0.04")) %>% 
  mutate(r2_label = ifelse(credible == "credible", r2_label, "")) %>% 
  ggplot(aes(x = xvar,
             y = Estimate,
             ymin = Q2.5,
             ymax = Q97.5,
             color=credible,
             label=r2_label)) +
  geom_hline(yintercept=0, size=0.9, alpha=.3)+
  #geom_point() +
  geom_pointrange(size=0.2, linewidth=0.7) +
  facet_wrap(~yvar, nrow=1) +
  theme_classic() +
  scale_x_discrete(labels=c(
    gamma_de = TeX(r'($\gamma_D$)'),
    gamma_ex = TeX(r'($\gamma_E$)'),
    lambda_de = TeX(r'($\lambda_D$)'),
    lambda_ex = TeX(r'($\lambda_E$)'),
    de_gap = "DE gap",
    sample_mean = "Median\nSamples"
  )) +
  labs(x = "",
       y= "Estimated Effects [±HDI]",
       color="") +
  scale_color_manual(values=c("black", "firebrick2", "darkviolet")) +
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle=90,
                                   vjust=c(0.5, 0.5, 0.5, 0.5, 0.5, 1), hjust=1),
        legend.position = "bottom",           # slightly below plot, left-aligned
        legend.justification = "left",          # anchor legend's top-left to the point
        legend.direction = "horizontal")+
  #geom_text(angle=90, size.unit = "pt", size=6,
  #          vjust = -1, color="firebrick2", parse=T) +
  guides(text = "none") +
  guides(color = "none") +
  theme(
    panel.spacing.x = unit(0, "pt"),
    panel.spacing.y = unit(0, "pt"),
    panel.border = element_blank(),
    strip.placement = "outside"
  ) +

  theme(
    strip.background = element_rect(fill="grey90", color=NA, linewidth = 1),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA)
  )

brms_plot







brms_plot_r2 <- data %>% mutate(xvar = fct_relevel(xvar, "lambda_de", "lambda_ex", "gamma_de", "gamma_ex", "de_gap", "sample_mean")) %>% 
  
  mutate(yvar = fct_relevel(yvar, "UIS_GESAMT", "UIS_A", "UIS_B", "UIS_C")) %>% 
  mutate(yvar = recode(yvar, UIS_GESAMT="IUS-18\nTotal Score",
                       UIS_A = "IUS-18\nImpaired Ability",
                       UIS_B = "IUS-18\nDistress",
                       UIS_C = "IUS-18\nVigilance")) %>%
  mutate(r2_label = paste0("R^2 == '", sprintf("%.3f", R2.Estimate), "'")) %>% 
  mutate(credible = recode(as.factor(credible), `TRUE` = "1", `FALSE` = "2")) %>% 
  mutate(credible = ifelse((xvar == "sample_mean") & (credible == "1"), "3", credible)) %>% 
  mutate(credible = factor(as.character(credible), levels=c("1", "2", "3"))) %>% 
  mutate(credible = recode(credible, `2` = "credible", `1` = "not credible", `3` = "credible with 𝘙² < 0.04")) %>% 
  mutate(r2_label = ifelse(credible == "credible", r2_label, "")) %>% 
  ggplot(aes(x = xvar,
             y = R2.Estimate,
             ymin = R2.Q2.5,
             ymax = R2.Q97.5,
             color=credible,
             label=r2_label)) +
  geom_hline(yintercept=0, size=0.9, alpha=.3)+
  geom_hline(yintercept=0.04, size=0.9, alpha=.3)+
  #geom_point() +
  geom_pointrange(size=0.2, linewidth=0.7) +
  facet_wrap(~yvar, nrow=1) +
  theme_classic() +
  scale_x_discrete(labels=c(
    gamma_de = TeX(r'($\gamma_D$)'),
    gamma_ex = TeX(r'($\gamma_E$)'),
    lambda_de = TeX(r'($\lambda_D$)'),
    lambda_ex = TeX(r'($\lambda_E$)'),
    de_gap = "DE gap",
    sample_mean = "Median\nSamples"
  )) +
  labs(x = "Model Parameters/Behavioral Indices",
       y= "𝘙² [±HDI]",
       color="") +
  scale_color_manual(values=c("black", "firebrick2", "darkviolet")) +
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle=90,
                                   vjust=c(0.5, 0.5, 0.5, 0.5, 0.5, 1), hjust=1),
        legend.position = "bottom",           # slightly below plot, left-aligned
        legend.justification = "left",          # anchor legend's top-left to the point
        legend.direction = "horizontal")+
  #geom_text(angle=90, size.unit = "pt", size=6,
  #          vjust = -1, color="firebrick2", parse=T) +
  guides(text="none") +
  theme(
    panel.spacing.x = unit(0, "pt"),
    panel.spacing.y = unit(0, "pt"),
    panel.border = element_blank(),
    strip.placement = "outside"
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  )

brms_plot_r2




library(patchwork)

brms_super_plot <- brms_plot / plot_spacer() /brms_plot_r2 +
  plot_layout(heights=c(1, -0.3, 1))
ggsave("brms_plot_super.png",plot=brms_super_plot, dpi=300, units="cm", width=17, height=12.5)










