




rm(list=ls())

DATA_FILE_TO_USE <- "data/prediction_check/probs_TK_parameters.csv"

library(tidyverse)
library(patchwork)
dat <- read.csv(DATA_FILE_TO_USE)
dat$group <- dat$group.x

# optim EV





# Description
exp_ev_max_1 <- dat %>% filter(cond=="desc") %>% 
  filter(!is.na(chosen_optim_ev)) %>% group_by(pro_id) %>% 
  summarize(m_behave = mean(chosen_optim_ev),
            sem_behave = sd(chosen_optim_ev)/sqrt(length(chosen_optim_ev)),
            m_model = mean(modprob_optim_ev),
            sem_model = sd(modprob_optim_ev)/sqrt(length(modprob_optim_ev)),
            group=dplyr::first(group)) %>% 
  arrange(m_behave) %>% 
  group_by(group) %>%
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%  pivot_longer(cols=c(m_behave, m_model)) %>% 
  mutate(name = recode(name, m_behave = "Behavioral Proportions",
                       m_model ="Model Predictions")) %>% 
  mutate(upper = ifelse(name=="Behavioral Proportions", value + sem_behave, value + sem_model),
         lower = ifelse(name=="Behavioral Proportions", value - sem_behave, value - sem_model),
         group = recode(group, control = "Control", patient = "Patient")) %>% 
  ggplot(aes(x = pro_id, color=name, y = value, ymax = upper, ymin = lower)) +
  geom_hline(yintercept=.5)+
  geom_point(size = 1, alpha=.5) +
  theme_classic() + ylim(0, 1) +
  labs(subtitle="Maximize\nExpected Value",
       x = "",
       y = "Proportion of Choices", color="") +
  scale_color_manual(values=c("black", "slateblue1")) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=12)) +
  guides(color="none") + facet_wrap(~group, nrow=2)



# Experie
exp_ev_max_2 <- dat %>%
  filter(cond == "exp", !is.na(chosen_optim_ev)) %>%
  group_by(pro_id) %>%
  summarise(
    m_behave   = mean(chosen_optim_ev),
    sem_behave = sd(chosen_optim_ev) / sqrt(n()),
    m_model    = mean(modprob_optim_ev),
    sem_model  = sd(modprob_optim_ev) / sqrt(n()),
    group      = dplyr::first(group),
    .groups    = "drop"
  ) %>%
  arrange(group, m_behave) %>%      # sort within group
  group_by(group) %>%               # re-index x within facet
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%
  pivot_longer(cols = c(m_behave, m_model)) %>%
  mutate(
    name = recode(name,
                  m_behave = "Behavioral Proportions",
                  m_model  = "Model Predictions"),
    upper = ifelse(name == "Behavioral Proportions", value + sem_behave, value + sem_model),
    lower = ifelse(name == "Behavioral Proportions", value - sem_behave, value - sem_model),
    group = recode(group, control="Control", patient = "Patient")
  ) %>%
  ggplot(aes(x = pro_id, y = value, color = name, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = .5) +
  geom_point(size = 1, alpha = .5) +
  theme_classic() + ylim(0, 1) +
  labs(
    subtitle = "",
    x = "Participant",
    y = "Proportions\nof Choices",
    color = ""
  ) +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 12)
  ) +
  guides(color = "none") +
  facet_wrap(~group, nrow=2)


#exp_ev_max_2



# SEP





exp_se_freq1 <- dat %>%
  filter(cond == "desc", !is.na(chosen_se_freq)) %>%
  group_by(pro_id) %>%
  summarise(
    m_behave   = mean(chosen_se_freq),
    sem_behave = sd(chosen_se_freq) / sqrt(n()),
    m_model    = mean(modpob_se_freq),
    sem_model  = sd(modpob_se_freq) / sqrt(n()),
    group      = dplyr::first(group),
    .groups    = "drop"
  ) %>%
  arrange(group, m_behave) %>%
  group_by(group) %>%
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%
  pivot_longer(cols = c(m_behave, m_model)) %>%
  mutate(
    name  = recode(name, m_behave = "Behavioral Proportions", m_model = "Model Predictions"),
    upper = ifelse(name == "Behavioral Proportions", value + sem_behave, value + sem_model),
    lower = ifelse(name == "Behavioral Proportions", value - sem_behave, value - sem_model),
    group = recode(group, control = "Control", patient = "Patient")
  ) %>%
  ggplot(aes(x = pro_id, color = name, y = value, ymax = upper, ymin = lower)) +
  geom_hline(yintercept = .5) +
  geom_point(size = 1, alpha = .5) +
  theme_classic() + ylim(0, 1) +
  labs(subtitle = "Minimize\nSide Effect Prob.", x = "", y = "", color = "") +
  scale_color_manual(values = c("black", "slateblue1")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    text = element_text(size = 12)
  ) +
  guides(color = "none") +
  facet_wrap(~group, nrow = 2)









exp_se_freq2 <- dat %>%
  filter(cond == "exp", !is.na(chosen_se_freq)) %>%
  group_by(pro_id) %>%
  summarise(
    m_behave   = mean(chosen_se_freq),
    sem_behave = sd(chosen_se_freq) / sqrt(n()),
    m_model    = mean(modpob_se_freq),
    sem_model  = sd(modpob_se_freq) / sqrt(n()),
    group      = dplyr::first(group),
    .groups    = "drop"
  ) %>%
  arrange(group, m_behave) %>%
  group_by(group) %>%
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%
  pivot_longer(cols = c(m_behave, m_model)) %>%
  mutate(
    name  = recode(name, m_behave = "Behavioral Proportions", m_model = "Model Predictions"),
    upper = ifelse(name == "Behavioral Proportions", value + sem_behave, value + sem_model),
    lower = ifelse(name == "Behavioral Proportions", value - sem_behave, value - sem_model),
    group = recode(group, control = "Control", patient = "Patient")
  ) %>%
  ggplot(aes(x = pro_id, color = name, y = value, ymax = upper, ymin = lower)) +
  geom_hline(yintercept = .5) +
  geom_point(size = 1, alpha = .5) +
  theme_classic() + ylim(0, 1) +
  labs(subtitle = "", x = "Participant", y = "Proportions of Choices", color = "") +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 12)
  ) +
  guides(color = "none") +
  facet_wrap(~group, nrow = 2)


exp_ev_be_max_1 <- dat %>%
  filter(cond == "desc", !is.na(chosen_optim_be_ev)) %>%
  group_by(pro_id) %>%
  summarise(
    m_behave   = mean(chosen_optim_be_ev),
    sem_behave = sd(chosen_optim_be_ev) / sqrt(n()),
    m_model    = mean(modprob_optim_be_ev),
    sem_model  = sd(modprob_optim_be_ev) / sqrt(n()),
    group      = dplyr::first(group),
    .groups    = "drop"
  ) %>%
  arrange(group, m_behave) %>%
  group_by(group) %>%
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%
  pivot_longer(cols = c(m_behave, m_model)) %>%
  mutate(
    name  = recode(name, m_behave = "Behavioral Proportions", m_model = "Model Predictions"),
    upper = ifelse(name == "Behavioral Proportions", value + sem_behave, value + sem_model),
    lower = ifelse(name == "Behavioral Proportions", value - sem_behave, value - sem_model),
    group = recode(group, control = "Control", patient = "Patient")
  ) %>%
  ggplot(aes(x = pro_id, color = name, y = value, ymax = upper, ymin = lower)) +
  geom_hline(yintercept = .5) +
  geom_point(size = 1, alpha = .5) +
  theme_classic() + ylim(0, 1) +
  labs(subtitle = "Maximize\nTreatment", x = "", y = "", color = "") +
  scale_color_manual(values = c("black", "slateblue1")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    text = element_text(size = 12)
  ) +
  guides(color = "none") +
  facet_wrap(~group, nrow = 2)



exp_ev_be_max_2 <- dat %>%
  filter(cond == "exp", !is.na(chosen_optim_be_ev)) %>%
  group_by(pro_id) %>%
  summarise(
    m_behave   = mean(chosen_optim_be_ev),
    sem_behave = sd(chosen_optim_be_ev) / sqrt(n()),
    m_model    = mean(modprob_optim_be_ev),
    sem_model  = sd(modprob_optim_be_ev) / sqrt(n()),
    group      = dplyr::first(group),
    .groups    = "drop"
  ) %>%
  arrange(group, m_behave) %>%
  group_by(group) %>%
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%
  pivot_longer(cols = c(m_behave, m_model)) %>%
  mutate(
    name  = recode(name, m_behave = "Behavioral Proportions", m_model = "Model Predictions"),
    upper = ifelse(name == "Behavioral Proportions", value + sem_behave, value + sem_model),
    lower = ifelse(name == "Behavioral Proportions", value - sem_behave, value - sem_model),
    group = recode(group, control = "Control", patient = "Patient")
  ) %>%
  ggplot(aes(x = pro_id, color = name, y = value, ymax = upper, ymin = lower)) +
  geom_hline(yintercept = .5) +
  geom_point(size = 1, alpha = .5) +
  theme_classic() + ylim(0, 1) +
  labs(subtitle = "", x = "Participant", y = "", color = "") +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 12)
  ) +
  guides(color = "none") +
  facet_wrap(~group, nrow = 2)


exp_se_raw_1 <- dat %>%
  filter(cond == "desc", !is.na(chosen_optim_se_raw)) %>%
  group_by(pro_id) %>%
  summarise(
    m_behave   = mean(chosen_optim_se_raw),
    sem_behave = sd(chosen_optim_se_raw) / sqrt(n()),
    m_model    = mean(modprob_optim_se_raw),
    sem_model  = sd(modprob_optim_se_raw) / sqrt(n()),
    group      = dplyr::first(group),
    .groups    = "drop"
  ) %>%
  arrange(group, m_behave) %>%
  group_by(group) %>%
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%
  pivot_longer(cols = c(m_behave, m_model)) %>%
  mutate(
    name  = recode(name, m_behave = "Behavioral Proportions", m_model = "Model Predictions"),
    upper = ifelse(name == "Behavioral Proportions", value + sem_behave, value + sem_model),
    lower = ifelse(name == "Behavioral Proportions", value - sem_behave, value - sem_model),
    group = recode(group, control = "Control", patient = "Patient")
  ) %>%
  ggplot(aes(x = pro_id, color = name, y = value, ymax = upper, ymin = lower)) +
  geom_hline(yintercept = .5) +
  geom_point(size = 1, alpha = .5) +
  theme_classic() + ylim(0, 1) +
  labs(subtitle = "Minimize\nSide Effect", x = "", y = "", color = "") +
  scale_color_manual(values = c("black", "slateblue1")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    text = element_text(size = 12)
  ) +
  guides(color = "none") +
  facet_wrap(~group, nrow = 2)

exp_se_raw_2 <- dat %>%
  filter(cond == "exp", !is.na(chosen_optim_se_raw)) %>%
  group_by(pro_id) %>%
  summarise(
    m_behave   = mean(chosen_optim_se_raw),
    sem_behave = sd(chosen_optim_se_raw) / sqrt(n()),
    m_model    = mean(modprob_optim_se_raw),
    sem_model  = sd(modprob_optim_se_raw) / sqrt(n()),
    group      = dplyr::first(group),
    .groups    = "drop"
  ) %>%
  arrange(group, m_behave) %>%
  group_by(group) %>%
  mutate(pro_id = factor(row_number())) %>%
  ungroup() %>%
  pivot_longer(cols = c(m_behave, m_model)) %>%
  mutate(
    name  = recode(name, m_behave = "Behavioral Proportions", m_model = "Model Predictions"),
    upper = ifelse(name == "Behavioral Proportions", value + sem_behave, value + sem_model),
    lower = ifelse(name == "Behavioral Proportions", value - sem_behave, value - sem_model),
    group = recode(group, control = "Control", patient = "Patient")
  ) %>%
  ggplot(aes(x = pro_id, color = name, y = value, ymax = upper, ymin = lower)) +
  geom_hline(yintercept = .5) +
  geom_point(size = 1, alpha = .5) +
  theme_classic() + ylim(0, 1) +
  labs(subtitle = "", x = "Participant", y = "", color = "") +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 12)
  ) +
  guides(color = "none") +
  facet_wrap(~group, nrow = 2)


fit_dat <- read.csv("data/fit/Prelec_parameters.csv")
plot_data <- fit_dat %>% select(iperf_de, iperf_ex)


ml1 = mean(plot_data$iperf_de)
legend_data <- data.frame(x = c(1, 1), y = c(99, 99), label = c("Description", "Experience"))

loo1 <- plot_data %>%
  arrange(iperf_de) %>%
  mutate(id = 1:length(iperf_de)) %>%
  ggplot(aes(x = id, y = iperf_de)) +
  geom_hline(yintercept = 0.5, color = "black") +
  geom_hline(yintercept = ml1, color = "grey", linewidth = 1.4, alpha = .5) +
  geom_point(aes(color = "Description"), size = 1, alpha = .5) +  # simulate color mapping
  geom_point(data = legend_data, aes(x = x, y = y, color = label)) +  # adds legend entries
  ylim(c(0, 1)) +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 12, color = "black"),
        plot.tag = element_text(size = 14, face = "bold")) +
  labs(x = "Participant",
       y = "LOO Accuracy",
       subtitle = "Approximate out of sample choice accuracy",
       tag = "A)") +
  scale_color_manual(name = "Condition", values = c("Description" = "slateblue1", "Experience" = "firebrick")) +
  theme(legend.position = "bottom") +  # move legend to bottom
  guides(color = "none")  +# enlarge legend points
  theme(legend.justification = c(0.3, 0))

loo1




ml2 = mean(plot_data$iperf_ex)

loo2 <- plot_data %>% 
  arrange(iperf_ex) %>% 
  mutate(id = 1:length(iperf_ex)) %>% 
  ggplot(aes(x = id, y = iperf_ex)) + 
  geom_hline(yintercept=0.5, color="black")+
  geom_hline(yintercept=ml2, color="grey", linewidth=1.4, alpha=.5) +
  geom_point( size = 1, alpha=.5,  color="firebrick") + 
  ylim(c(0,1)) +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),) + 
  theme(text=element_text(size=12, color="black")) +
  labs(x="Participant", y = "") 
loo2





rec_data = read.csv2('fit/parameter_recovery_medians.csv')

#
library(ggplot2)
library(patchwork)
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
  filter(parameter != "theta") %>% 
  mutate(condition = recode(condition, de = "Description", ex="Experience"),
         group = recode(group, pat = "Patient", con = "Control"),
         parameter = recode(parameter, gamma="𝛾", lambda="𝜆")) %>% 
  ggplot(aes(x=gen, y=rec, color=condition, shape=group)) + 
  geom_line(stat="smooth", method="lm", alpha=.5) +
  geom_point(alpha=.4) + 
  facet_wrap(~parameter, scales="free") +
  theme_classic() +
  theme(text=element_text(size=12)) +
  scale_color_manual(values=c("slateblue1", "firebrick")) +
  labs(x = "Generated", y = "Recovered", tag="B)") +
  theme(plot.tag = element_text(size = 14, face = "bold")) +
  guides(color="none", shape="none")


rec_pat <- rec_data %>% filter(group=="pat")

cors <- c(
  cor(rec_pat$lambda_de_gen, rec_pat$lambda_de_rec),
  cor(rec_pat$lambda_ex_gen, rec_pat$lambda_ex_rec),
  cor(rec_pat$gamma_de_gen, rec_pat$gamma_de_rec),
  cor(rec_pat$gamma_ex_gen, rec_pat$gamma_ex_rec)
)

rec_con <- rec_data %>% filter(group=="con")

cors <- c(
  cors,
  cor(rec_con$lambda_de_gen, rec_con$lambda_de_rec),
  cor(rec_con$lambda_ex_gen, rec_con$lambda_ex_rec),
  cor(rec_con$gamma_de_gen, rec_con$gamma_de_rec),
  cor(rec_con$gamma_ex_gen, rec_con$gamma_ex_rec)
)
quantile(cors)

legend_plot <- ggplot(
  data.frame(
    cond = c("Description", "Experience"),
    x = 0, y = 0
  ),
  aes(x, y, color = cond)
) +
  geom_point(size = 3, alpha=1) +
  scale_color_manual(
    values = c("Description" = "slateblue1",
               "Experience"  = "firebrick")
  ) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        text=element_text(size=12),
        legend.text=element_text(size=10),
        legend.key.size=unit(1, "cm")) +
  guides(color=guide_legend(nrow=1, override.aes=list(size=2, linewidth=2)),
         linetype=guide_legend(nrow=1, override.aes=list(size=2, linewidth=2))) +
  theme(plot.margin=margin(-15,0,0,0),
        legend.margin=margin(-15,0,0,0))


dummy_legend <- data.frame(
  label1 = c("Description", "Experience"),
  label2 = factor(c("Controls", "Patients"),
                  levels = c("Patients", "Controls")),
  x = c(0, 0),
  y = c(0, 0)
) %>%
  ggplot(aes(color = label1, shape = label2, x = x, y = y)) +
  geom_point(alpha = 0) +
  theme_void() +
  scale_color_manual(values = c("Description" = "slateblue1",
                                "Experience" = "firebrick")) +
  scale_shape_manual(values = c("Controls" = 16,
                                "Patients" = 17)) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",        # 👈 KEY: stack guides
    legend.title = element_blank(),
    text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, "cm"),
    plot.margin = margin(-15, 0, 0, 0),
    legend.margin = margin(-15, 0, 0, 0)
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(size = 3, alpha = 1)
    ),
    shape = guide_legend(
      nrow = 1,                     # 👈 also horizontal
      override.aes = list(size = 3, alpha = 1)
    )
  )


negative_height <- 0
just_label <- 4.75
expand_fc = 2


maximize_expected_value <- wrap_plots(
  exp_ev_max_1 + scale_x_discrete(expand=expansion(add = expand_fc)) +
    labs(tag = "") +
    theme(
      plot.tag = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.margin = margin(2, 5, 1, 5)
    ) + theme(axis.title.y = element_text(hjust =just_label)), plot_spacer(),
  exp_ev_max_2 +scale_x_discrete(expand=expansion(add = expand_fc))+
    labs(y = NULL) +
    theme(
      plot.margin = margin(1, 5, 2, 5)
    ),
  ncol = 1,
  heights = c(1,negative_height, 1)
)

minimize_side_effect <- ((exp_se_raw_1 +scale_x_discrete(expand=expansion(add = expand_fc)) +
                           theme(plot.margin = margin(5, 5, 0, 5))) / plot_spacer() /
  (exp_se_raw_2 +scale_x_discrete(expand=expansion(add = expand_fc)) +
     theme(plot.margin = margin(0, 5, 0, 5))) ) + plot_layout(heights=c(1,negative_height, 1))

minimize_side_effect_p <- ((exp_se_freq1 +scale_x_discrete(expand=expansion(add = expand_fc)) +
                             theme(plot.margin = margin(5, 5, 0, 5))) /
                             plot_spacer() /
  (exp_se_freq2 +scale_x_discrete(expand=expansion(add = expand_fc)) +
     theme(plot.margin = margin(0, 5, 0, 5))) ) + plot_layout(heights=c(1,negative_height, 1))


maximize_treatment <-( (exp_ev_be_max_1 +scale_x_discrete(expand=expansion(add = expand_fc)) +
                         theme(plot.margin = margin(5, 5, 0, 5))) / plot_spacer() /
  (exp_ev_be_max_2 +scale_x_discrete(expand=expansion(add = expand_fc)) +
     theme(plot.margin = margin(0, 5, 0, 5))) ) + plot_layout(heights=c(1,negative_height, 1))

# Make the final layout
top_row <- (maximize_expected_value | 
  minimize_side_effect |
  minimize_side_effect_p |
  maximize_treatment) +
  plot_layout(ncol = 4, widths = c(1, 1, 1, 1)) #& 
  #theme(plot.margin = margin(2, 5, 2, 5))  # tight margins between plots

pp_check_plot <- (top_row / legend_plot) +
  plot_layout(heights = c(1, 0.0001))
pp_check_plot



ggsave("pp_check.png", pp_check_plot, dpi=300, units="cm", width=16, height=13)




rec_loo_plot <- ((loo1 | loo2) /  rec_plot / dummy_legend) + plot_layout(heights=c(1,1,0.1))



ggsave("rec_loo.png", dpi=300, plot=rec_loo_plot, units="cm", width=12, height=10)

