#script for mediation analysis of the effect of pressure and glove factors on information retention via mental effort

# Install required packages if not already installed
required_packages <- c("brms", "tidyverse", "tidybayes", "posterior")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load libraries
library(brms)
library(tidyverse)   # includes readr, dplyr, ggplot2, tidyr, tibble
library(tidybayes)
library(posterior)

# Read the data
df_ME <- read_csv("Data.csv") %>% 
  mutate(
    pressure = factor(if_else(Group %in% c(1,2), 1, 0), labels = c("NoPressure","Pressure")),
    gloves   = factor(if_else(Group %in% c(1,3), 1, 0), labels = c("NoGloves","Gloves"))
  )

#rename Dual Task into reactTime
df_ME <- df_ME %>% 
  rename(
    reactTime = `Dual Task (in seconds)`
  )

# Inspect the raw RT distribution
summary(df_ME$reactTime)
#view(df_ME)

#Prior‐only predictive check with lognormal likelihood (prior already adjusted here)
prior_only_lognorm <- brm(
  formula      = reactTime ~ pressure * gloves,
  data         = df_ME,
  family       = lognormal(), 
  prior = c(
    prior(normal(log(5.867), 0.3), class = "Intercept"),
    prior(normal(0, 0.3), class = "b"),
    prior(normal(0, 0.2), class = "sigma", lb = 0)
  ),
  sample_prior = "only",
  chains       = 2, iter = 2000, refresh = 0
)

#visualise prior
disPlotprior_ME<-pp_check(prior_only_lognorm, type = "stat", stat = "mean", resp = "reactTime")
#Add a theme layer to increase font size
disPlotprior_ME + theme(
  plot.title = element_text(size = 20),      # Main plot title
  axis.title = element_text(size = 18),      # Axis titles (e.g., "density")
  axis.text = element_text(size = 16),       # Axis labels (the numbers)
  legend.title = element_text(size = 20),    # Legend title
  legend.text = element_text(size = 18)      # Legend labels
)


#define priors. We center the intercept on the log of the empirical mean of RT reaction time
mu_rt <- mean(df_ME$reactTime, na.rm = TRUE)
priors_lognorm <- c(
  prior(normal(log(5.867), 0.3), class = "Intercept",resp = "reactTime"),  # log-RT intercept
  prior(normal(0, 0.3), class = "b",resp = "reactTime"), # slopes on log-RT
  prior(normal(0, 0.2), class = "sigma", lb = 0,resp = "reactTime"),
  
  prior(normal(1.5, 1), class = "Intercept", resp = "PostResults"),
  prior(normal(0,   1), class = "b",resp = "PostResults")
)

#define two linked regression models
m1 <- bf(
  reactTime ~ pressure + gloves + pressure:gloves,
  family = lognormal()
)

# Information retention (binomial outcome, predicted by both factors and mental effort)
m2 <- bf(
  Post_Results | trials(10) ~ pressure + gloves + pressure:gloves +  log(reactTime),
  family = binomial(link = "logit")
)

# default_prior(
#   m1 + m2 + set_rescor(FALSE),
#   data   = df_ME,
#   family = list(lognormal(), binomial("logit"))
# )

#fit the multivariate model
fit_mv_ME <- brm(
  m1 + m2 + set_rescor(FALSE),  # no residual correlation between outcomes
  data = df_ME,
  prior = priors_lognorm,
  family = list(lognormal(), binomial("logit")),
  chains = 4, iter = 4000, cores = 4,
  sample_prior = "no"  # use real data
)
summary(fit_mv_ME)

#visualise the posterior distribution
disPlotpost_ME<-pp_check(fit_mv_ME, type = "stat", stat = "mean", resp = "reactTime")
#Add a theme layer to increase font size
disPlotpost_ME + theme(
  plot.title = element_text(size = 20),      # Main plot title
  axis.title = element_text(size = 18),      # Axis titles (e.g., "density")
  axis.text = element_text(size = 16),       # Axis labels (the numbers)
  legend.title = element_text(size = 20),    # Legend title
  legend.text = element_text(size = 18)      # Legend labels
)

#convert the model output into a tidy draws_df format
posterior_ME <- as_draws_df(fit_mv_ME)
#view(posterior_ME)


# ##pressure##
# posterior_ME %>%
#   mutate(
#     indirect = b_reactTime_pressurePressure * b_PostResults_logreactTime,
#     direct = b_PostResults_pressurePressure,
#     total = indirect + direct
#   ) %>%
#   pivot_longer(cols = c(direct, indirect, total),
#                names_to = "Effect", values_to = "Estimate") %>%
#   ggplot(aes(x = Estimate, y = Effect)) +
#   stat_halfeye(.width = c(0.66, 0.95), fill = "lightblue") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   labs(
#     title = "Posterior Distributions of Mediation Effects",
#     subtitle = "Direct, Indirect, and Total Effects of Pressure",
#     x = "Effect size (log-odds)",
#     y = NULL
#   ) +
#   theme_minimal()

# store the effects in a tibble
indirect_pressure_ME <- posterior_ME$b_PostResults_logreactTime * posterior_ME$b_reactTime_pressurePressure
direct_pressure_ME   <- posterior_ME$b_PostResults_pressurePressure
total_pressure_ME    <- indirect_pressure_ME + direct_pressure_ME

posterior_ME_summary_pressure <- tibble(
  indirect = indirect_pressure_ME,
  direct   = direct_pressure_ME,
  total    = total_pressure_ME
)

posterior_ME_summary_pressure %>% 
  summarise(across(everything(), list(mean = mean, lower = ~quantile(.x, 0.025), upper = ~quantile(.x, 0.975))))

view(posterior_ME_summary_pressure)

#get % for mediation
mean(posterior_ME_summary_pressure$indirect < 0)  # For pressure


# ###Gloves effect###
# posterior_ME %>%
#   mutate(
#     indirect = b_reactTime_glovesGloves * b_PostResults_logreactTime,
#     direct = b_PostResults_glovesGloves,
#     total = indirect + direct
#   ) %>%
#   pivot_longer(cols = c(direct, indirect, total),
#                names_to = "Effect", values_to = "Estimate") %>%
#   ggplot(aes(x = Estimate, y = Effect)) +
#   stat_halfeye(.width = c(0.66, 0.95), fill = "lightblue") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   labs(
#     title = "Posterior Distributions of Mediation Effects",
#     subtitle = "Direct, Indirect, and Total Effects of Pressure",
#     x = "Effect size (log-odds)",
#     y = NULL
#   ) +
#   theme_minimal()

# store the effects in a tibble
indirect_gloves_ME <- posterior_ME$b_PostResults_logreactTime * posterior_ME$b_reactTime_glovesGloves
direct_gloves_ME   <- posterior_ME$b_PostResults_glovesGloves
total_gloves_ME    <- indirect_gloves_ME + direct_gloves_ME

posterior_ME_summary_gloves <- tibble(
  indirect = indirect_gloves_ME,
  direct   = direct_gloves_ME,
  total    = total_gloves_ME
)

posterior_ME_summary_gloves %>% 
  summarise(across(everything(), list(mean = mean, lower = ~quantile(.x, 0.025), upper = ~quantile(.x, 0.975))))

#get % for mediation
mean(posterior_ME_summary_gloves$indirect < 0)  # For pressure

