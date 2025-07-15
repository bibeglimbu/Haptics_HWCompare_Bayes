#script for mediation analysis of the effect of pressure and glove factors on information retention via perceived workload

# Install required packages if not already installed
required_packages <- c("brms", "tidyverse", "tidybayes", "posterior")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load libraries
library(brms)
library(tidyverse)
library(tidybayes)
library(posterior)

#read the data
df_PW <- read.csv("Data.csv")

#View(df)

df_PW$Group <- factor(df_PW$Group)
df_PW$Group <- relevel(df_PW$Group, ref = "4")
levels(df_PW$Group)


df_PW <- read_csv("Data.csv") %>% 
  mutate(
    pressure = factor(
      if_else(Group %in% c(1, 2), 1, 0),
      labels = c("NoPressure", "Pressure")
    ),
    gloves = factor(
      if_else(Group %in% c(1, 3), 1, 0),
      labels = c("NoGloves", "Gloves")
    )
  )
contrasts(df_PW$pressure)
contrasts(df_PW$gloves)

#Normalise Results column to Results_Rescaled to fit a beta model
df_PW$Result_rescaled <- df_PW$Result / 100

#Estimate a multivariate prior for the Perceived workload
fit_beta_prior <- brm(
  Result_rescaled ~ pressure * gloves,
  data = df_PW,
  family = Beta(),
  prior = c(
    prior(normal(0.4, 0.3), class = "Intercept"),
    prior(normal(0, 0.3), class = "b"),
    prior(gamma(30, 2), class = "phi")
  ),
  sample_prior = "only",
  chains = 4, iter = 4000, cores = 4
)

disPlotpri<-pp_check(fit_beta_prior, type = "dens_overlay", ndraws = 100)
#Add a theme layer to increase font size
disPlotpri + theme(
  plot.title = element_text(size = 20),      # Main plot title
  axis.title = element_text(size = 18),      # Axis titles (e.g., "density")
  axis.text = element_text(size = 16),       # Axis labels (the numbers)
  legend.title = element_text(size = 20),    # Legend title
  legend.text = element_text(size = 18)      # Legend labels
)

#We will simply reuse the prior from Bayes_GroupEffects for immediate recall
#define two linked regression models
m1 <- bf(
  Result_rescaled ~ pressure + gloves + pressure:gloves,
  family = Beta()
)

# Information retention (binomial outcome, predicted by both factors and perceived workload)
m2 <- bf(
  Post_Results | trials(10) ~ pressure + gloves + pressure:gloves + Result_rescaled,
  family = binomial(link = "logit")
)

priors <- c(
  # Priors for Result_rescaled (Beta model)
  prior(normal(0.4, 0.3), class = "Intercept", resp = "Resultrescaled"),
  prior(normal(0, 0.3), class = "b", resp = "Resultrescaled"),
  prior(gamma(30, 2), class = "phi", resp = "Resultrescaled"),
  
  # Priors for Post_Results (Binomial model)
  prior(normal(1.5, 1), class = "Intercept", resp = "PostResults"),
  prior(normal(0, 1), class = "b", resp = "PostResults")
)

#fit the multivariate regression model
fit_mv_PW <- brm(
  m1 + m2 + set_rescor(FALSE),  # no residual correlation between outcomes
  data = df_PW,
  prior = priors,
  chains = 4, iter = 4000, cores = 4,
  sample_prior = "no",  # now use real data
  set.seed(1234)
)
#print summary of the posterior 
summary(fit_mv_PW)

#visualise Post
disPlotpost<-pp_check(fit_mv_PW, type = "dens_overlay", resp = "Resultrescaled", ndraws = 100)
#Add a theme layer to increase font size
disPlotpost + theme(
  plot.title = element_text(size = 20),      # Main plot title
  axis.title = element_text(size = 18),      # Axis titles (e.g., "density")
  axis.text = element_text(size = 16),       # Axis labels (the numbers)
  legend.title = element_text(size = 20),    # Legend title
  legend.text = element_text(size = 18)      # Legend labels
)


#convert the model output into a tidy draws_df format
posterior_PW <- as_draws_df(fit_mv_PW)
#view(posterior_PW)

# ###Pressure effect###
#Visualise mediation effects (uncomment if you need to use)
# posterior_PW %>%
#   mutate(
#     indirect = b_Resultrescaled_pressurePressure * b_PostResults_Result_rescaled,
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
indirect_pressure_PW <- posterior_PW$b_PostResults_Result_rescaled * posterior_PW$b_Resultrescaled_pressurePressure
direct_pressure_PW   <- posterior_PW$b_PostResults_pressurePressure
total_pressure_PW    <- indirect_pressure + direct_pressure

posterior_PW_summary_pressure <- tibble(
  indirect = indirect_pressure_PW,
  direct   = direct_pressure_PW,
  total    = total_pressure_PW
)

posterior_PW_summary_pressure %>% 
  summarise(across(everything(), list(mean = mean, lower = ~quantile(.x, 0.025), upper = ~quantile(.x, 0.975))))
#get % for mediation
mean(posterior_PW_summary_pressure$indirect < 0)  # For pressure

###Gloves effect###
# posterior_PW %>%
#   mutate(
#     indirect = b_Resultrescaled_glovesGloves * b_PostResults_Result_rescaled,
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
indirect_gloves_PW <- posterior_PW$b_PostResults_Result_rescaled * posterior_PW$b_Resultrescaled_glovesGloves
direct_gloves_PW   <- posterior_PW$b_PostResults_glovesGloves
total_gloves_PW    <- indirect_gloves_PW + direct_gloves_PW

posterior_PW_summary_gloves <- tibble(
  indirect = indirect_gloves_PW,
  direct   = direct_gloves_PW,
  total    = total_gloves_PW
)

posterior_PW_summary_gloves %>% 
  summarise(across(everything(), list(mean = mean, lower = ~quantile(.x, 0.025), upper = ~quantile(.x, 0.975))))
#get % for mediation
mean(posterior_PW_summary_gloves$indirect < 0)  # For pressure
