#Script for analysing the effects of the group on information retention directly (overlaps with mediation analysis to an extent)

# Install required packages if not already installed
required_packages <- c("brms", "tidyverse", "tidybayes", "ggplot2", "posterior")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load the libraries
library(brms)
library(tidyverse)  # Includes readr, dplyr, tidyr, ggplot2 etc.
library(tidybayes)
library(ggplot2)
library(posterior)

#Read data
df_GE <- read.csv("Data.csv")
#View(df_GE)

#Set factors
df_GE$Group <- factor(df_GE$Group)
df_GE$Group <- relevel(df_GE$Group, ref = "4")
levels(df_GE$Group)

# Simulated test for prior (uncomment if you want to test it again. The prior set here has already been adjusted)
prior_model_GE <- brm(
  Post_Results | trials(10) ~ Group,
  data = df_GE,
  family = binomial(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(1.5, 1.0), class = "Intercept")
  ),
  sample_prior = "only",
  chains = 4,
  iter = 4000,
  cores = 4
)

# Visualise the prior with posterior prediction check
disPlotprior<-pp_check(prior_model_GE, type = "bars", ndraws = 100)
#Add a theme layer to increase font size
disPlotprior + theme(
  plot.title = element_text(size = 20),      # Main plot title
  axis.title = element_text(size = 18),      # Axis titles (e.g., "density")
  axis.text = element_text(size = 16),       # Axis labels (the numbers)
  legend.title = element_text(size = 20),    # Legend title
  legend.text = element_text(size = 18)      # Legend labels
)

#fit a posterior model
posterior_model_GE <- brm(
  Post_Results | trials(10) ~ Group,
  data = df_GE,
  family = binomial(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    # Tuned intercept prior reflecting ~82% baseline recall
    prior(normal(1.5, 1), class = "Intercept")
  ),
  chains = 4,
  iter = 4000,
  cores = 4
)

#print summary of the posterior 
summary(posterior_model_GE)

#visualise Post
disPlotpost<-pp_check(posterior_model_GE, type = "bars", ndraws = 100)
#Add a theme layer to increase font size
disPlotpost + theme(
  plot.title = element_text(size = 20),      # Main plot title
  axis.title = element_text(size = 18),      # Axis titles (e.g., "density")
  axis.text = element_text(size = 16),       # Axis labels (the numbers)
  legend.title = element_text(size = 20),    # Legend title
  legend.text = element_text(size = 18)      # Legend labels
)


# print posterior distributions of groups
posterior_model_GE %>%
  spread_draws(b_Intercept, b_Group1, b_Group2, b_Group3) %>%
  mutate(
    p_Group1 = plogis(b_Intercept + b_Group1),
    p_Group2 = plogis(b_Intercept + b_Group2),
    p_Group3 = plogis(b_Intercept + b_Group3),
    p_Group4 = plogis(b_Intercept)
  ) %>%
  select(starts_with("p_")) %>%
  pivot_longer(cols = everything(), names_to = "Group", values_to = "Probability") %>%
  ggplot(aes(x = Probability, y = Group)) +
  stat_halfeye(.width = c(0.66, 0.95)) +
  labs(
    title = "Posterior Distributions of Group Probabilities",
    x = "Probability",
    y = NULL
  ) +
  theme_minimal(base_size = 16) +  # base font size for theme_minimal
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )


#convert the model output into a tidy draws_df format
posterior_GE <- as_draws_df(posterior_model_GE)

#calculating the posterior probability that the regression coefficient for Groups is less than 0
hypothesis(model, "Group1 < 0")

hypothesis(model, "Group2 < 0")

hypothesis(model, "Group3 < 0")

