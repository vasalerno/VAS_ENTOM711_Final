## CODE FOR VAS ENTOM 711 Final Project ##
# Coding assisted and troubleshooting by ChatGPT, annotations provided by VAS

#Load packages and check if packages are already installed
packages <- c(
  "tidyverse", "MASS", "brms", "lme4",
  "psych", "GGally", "corrplot", "DHARMa"
)

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p)
}

library(MASS)
library(tidyverse)
library(brms)
library(lme4)
library(psych)
library(GGally)
library(corrplot)
library(DHARMa)
library(ggrepel)

#Call in a randomized, but repeatable dataset base
set.seed(123)

#DATA SIMULATION - based on VAS actual osmiaCAM research


#set number of bees
n_bees <- 40 
#set time points per bee (times measured across research timeline)
n_time <- 20
#set amount of units
n_units <- 4

#create dataset of all set data simulations
data <- expand.grid(
  bee_id = factor(1:n_bees),
  time = 1:n_time
)

#assign bees to random units (this helps to make the simualted dataset random/unknown like it would be in reality)
bee_unit <- sample(1:n_units, n_bees, replace = TRUE)
#rematch bee-unt to dataset
data$unit <- factor(bee_unit[as.numeric(data$bee_id)])

#simulate temperature (sets base temp (15) and then has variation and random noise; also created different values for each unit to add realistic aspect)
data$temp <- 15 +
  10 * sin(data$time / 3) +
  rnorm(nrow(data), 0, 2) +
  as.numeric(data$unit)

View(data)

#create an individual effect so that each bee has its own "behavior"
bee_effect <- rnorm(n_bees, 0, 0.5)

#create latent behavior so that the dataset better represents hidden behaviors that would be at play in the actual experiment
latent_behavior <- 0.3 * data$temp +
  bee_effect[as.numeric(data$bee_id)] +
  rnorm(nrow(data), 0, 1)

#Activity is what I want to monitor, so this creates activity counts per bee also accounting for the latent behavior factor
data$activity <- rnbinom(nrow(data), mu = exp(1 + latent_behavior), size = 1)

#Sets provisioning behavior accounting for latent behavior and random noise
data$provisioning <- 5 + latent_behavior + rnorm(nrow(data), 0, 1)

#Sets nesting success as a binary value and then created a probabilty that it happens based on latent_behavior
prob_success <- plogis(-1 + 0.5 * latent_behavior)
data$nest_success <- rbinom(nrow(data), 1, prob_success)

#created a behavioral data set that uses the three behaviors I want to model and standardizes them
behav_data <- data %>%
  dplyr::select(activity, provisioning, nest_success) %>%
  scale() %>%
  as.data.frame()

# PRINCIPLE COMPONENT ANALYSIS (PCA)

# PCA to explore dimensionality
pca <- prcomp(
  behav_data,
  center = TRUE,
  scale. = TRUE
)

# Examine variance explained
summary(pca)

pca_scores <- as.data.frame(pca$x)
pca_scores$temp <- data$temp
pca_scores$unit <- data$unit

pca_loadings <- as.data.frame(pca$rotation)
pca_loadings$behavior <- rownames(pca_loadings)

# scale loadings for visuals
scale_factor <- min(
  diff(range(pca_scores$PC1)),
  diff(range(pca_scores$PC2))
)

pca_loadings$PC1 <- pca_loadings$PC1 * scale_factor
pca_loadings$PC2 <- pca_loadings$PC2 * scale_factor


plot_pca <- ggplot() +
  
  # PCA SCORES (observations)
  geom_point(
    data = pca_scores,
    aes(x = PC1, y = PC2, color = temp),
    alpha = 0.6,
    size = 2
  ) +
  
  scale_color_gradient(low = "steelblue", high = "red") +
  
  # PCA LOADINGS (variables)
  geom_segment(
    data = pca_loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 1,
    color = "black"
  ) +
  
  geom_text_repel(
    data = pca_loadings,
    aes(x = PC1, y = PC2, label = behavior),
    size = 4,
    color = "black"
  ) +
  
  theme_minimal(base_size = 13) +
  labs(
    title = "PCA Plot of Behaviors of Interest",
    x = "PC1",
    y = "PC2",
    color = "Temperature (C)"
  )

plot_pca

# Loadings
pca_loadings1 <- data.frame(
  behavior = rownames(pca$rotation),
  PC1 = pca$rotation[,1]
)

pca_loadings_plot <- ggplot(pca_loadings1, aes(x = behavior, y = PC1)) +
  geom_col() +
  labs(x = "Behavior", y = "PC1 Loading") +
  theme_minimal()

pca_loadings_plot

# Variance in PCA
pca_var <- data.frame(
  PC = paste0("PC", 1:length(pca$sdev)),
  Variance = (pca$sdev^2) / sum(pca$sdev^2)
)

pca_var_plot <- ggplot(pca_var, aes(x = PC, y = Variance)) +
  geom_col() +
  labs(y = "Proportion of Variance Explained") +
  theme_minimal()

pca_var_plot


# FACTOR ANALYSIS

#created a behavioral data set that uses the three behaviors I want to model and standardizes them
behav_data <- data %>%
  dplyr::select(activity, provisioning, nest_success) %>%
  scale() %>%
  as.data.frame()

# Run Exploratory Factor Analysis (EFA)

efa <- psych::fa(
  behav_data,
  nfactors = 1,
  rotate = "none",
  scores = "regression"
)

loadings_df <- data.frame(
  behavior = rownames(efa$loadings),
  loading = as.numeric(efa$loadings[,1])
)

efa_loadings <- data.frame(
  behavior = rownames(efa$loadings),
  loading = as.numeric(efa$loadings[,1])
)

efa_loadings_plot<- ggplot(efa_loadings, aes(x = behavior, y = loading)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = behavior, xend = behavior, y = 0, yend = loading)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = round(loading, 2))) +
  theme_minimal() +
  labs(title = "EFA Loadings: Latent Behavioral Structure",
       x = "Behavior",
       y = "Factor Loading")

efa_loadings_plot

#reformat the data to look at observed behavior and its associated latent score
data_long <- data %>%
  dplyr::select(activity, provisioning, nest_success, behavior_factor) %>%
  tidyr::pivot_longer(cols = c(activity, provisioning, nest_success))

behav_latent_plot <- ggplot(data_long, aes(x = behavior_factor, y = value)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  facet_wrap(~name, scales = "free_y") +
  theme_minimal() +
  labs(title = "EFA Factor vs Observed Behaviors",
       x = "Latent Behavior Factor",
       y = "Behavior Value")

behav_latent_plot

#VISUALIZING THE DATA

# Behavioral Relationships
GGally::ggpairs(
  data[, c("activity", "provisioning", "nest_success")]
)

corr <- cor(data[, c("activity", "provisioning", "nest_success")])
corrplot::corrplot(corr, method = "color")


