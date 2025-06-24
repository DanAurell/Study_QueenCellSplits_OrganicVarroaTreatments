#### Setup ----

library(MASS)
library(tidyverse)
library(lme4)
library(emmeans)
library(ggpubr) # For combo plot
library(scales) # To decide on color scheme

#### Read in data ----
an_dat <- read.csv("2 analyses for split-treat/curated data/an_dat_up.csv")


#### Prep data for analysis ----


dat_devel <- an_dat %>% 
  rename(
    fob_44 = ch44_fob,
    fob_58 = ch58_fob,
    fob_72 = ch72_fob,
    fob_77 = frames_bees,
    fob_105 = post_quick_fob
  ) %>% 
  mutate(
    fob_77 = if_else(weak == 1, 0, fob_77),# Set FOB to 0 for weak colonies
    trt_label = fct_relevel(trt_label, "Control", "Apivar", "Amitraz EC", "OA Dribble", "5x OA Dribble", "OA Vapor", "HopGuard")
  )

#### Prep data for d77 FOB analysis ----
dat_devel1 <- dat_devel %>% 
  filter(established == 1,
         QL == 0,
         DL == 0, 
         efb == 0
  ) # Still includes 5 "weak" colonies

# Data prep to add sample size
dat_devel1_n <- dat_devel1 %>% 
  group_by(trt_label) %>% 
  summarise(n = n())

dat_devel1_n <- dat_devel1_n %>% 
  mutate(
    ssize = paste0("n=",n)
  )


#### Data prep for weight at honey harvest analysis ----
dat_devel2 <- dat_devel %>% 
  filter(established == 1,
         QL == 0,
         DL == 0, 
         efb == 0,
         weak == 0
  ) %>% 
  filter(is.na(post_issues))
# Excludes 5 "weak" colonies and 10 that had issues on day 103

dat_devel2_n <- dat_devel2 %>% 
  group_by(trt_label) %>% 
  summarise(n = n())

dat_devel2_n <- dat_devel2_n %>% 
  mutate(
    ssize = paste0("n=",n)
  )


#### Analyze day 77 Colony strength ----

m.bees1 <- lm(fob_77~trt_label, data = dat_devel1) 
m.bees1.0 <- lm(fob_77~1, data = dat_devel1) 
anova(m.bees1.0, m.bees1)
# No improvement to include a treatment effect (F 6,121 = 0.0784, p=0.9981)


# Send emmeans predictions to its own strength item
strength1 <- emmip(m.bees1, ~trt_label, type = "response", CIs = TRUE, plotit = FALSE)
# The results were so similar that I am very happy sticking with the simpler model (no random effect)

emplot_strength <- emmip(m.bees1, ~trt_label, type = "response", CIs = TRUE, plotit = FALSE)
emplot_strength <- emplot_strength %>% 
  rename(
    fob_77 = yvar
  )


emplot_strength <- emplot_strength %>% 
  mutate(
    trt_label = fct_relevel(trt_label, "Control", "Apivar", "Amitraz EC", "OA Dribble", "5x OA Dribble", "OA Vapor", "HopGuard")
  )

#### Analyze data for colony weight post honey production ----

m.weight2 <- lmer(post_kgs_C ~ trt_label + (1 | new_yard), data = dat_devel2)
m.weight2.0 <- lmer(post_kgs_C ~ 1 + (1 | new_yard), data = dat_devel2)
anova(m.weight2, m.weight2.0) # No effect of treatment (ChiSq 6 = 3.0146, p=0.807)

# Graph predictions based on m2 to estimate effect size

weight2 <- emmip(m.weight2, ~trt_label, type = "response", CIs = TRUE, plotit = FALSE)

emplot_weight <- emmip(m.weight2, ~trt_label, type = "response", CIs = TRUE, plotit = FALSE)
emplot_weight <- emplot_weight %>% 
  rename(
    post_kgs_C = yvar
  )


emplot_weight <- emplot_weight %>% 
  mutate(
    trt_label = fct_relevel(trt_label, "Control", "Apivar", "Amitraz EC", "OA Dribble", "5x OA Dribble", "OA Vapor", "HopGuard")
  )


#### Plot day 77 colony strength ----

# Color and treatment names (for horizontal layout)

paneA_color <- ggplot() +
  geom_boxplot(data = dat_devel1, aes(x = trt_label, y = fob_77, fill = trt_type), outlier.shape = NA, color = "black") +
  scale_colour_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "#F8766D", "#619CFF"), 
    aesthetics = "fill"
  ) +  
  geom_point(aes(x = trt_label, y = fob_77), shape=20, size=10, color="black", fill="black", data = emplot_strength) +
  geom_text(data = dat_devel1_n, aes(x = trt_label, y = 0, label = ssize), size = 5) +
  #xlab("Treatment") +
  ylab("Colony strength\n(frames of bees)") +
  theme_classic( base_size = 30) 

# No color
paneA_grey <- ggplot() +
  geom_boxplot(data = dat_devel1, aes(x = trt_label, y = fob_77, fill = trt_type), outlier.shape = NA, color = "black") +
  scale_colour_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "grey85", "grey85"), 
    aesthetics = "fill"
  ) +  
  geom_point(aes(x = trt_label, y = fob_77), shape=20, size=10, color="black", fill="black", data = emplot_strength) +
  geom_text(data = dat_devel1_n, aes(x = trt_label, y = 0, label = ssize), size = 5) +
  #xlab("Treatment") +
  ylab("Colony strength\n(frames of bees)") +
  theme_classic( base_size = 30) 



## Color plot with labels
paneA_color +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(), legend.position = "top")

# ggsave("2 analyses for split-treat/outputs/Fig4 labels 2025-04-06.png", width = 8.625, height = 8, units = "in")

## Color plot, no labels
paneA_color <- paneA_color +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(), legend.position = "none")


## Greyscale plot, no labels
paneA_grey <- paneA_grey +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(), legend.position = "none")


#### Plot post-harvest colony weight ----

# Color

paneB_color <- ggplot() +
  geom_boxplot(aes(x = trt_label, y = post_kgs_C, fill = trt_type), outlier.shape = NA, color = "black", data = dat_devel2) +
  geom_point(aes(x = trt_label, y = post_kgs_C), shape=20, size=10, color="black", fill="black", data = emplot_weight) +
  geom_text(data = dat_devel2_n, aes(x = trt_label, y = 40, label = ssize), size = 5) +
  scale_colour_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "#F8766D", "#619CFF"), 
    aesthetics = "fill"
  ) +  
  scale_y_continuous(breaks = c(40,60,80,100), limits = c(40,115)) +
  #xlab("Treatment") +
  ylab("\nHive weight (kg)") +
  theme_classic( base_size = 30) +
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank(), legend.position = "none")


# Greyscale

paneB_grey <- ggplot() +
  geom_boxplot(aes(x = trt_label, y = post_kgs_C, fill = trt_type), outlier.shape = NA, color = "black", data = dat_devel2) +
  geom_point(aes(x = trt_label, y = post_kgs_C), shape=20, size=10, color="black", fill="black", data = emplot_weight) +
  geom_text(data = dat_devel2_n, aes(x = trt_label, y = 40, label = ssize), size = 5) +
  scale_colour_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "grey85", "grey85"), 
    aesthetics = "fill"
  ) +  
  scale_y_continuous(breaks = c(40,60,80,100), limits = c(40,115)) +
  #xlab("Treatment") +
  ylab("\nHive weight (kg)") +
  theme_classic( base_size = 30) +
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank(), legend.position = "none")


#### Combine plot panes ----

# Color

ggarrange(paneA_color, paneB_color, ncol = 2, nrow = 1, labels = "AUTO", font.label = list(size = 30))

# ggsave("2 analyses for split-treat/outputs/Fig4_strength_weight_color 2025-04-06.png", width = 17.25, height = 8, units = "in")
# ggsave("2 analyses for split-treat/outputs/Fig4_strength_weight_color 2025-04-06.tiff", width = 17.25, height = 8, units = "in")

# Greyscale

ggarrange(paneA_grey, paneB_grey, ncol = 2, nrow = 1, labels = "AUTO", font.label = list(size = 30))

# ggsave("2 analyses for split-treat/outputs/Fig4_strength_weight_grey 2025-04-06.png", width = 17.25, height = 8, units = "in")
# ggsave("2 analyses for split-treat/outputs/Fig4_strength_weight_grey 2025-04-06.tiff", width = 17.25, height = 8, units = "in")


#### A final exploratory analysis ----

# Estimate how much more honey was produced in the effective treatments 
# than in Control and OA Vapor

unique(dat_devel2$trt_label)


dat_devel2 <- dat_devel2 %>% 
  mutate(
    effYN = case_when(
      trt_label == "Control" ~ "Not Effective",
      trt_label == "OA Vapor" ~ "Not Effective",
      trt_label == "Amitraz EC" ~ "Effective",
      trt_label == "OA Dribble" ~ "Effective",
      trt_label == "HopGuard" ~ "Effective",
      trt_label == "5x OA Dribble" ~ "Effective",
      trt_label == "Apivar" ~ "Effective"
    )
  )

m.weight.eff <- lmer(post_kgs_C ~ effYN + (1 | new_yard), data = dat_devel2)
m.weight.eff.0 <- lmer(post_kgs_C ~ 1 + (1 | new_yard), data = dat_devel2)
anova(m.weight.eff, m.weight.eff.0) # Effect of treatment was not significant (ChiSq 1 = 1.8299, p=0.1761)
summary(m.weight.eff) 
# Est gain 4.028 kg * 2.205 = 8.88 lb * $2.69 perlb = $23.89 from treating
# 2024 honey price at https://release.nass.usda.gov/reports/hony0325.pdf

confint(m.weight.eff) 
  # From a loss of 1.822262 to a gain of 9.797710 kg
  # From a loss of 4.018088 lbs to a gain of 21.60395 lbs

