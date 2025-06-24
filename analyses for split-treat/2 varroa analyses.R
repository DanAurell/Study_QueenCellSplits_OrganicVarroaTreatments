#### Setup ----


library(scales)
library(MASS)
library(tidyverse)
library(emmeans)
library(broom)
library(multcomp)
library(lubridate)
library(mgcv)
library(ggpubr) # For combo plot
library(gt)
library(lme4)


#### Read in data ----
an_dat <- read.csv("analyses for split-treat/curated data/an_dat_up.csv")

# For Varroa analyses
# Exclude any colonies that had EFB or queen issues or were removed due to weakness by day 77

dat <- an_dat %>% 
  filter(issues == 0) %>% 
  mutate(ch2_wash_bees = ch2_wash_workers + ch2_wash_drones + ch2_wash_queens,
         ch77_wash_bees = ch77_wash_workers + ch77_wash_drones + ch77_wash_queens,
         chpost_wash_bees = chpost_wash_workers + chpost_wash_drones + chpost_wash_queens
  )

# Impute number of bees in sample for two colonies on day 77
mean(an_dat$ch77_wash_bees, na.rm=T) # 362

dat <- dat %>% 
  mutate(
    ch77_wash_bees = ifelse(colony_num %in% c("22-078", "22-108"), 362, ch77_wash_bees)
  )

dat <- dat %>% 
  mutate(
    ch2_perc_var = 100*ch2_var_total/ch2_wash_bees,
    ch77_perc_var = 100*ch77_var_total/ch77_wash_bees,
    chpost_perc_var = 100*chpost_var_total/chpost_wash_bees
  )

# 117 of 177 colonies still in trial on d77: 66% of colonies.


dat <- dat %>%
  mutate(trt_label = fct_relevel(trt_label, "Control", "Apivar", "Amitraz EC", "OA Dribble",  "5x OA Dribble","OA Vapor", "HopGuard"))

# calculate the number of days between day 2 sampling and post sampling
dat <- dat %>% 
  mutate(
    num_days = mdy(chpost_date) - mdy(ch2_date),
    days_numeric = as.numeric(num_days)
  )

unique(dat$days_numeric) 
# 95 to 113 days between day 2 and "Harvest" mites samples
mean(unique(dat$days_numeric)) # Mean 102.75 ; round to 103 for predictions;
# The mean difference is 103 days between day 2 and "harvest", 
# that means the predictions for "harvest" are reasonable to make for "day 105"


#### Assess infestation of parent colonies ----

an_dat %>% 
  summarise(meanvar = mean(varroa_5cup),
            minvar =  min(varroa_5cup),
            maxvar =  max(varroa_5cup)
  )

#### Infestation of stratum 1 per yard ----

an_dat %>% 
  filter(stratum == 1) %>% 
  group_by(new_yard) %>% 
  summarise(minvar =  min(varroa_5cup),
            maxvar =  max(varroa_5cup)
  )

#### Amitraz resistance bioassay ----

# At Wire, the tests went over the three-hour time - don't use these

dat_amres <- an_dat %>% 
  filter(new_yard != "Wire") %>% 
  filter(!is.na(ch2_var_tray), !is.na(ch2_var_wash)) %>% 
  select(colony_num, new_yard, parent_id, ch2_var_tray, ch2_var_wash, ch2_var_total)

dat_amres <- dat_amres %>% 
  mutate(resistance = 100*ch2_var_wash/ch2_var_total)

mean(dat_amres$resistance) # Mean calculated at the colony level


#### Mean Varroa infestations on day 2, 77, harvest ----

dat %>% 
  summarise(mean = mean(ch2_perc_var))

dat %>% 
  summarise(mean = mean(ch77_perc_var))

dat %>% 
  summarise(mean = mean(chpost_perc_var))

dat %>%
  group_by(trt_label) %>% 
  summarise(mean2 = mean(ch2_perc_var),
            mean77 = mean(ch77_perc_var),
            meanpost = mean(chpost_perc_var)
  )


# Greatest numerical difference between 5x OA dribble and HopGuard

# 2.05/0.89 = 2.3 times as many mites in the HopGuard group as in the 5x OA dribble. 
# This justifies correcting for it statistically, 
# whether pre-treatment differences are statistically significant or not


#### Make long version of Varroa data ----

# Rename to allow pivoting using underscore as separator

dat <- dat %>% 
  rename(ch2_var.total = ch2_var_total,
         ch2_wash.bees = ch2_wash_bees,
         ch2_perc.var = ch2_perc_var,
         ch77_var.total = ch77_var_total,
         ch77_wash.bees = ch77_wash_bees,
         ch77_perc.var = ch77_perc_var,
         
         chpost_var.total = chpost_var_total,
         chpost_wash.bees = chpost_wash_bees,
         chpost_perc.var = chpost_perc_var
  )


dat_piv <- dat %>% 
  pivot_longer(c(ch2_var.total, ch2_wash.bees, ch2_perc.var, 
                 ch77_var.total, ch77_wash.bees, ch77_perc.var,
                 chpost_var.total, chpost_wash.bees, chpost_perc.var
  ),
  names_to = c("day", ".value"), 
  names_sep = "_"
  )

dat_piv <- dat_piv %>% mutate(post_issues = ifelse(is.na(post_issues), 0, post_issues))

dat_piv <- dat_piv %>% 
  filter(
    !((post_issues == 1) & (day == "chpost"))
  )
# excluded post data for 10 colonies that had issues, 
# just for post honey production time point


dat_piv <- dat_piv %>%
  mutate(trt_label = fct_relevel(trt_label, "Control", "Apivar", "Amitraz EC", "OA Dribble",  "5x OA Dribble","OA Vapor", "HopGuard"))

# Rename "day" values to 2, 77, Harvest

dat_piv <- dat_piv %>% 
  mutate(day = case_when(
    day == "ch2" ~ "2", 
    day == "ch77" ~ "77", 
    day == "chpost" ~ "Harvest"
  ))
           


#### Determine how many were above/below 3% threshold for each treatment each day


dat_piv <- dat_piv %>%
  mutate(above = if_else(perc.var >= 3, 1, 0))


# This above/below information is now saved in a df
above_df <- dat_piv %>% 
  group_by(day, trt_label) %>% 
  summarize(n=n(),
            above_num = sum(above),
            above_perc = 100*above_num/n,
            below_num = n-above_num,
            below_perc = 100-above_perc
  ) %>% 
  select(day, trt_label, below_perc, below_num, n, everything())


#### Sample size per date ----

ssize <- dat_piv %>% 
  group_by(trt_label, day) %>% 
  summarize(n=n())

ssize2 <- ssize
ssize2$day[ssize2$day == "Harvest"] <- 105
ssize2$day <- as.numeric(ssize2$day)

#### Supp Fig S2 - Graph Varroa levels on day 2, 77, post ----

dat_piv %>% 
  ggplot(aes(x=trt_label, y = perc.var, fill = trt_label)) +
  geom_boxplot() +
  geom_text(data = ssize, aes(x = trt_label, y = -1, label = n), size = 5) +
  geom_hline(aes(yintercept = 0), linewidth = 1) +
  ylab("Varroa infestation of\nadult bees (%)") +
  facet_wrap(~day) +
  theme_classic(base_size = 30) +  #base_size = 30
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank(), legend.position = "none")

ggsave("analyses for split-treat/outputs/SuppFig_S2_varroa_daily 2025-04-04.png", width = 12, height = 14, units = "in")


#### Assess pre-treatment differences statistically ----

# Included only non-issues colonies in reported analysis

m.day2 <- glm.nb(ch2_var.total~trt_label, data = dat)
summary(m.day2)

m.day2.0 <- glm.nb(ch2_var.total~1, data = dat)
summary(m.day2.0)

anova(m.day2, m.day2.0) 
# Results: Pre-treatment differences were not significant (ChiSq 6 = 4.023; P=0.674)

## Qualitatively same result if include all colonies
m.day2 <- glm.nb(ch2_var_total~trt_label, data = an_dat)
summary(m.day2)

m.day2.0 <- glm.nb(ch2_var_total~1, data = an_dat)
summary(m.day2.0)

anova(m.day2, m.day2.0) # Pre-treatment differences were not significant (P=0.84)


#### ANCOVA d77 Varroa infestation ----

# I considered including cohort (using column new_yard) as a random effect
  # To increase power to detect treatment effects
  # But often gave singular fits, and estimates are so similar without random effect
  # Decision: not to use random effect of cohort/new_yard

# Correcting for number of bees in sample by including it as an offset

m3 <- glm.nb(ch77_var.total ~ 
               trt_label + 
               ch2_perc.var +
               offset(log(ch77_wash.bees)), 
             data=dat)
summary(m3)

m0 <- glm.nb(ch77_var.total ~ 
               # trt_label + 
               ch2_perc.var + 
               offset(log(ch77_wash.bees)), 
             data=dat)
summary(m0)

anova(m3, m0) 
# Results: Treatment significantly influenced Varroa load on day 77 
# (ChiSq 6 = 37.27767, P=1.55447e-06)

#### Post Honey Production Dataset "post" ----

dat_post <- 
  dat %>% 
  filter(is.na(post_issues))




#### ANCOVA day 105 Varroa infestation ----

# Run models for post-honey Varroa
# The 4 experimental blocks had a staggered start, 
# but ended at nearly the same date (due to the synchronized honey harvest )
# To account for the different length of time that each yard was in the experiment, 
# I include the number of days as a predictor, 
# and make predictions normalized to ~103 days of interval between pre and post

m4 <- glm.nb(chpost_var.total ~ 
               trt_label + 
               ch2_perc.var + 
               days_numeric +
               offset(log(chpost_wash.bees)), 
             data=dat_post)

summary(m4)


m4.0 <- glm.nb(chpost_var.total ~ 
                 trt_label + 
                 ch2_perc.var + 
                 # days_numeric +
                 offset(log(chpost_wash.bees)), 
               data=dat_post)

summary(m4.0)

anova(m4, m4.0)
# Effect of days was significant p = 0.02; go with m4

m4.00 <- glm.nb(chpost_var.total ~ 
                  # trt_label + 
                  ch2_perc.var + 
                  days_numeric +
                  offset(log(chpost_wash.bees)), 
                data=dat_post)

summary(m4.00)
anova(m4, m4.00)
# Effect of treatment was significant on day 105 as well (ChiSq 6 = 26.6724; P=0.0001667661) 



#### Day 77 contrasts ----

K <- contrMat(table(dat$trt_label), type="Tukey") # all contrasts
K_ctl <- rbind(K[c(1:6),]) # Curious what happens if I just compare each treatment to control

# For reporting results, all contrasts, multiple comparison corrected
summary(glht(m3, linfct = mcp(trt_label = K)),
        test = adjusted("single-step"))

# Curious how much influence the p-value correction has
summary(glht(m3, linfct = mcp(trt_label = K)),
        test = adjusted("none"))

# I want to be conservative and make sure I'm not unnecessarily writing off OA vapor. 
# If I test each treatment against control, is OA still stat similar to control?
# Yes, OA is still not distinguishable from control.
summary(glht(m3, linfct = mcp(trt_label = K_ctl)),
        test = adjusted("single-step"))

# On day 77, all treatments except OA vapor had lower Varroa than Controls 

# Amitraz EC significantly lower than untreated Controls 
  # (z = -4.83; P < 0.001)
# Apivar lower (z = -3.95; P = 0.002)

# HopGuard significantly lower than untreated Controls (z = -4.15; P < 0.001)
# Also 5x OA dribble (z = -4.435; P < 0.001)
# Also OA dribble (z = -3.518; P = 0.008)
# Control vs. OA vapor was not significant (z = -1.49; P = 0.748)

# OA vapor had a higher Varroa infestation rate than Amitraz EC (z = 3.44; P = 0.010) 


# A - Control
# BC - Apivar
# C - AEC
# BC - OA dribble
# BC- 5x OA dribble
# AB - OA Vapor
# BC - HopGuard

cont77 <- summary(glht(m3, linfct = mcp(trt_label = K)),
                     test = adjusted("single-step"))

cont77 <- tidy(cont77)

cont77$adj.p.value <- sprintf("%.3f", cont77$adj.p.value)
cont77$statistic <- sprintf("%.2f", cont77$statistic)
# After the fact manually changed 0.000 to < 0.001



## Make the table of contrasts

cont77 <- cont77 %>% 
  select(contrast, statistic, adj.p.value)

cont_tbl <- gt(cont77)

cont_rtf <- cont_tbl %>%
  as_rtf()

my_conn <- file("analyses for split-treat/outputs/SuppTable_S3_2025-04-05.RTF", "w")
writeLines(cont_rtf, my_conn)
close(my_conn)







#### Day 105 contrasts ----

# For interest, but not reported in the paper

K <- contrMat(table(dat_post$trt_label), type="Tukey")

summary(glht(m4, linfct = mcp(trt_label = K)),
        test = adjusted("single-step"))



# On Day 105,

# AEC different from control
# OA Vapor different from AEC
# 5x OA dribble different from control
# OA Vapor different from 5x OA dribble
# OA Vapor different from OA dribble


# A - Control - 6
# ABC - Apivar 2.3
# C - AEC - 1.5
# BC - OA dribble - 2.2
# C - 5x OA Dribble - 1.9
# A - OA Vapor - 6.7
# ABC - HopGuard - 2.3






#### Day 77 Varroa predictions ----

pred_77 <- emmip(m3, ~ trt_label, 
                 type = "response", 
                 CIs = TRUE, 
                 plotit = FALSE, 
                 offset = log(100))
# this is based on the non-included variable (Varroa infestation on day 2) being set to its mean... 
# Predictions are now set to the value of offset = 100... (mites per 100 bees)

pred_77$day <- 77

# Check this against the actual means per group - Matches closely
dat %>% 
  group_by(trt_label) %>% 
  summarise(varroa = mean(ch77_perc.var))


#### Day 105 Varroa predictions ----

# Make sure these predictions are made based on the same value of the covariate
# To make sure that the predictions are not shifted up or down because 10 fewer colonies were used

mean(dat$ch2_perc.var, na.rm = TRUE) # 1.42 Varroa per 100 bees

pred_post <- emmip(m4, ~ trt_label, 
                   type = "response", 
                   CIs = TRUE, 
                   plotit = FALSE, 
                   offset = log(100), 
                   at = list(ch2_perc.var = 1.42)) # Predict based on 1.42 on ch2

pred_post$day <- 105

# Combine day 2, 77, 105 into one plottable data frame

pred_pre <- data.frame(pred_77$trt_label, yvar = c(rep(1.42, 7)), day = c(rep(2, 7)))
pred_pre <- pred_pre %>% 
  rename(trt_label = pred_77.trt_label)

pred <- bind_rows(pred_pre, pred_77, pred_post)

# Combine predictions and sample size for graphing

pred <- left_join(pred, ssize2) # Data for 117 colonies fed into Day 2

pred <- pred %>% 
  select(day, trt_label, n, yvar, LCL, UCL)

pred <- pred %>% 
  mutate(
    # trt_label_n = paste0(trt_label, " (n=", n, ")")     # Adding to legend
    ssize = paste0("n=",n)
  )

pred$order <- c(rep(1:7, 3))

# Add a trt_type column
pred <- pred %>% 
  mutate(
    trt_type = case_when(
      trt_label == "Control" ~ "Control",
      trt_label == "Apivar" ~ "Synthetic",
      trt_label == "Amitraz EC" ~ "Synthetic",
      trt_label == "OA Dribble" ~ "Organic",
      trt_label == "5x OA Dribble" ~ "Organic",
      trt_label == "OA Vapor" ~ "Organic",
      trt_label == "HopGuard" ~ "Organic"
    )
  )

# Add CLD to pred table
pred$cld <- c(rep("X", 7), c("A", "BC", "C", "BC", "BC", "AB", "BC"), c("A", "ABC", "C", "BC", "C", "A", "ABC"))


# Combine predictions and sample size for tables

pred2 <- left_join(pred, ssize2)

pred2 <- pred2 %>% 
  select(day, trt_label, n, yvar, LCL, UCL, -ssize)

pred2$yvar <- sprintf("%.2f", pred2$yvar)
pred2$LCL <- sprintf("%.2f", pred2$LCL)
pred2$UCL <- sprintf("%.2f", pred2$UCL)

pred2 <- pred2 %>% 
  mutate(
    Estimate = paste0(yvar, " [", LCL, ", ", UCL, "]")
  ) %>% 
  select(-c(yvar, LCL, UCL))

pred3<- pred2 %>% 
  group_by(day) %>% 
  arrange(trt_label, .by_group = TRUE) %>% 
  ungroup()



pred_tbl <- gt(pred3)

pred_rtf <- pred_tbl %>%
  as_rtf()

my_conn <- file("analyses for split-treat/outputs/SuppTable_S2_2025-04-05.RTF", "w")
writeLines(pred_rtf, my_conn)
close(my_conn)

#### Day 77 Efficacy ----

tidy(m3)
summary(m3)

eff77 <- tidy(m3)

eff77_CI <- as.data.frame(confint(m3))
eff77 <- cbind(eff77, eff77_CI)
rm(eff77_CI)
eff77 <- eff77 %>% 
  rename(LCL = "2.5 %", UCL = "97.5 %")

eff77 <- eff77 %>% 
  mutate(times_as = exp(estimate),
         efficacy = 100*(1-times_as),
         
         times_as_UCL = exp(UCL),
         efficacy_LCL = 100*(1-times_as_UCL),
         
         times_as_LCL = exp(LCL),
         efficacy_UCL = 100*(1-times_as_LCL)
  )

eff77 <- eff77[2:7,]

eff77$day <- 77


#### Day 105 Efficacy ----

effpost <- tidy(m4)
effpost_CI <- as.data.frame(confint(m4))

effpost <- cbind(effpost, effpost_CI)
rm(effpost_CI)
effpost <- effpost %>% 
  rename(LCL = "2.5 %", UCL = "97.5 %")

effpost <- effpost %>% 
  mutate(times_as = exp(estimate),
         efficacy = 100*(1-times_as),
         
         times_as_UCL = exp(UCL),
         efficacy_LCL = 100*(1-times_as_UCL),
         
         times_as_LCL = exp(LCL),
         efficacy_UCL = 100*(1-times_as_LCL)
  )

effpost <- effpost[2:7,]

effpost$day <- 105


#### Combine day 77 and day 105 efficacy ----

eff_all <- rbind(eff77, effpost)

eff_all <- eff_all %>% 
  mutate(
    trt_label = case_when(term == "trt_labelApivar" ~ "Apivar",
                          term == "trt_labelOA Dribble" ~ "OA Dribble", 
                          term == "trt_label5x OA Dribble" ~ "5x OA Dribble",
                          term == "trt_labelOA Vapor" ~ "OA Vapor",
                          term == "trt_labelHopGuard" ~ "HopGuard",
                          term == "trt_labelAmitraz EC" ~ "Amitraz EC", 
                          TRUE ~ "NA"
    )
  )

eff_all <- eff_all %>% 
  mutate(
    trt_label = fct_relevel(trt_label, "Apivar", "Amitraz EC", "OA Dribble", "5x OA Dribble", "OA Vapor", "HopGuard")
  )

# Add a trt_type column
eff_all <- eff_all %>% 
  mutate(
    trt_type = case_when(
      trt_label == "Apivar" ~ "Synthetic",
      trt_label == "OA Dribble" ~ "Organic",
      trt_label == "5x OA Dribble" ~ "Organic",
      trt_label == "OA Vapor" ~ "Organic",
      trt_label == "HopGuard" ~ "Organic",
      trt_label == "Amitraz EC" ~ "Synthetic"
    )
  )

# Add sample size for the days
eff_all <- left_join(eff_all, ssize2)

# Select the columns of interest to put in table
eff_all2 <- eff_all %>% 
  select(day, trt_label, n, efficacy, efficacy_LCL, efficacy_UCL)

eff_all2$efficacy <- sprintf("%.2f", eff_all2$efficacy)
eff_all2$efficacy_LCL <- sprintf("%.2f", eff_all2$efficacy_LCL)
eff_all2$efficacy_UCL <- sprintf("%.2f", eff_all2$efficacy_UCL)


eff_all2 <- eff_all2 %>% 
  mutate(
    Efficacy = paste0(efficacy, " [", efficacy_LCL, ", ", efficacy_UCL, "]")
  ) %>% 
  select(-c(efficacy, efficacy_LCL, efficacy_UCL))

eff_all2<- eff_all2 %>% 
  group_by(day) %>% 
  arrange(trt_label, .by_group = TRUE) %>% 
  ungroup()




eff_all_tbl <- gt(eff_all2)

eff_all_rtf <- eff_all_tbl %>%
  as_rtf()

my_conn <- file("analyses for split-treat/outputs/SuppTable_S4_2025-04-05.RTF", "w")
writeLines(eff_all_rtf, my_conn)
close(my_conn)

#### Make an efficacy object for Day 77 plotting ----

eff77_out <- eff_all %>% 
  filter(day == 77)


eff77_out <- eff77_out %>% 
  mutate(
    ssize = paste0("n=",n)
  )

# Add significance to eff77_3b object
eff77_out$sig <- c("*", "*", "*", "*", "N.S.", "*")




#### Fig 1 - Varroa infestation predictions bar plot ----

### Varroa predictions Day 2

Fig1_PaneA <- pred_pre %>% 
  filter(trt_label == "Control") %>% 
  ggplot(mapping = aes(y = yvar, x = "All treatments")) + # fill = trt_type # All treatments \n(n=117)
  geom_hline(yintercept = 3, linetype = 2, colour = "red") +
  geom_bar(stat = "identity", color = "black", fill = "grey85") +
  # scale_color_manual() +        # Put in control / synthetic / organic
  geom_text(aes(label = "n=117", y = -0.5), size = 5) + # 
  ylab("Varroa infestation of\nadult bees (%)") +
  ggtitle("Day 2") +
  scale_y_continuous(breaks = c(0,3,6,9, 12), limits = c(-0.5,12)) +
  theme_classic(base_size = 30) + # base_size = 30
  theme(axis.text.x=element_text(angle=45,hjust=1),  
        axis.title.x=element_blank(), 
        axis.title.y=element_text(size = 27), 
        legend.position = "none",
        plot.title = element_text(size = 25, face = "bold")
  )

Fig1_PaneA_grey <- pred_pre %>% 
  filter(trt_label == "Control") %>% 
  ggplot(mapping = aes(y = yvar, x = "All treatments")) + # fill = trt_type # All treatments \n(n=117)
  geom_hline(yintercept = 3, linetype = 2, colour = "grey") +
  geom_bar(stat = "identity", color = "black", fill = "grey85") +
  # scale_color_manual() +        # Put in control / synthetic / organic
  geom_text(aes(label = "n=117", y = -0.5), size = 5) + # 
  ylab("Varroa infestation of\nadult bees (%)") +
  ggtitle("Day 2") +
  scale_y_continuous(breaks = c(0,3,6,9, 12), limits = c(-0.5,12)) +
  theme_classic(base_size = 30) + # base_size = 30
  theme(axis.text.x=element_text(angle=45,hjust=1),  
        axis.title.x=element_blank(), 
        axis.title.y=element_text(size = 27), 
        legend.position = "none",
        plot.title = element_text(size = 25, face = "bold")
  )


### Varroa predictions Day 77

# mutate(trt_label_n = fct_reorder(trt_label, as.integer(order)))
  
Fig1_PaneB <- pred %>% 
  filter(day == 77) %>%
  ggplot(mapping = aes(y = yvar, x = trt_label, fill = trt_type)) + # fill = trt_type
  geom_hline(yintercept = 3, linetype = 2, colour = "red") +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0.5) +
  geom_text(aes(label = cld, y = UCL + 0.5), size = 5) +
  geom_text(aes(label = ssize, y = -0.5), size = 5) +
  scale_color_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "#F8766D", "#619CFF"), 
    aesthetics = "fill"
  ) +
  ylab("Varroa infestation of\nadult bees (%)") +
  ggtitle("Day 77") +
  scale_y_continuous(breaks = c(0,3,6,9,12), limits = c(-0.5,12)) +
  theme_classic(base_size = 30) + # base_size = 30
  theme(axis.text.x=element_text(angle=45,hjust=1),  
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 25, face = "bold")
  )


Fig1_PaneB_grey <- pred %>% 
  filter(day == 77) %>%
  ggplot(mapping = aes(y = yvar, x = trt_label, fill = trt_type)) + # fill = trt_type
  geom_hline(yintercept = 3, linetype = 2, colour = "grey") +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0.5) +
  geom_text(aes(label = cld, y = UCL + 0.5), size = 5) +
  geom_text(aes(label = ssize, y = -0.5), size = 5) +
  scale_color_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "grey85", "grey85"), 
    aesthetics = "fill"
  ) +
  ylab("Varroa infestation of\nadult bees (%)") +
  ggtitle("Day 77") +
  scale_y_continuous(breaks = c(0,3,6,9,12), limits = c(-0.5,12)) +
  theme_classic(base_size = 30) + # base_size = 30
  theme(axis.text.x=element_text(angle=45,hjust=1),  
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 25, face = "bold")
  )


### Varroa predictions Day 105

# mutate(trt_label_n = fct_reorder(trt_label, as.integer(order)))
  
  
Fig1_PaneC <- pred %>%
  filter(day == 105) %>%
  ggplot(mapping = aes(y = yvar, x = trt_label, fill = trt_type)) + # fill = trt_type
  geom_hline(yintercept = 3, linetype = 2, colour = "red") +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0.5) +
  geom_text(aes(label = cld, y = UCL + 0.5), size = 5) +
  geom_text(aes(label = ssize, y = -0.5), size = 5) +
  scale_color_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "#F8766D", "#619CFF"), 
    aesthetics = "fill"
  ) +
  ylab("Varroa infestation of\nadult bees (%)") +
  ggtitle("Day 105") +
  scale_y_continuous(breaks = c(0,3,6,9,12), limits = c(-0.5,12)) +
  theme_classic(base_size = 30) + # base_size = 30
  theme(
    axis.text.x=element_text(angle=45,hjust=1),  
    axis.title.x=element_blank(), 
    axis.title.y=element_blank(), 
    legend.position = "none",
    plot.title = element_text(size = 25, face = "bold")
  )


Fig1_PaneC_grey <- pred %>%
  filter(day == 105) %>%
  ggplot(mapping = aes(y = yvar, x = trt_label, fill = trt_type)) + # fill = trt_type
  geom_hline(yintercept = 3, linetype = 2, colour = "grey") +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0.5) +
  geom_text(aes(label = cld, y = UCL + 0.5), size = 5) +
  geom_text(aes(label = ssize, y = -0.5), size = 5) +
  scale_color_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("white", "grey85", "grey85"), 
    aesthetics = "fill"
  ) +
  ylab("Varroa infestation of\nadult bees (%)") +
  ggtitle("Day 105") +
  scale_y_continuous(breaks = c(0,3,6,9,12), limits = c(-0.5,12)) +
  theme_classic(base_size = 30) + # base_size = 30
  theme(
    axis.text.x=element_text(angle=45,hjust=1),  
    axis.title.x=element_blank(), 
    axis.title.y=element_blank(), 
    legend.position = "none",
    plot.title = element_text(size = 25, face = "bold")
  )

# Stitch together color figure

ggarrange(Fig1_PaneA, Fig1_PaneB, Fig1_PaneC, ncol = 3, nrow = 1, heights = c(1, 1, 1), widths = c(0.43, 1, 1), labels = NULL, font.label = list(size = 30))

ggsave("analyses for split-treat/outputs/Fig1_varroa_pred_2025-04-05.png", width = 17.25, height = 8, units = "in")
ggsave("analyses for split-treat/outputs/Fig1_varroa_pred_2025-04-05.tiff", width = 17.25, height = 8, units = "in")


# Stitch together greyscale figure

ggarrange(Fig1_PaneA_grey, Fig1_PaneB_grey, Fig1_PaneC_grey, ncol = 3, nrow = 1, heights = c(1, 1, 1), widths = c(0.43, 1, 1), labels = NULL, font.label = list(size = 30))

ggsave("analyses for split-treat/outputs/Fig1_varroa_pred_grey_2025-04-05.png", width = 17.25, height = 8, units = "in")
ggsave("analyses for split-treat/outputs/Fig1_varroa_pred_grey_2025-04-05.tiff", width = 17.25, height = 8, units = "in")



#### Fig 2 - Day 77 efficacy ----


# Applying this to efficacy on day 77

eff77_out %>%
  ggplot(aes(y = efficacy, x = trt_label, fill = trt_type)) + 
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=efficacy_LCL, ymax=efficacy_UCL), width=0.5) +
  geom_text(aes(label = sig, y = efficacy_UCL + 8), size = 10) +
  geom_text(aes(label = ssize, y = -10), size = 8) +
  scale_colour_manual(
    limits = c("Control", "Synthetic", "Organic"),
    labels = c("Control", "Synthetic", "Organic"),
    values = c("black", "#F8766D", "#619CFF"), 
    aesthetics = "fill"
  ) +  
  ylab("Efficacy (%)") +
  scale_y_continuous(breaks = c(-20, 0,20,40,60,80,100)) +
  theme_classic(base_size = 35) +
  
  theme(axis.text.x=element_text(angle=45,hjust=1),  axis.title.x=element_blank(), legend.position = "none")

# Save color version

ggsave("analyses for split-treat/outputs/Fig2_efficacy_2025-04-05.png", width = 12, height = 8.625, units = "in")
ggsave("analyses for split-treat/outputs/Fig2_efficacy_2025-04-05.tiff", width = 12, height = 8.625, units = "in")




# Greyscale version

eff77_out %>%
  ggplot(aes(y = efficacy, x = trt_label)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey85") +
  geom_errorbar(aes(ymin=efficacy_LCL, ymax=efficacy_UCL), width=0.5) +
  geom_text(aes(label = sig, y = efficacy_UCL + 8), size = 10) +
  geom_text(aes(label = ssize, y = -10), size = 8) +
  ylab("Efficacy (%)") +
  scale_y_continuous(breaks = c(-20, 0,20,40,60,80,100)) +
  theme_classic(base_size = 35) +
  
  theme(axis.text.x=element_text(angle=45,hjust=1),  axis.title.x=element_blank(), legend.position = "none")



# Save greyscale version

ggsave("analyses for split-treat/outputs/Fig2_efficacy_grey_2025-04-05.png", width = 12, height = 8.625, units = "in")
ggsave("analyses for split-treat/outputs/Fig2_efficacy_grey_2025-04-05.tiff", width = 12, height = 8.625, units = "in")



#### BACI analysis ----

# Suggested by a reviewer

# Data prep for BACI

# Reverse renaming

dat <- dat %>% 
  rename(ch2_var_total = ch2_var.total,
        ch2_wash_bees = ch2_wash.bees,
        ch2_perc_var = ch2_perc.var,
        ch77_var_total = ch77_var.total,
        ch77_wash_bees = ch77_wash.bees,
        ch77_perc_var = ch77_perc.var,
         
        chpost_var_total = chpost_var.total,
        chpost_wash_bees = chpost_wash.bees,
        chpost_perc_var = chpost_perc.var
  )

dat_slim <- dat %>% 
  select(colony_num,
         trt_label,
         ch2_var_total,
         ch2_wash_bees,
         ch2_perc_var,
         ch77_var_total,
         ch77_wash_bees,
         ch77_perc_var,
         issues) %>% 
  mutate(
    delta = ch77_perc_var - ch2_perc_var
  )
# dat excludes colonies with issues, so all the colonies have paired data for d2 and d77

hist(dat_slim$delta)




dat_piv_slim <- dat_piv %>% 
  select(colony_num,
         trt_label,
         day,
         var.total,
         wash.bees,
         perc.var,
         issues)


dat_piv_slim077 <- dat_piv_slim %>% 
  filter(day %in% c("2", "77"))


## Statistically compare delta between treatments

# "model BACI 1"
m.b1 <- lm(delta ~ trt_label, data = dat_slim)

TukeyHSD(aov(m.b1))
# Results: All treatments are different from the control, except OA Vapor
# No other treatments differ

# Qualititatively this is very similar to efficacy comparisons obtained previously.
# With other method, the only additional significant comparison: OA Vapor was worse than Amitraz EC


hist(resid(m.b1)) # Not normal residuals - but also not terrible


## Next step, fit negative binomial model, then compare differences between differences

# Modeled this on m3, which was used for main analysis


m.b2 <- glmer.nb(var.total ~ trt_label + day + trt_label:day +
                   (1|colony_num) +
                   offset(log(wash.bees)),
                 data = dat_piv_slim077
)
summary(m.b2)

m.b2.0 <- glmer.nb(var.total ~ trt_label + day +
                     (1|colony_num) +
                     offset(log(wash.bees)),
                   data = dat_piv_slim077
)
anova(m.b2, m.b2.0)


# Investigate convergence warning of model
# Based on https://rdrr.io/cran/lme4/man/allFit.html

## show available methods
allFit(show.meth.tab=TRUE) 
gm_all <- allFit(m.b2)
ss <- summary(gm_all)
ss$which.OK            ## logical vector: which optimizers worked?
## the other components only contain values for the optimizers that worked
ss$llik                ## vector of log-likelihoods
ss$fixef               ## table of fixed effects
ss$sdcor               ## table of random effect SDs and correlations
ss$theta               ## table of random effects parameters, Cholesky scale

# Set up the comparisons
emm3 <- emmeans(m.b2, ~ day * trt_label) 

# Perform the "contrast of contrasts" procedure
pairs(pairs(emm3, simple = "day"), simple = "trt_label") # Finally figured it out!!!

# I confirmed mathematically that these "contrasts of contrasts" line up with 
# how efficacy is calculated with the Henderson-Tilton formula

df <- tidy(pairs(pairs(emm3, simple = "day"), simple = "trt_label"))
df_CI <- confint(pairs(pairs(emm3, simple = "day"), simple = "trt_label"))

LCL <- df_CI$asymp.LCL
UCL <- df_CI$asymp.UCL

df2 <- data.frame(df, UCL, LCL)

names(df2)

df3 <- df2 %>% 
  mutate(times_as = exp(estimate),
         efficacy = 100*(1-times_as),
         
         times_as_UCL = exp(UCL),
         efficacy_LCL = 100*(1-times_as_UCL),
         
         times_as_LCL = exp(LCL),
         efficacy_UCL = 100*(1-times_as_LCL)
  )


# Filter based on grepl to include "Control"
df3 <- df3 %>% 
  select(contrast1, efficacy, efficacy_LCL, efficacy_UCL, everything()) %>% 
  filter(grepl("Control", contrast1))

# The estimates align well, but the confidence intervals are wider with the BACI analysis










