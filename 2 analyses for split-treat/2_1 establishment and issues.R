library(MASS)
library(tidyverse)
# library(stats)
# library(lme4)



an_dat <- read.csv("2 analyses for split-treat/curated data/an_dat_up.csv")


#### Sum of issue types and issues per treatment ---- 

sum(an_dat$failed) #23
sum(an_dat$established) #154

# 154/177 = 87.0 %

sum(an_dat$QL) # 9
sum(an_dat$DL) # 4
sum(an_dat$efb) # 19
sum(an_dat$weak) # 5

# 23+9+4+19+5 = 60

sum(an_dat$issues) # 60

# Most constrained dataset includes 117 colonies of the initial 177



#### Establishment rate by treatment ----

established_summary <- an_dat %>% 
  group_by(trt_label) %>% 
  summarize(n = n(),
            established_num = sum(established),
            established_p = established_num/n
  )

established_summary <- established_summary %>% 
  mutate(failed_num = n - established_num,
         failed_p = 1-established_p)

sum(established_summary$n)



#### Analyze establishment rate by treatment - GLM framework ----

m1 <- glm(established ~ trt_label, family = binomial, data = an_dat)
m0 <- glm(established ~ 1, family = binomial, data = an_dat)

summary(m1) 
summary(aov(m1)) # To generate an overall p value
# Results: 
# No significant difference in establishment rate across groups (F 6, 170 = 1.24, P=0.288)



#### Analyze issue rate by treatment - GLM framework ----



m3 <- glm(issues ~ trt_label, family = binomial, data = an_dat)

summary(aov(m3)) # To generate an overall p value
# No significant difference in occurrence of issues between treatments (F 6, 170 = 0.794, P=0.576)



#### Plot establishment rate by treatment ----


established_summary_piv <- established_summary %>% 
  pivot_longer(c(established_p, failed_p), names_to = "outcome", values_to = "proportion") 

established_summary_piv <- established_summary_piv %>% 
  mutate(
    outcome_label = if_else(outcome == "established_p", "Established",
                            if_else(outcome == "failed_p", "Did not establish", NA)),
    percent = 100*proportion
  )

# Relevel factors for plotting
established_summary_piv <- established_summary_piv %>%
  mutate(trt_label = fct_relevel(trt_label, "Control", "Apivar", "Amitraz EC", "OA Dribble", "5x OA Dribble", "OA Vapor", "HopGuard")
  )

established_summary_piv <- established_summary_piv %>%
  mutate(outcome_label = fct_relevel(outcome_label, "Did not establish", "Established")
  )


# Plot establish vs. not
established_summary_piv %>% 
  ggplot(aes(x = trt_label, y = percent, fill = outcome_label)) +
  geom_col(position = "stack", color = "black") +
  scale_colour_manual(
    limits = c("Did not establish", "Established"),
    labels = c("Did not establish", "Established"),
    values = c("white", "black"), 
    aesthetics = "fill"
  ) +
  scale_y_continuous(limits = c(0,100)) + 
  ylab("Percent of colonies") +
  theme_classic(base_size = 30) + 
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank(), legend.title = element_blank())


# ggsave("2 analyses for split-treat/outputs/Fig3_established_p_2025-04-05.png", width = 12, height = 8.625, units = "in")

# ggsave("2 analyses for split-treat/outputs/Fig3_established_p_2025-04-05.tiff", width = 12, height = 8.625, units = "in")





#### Plot colony issues ----


# Issues

issues_summary <- an_dat %>% 
  group_by(trt_label) %>% 
  summarize(n = n(),
            
            efb_num = sum(efb),
            efb_p = efb_num/n,
            
            weak_num = sum(weak),
            weak_p = weak_num/n,
            
            failed_num = sum(failed),
            failed_p = failed_num/n,
            
            QL_num = sum(QL),
            QL_p = QL_num/n,
            
            DL_num = sum(DL),
            DL_p = DL_num/n
            
  )







sum(an_dat$issues) #60 colonies had any issue
# i.e. 117 had no issue
# 117/177 i.e., 66.1 % of colonies had no issue

# Check if this lines up with the number of colonies for which I have Varroa data on day 77. It does. 



### Prep data for colony issues graphing


issues_summary_pivoted <- issues_summary %>% 
  pivot_longer(c(efb_p, weak_p, failed_p, QL_p, DL_p), names_to = "reason", values_to = "proportion") 


issues_summary_pivoted <- issues_summary_pivoted %>% 
  mutate(
    reason_label = if_else(reason == "efb_p", "EFB",
                           if_else(reason == "weak_p", "Weak",
                                   if_else(reason == "failed_p", "Did not establish",
                                           if_else(reason == "DL_p", "Drone layer",
                                                   if_else(reason == "QL_p", "Queenless", "NA"
                                                   )
                                                   
                                                   
                                           )))))

# Relevel factors for plotting

issues_summary_pivoted <- issues_summary_pivoted %>%
  mutate(trt_label = fct_relevel(trt_label, "Control", "Apivar", "Amitraz EC", "OA Dribble", "5x OA Dribble", "OA Vapor", "HopGuard"),
         reason_label = fct_relevel(reason_label,
                                    "Weak",
                                    "EFB",
                                    "Queenless",
                                    "Drone layer", 
                                    "Did not establish"
         )
  )







### Plotting colony issues


# Graph with queen issues, "weak" and "EFB"

issues_summary_pivoted %>% 
  ggplot(aes(x = trt_label, y = proportion, fill = reason_label)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(limits = c(0,1)) + 
  ylab("Proportion of colonies") +
  theme_classic(base_size = 30) + 
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank(), legend.title = element_blank())

# ggsave("2 analyses for split-treat/outputs/FigS3_CoD_2025-04-05.png", width = 12, height = 8.625, units = "in")
# ggsave("2 analyses for split-treat/outputs/FigS3_CoD_2025-04-05.tiff", width = 12, height = 8.625, units = "in")

