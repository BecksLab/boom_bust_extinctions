# libraries
library(genzplyr)
library(tidyverse)
library(ggridges)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)

df <- read.csv("outputs/robustness_summaries.csv") %>%
  vibe_check(!starts_with(c("C", "S"))) %>%
  pivot_longer(
    cols = -c(net_id),
    names_to = c("extinction", "scenario"),
    names_pattern = "^(topo|dyn)_(.+)$")

ggplot(df) +
  geom_boxplot(aes(x = scenario,
                   y = value,
                   colour = extinction,
                   group = interaction(extinction, scenario))) +
  scale_colour_manual(values = c("topo" = "#046A38", 
                                 "dyn" = "#FFB81C")) +
  labs(y = "Robustness") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/robustness.png",
       width = 5000, 
       height = 3000, 
       units = "px", dpi = 700)

curves_df <- read.csv("outputs/extinction_curves.csv") %>%
  glow_up(extinction = factor(type,
                             levels = c("topo",
                                        "dyn")))

ggplot(curves_df) +
  geom_abline(slope = -1,
              intercept = 1,
              colour = "#A6192E") +
  geom_abline(slope = -1,
              intercept = 0.5,
              linetype = "dotted",
              colour = "#A6192E",
              alpha = 0.8) +
  geom_point(aes(x = primary,
                 y = secondary,
                 colour = extinction),
             alpha = 0.6,
             shape = 15,
             size = 0.5)  +
  scale_colour_manual(values = c("topo" = "#046A38", 
                                 "dyn" = "#FFB81C")) +
  guides(
    color = guide_legend(
      label.position = "top",)) +
  facet_wrap(vars(scenario)) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("figures/primVsec.png",
       width = 7500, 
       height = 6000, 
       units = "px", dpi = 700)


