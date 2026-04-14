# libraries
library(genzplyr)
library(tidyverse)

df <- read.csv("outputs/robustness_summaries.csv") %>%
  pivot_longer(
    cols = -c(net_id, S, C),
    names_to = c(".value", "scenario"),
    names_pattern = "^(topo|dyn)_(.*)$")

ggplot(df) +
  geom_abline(slope = 1) +
  geom_point(aes(x = topo,
                 y = dyn),
             alpha = 0.7) +
  facet_wrap(vars(scenario)) +
  xlim(0,0.5) +
  theme_classic()

ggplot(df %>%
         pivot_longer(-c(net_id, S, C, scenario),
                      names_to = "extinction")) +
  geom_boxplot(aes(y = scenario,
                   x = value,
                   colour = extinction)) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

curves_df <- read.csv("outputs/extinction_curves.csv")

ggplot(curves_df) +
  geom_abline(slope = -1,
              intercept = 1) +
  geom_abline(slope = -1,
              intercept = 0.5,
              linetype = "dashed") +
  geom_point(aes(x = primary,
                 y = secondary,
                 colour = type),
             alpha = 0.2,
             shape = 15,
             size = 0.5) +
  facet_wrap(vars(scenario)) +
  ylim(0, 1) + 
  theme_classic()
