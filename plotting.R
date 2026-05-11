# libraries
library(genzplyr)
library(tidyverse)

df <- read.csv("outputs/robustness_summaries.csv") %>%
  pivot_longer(
    cols = -c(net_id, S, C),
    names_to = c(".value", "scenario"),
    names_pattern = "^(topo|dyn)_(.*)$")

ggplot(df) +
  geom_point(aes(x = topo,
                 y = dyn),
             alpha = 0.6,
             colour = "#EAAA00") +
  geom_abline(slope = 1,
              colour = "#A6192E") +
  facet_wrap(vars(scenario)) +
  xlim(0,0.5) +
  theme_classic()

ggplot(df %>%
         pivot_longer(-c(net_id, S, C, scenario),
                      names_to = "extinction")) +
  geom_boxplot(aes(y = scenario,
                   x = value,
                   colour = extinction)) +
  scale_colour_manual(values = c("topo" = "#046A38", "dyn" = "#FFB81C")) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

curves_df <- read.csv("outputs/extinction_curves.csv")

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
                 colour = type),
             alpha = 0.3,
             shape = 15,
             size = 0.5) +
  scale_colour_manual(values = c("topo" = "#046A38", "dyn" = "#FFB81C")) +
  guides(
    color = guide_legend(
      label.position = "top",
      override.aes = list(shape = 15, size = 5, alpha = 1)
    )
  ) +
  facet_wrap(vars(scenario)) +
  ylim(0, 1) + 
  theme_classic()

df %>%
  squad_up(scenario) %>%
  no_cap(mean_topo = mean(topo, na.rm = TRUE),
         mean_dynam = mean(dyn, na.rm = TRUE),
         mean_C = mean(C),
         mean_S = mean(S))

summary(df)

spp_df <- read.csv("outputs/species_metadata.csv")
