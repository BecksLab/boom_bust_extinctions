# libraries
library(genzplyr)
library(tidyverse)

df <- read.csv("outputs/paleo_robustness_summaries.csv") %>%
  pivot_longer(
    cols = -c(net_id, S_final, C_final, C_initial, S_initial, net_type),
    names_to = c(".value", "scenario"),
    names_pattern = "^(topo|dyn)_(.*)$")

ggplot(df) +
  geom_point(aes(x = topo,
                 y = dyn),
             alpha = 0.6,
             colour = "#EAAA00") +
  geom_abline(slope = 1,
              colour = "#A6192E") +
  facet_grid(cols = vars(scenario), 
             rows = vars(net_type)) +
  xlim(0,0.5) +
  ylim(0,0.5) +
  theme_classic()

ggplot(df %>%
         pivot_longer(-c(net_id, S_final, C_final, C_initial, S_initial, scenario, net_type),
                      names_to = "extinction")) +
  geom_boxplot(aes(y = net_type,
                   x = value,
                   colour = extinction)) +
  scale_colour_manual(values = c("topo" = "#046A38", "dyn" = "#FFB81C")) +
  coord_flip()+
  facet_wrap(vars(scenario)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df) +
  geom_boxplot(aes(y = S_final,
                   x = net_type),
               colour = "#EAAA00") +
  theme_classic()

curves_df <- read.csv("outputs/paleo_extinction_curves.csv") %>%
  glow_up(
    net_type = str_remove(type, "^[^_]+_"),
    type = str_extract(type, "^[^_]+")
  )

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
             alpha = 0.6,
             shape = 15,
             size = 0.5) +
  scale_colour_manual(values = c("topo" = "#046A38", "dyn" = "#FFB81C")) +
  guides(
    color = guide_legend(
      label.position = "top",
      override.aes = list(shape = 15, size = 5, alpha = 1)
    )
  ) +
  facet_grid(cols = vars(scenario), 
             rows = vars(net_type)) +
  ylim(0, 1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/paleo_primVsec.png",
       width = 8000, 
       height = 4000, 
       units = "px", dpi = 700)

df %>%
  squad_up(scenario, net_type) %>%
  no_cap(mean_topo = mean(topo, na.rm = TRUE),
         mean_dyn = mean(dyn, na.rm = TRUE),
         mean_C = mean(C_final),
         mean_S = mean(S_final))

summary(df)

spp_df <- read.csv("outputs/paleo_species_metadata.csv")

spp_df %>%
  squad_up(net_type, net_id, tropic_class) %>%
  tally() %>%
  pivot_wider(names_from = tropic_class,
              values_from = n) %>%
  glow_up(S = rowSums(pick(basal:intermediate), na.rm = TRUE),
          basal = basal/S,
          top = top/S,
          intermediate = intermediate/S)


df %>%
  vibe_check(-c(S_final, C_final, C_initial, S_initial)) %>%
  pivot_longer(cols = c(topo, dyn),
               names_to = "mechanism",
               values_to = "robustness") %>%
  pivot_wider(names_from = net_type,
              values_from = robustness) %>%
  ggplot() +
  geom_abline(slope = 1,
              intercept = 0,
              colour = "#A6192E") +
  geom_point(aes(x = niche,
                 y = down,
                 colour = mechanism)) +
  scale_colour_manual(values = c("topo" = "#046A38", "dyn" = "#FFB81C"))  +
  facet_grid(cols = vars(scenario),
             rows = vars(mechanism)) +
  theme_bw()

read.csv("outputs/extinction_curves.csv") %>%
  ggplot() +
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
             alpha = 0.6,
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
  theme_bw()
  
                 
                 
             