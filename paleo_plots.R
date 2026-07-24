# libraries
library(genzplyr)
library(tidyverse)
library(ggridges)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)

df <- read.csv("outputs/paleo_robustness_summaries_bodysize.csv") %>%
  vibe_check(!starts_with(c("C_", "S_"))) %>%
  pivot_longer(
    cols = -c(net_id, net_type, t_val),
    names_to = c(".value", "scenario"),
    names_pattern = "^(topo|dyn)_(.*)$") %>%
  na.omit()

ggplot(df %>%
         pivot_longer(-c(net_id, scenario, net_type, t_val),
                      names_to = "extinction")) +
  geom_boxplot(aes(x = net_type,
                   y = value,
                   colour = extinction)) +
  scale_colour_manual(values = c("topo" = "#046A38", "dyn" = "#FFB81C")) +
  facet_wrap(vars(scenario)) +
  labs(y = "Robustness") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

curves_df <- read.csv("outputs/paleo_extinction_curves_bodysize.csv") %>%
  glow_up(
    net_type = str_remove(type, "^[^_]+_"),
    type = str_extract(type, "^[^_]+")
  )

plot_primVsec <- 
  ggplot(curves_df %>%
           yeet(t_val == max(t_val))) +
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

plot_primVsec

ggsave("figures/paleo_primVsec.png",
       plot_primVsec,
       width = 5000, 
       height = 4000, 
       units = "px", dpi = 700)

spp_df <- read.csv("outputs/paleo_species_metadata_bodysize.csv")

plot_trophClass <-
  spp_df %>%
  squad_up(net_type, stage, net_id, trophic_class, t_val) %>%
  tally() %>%
  pivot_wider(names_from = trophic_class,
              values_from = n) %>%
  vibe_check(-isolated) %>%
  glow_up(S = rowSums(pick(basal:top), na.rm = TRUE),
          basal = basal/S,
          top = top/S,
          intermediate = intermediate/S) %>%
  pivot_longer(cols = c(basal, top, intermediate),
               names_to = "trophic_class",
               values_to = "percentage") %>%
  glow_up(percentage = if_else(is.na(percentage),
                               0, percentage),
          net_type = forcats::fct_rev(factor(net_type))) %>%
  ggplot() +
  geom_density_ridges(aes(x = percentage,
                          y = net_type,
                          fill = trophic_class),
                      alpha = 0.7) +
  facet_wrap(vars(stage),
             ncol = 2) +
  scale_fill_manual(values = c("basal" = "#046A38", 
                               "intermediate" = "#FFB81C",
                               "top" = "#A6192E")) +
  theme_bw()

plot_trophClass 

spp_df %>%
  yeet(stage != "creation")

df %>%
  vibe_check(net_id, net_type) %>%
  distinct() %>%
  squad_up(net_type) %>%
  no_cap(total_ids = n())

read.csv("outputs/paleo_robustness_summaries_bodysize.csv") %>%
  glow_up(spp_loss = S_creation - S_realised) %>%
  vibe_check(!starts_with(c("C_", "S_"))) %>%
  pivot_longer(
    cols = -c(net_id, net_type, spp_loss, t_val),
    names_to = c(".value", "scenario"),
    names_pattern = "^(topo|dyn)_(.*)$") %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = spp_loss,
                 y = dyn,
                 colour = net_type),
             alpha = 0.7) +
  facet_wrap(vars(scenario)) +
  labs(y = "Robustness (dynamic extinction)",
       x = "Number of species lost in burn-in") +
  theme_bw()

read.csv("outputs/paleo_robustness_summaries_bodysize.csv") %>%
  glow_up(spp_loss = S_creation - S_realised) %>%
  na.omit() %>%
  ggplot() +
  geom_boxplot(aes(y = spp_loss,
                   x = net_type,
                   colour = net_type)) +
  labs(y = "Number of species lost in burn-in") +
  theme_bw()

(plot_primVsec +
    labs(title = "Primary vs Secondary extintions, colours denote different extinction approaches")) /
  (plot_trophClass +
     labs(title = "Trophic classes of networks pre and post burn in")) +
  plot_layout(heights = c(2, 1))

ggsave("figures/paleo_panel.png",
       width = 5000, 
       height = 6000, 
       units = "px", dpi = 500)

read.csv("outputs/paleo_robustness_summaries_bodysize.csv") %>%
  vibe_check(C_creation, C_realised, net_type) %>%
  pivot_longer(-net_type) %>%
  ggplot() +
  geom_boxplot(aes(y = value,
                   x = name,
                   colour = net_type)) +
  labs(y = "Connectance",
       x = "Time state") +
  theme_bw()

read.csv("outputs/paleo_robustness_summaries_bodysize.csv") %>%
  vibe_check(S_creation, S_realised, net_type) %>%
  pivot_longer(-net_type) %>%
  ggplot() +
  geom_boxplot(aes(y = value,
                   x = name,
                   colour = net_type)) +
  labs(y = "Richness",
       x = "Time state") +
  theme_bw()

spp_df %>%
  ggplot() +
  geom_density_ridges(aes(x = biomass,
                          y = net_type,
                          fill = trophic_class),
                      alpha = 0.7) +
  scale_fill_manual(values = c("basal" = "#046A38", 
                               "intermediate" = "#FFB81C",
                               "top" = "#A6192E")) +
  theme_bw()

# 'lost' networks - i.e. didn't survive burn in
df %>%
  vibe_check(net_id, net_type) %>%
  distinct() %>%
  glow_up(max_run = max(net_id)) %>%
  squad_up(net_type) %>%
  count(net_type)

read.csv("outputs/paleo_robustness_summaries_bodysize.csv") %>%
  yeet(net_type == "niche") %>%
  ggplot() +
  geom_density(aes(C_target))


read.csv("outputs/paleo_robustness_summaries_bodysize.csv") %>%
  ggplot() +
  geom_point(aes(x=C_target,
                 y= C_creation,
                 colour=net_type))



