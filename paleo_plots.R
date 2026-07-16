# libraries
library(genzplyr)
library(tidyverse)
library(ggridges)
library(patchwork)

df <- read.csv("outputs/paleo_robustness_summaries.csv") %>%
  vibe_check(!starts_with(c("C_", "S_"))) %>%
  pivot_longer(
    cols = -c(net_id, net_type, t_val),
    names_to = c(".value", "scenario"),
    names_pattern = "^(topo|dyn)_(.*)$")

ggplot(df) +
  geom_point(aes(x = topo,
                 y = dyn,
                 shape = factor(t_val)),
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
         pivot_longer(-c(net_id, scenario, net_type, t_val),
                      names_to = "extinction") %>%
         yeet(extinction == "dyn")) +
  geom_boxplot(aes(x = net_type,
                   y = value,
                   colour = factor(t_val))) +
  scale_colour_manual(values = c("50" = "#046A38", "500" = "#FFB81C")) +
  facet_wrap(vars(scenario)) +
  labs(y = "Robustness") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

curves_df <- read.csv("outputs/paleo_extinction_curves.csv") %>%
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

spp_df <- read.csv("outputs/paleo_species_metadata.csv")

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
  facet_grid(cols = vars(t_val),
             rows = vars(stage)) +
  scale_fill_manual(values = c("basal" = "#046A38", 
                               "intermediate" = "#FFB81C",
                               "top" = "#A6192E")) +
  theme_bw()

plot_trophClass 

df %>%
  vibe_check(net_id, net_type) %>%
  distinct() %>%
  squad_up(net_type) %>%
  tally()

read.csv("outputs/paleo_robustness_summaries.csv") %>%
  glow_up(spp_loss = S_creation - S_realised) %>%
  vibe_check(!starts_with(c("C_", "S_"))) %>%
  pivot_longer(
    cols = -c(net_id, net_type, spp_loss, t_val),
    names_to = c(".value", "scenario"),
    names_pattern = "^(topo|dyn)_(.*)$") %>%
  ggplot() +
  geom_point(aes(x = spp_loss,
                 y = dyn,
                 colour = net_type,
                 size = t_val),
             alpha = 0.7) +
  facet_wrap(vars(scenario)) +
  labs(y = "Robustness (dynamic extinction)",
       x = "Number of species lost in burn-in") +
  scale_colour_manual(values = c("#A4A6DE", "#35365E", "#DDCBA4", "#897D63", "pink")) +
  theme_bw()

(plot_primVsec +
  labs(title = "Primary vs Secondary extintions, colours donate different extinction approaches")) /
  (plot_trophClass +
  labs(title = "Trophic classes of networks pre and post burn in. Burn-in/extinction time length is shown as columns")) +
  plot_layout(heights = c(2, 1))

ggsave("figures/paleo_panel.png",
       width = 5000, 
       height = 6000, 
       units = "px", dpi = 500)
