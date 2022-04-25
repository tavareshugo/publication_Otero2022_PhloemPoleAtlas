library(tidyverse)
library(readxl)
library(janitor)
library(patchwork)
library(lme4)

# Read data ---------------------------------------------------------------

read_metabolites <- function(x){
  # read info first rows
  sample_info <- read_excel(x) %>%
    slice(1:4) %>%
    select(-1) %>%
    mutate(field = c("condition", "dps", "genotype", "tissue")) %>%
    pivot_longer(-field, names_to = "sample", values_to = "value") %>%
    pivot_wider(names_from = "field", values_from = "value")

  # read peak area and add column names
  peak_area <- read_excel(x) %>%
    slice(-(1:4)) %>%
    rename(`metabolite` = `...1`) %>%
    pivot_longer(-metabolite, names_to = "sample", values_to = "peak_area") %>%
    mutate(peak_area = as.numeric(peak_area))

  # output joined tables
  peak_area %>% full_join(sample_info, by = "sample")
}

# read all files
met <- list.files("data/metabolites/",
                  pattern = ".xlsx", recursive = TRUE, full.names = TRUE) %>%
  map_dfr(read_metabolites)

# bit more tidying
met <- met %>%
  mutate(dps = str_remove(dps, "dps"),
         tissue = case_when(tissue == "L" ~ "Leaf",
                            tissue == "R" ~ "Root",
                            TRUE ~ NA_character_))

# remove data that cannot be interpreted due to signal saturation
met <- met %>%
  mutate(peak_area = ifelse(condition == "sucrose" & metabolite == "Sucrose" & tissue == "Root",
                            NA, peak_area))

# make two separate tables for conditions
met_nosuc <- met %>% filter(condition == "no sucrose") %>% drop_na(peak_area)
met_suc <- met %>% filter(condition == "sucrose") %>% drop_na(peak_area)

# mean-variance relationship
# log-scale is probably more adequate to model these skewed data
met %>%
  group_by(tissue, metabolite, condition, genotype, dps) %>%
  summarise(mean = mean(peak_area),
            sd = sd(peak_area),
            var = var(peak_area)) %>%
  ggplot(aes(mean, var)) +
  geom_point(aes(colour = genotype)) +
  scale_x_log10() + scale_y_log10()


# Fit full model ------------------------------------------------

# fit hierarchical model accounting for sample variation
fit <- lmer(log2(peak_area) ~ genotype*dps*metabolite*tissue*condition + (1|sample),
            data = met)

# residuals look good - choice to model on the log-scale definitely makes sense
# differences therefore represent log-fold-change
plot(fit)

# marginal mean difference per metabolite, day and tissue
fit_contrasts <- fit %>%
  emmeans::emmeans(~ genotype | metabolite + dps + tissue + condition) %>%
  emmeans::contrast("pairwise") %>%
  as_tibble() %>%
  mutate(padj = p.adjust(p.value, "fdr"))

# metabolites with some evidence of difference in each tissue
fit_contrasts %>%
  filter(padj < 0.05) %>%
  distinct(tissue, metabolite, condition) %>%
  count(tissue, condition)


# Visualise ---------------------------------------------------------------

# plot of means
met %>%
  filter(metabolite == "Sucrose") %>%
  ggplot(aes(dps, log2(peak_area))) +
  geom_pointrange(aes(colour = genotype),
                  stat = "summary", fun.data = mean_cl_boot) +
  geom_line(stat = "summary", fun = mean, aes(group = genotype, colour = genotype)) +
  facet_grid(condition ~ tissue) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) +
  labs(x = "Days post sowing", y = "Peak area", colour = "Genotype")

# plot for comparison between mutant and WT - no sucrose
p1_nosuc <- fit_contrasts %>%
  filter(condition == "no sucrose") %>%
  mutate(highlight = ifelse(padj < 0.05, estimate, NA),
         metabolite = str_replace(metabolite, " ", "\n")) %>%
  ggplot(aes(dps, estimate, ymin = estimate - 2*SE, ymax = estimate + 2*SE)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = tissue, group = tissue),
            position = position_dodge(width = 0.7)) +
  geom_linerange(aes(colour = tissue),
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(y = highlight, group = tissue), shape = "asterisk",
             position = position_dodge(width = 0.7)) +
  facet_wrap( ~ metabolite, nrow = 4) +
  scale_colour_manual(values = c("Root" = "#7570b3", "Leaf" = "#1b9e77")) +
  labs(x = "Days post-sowing", y = "log2(3papl/WT)") +
  theme_classic(base_size = 14)

# visualise some sugars
p2_nosuc <- fit %>%
  emmeans::emmeans(~ genotype | metabolite + dps + tissue + condition) %>%
  as_tibble() %>%
  filter(metabolite %in% c("Sucrose", "Fructose", "Glucose") & condition == "no sucrose") %>%
  ggplot(aes(dps, emmean)) +
  geom_line(aes(group = genotype, colour = genotype),
            size = 1) +
  geom_linerange(aes(colour = genotype, ymin = asymp.LCL, ymax = asymp.UCL),
                 alpha = 0.8, size = 2) +
  geom_point(data = met %>% filter(metabolite %in% c("Sucrose", "Fructose", "Glucose") & condition == "no sucrose"),
             aes(dps, log2(peak_area), colour = genotype),
             position = position_dodge(0.7),
             alpha = .8) +
  facet_grid(metabolite ~ tissue) +
  ggthemes::scale_colour_tableau() +
  theme_classic(base_size = 14) +
  labs(x = "Days post sowing", y = "log2(peak area)", colour = "Genotype")


# plot for comparison between mutant and WT - sucrose
p1_suc <- fit_contrasts %>%
  filter(condition == "sucrose") %>%
  mutate(highlight = ifelse(padj < 0.05, estimate, NA),
         metabolite = str_replace(metabolite, " ", "\n")) %>%
  ggplot(aes(dps, estimate, ymin = estimate - 2*SE, ymax = estimate + 2*SE)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = tissue, group = tissue),
            position = position_dodge(width = 0.7)) +
  geom_linerange(aes(colour = tissue),
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(y = highlight, group = tissue), shape = "asterisk",
             position = position_dodge(width = 0.7)) +
  facet_wrap( ~ metabolite, nrow = 4) +
  scale_colour_manual(values = c("Root" = "#7570b3", "Leaf" = "#1b9e77")) +
  labs(x = "Days post-sowing", y = "log2(3papl/WT)") +
  theme_classic(base_size = 14)

# confidence intervals for sucrose
p2_suc <- fit %>%
  emmeans::emmeans(~ genotype | metabolite + dps + tissue + condition) %>%
  as_tibble() %>%
  filter(metabolite %in% c("Sucrose", "Fructose", "Glucose") & condition == "sucrose") %>%
  ggplot(aes(dps, emmean)) +
  geom_line(aes(group = genotype, colour = genotype),
            size = 1) +
  geom_linerange(aes(colour = genotype, ymin = asymp.LCL, ymax = asymp.UCL),
                 alpha = 0.8, size = 2) +
  geom_point(data = met %>% filter(metabolite %in% c("Sucrose", "Fructose", "Glucose") & condition == "sucrose"),
             aes(dps, log2(peak_area), colour = genotype),
             position = position_dodge(0.7),
             alpha = .8) +
  facet_grid(metabolite ~ tissue) +
  ggthemes::scale_colour_tableau() +
  theme_classic(base_size = 14) +
  labs(x = "Days post sowing", y = "log2(peak area)", colour = "Genotype")

# save PDFs
pdf("documents/pdf for figures/metabolomics/metabolites_no_sucrose.pdf", width = 7.5, height = 6)
p1_nosuc
p2_nosuc
dev.off()

pdf("documents/pdf for figures/metabolomics/metabolites_sucrose.pdf", width = 7.5, height = 6)
p1_suc
p2_suc
dev.off()


# Fit model - no sucrose ------------------------------------------------

# fit hierarchical model accounting for sample variation
fit <- lmer(log2(peak_area) ~ genotype*dps*metabolite*tissue + (1|sample),
            data = met_nosuc)

# residuals look good - choice to model on the log-scale definitely makes sense
# differences therefore represent log-fold-change
plot(fit)

# marginal mean difference per metabolite, day and tissue
fit_contrasts <- fit %>%
  emmeans::emmeans(~ genotype | metabolite + dps + tissue) %>%
  emmeans::contrast("pairwise") %>%
  as_tibble() %>%
  mutate(padj = p.adjust(p.value, "fdr"))

# plot of means
met_nosuc %>%
  filter(metabolite == "Sucrose") %>%
  ggplot(aes(dps, log2(peak_area))) +
  geom_pointrange(aes(colour = genotype),
                  stat = "summary", fun.data = mean_cl_boot) +
  geom_line(stat = "summary", fun = mean, aes(group = genotype, colour = genotype)) +
  facet_grid( ~ tissue) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) +
  labs(x = "Days post sowing", y = "Peak area", colour = "Genotype")

fit %>%
  emmeans::emmeans(~ genotype | metabolite + dps + tissue) %>%
  as_tibble() %>%
  filter(metabolite == "Sucrose") %>%
  ggplot(aes(dps, emmean)) +
  geom_line(aes(group = genotype, colour = genotype),
            size = 1) +
  geom_linerange(aes(colour = genotype, ymin = asymp.LCL, ymax = asymp.UCL),
                 alpha = 0.8, size = 3) +
  geom_point(data = met_nosuc %>% filter(metabolite == "Sucrose"),
             aes(dps, log2(peak_area), colour = genotype),
             position = position_dodge(0.7),
             alpha = .8) +
  facet_grid( ~ tissue) +
  ggthemes::scale_colour_tableau() +
  theme_classic(base_size = 14) +
  labs(x = "Days post sowing", y = "log2(peak area)", colour = "Genotype")

# plot for comparison between mutant and WT
fit_contrasts %>%
  mutate(highlight = ifelse(padj < 0.05, estimate, NA)) %>%
  ggplot(aes(dps, estimate, ymin = estimate - 2*SE, ymax = estimate + 2*SE)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = tissue, group = tissue),
            position = position_dodge(width = 0.7)) +
  geom_linerange(aes(colour = tissue),
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(y = highlight, group = tissue), shape = "asterisk",
             position = position_dodge(width = 0.7)) +
  facet_wrap( ~ metabolite, nrow = 4) +
  scale_colour_manual(values = c("Root" = "#7570b3", "Leaf" = "#1b9e77")) +
  labs(x = "Days post-sowing", y = "log2(3papl/WT)") +
  theme_classic(base_size = 14)

# metabolites with some evidence of difference in each tissue
fit_contrasts %>%
  filter(padj < 0.05) %>%
  distinct(tissue, metabolite)


# Fit model - sucrose ------------------------------------------------

# fit hierarchical model accounting for sample variation
fit <- lmer(log2(peak_area) ~ genotype*dps*metabolite*tissue + (1|sample),
            data = met_suc)

# residuals look good - choice to model on the log-scale definitely makes sense
# differences therefore represent log-fold-change
plot(fit)

# marginal mean difference per metabolite, day and tissue
fit_contrasts <- fit %>%
  emmeans::emmeans(~ genotype | metabolite + dps + tissue) %>%
  emmeans::contrast("pairwise") %>%
  as_tibble() %>%
  mutate(padj = p.adjust(p.value, "fdr"))

# plot of means
met_suc %>%
  filter(metabolite == "Sucrose") %>%
  ggplot(aes(dps, log2(peak_area))) +
  geom_pointrange(aes(colour = genotype),
                  stat = "summary", fun.data = mean_cl_boot) +
  geom_line(stat = "summary", fun = mean, aes(group = genotype, colour = genotype)) +
  facet_grid( ~ tissue) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) +
  labs(x = "Days post sowing", y = "Peak area", colour = "Genotype")

fit %>%
  emmeans::emmeans(~ genotype | metabolite + dps + tissue) %>%
  as_tibble() %>%
  filter(metabolite == "Sucrose") %>%
  ggplot(aes(dps, emmean)) +
  geom_line(aes(group = genotype, colour = genotype),
            size = 1) +
  geom_linerange(aes(colour = genotype, ymin = asymp.LCL, ymax = asymp.UCL),
                 alpha = 0.8, size = 3) +
  geom_point(data = met_suc %>% filter(metabolite == "Sucrose"),
             aes(dps, log2(peak_area), colour = genotype),
             position = position_dodge(0.7),
             alpha = .8) +
  facet_grid( ~ tissue) +
  ggthemes::scale_colour_tableau() +
  theme_classic(base_size = 14) +
  labs(x = "Days post sowing", y = "log2(peak area)", colour = "Genotype")

# plot for comparison between mutant and WT
fit_contrasts %>%
  mutate(highlight = ifelse(padj < 0.05, estimate, NA)) %>%
  ggplot(aes(dps, estimate, ymin = estimate - 2*SE, ymax = estimate + 2*SE)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = tissue, group = tissue),
            position = position_dodge(width = 0.7)) +
  geom_linerange(aes(colour = tissue),
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(y = highlight, group = tissue), shape = "asterisk",
             position = position_dodge(width = 0.7)) +
  facet_wrap( ~ metabolite, nrow = 4) +
  scale_colour_manual(values = c("Root" = "#7570b3", "Leaf" = "#1b9e77")) +
  labs(x = "Days post-sowing", y = "log2(3papl/WT)") +
  theme_classic(base_size = 14)

# metabolites with some evidence of difference in each tissue
fit_contrasts %>%
  filter(padj < 0.05) %>%
  distinct(tissue, metabolite)






# DEPRECATED --------------------------------------------------------------

# separate testing for each group (lower statistical power)
tests <- met %>%
  filter(condition == "no sucrose") %>%
  mutate(value = log10(peak_area),
         dps = factor(dps)) %>%
  group_nest(condition, tissue, metabolite) %>%
  mutate(test = map(data, ~ lm(value ~ genotype*dps, data = .x)))

# residuals
tests %>%
  mutate(resid = map(test, resid)) %>%
  unnest(resid) %>%
  ggplot(aes(resid)) +
  geom_density(aes(group = interaction(condition, tissue, metabolite)))

# emmeans estimates
tests <- tests %>%
  mutate(contrast = map(test, ~ .x %>%
                          emmeans::emmeans(~ genotype | dps) %>%
                          emmeans::contrast("pairwise") %>%
                          as_tibble()))

