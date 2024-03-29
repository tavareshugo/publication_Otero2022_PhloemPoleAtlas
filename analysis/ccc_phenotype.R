library(tidyverse)
library(patchwork)
theme_set(theme_minimal(base_size = 14))

ccc <- read_csv("data/raw/triple_mutant_phenotype/compilation_no sucrose_exp_id_SO.csv")

# clean data
ccc <- ccc %>%
  mutate(days = as.numeric(str_remove(days, "dps")),
         comp = case_when(str_detect(genotype, "_comp_") ~ "Complemented",
                          genotype == "wt" ~ "WT",
                          genotype == "ccc" ~ "ccc",
                          TRUE ~ NA_character_),
         genotype_tidy = case_when(genotype == "ccc" ~ "triple",
                                   genotype == "wt" ~ "WT",
                                   TRUE ~ genotype)) %>%
  mutate(genotype_tidy = genotype_tidy %>% str_remove("ccc_comp_") %>% str_replace("_", "-"))

# experimental design is all over the place...
ccc %>%
  count(stock, genotype, experiment_id, days) %>%
  ggplot(aes(stock, genotype)) +
  geom_label(aes(label = n, fill = factor(days))) +
  facet_wrap(~ experiment_id) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Days")

# keep only samples measured at 6 days post sowing
ccc <- ccc %>% filter(days == 6)


# Analysis ----------------------------------------------------------------

# average of each genotype
ccc %>%
  filter(genotype %in% c("ccc", "wt")) %>%
  mutate(stock = fct_relevel(stock, "so488", "original")) %>%
  ggplot(aes(stock, root_length)) +
  geom_pointrange(stat = "summary", fun.data = mean_cl_boot,
                  aes(group = interaction(stock, experiment_id, genotype),
                      colour = experiment_id),
                  position = position_dodge(width = 0.7)) +
  facet_grid(~ genotype, scales = "free_x", space = "free_x") +
  ggthemes::scale_colour_tableau("Classic Green-Orange 12") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Seed stock", y = "Root Length (mm)", colour = "Experiment")

ccc %>%
  filter(genotype %in% c("ccc", "wt")) %>%
  mutate(strata = paste(stock, experiment_id, genotype, sep = ".")) %>%
  mutate(strata = fct_reorder(strata, root_length)) %>%
  ggplot(aes(strata, root_length)) +
  geom_pointrange(stat = "summary", fun.data = mean_cl_boot,
                  aes(colour = stock)) +
  facet_grid(~ genotype, scales = "free_x", space = "free_x") +
  ggthemes::scale_colour_tableau("Classic Green-Orange 12") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "Root Length (mm)", colour = "Seed Stock")


# all genotypes
ccc %>%
  mutate(stock = fct_relevel(stock, "so488", "original"),
         genotype = fct_relevel(genotype, "wt", "ccc")) %>%
  ggplot(aes(stock, root_length)) +
  geom_pointrange(stat = "summary", fun.data = mean_cl_boot,
                  aes(group = interaction(stock, experiment_id, genotype),
                      colour = experiment_id),
                  position = position_dodge(width = 0.7)) +
  # facet_grid(~ genotype, scales = "free_x", space = "free_x") +
  facet_wrap(~ genotype, scales = "free_x") +
  ggthemes::scale_colour_tableau("Color Blind") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Root Length (mm)", colour = "Line")


# Mixed Model -------------------------------------------------------------
# this is not ideal because of bimodal distributions

library(lme4)

# fit model
ccc$genotype <- fct_relevel(ccc$genotype, "wt", "ccc")
fit <- lmer(log(root_length) ~ genotype + (1|experiment_id:stock:genotype),
            data = ccc)
summary(fit)

marginal_means <- emmeans::emmeans(fit, ~ genotype, type = "response") %>%
  as_tibble()

marginal_means %>%
  mutate(genotype = fct_reorder(genotype, response)) %>%
  ggplot(aes(genotype)) +
  geom_pointrange(aes(y = response, ymin = lower.CL, ymax = upper.CL)) +
  labs(x = "Genotype", y = "Root Length (mm)") +
  coord_flip()


# Bootstrap differences ---------------------------------------------------
library(rsample)

# bootstrapped average of each genotype
resampled <- ccc %>%
  mutate(strata = paste(stock, experiment_id, genotype, sep = ".")) %>%
  select(genotype, root_length, strata) %>%
  bootstraps(times = 500) %>%
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) %>%
      group_by(strata, genotype) %>%
      summarise(mean_length = median(root_length)) %>%
      group_by(genotype) %>%
      summarise(mean_length = median(mean_length)) %>%
      ungroup() %>%
      mutate(mean_length = exp(mean_length))
  }))

resampled %>%
  unnest(boot) %>%
  group_by(genotype) %>%
  summarise(mid = median(mean_length),
            up = quantile(mean_length, 0.975),
            lo = quantile(mean_length, 0.025)) %>%
  ungroup() %>%
  mutate(genotype = fct_reorder(genotype, mid)) %>%
  ggplot(aes(genotype, mid)) +
  # geom_pointrange(stat = "summary", fun.data = median_hilow) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  labs(x = "Genotype", y = "Root Length (mm)") +
  coord_flip()


# Paper Figures -------------------------------------------------------------

# average per genotype
p1 <- resampled %>%
  unnest(boot) %>%
  group_by(genotype) %>%
  summarise(mid = median(mean_length),
            up = quantile(mean_length, 0.975),
            lo = quantile(mean_length, 0.025)) %>%
  ungroup() %>%
  mutate(genotype = fct_reorder(genotype, mid)) %>%
  ggplot(aes(genotype, mid)) +
  # geom_pointrange(stat = "summary", fun.data = median_hilow) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  labs(x = "Genotype", y = "Root Length (mm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# variation across experiments
p2 <- ccc %>%
  filter(genotype %in% c("ccc", "wt")) %>%
  mutate(stock = fct_relevel(stock, "so488", "original")) %>%
  ggplot(aes(stock, root_length)) +
  geom_pointrange(stat = "summary", fun.data = mean_cl_boot,
                  aes(group = interaction(stock, experiment_id, genotype),
                      colour = experiment_id),
                  position = position_dodge(width = 0.7)) +
  facet_grid( ~ genotype, scales = "free_x", space = "free_x") +
  ggthemes::scale_colour_tableau("Classic Green-Orange 12") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Seed stock", y = "Root Length (mm)", colour = "Experiment")

p3 <- ccc %>%
  filter(experiment_id %in% c("tech_rep_May", "techrep_2021_04_22") &
           genotype %in% c("ccc", "wt")) %>%
  group_by(genotype) %>% filter(n_distinct(experiment_id) == 2) %>% ungroup() %>%
  ggplot(aes(paste(genotype, stock), root_length)) +
  geom_violin(aes(fill = genotype), scale = "width") +
  ggbeeswarm::geom_quasirandom(aes(group = genotype), dodge.width = 1,
                               size = 0.5) +
  facet_grid(~ experiment_id, scales = "free_x", space = "free_x") +
  ggthemes::scale_fill_tableau("Tableau 10") +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "Root Length (mm)", colour = "Genotype")

pdf("documents/pdf for figures/ccc_phenotype.pdf", width = 7.5, height = 8.5)
(p1 / p2 / p3) + plot_annotation(tag_levels = "A")
dev.off()



# CRISPR mutant -----------------------------------------------------------

crispr <- read_csv("data/raw/triple_mutant_phenotype/crispr_combined.csv")

# bootstrapped average of each genotype
resampled <- crispr %>%
  mutate(strata = paste(condition, experiment, genotype, sep = ".")) %>%
  select(genotype, condition, root_length, strata) %>%
  bootstraps(times = 500) %>%
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) %>%
      group_by(strata, genotype, condition) %>%
      summarise(mean_length = mean(root_length), .groups = "drop") %>%
      group_by(genotype, condition) %>%
      summarise(mean_length = mean(mean_length), .groups = "drop") %>%
      ungroup() %>%
      mutate(mean_length = exp(mean_length))
  }))


p1 <- crispr %>%
  ggplot(aes(genotype, root_length)) +
  geom_violin(scale = "width") +
  ggbeeswarm::geom_quasirandom(alpha = 0.5) +
  facet_grid(experiment ~ condition) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# average per genotype/condition
p2 <- resampled %>%
  unnest(boot) %>%
  group_by(genotype, condition) %>%
  summarise(mid = median(mean_length),
            up = quantile(mean_length, 0.975),
            lo = quantile(mean_length, 0.025)) %>%
  ungroup() %>%
  # mutate(genotype = fct_relevel(genotype, "wt", "wt_so679")) %>%
  ggplot(aes(genotype, mid)) +
  # geom_pointrange(stat = "summary", fun.data = median_hilow) +
  geom_pointrange(aes(ymin = lo, ymax = up),
                  position = position_dodge(width = 0.5)) +
  facet_grid(~ condition) +
  labs(x = "Genotype", y = "Root Length (mm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 + ggtitle("Raw data distributions") +
  p2 + ggtitle("Overall average per genotype (with 95% CI)") +
  plot_layout(nrow = 2)

library(lme4)

# fit model
fit <- lmer(log(root_length) ~ genotype + (1|experiment),
            data = crispr %>% filter(condition == "0suc"))
summary(fit)

marginal_means <- emmeans::emmeans(fit, ~ genotype, type = "response") %>%
  as_tibble()

marginal_means %>%
  mutate(genotype = fct_reorder(genotype, response)) %>%
  ggplot(aes(genotype)) +
  geom_pointrange(aes(y = response, ymin = lower.CL, ymax = upper.CL)) +
  labs(x = "Genotype", y = "Root Length (mm)") +
  coord_flip()


# Both datasets ------------

crispr <- read_csv("data/raw/triple_mutant_phenotype/crispr_combined.csv")

# tidy to match ccc data
crispr <- crispr %>%
  mutate(genotype = ifelse(str_detect(genotype, "crispr"), str_replace(genotype, "_", "-"), genotype)) %>%
  separate(genotype, c("genotype", "stock"), sep = "_") %>%
  # mutate(stock = ifelse(is.na(stock), "new_stock", stock)) %>%
  mutate(stock = case_when(is.na(stock) & genotype == "ccc" ~ "unknown",
                           stock == "shared" & genotype == "ccc" ~ "shared_march21",
                           is.na(stock) & genotype == "wt" ~ "so656",
                           is.na(stock) ~ "new_stock",
                           TRUE ~ stock)) %>%
  select(genotype, root_length, stock, experiment_id = experiment, everything())

# bind the two datasets together
ccc_crispr <- crispr %>% filter(condition == "0suc") %>% select(-condition) %>%
  bind_rows(ccc %>% select(genotype, root_length, stock, experiment_id)) %>%
  # remove failed crispr mutant
  filter(genotype != "crispr_63")

# bootstrapped average of each genotype
resampled <- ccc_crispr %>%
  mutate(strata = paste(stock, experiment_id, genotype, sep = ".")) %>%
  select(genotype, root_length, strata) %>%
  bootstraps(times = 500) %>%
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) %>%
      group_by(strata, genotype) %>%
      summarise(mean_length = median(root_length)) %>%
      group_by(genotype) %>%
      summarise(mean_length = median(mean_length)) %>%
      ungroup() %>%
      mutate(mean_length = exp(mean_length))
  }))

# average per genotype
p1 <- resampled %>%
  unnest(boot) %>%
  group_by(genotype) %>%
  summarise(mid = median(mean_length),
            up = quantile(mean_length, 0.975),
            lo = quantile(mean_length, 0.025)) %>%
  ungroup() %>%
  mutate(genotype = fct_reorder(genotype, mid)) %>%
  ggplot(aes(genotype, mid)) +
  # geom_pointrange(stat = "summary", fun.data = median_hilow) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  labs(x = "Genotype", y = "Root Length (mm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# variation across experiments
p2 <- ccc_crispr %>%
  filter(genotype %in% c("ccc", "wt", "cdf2", "cog1-6cdf4", "crispr-61", "crispr-66")) %>%
  mutate(genotype = factor(genotype, levels = c("ccc",
                                             "crispr-66",
                                             "wt",
                                             "crispr-61",
                                             "cog1-6cdf4",
                                             "cdf2"))) %>%
  mutate(stock = fct_relevel(stock, "so488", "original")) %>%
  ggplot(aes(stock, root_length)) +
  geom_pointrange(stat = "summary", fun.data = mean_cl_boot,
                  aes(group = interaction(stock, experiment_id, genotype),
                      colour = experiment_id),
                  position = position_dodge(width = 0.7)) +
  facet_grid( ~ genotype, scales = "free_x", space = "free_x") +
  ggthemes::scale_colour_tableau("Classic Green-Orange 12") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Seed stock", y = "Root Length (mm)", colour = "Experiment")

pdf("documents/pdf for figures/ccc_crispr_phenotype.pdf", width = 10, height = 8.5)
(p1 / p2) + plot_annotation(tag_levels = "A")
dev.off()


# Repeat experiment -------------

ccc_repeat <- read_csv("data/raw/triple_mutant_phenotype/phloem_mutants_Rfile.csv")

resampled <- ccc_repeat %>%
  # filter(sucrose == "zero") %>%
  mutate(strata = paste(sucrose, replicate, genotype, sep = ".")) %>%
  select(genotype, sucrose, root_length, strata) %>%
  bootstraps(times = 500) %>%
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) %>%
      group_by(strata, genotype, sucrose) %>%
      summarise(mean_length = median(root_length)) %>%
      group_by(genotype, sucrose) %>%
      summarise(mean_length = median(mean_length)) %>%
      ungroup() %>%
      mutate(mean_length = exp(mean_length))
  }))

# average per genotype
p1 <- resampled %>%
  unnest(boot) %>%
  group_by(genotype, sucrose) %>%
  summarise(mid = median(mean_length),
            up = quantile(mean_length, 0.975),
            lo = quantile(mean_length, 0.025)) %>%
  ungroup() %>%
  mutate(genotype = factor(genotype, levels = c("wt",
                                                "3papl",
                                                "3papl-2",
                                                "apl",
                                                "p1p2",
                                                "6se")),
         sucrose = factor(sucrose, levels = c("sucrose", "zero"))) %>%
  ggplot(aes(genotype, mid, colour = sucrose)) +
  # geom_pointrange(stat = "summary", fun.data = median_hilow) +
  geom_pointrange(aes(ymin = lo, ymax = up),
                  position = position_dodge(0.5)) +
  labs(x = "Genotype", y = "Root Length (mm)", colour = "Sucrose") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim  = c(0, 3))

# variation across experiments
p2 <- ccc_repeat %>%
  mutate(genotype = factor(genotype, levels = c("wt",
                                                "3papl",
                                                "3papl-2",
                                                "apl",
                                                "p1p2",
                                                "6se")),
         sucrose = factor(sucrose, levels = c("sucrose", "zero"))) %>%
  ggplot(aes(replicate, root_length)) +
  geom_pointrange(stat = "summary", fun.data = mean_cl_boot,
                  aes(group = interaction(replicate, genotype, sucrose),
                      shape = replicate, colour = sucrose),
                  position = position_dodge(width = 0.5)) +
  facet_grid( ~ genotype, scales = "free_x", space = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Replicate", y = "Root Length (mm)", colour = "Sucrose") +
  coord_cartesian(ylim  = c(0, 3))

pdf("documents/pdf for figures/ccc_crispr_phenotype_repeat.pdf", width = 10, height = 8.5)
(p1 / p2) + plot_annotation(tag_levels = "A")
dev.off()


# Deprecated ---------------
library(rsample)

# This analysis is wrong because genotypes are not actually paired in a straightforward way
# leaving it here for reference, but this is deprecated

# assume that stock "so488" is the "original" stock - wrong assumption!
# retain only stocks that have several genotypes grown at the same time
ccc_filter <- ccc %>%
  mutate(stock = str_replace(stock, "so488", "original")) %>%
  filter(stock %in% c("febr_21", "flowering", "original")) %>%
  group_by(stock, experiment_id) %>%
  filter(any(genotype %in% "ccc") & any(genotype %in% "wt")) %>%
  ungroup()

# bootstrapped difference between each genotype and WT
resampled <- ccc_filter %>%
  mutate(strata = paste(stock, experiment_id, sep = ".")) %>%
  select(genotype_tidy, root_length, strata) %>%
  bootstraps(times = 500, strata = strata) %>%
  mutate(dif = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      group_by(strata, genotype_tidy) %>%
      summarise(median_length = median(root_length)) %>%
      group_by(strata) %>%
      mutate(dif = median_length - median_length[genotype_tidy == "WT"]) %>%
      ungroup() %>%
      filter(genotype_tidy != "WT")
  }))

# plot resampled differences
p1 <- resampled %>%
  unnest(dif) %>%
  separate(strata, c("stock", "experiment"), sep = "\\.", remove = FALSE) %>%
  mutate(genotype_tidy = str_replace(genotype_tidy, "papl1papl2cdf2", "triple")) %>%
  ggplot(aes(genotype_tidy, dif)) +
  geom_pointrange(stat = "summary", fun.data = median_hilow,
                  aes(group = strata, colour = stock),
                  position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0) +
  labs(x = "Genotype", y = "Median Difference to WT") +
  facet_wrap(~ experiment, nrow = 2) +
  scale_colour_brewer(palette = "Dark2")

# plot distributions
p2 <- ccc_filter %>%
  mutate(genotype_tidy = str_replace(genotype_tidy, "papl1papl2cdf2", "triple")) %>%
  ggplot(aes(x = genotype_tidy, y = root_length)) +
  geom_violin(aes(colour = stock), fill = "lightgrey", scale = "width") +
  ggbeeswarm::geom_quasirandom(aes(colour = stock),
                               dodge.width = 0.9, size = 0.5) +
  geom_linerange(stat = "summary", fun.data = mean_cl_boot,
                 aes(group = stock),
                 position = position_dodge(0.9),
                 size = 1) +
  facet_wrap(~ experiment_id, nrow = 2) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Genotype", y = "Root Length (mm)") +
  theme(legend.position = "none")

library(patchwork)
(p2 / p1) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme_classic() & scale_x_discrete(guide = guide_axis(n.dodge = 2))


# mixed model
library(lme4)

ccc_filter$genotype_tidy <- fct_relevel(ccc_filter$genotype_tidy, "WT", "papl1papl2cdf2")
fit <- lmer(log(root_length) ~ genotype_tidy + (genotype_tidy|experiment_id:stock), data = ccc_filter)
summary(fit)

emmeans::emmeans(fit, ~ genotype_tidy, type = "response") %>%
  as_tibble() %>%
  ggplot(aes(genotype_tidy)) +
  geom_linerange(aes(ymin = lower.CL, ymax = upper.CL)) +
  labs(x = "Genotype", y = "Relative Difference to WT") +
  scale_y_log10() +
  geom_hline(yintercept = 1)
