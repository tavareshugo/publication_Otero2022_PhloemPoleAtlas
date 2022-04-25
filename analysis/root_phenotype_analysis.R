library(tidyverse)
library(readxl)
library(rsample)

set.seed(20220425) # for reproducibility

# Figure S12a ----

# read data
s12a <- read_csv("data/raw/triple_mutant_phenotype/Re _Your_manuscript,_NPLANTS-210711374C/FigS12a.csv")

# tidy
s12a <- s12a |>
  rename(condition = characteristic) |>
  mutate(condition = case_when(condition == "0suc" ~ "0% sucrose",
                               condition == "suc05" ~ "0.5% sucrose",
                               condition == "suc1%" ~ "1% sucrose",
                               condition == "suc5" ~ "5% sucrose",
                               condition == "24hours" ~ "24 hours",
                               TRUE ~ NA_character_)) |>
  mutate(condition = factor(condition,
                            levels = c("1% sucrose",
                                       "0.5% sucrose",
                                       "0% sucrose",
                                       "24 hours")))

# plot in the figure
s12a |>
  ggplot(aes(condition, root_length)) +
  geom_violin(aes(fill = genotype), scale = "width") +
  ggbeeswarm::geom_quasirandom(aes(group = genotype),
                               dodge.width = 1) +
  theme_classic()

# contrast WT vs Mutant within each treatment (bootstrap)
s12a_contrast <- s12a |>
  # mutate(strata = interaction(genotype, condition)) |>
  # bootstraps(times = 500, strata = strata) |>
  bootstraps(times = 500) |>
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) |>
      group_by(genotype, condition) |>
      summarise(mean_length = mean(root_length), .groups = "drop") |>
      pivot_wider(names_from = genotype, values_from = mean_length) |>
      mutate(dif = wt - ccc)
  }))

# plot bootstrapped differences and calculate bootstrap p-value
s12a_contrast |>
  unnest(boot) |>
  ggplot(aes(dif)) +
  geom_density(aes(colour = condition)) +
  geom_vline(xintercept = 0)

s12a_contrast_summary <- s12a_contrast |>
  unnest(boot) |>
  group_by(condition) |>
  summarise(mid = median(dif),
            lo = quantile(dif, 0.025),
            hi = quantile(dif, 0.975),
            pval = (sum(abs(dif) < abs(dif - mean(dif))) + 1)/(n() + 1),
            .groups = "drop") |>
  # convert back to wt/mutant ratio for more interpretable values
  mutate(mid = exp(mid),
         lo = exp(lo),
         hi = exp(hi)) |>
  mutate(padj = p.adjust(pval, method = "bonferroni"))

# write to temporary file to compile into supplementary spreadsheet
write_csv(s12a_contrast_summary, "~/temp/figS12a.csv")


# Fig S12b -----

# read data
s12b <- read_csv("data/raw/triple_mutant_phenotype/Re _Your_manuscript,_NPLANTS-210711374C/FigS12b.csv")

# tidy
s12b <- s12b |>
  mutate(condition = as.factor(as.numeric(str_remove(days, "days")))) |>
  select(-days)

# contrast WT vs Mutant within each treatment (bootstrap)
s12b_contrast <- s12b |>
  filter(condition == 6) |>
  # mutate(strata = interaction(genotype, condition)) |>
  # bootstraps(times = 500, strata = strata) |>
  bootstraps(times = 500) |>
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) |>
      group_by(genotype, condition) |>
      summarise(mean_length = mean(root_length), .groups = "drop") |>
      pivot_wider(names_from = genotype, values_from = mean_length) |>
      pivot_longer(cols = c(-wt, -condition),
                   names_to = "genotype",
                   values_to = "root_length") |>
      mutate(dif = wt - root_length)
  }))

# plot bootstrapped differences and calculate bootstrap p-value
s12b_contrast |>
  unnest(boot) |>
  ggplot(aes(dif)) +
  geom_density(aes(colour = genotype)) +
  geom_vline(xintercept = 0)

s12b_contrast_summary <- s12b_contrast |>
  unnest(boot) |>
  group_by(genotype) |>
  summarise(mid = median(dif),
            lo = quantile(dif, 0.025),
            hi = quantile(dif, 0.975),
            pval = (sum(abs(dif) < abs(dif - mean(dif))) + 1)/(n() + 1),
            .groups = "drop") |>
  # convert back to wt/mutant ratio for more interpretable values
  mutate(mid = exp(mid),
         lo = exp(lo),
         hi = exp(hi)) |>
  mutate(padj = p.adjust(pval, method = "bonferroni"))

# write to temporary file to compile into supplementary spreadsheet
write_csv(s12b_contrast_summary, "~/temp/figS12b.csv")


# Fig S12e ----

s12e <- read_csv("data/raw/triple_mutant_phenotype/Re _Your_manuscript,_NPLANTS-210711374C/FigS12e_phloem_mutants_Rfile.csv")

s12e <- s12e |> rename(condition = sucrose)

# contrast WT vs Mutant within each treatment (bootstrap)
s12e_contrast <- s12e |>
  # mutate(strata = interaction(genotype, condition)) |>
  # bootstraps(times = 500, strata = strata) |>
  bootstraps(times = 500) |>
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) |>
      group_by(genotype, condition, replicate) |>
      summarise(mean_length = mean(root_length), .groups = "drop") |>
      group_by(genotype, condition) |>
      summarise(mean_length = mean(mean_length), .groups = "drop") |>
      pivot_wider(names_from = genotype, values_from = mean_length) |>
      pivot_longer(cols = c(-wt, -condition),
                   names_to = "genotype",
                   values_to = "root_length") |>
      mutate(dif = wt - root_length)
  }))


# plot bootstrapped differences and calculate bootstrap p-value
s12e_contrast |>
  unnest(boot) |>
  ggplot(aes(dif)) +
  geom_density(aes(colour = genotype)) +
  geom_vline(xintercept = 0) +
  facet_grid(~ condition)

s12e_contrast_summary <- s12e_contrast |>
  unnest(boot) |>
  group_by(genotype, condition) |>
  summarise(mid = median(dif),
            lo = quantile(dif, 0.025),
            hi = quantile(dif, 0.975),
            pval = (sum(abs(dif) < abs(dif - mean(dif))) + 1)/(n() + 1),
            .groups = "drop") |>
  # convert back to wt/mutant ratio for more interpretable values
  mutate(mid = exp(mid),
         lo = exp(lo),
         hi = exp(hi)) |>
  mutate(padj = p.adjust(pval, method = "bonferroni")) |>
  arrange(condition, genotype)

# write to temporary file to compile into supplementary spreadsheet
write_csv(s12e_contrast_summary, "~/temp/figS12e.csv")


# Fig 4G ----

fig4g <- read_csv("data/raw/triple_mutant_phenotype/compilation_no sucrose_exp_id_SO.csv")

# clean data
fig4g <- fig4g %>%
  mutate(days = as.numeric(str_remove(days, "dps")),
         comp = case_when(str_detect(genotype, "_comp_") ~ "Complemented",
                          genotype == "wt" ~ "WT",
                          genotype == "ccc" ~ "ccc",
                          TRUE ~ NA_character_),
         genotype_tidy = ifelse(genotype == "ccc", "3papl", genotype)) %>%
  mutate(genotype_tidy = genotype_tidy %>% str_remove("ccc_comp_") %>% str_replace("_", "-"))

# keep only samples measured at 6 days post sowing
# other timepoints have too much missing data (most genotypes not assayed)
fig4g <- fig4g %>% filter(days == 6)

# contrast WT vs Mutant within each treatment (bootstrap)
fig4g_contrast <- fig4g |>
  select(genotype = genotype_tidy, root_length, stock, experiment_id) |>
  # mutate(strata = interaction(genotype)) |>
  # bootstraps(times = 500, strata = strata) |>
  bootstraps(times = 500) |>
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) |>
      # average per stock and experimental batch
      group_by(genotype, stock, experiment_id) |>
      summarise(mean_length = mean(root_length), .groups = "drop") |>
      # average-of-averages (per stock)
      group_by(stock, genotype) |>
      summarise(mean_length = mean(mean_length), .groups = "drop") |>
      # average-of-averages (per genotype)
      group_by(genotype) |>
      summarise(mean_length = mean(mean_length), .groups = "drop") |>
      pivot_wider(names_from = genotype, values_from = mean_length) |>
      pivot_longer(cols = c(-wt),
                   names_to = "genotype",
                   values_to = "root_length") |>
      mutate(dif = wt - root_length)
  }))

# plot bootstrapped differences and calculate bootstrap p-value
fig4g_contrast |>
  unnest(boot) |>
  ggplot(aes(exp(dif), genotype)) +
  ggridges::geom_density_ridges(alpha = 0.5) +
  geom_vline(xintercept = 1) +
  labs(x = "Relative Difference (WT/others)", y = "Genotype Compared to WT")

fig4g_contrast_summary <- fig4g_contrast |>
  unnest(boot) |>
  group_by(genotype) |>
  summarise(mid = median(dif),
            lo = quantile(dif, 0.025),
            hi = quantile(dif, 0.975),
            pval = (sum(abs(dif) < abs(dif - mean(dif))) + 1)/(n() + 1),
            .groups = "drop") |>
  # convert back to wt/mutant ratio for more interpretable values
  mutate(mid = exp(mid),
         lo = exp(lo),
         hi = exp(hi)) |>
  mutate(padj = p.adjust(pval, method = "bonferroni")) |>
  arrange(genotype)

# write to temporary file to compile into supplementary spreadsheet
write_csv(fig4g_contrast_summary, "~/temp/fig4g.csv")


# Figure 4h, S12c & S12g ----

# first 3papl experiment batch
fig4h <- read_excel("data/raw/triple_mutant_phenotype/Re _Your_manuscript,_NPLANTS-210711374C/ccc_fig4h_wt_figS12g.xlsx", na = "NA") |>
  mutate(batch = "experiment1")

# second experiment batch
figS12c <- read_excel("data/raw/triple_mutant_phenotype/Re _Your_manuscript,_NPLANTS-210711374C/FigS12c.xlsx", na = "NA") |>
  mutate(batch = "experiment2")

fig4hS12 <- bind_rows(fig4h, figS12c) |>
  mutate(day_transfer = day_transfer |> str_remove("^day") |> str_remove("^d")) |>
  mutate(day_transfer = ifelse(is.na(day_transfer), "control", day_transfer)) |>
  mutate(genotype = str_remove(genotype, "_.*")) |>
  filter(!is.na(sucrose))


# contrast treated vs control within each day and for each genotype (bootstrap)
fig4hS12_contrast <- fig4hS12 |>
  # mutate(strata = interaction(genotype, sucrose, day_transfer)) |>
  # bootstraps(times = 500, strata = strata) |>
  bootstraps(times = 500) |>
  mutate(boot = map(splits, function(x){
    dat <- as_tibble(x)
    dat %>%
      mutate(root_length = log(root_length)) |>
      # average per treatment and experimental batch
      group_by(genotype, sucrose, day_transfer, batch) |>
      summarise(mean_length = mean(root_length), .groups = "drop") |>
      # average-of-averages (per treatment combination)
      group_by(genotype, sucrose, day_transfer) |>
      summarise(mean_length = mean(mean_length), .groups = "drop") |>
      # contrast
      pivot_wider(names_from = sucrose, values_from = mean_length) |>
      mutate(dif = yes - no)
  }))

# figure similar to the one in the paper
fig4hS12 |>
  ggplot(aes(sucrose, root_length)) +
  ggbeeswarm::geom_quasirandom() +
  facet_grid(genotype + batch ~ day_transfer)

# plot bootstrapped differences and calculate bootstrap p-value
fig4hS12_contrast |>
  unnest(boot) |>
  ggplot(aes(exp(dif), day_transfer)) +
  ggridges::geom_density_ridges(alpha = 0.5) +
  geom_vline(xintercept = 1) +
  facet_grid(genotype ~ .) +
  labs(x = "Relative Difference (sucrose / no sucrose)", y = "Transfer day") +
  scale_x_continuous(trans = "log2")

fig4hS12_contrast_summary <- fig4hS12_contrast |>
  unnest(boot) |>
  group_by(genotype, day_transfer) |>
  summarise(mid = median(dif),
            lo = quantile(dif, 0.025),
            hi = quantile(dif, 0.975),
            pval = (sum(abs(dif) < abs(dif - mean(dif))) + 1)/(n() + 1),
            .groups = "drop") |>
  # convert back to wt/mutant ratio for more interpretable values
  mutate(mid = exp(mid),
         lo = exp(lo),
         hi = exp(hi)) |>
  mutate(padj = p.adjust(pval, method = "bonferroni")) |>
  arrange(genotype, day_transfer)

fig4hS12_contrast_summary <- fig4hS12_contrast |>
  unnest(boot) |>
  group_by(genotype, day_transfer) |>
  summarise(mid = median(dif),
            lo = quantile(dif, 0.025),
            hi = quantile(dif, 0.975),
            pval = (sum(abs(dif) < abs(dif - mean(dif))) + 1)/(n() + 1),
            .groups = "drop") |>
  # convert back to wt/mutant ratio for more interpretable values
  mutate(mid = exp(mid),
         lo = exp(lo),
         hi = exp(hi)) |>
  mutate(padj = p.adjust(pval, method = "bonferroni")) |>
  arrange(genotype, day_transfer)


# write to temporary file to compile into supplementary spreadsheet
write_csv(fig4hS12_contrast_summary, "~/temp/fig4hS12.csv")
