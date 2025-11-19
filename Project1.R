data <- read.csv('~/AdvancedEcologyWorkshop/group2.csv')
library(tidyverse)
library('dplyr')
library('DescTools')
library('vegan')
library(ggplot2)
#Section 1
##1 – Use a discrete density-dependent and -independent population model to assess the population growth 
##dynamics of any three species within the 10 most abundant species. What can you infer from both models 
##(Modeli | Speciesi)? For each species, what is the most proper model to assess their growth? Why?
species_abundance <- data |>
  group_by(species) |>
  summarise(total_abundance = sum(n, na.rm = T)) |>
  arrange(desc(total_abundance)) |>
  slice(1:10)
top3species <- species_abundance$species[1:3]
top3species

data_over_time <- data |>
  filter(species %in% top3species) |>
  group_by(species, date) |>
  summarise(N = sum(n, na.rm = T)) |>
  arrange(species, date) |>
  ungroup()
data_over_time$date <- as.Date(data_over_time$date, format = "%m/%d/%Y")
ggplot(data_over_time, aes(x = date, y = N, color = species)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(title = "Abundance of Top 3 Species Over Time",
       x = "Date",
       y = "Abundance",
       color = "Species")

###Density-independent
data_w_lambda <- data_over_time |>
  arrange(species, date) |>
  mutate(lambda = lead(N)/N)
avg_lambda <- data_w_lambda |>
  group_by(species) |>
  summarise(
    mean_lambda = Gmean(lambda, na.rm = T),
    mean_sd = Gsd(lambda, na.rm = T)
  )
ggplot(data_w_lambda |> filter(!is.na(lambda)),
       aes(x = date, y = lambda, color = species)) +
  geom_point() +
  geom_line() +
  labs(title = "Lambda Over Time", x = "Date", y = "Lambda", color = "Species") +
  theme_bw()
trends <- avg_lambda |>
  mutate(pop_trends = case_when(
    mean_lambda > 1.0 ~ "Growing",
    mean_lambda == 1.0 ~ "Satble",
    mean_lambda < 1.0 ~ "Declining"))

years <- 30
pop <- data_w_lambda |>
  group_by(species) |>
  slice_max(date, n = 1) |>
  select(species, N_0 = N, t_0 = date) |>
  left_join(avg_lambda, by = "species")
sim <- tibble()
for(i in 1:nrow(pop)){
  Species <- pop$species[i]
  N0 <- pop$N_0[i]
  lambda <- pop$mean_lambda[i]
  t0 <- pop$t_0[i]
  N <- numeric(years + 1)
  N[1] <- N0
  for(t in 1:years) {
    N[t + 1] <- lambda * N[t]
  }
  projected_pop <- tibble(
    species = Species,
    year = 0:years,
    projection = N,
    t_i = t0 + years(0:years)
  )
  sim <- bind_rows(sim, projected_pop)
}
ggplot(sim, aes(x = year, y = projection, color = species)) +
  geom_point() +
  geom_line() +
  labs(x = "Years", y = "Projected Abundance", color = "Species") +
  theme_bw()

###Density-dependent
years <- 30
pop2 <- data_w_lambda |>
  group_by(species) |>
  slice_max(date, n = 1) |>
  select(species, N_0 = N, t_0 = date) |>
  left_join(avg_lambda, by = "species") |>
  mutate(r = log(mean_lambda))
K_estimate <- data_over_time |>
  group_by(species) |>
  summarise(K = mean(N, na.rm = T))
pop2 <- pop2 |>
  left_join(K_estimate, by = "species")
sim2 <- tibble()
for(i in 1:nrow(pop2)) {
  Species <- pop2$species[i]
  N0 <- pop2$N_0[i]
  r <- pop2$r[i]
  K <- pop2$K[i]
  t0 <- pop2$t_0[i]
  N <- numeric(years + 1)
  N[1] <- N0
  for(t in 1:years) {
    N[t + 1] <- N[t] + r * N[t] * (1 - N[t] / K)
    if (N[t + 1] < 0) {
      N[t + 1] <- 0
    }
  }
  projected_pop2 <- tibble(
    species = Species,
    year = 0:years,
    projection = N,
    t_i = t0 + years(0:years),
    K = K,
    r = r
  )
  sim2 <- bind_rows(sim2, projected_pop2)
}
ggplot(sim2, aes(x = year, y = projection, color = species)) +
  geom_line() +
  geom_hline(data = sim2 |> distinct(species, K),
             aes(yintercept = K, color = species),
             linetype = "dashed", alpha = 0.5) +
  labs(
    x = "Years",
    y = "Projected Abundance",
    color = "Species"
  ) +
  theme_bw()

##2 – Based on the selected models, how does lambda change across time?

##3 – Starting with the average abundance in 2000, perform a simulation (N = 100 replicates) of population 
##growth using the growth variability based on the abundance values from start to 1999 (i.e., use 100 samples 
##of lambda from 1980 to 1999’s frequency distribution to project population beyond 2000). How many times 
##are the real abundance values intercepted by the uncertainty bounds of the simulation (i.e., values within 
##the 95% confidence intervals)? Briefly discussed the performance of the models and potential reasons for 
##the results.

#Section 2
##1 – Does the sampling assigned to your group reach the asymptote of the species accumulation curve? Discuss results.
Rare <- data |> 
  group_by(date) |> 
  ungroup()
sp_rare  <- unique(Rare$species)
sp_rare_tib = tibble(n_samp = 1:length(sp_rare), n_spp = NA)
for (i in 1:length(sp_rare)){
  # sample ID to include
  samp = sp_rare[1:i]
  # include only sample numbers 
  d = Rare |> 
    filter(species %in% samp,
           n > 0)
  sp_rare_tib$n_spp[i] = length(unique(d$species))
}
ggplot(sp_rare_tib, aes(n_samp, n_spp))+
  geom_line(linewidth = 1)+
  labs(x = 'Number of Samples',
       y = 'Number of Species')+
  theme_bw()
###Need help from Santos - why is this not working?

##2 - Using the species that composed 75% of the abundance of species within the nektonic community, assess 
##the annual trends of alpha and beta diversity, and the inverse Simpson diversity index.
data2 <- data |>
  group_by(species) |>
  summarise(sum(n))
sum(data2$`sum(n)`)
16745*0.75
###Black drum + Red drum + Atlantic Croaker + Blue Crab + Striped Mullet

sp75 <- c('Black Drum', 'Red Drum', 'Atlantic Croaker', 'Blue Crab', 'Striped Mullet')
###Alpha: Shannon, Simpson. Beta: Bray-Curtis, Jaccard
###Per year
library(dplyr)
library(lubridate)

data <- data %>%
  mutate(date = mdy(date))
data <- data %>%
  mutate(year = year(date))
annual_abundance <- data %>%
  group_by(species, year) %>%
  summarise(total_abundance = sum(n, na.rm = TRUE), .groups = "drop")
annual_abundance <- annual_abundance |>
  filter(species %in% sp75)

df = unique(annual_abundance[c("year")])
df$H = NA
annual_abundance <- annual_abundance %>%
  rename(abundance = total_abundance)
for (i in 1:nrow(df)){
  d = annual_abundance |> filter(year == df$year[i])
  d = d |> 
    count(species,wt = abundance) |> 
    mutate(pi = n/sum(n),
           ln_pi = log(pi),
           p_ln_pi = pi*ln_pi)
  
  df$H[i] = -sum(d$p_ln_pi)
}
df |>
  summarise(mean_H = mean(H, na.rm = TRUE),
            sd_H = sd(H, na.rm = TRUE))
ggplot(df, aes(year, H)) +
  geom_line(linewidth = 1)+ 
  geom_point() +
  labs(x = 'Date',
       y = 'Shannon Diversity')+
  theme_bw()

df$d <- NA
for (i in 1:nrow(df)){
  d = annual_abundance |> filter(year == df$year[i])
  d = d |> count(species,wt = abundance) |> 
    mutate(pi = n/sum(n))
  
  df$d[i] = sum(d$pi^2)
}

df$even <- NA
df$inv <- NA
for(i in 1:nrow(df)){
df$even[i] = 1 - df$d[i]
df$inv[i] = 1/df$d[i]}

ggplot(df, aes(year, d)) +
  geom_line(linewidth = 1)+ 
  geom_point() +
  labs(x = 'Date',
       y = 'Simpson Diversity')+
  theme_bw()

ggplot(df, aes(year, inv)) +
  geom_line(linewidth = 1)+ 
  geom_point() +
  labs(x = 'Date',
       y = 'Inverse Simpson Diversity')+
  theme_bw()

mean(df$inv)

##3 – Using a community trajectory analysis, assess the annual stability of the fish community across the 
##sites from the assigned basin and gear type. Describe and discuss the main results.