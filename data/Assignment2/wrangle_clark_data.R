
    dat_OA <- read_csv("./data/Assignment2/OA_activitydat_20190302_BIOL3207.csv")
clark_meta <- read_csv("./data/Assignment2/clark_paper_data.csv")

summary_data <-dat_OA %>% group_by(species, treatment) %>%
              summarise(mean = mean(activity, na.rm = TRUE),
                        sd = sd(activity, na.rm = TRUE),
                        n = length(unique(animal_id))) %>%
              rename(Species = "species")

total <- cbind(clark_meta, summary_data)

final <- pivot_wider(total, names_from = treatment,
                     names_glue = "{treatment}_{.value}",
                     values_from = c("mean", "sd", "n"))

meta_data_full <- read_csv("./data/Assignment2/ocean_meta_data.csv")

dim(meta_data_full)
dim(final)

## Do some renaming of colnames so they match meta-Data_full
final2 <- final %>% rename("oa.mean" = CO2_mean,
                            "oa.sd" = CO2_sd,
                            "oa.n" = CO2_n,
                            "ctrl.mean" = control_mean,
                            "ctrl.sd" = control_sd,
                            "ctrl.n" = control_n)

# Reorder col names based on names in meta_data_full
final2 <- final2[names(meta_data_full)]

# Check columns are in same order
colnames(meta_data_full) == colnames(final2)

# Bind teh two dataframes
full_final <- rbind(meta_data_full, final2)

