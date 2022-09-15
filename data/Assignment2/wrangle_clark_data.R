
    dat_OA <- read_csv("./data/Assignment2/OA_activitydat_20190302_BIOL3207.csv")
clark_meta <- read_csv("./data/Assignment2/clark_paper_data.csv")

summary_data <-dat_OA %>% group_by(species, treatment) %>% summarise(mean = mean(activity, na.rm = TRUE), sd = sd(activity, na.rm = TRUE), n = length(unique(animal_id)))

total <- cbind(clark_meta, summary_data)

final <- pivot_wider(total, names_from = treatment,
                     names_glue = "{treatment}_{.value}",
                     values_from = c("mean", "sd", "n"))

meta_data_full <- read_csv("./data/Assignment2/ocean_meta_data.csv")

dim(meta_data_full)
dim(final)


