

clark_meta <- read_csv("./data/Assignment2/clark_paper_data.csv")

total <- cbind(clark_meta, summary_data)

final <- pivot_wider(total, names_from = Treatment,
                     names_glue = "{Treatment}_{.value}",
                     values_from = c("Mean Activity", "Standard Error", "N_inds", "N"))
