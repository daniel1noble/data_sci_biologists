
# Load tidyverse
library(tidyverse); library(ggforce)

# Load data
data <- read_csv("./data/clark_et-al_2020/OA_activitydat_20190302_BIOL3207.csv")

# Look at the data
head(data)
str(data)

## Data wrangling; get rid of NA's and remove unecssary columns like rownames from file, comments and loation
data_clean <- data %>% filter(!is.na(activity)) %>% select(-c("...1", "comment", "loc"))
str(data_clean)

# Summarise the mean and SE, n for each speceis by treatment combination which we can add to the plot!
summary_data <- data_clean %>% group_by(species, treatment) %>% summarise(x = mean(activity), se=sd(activity) / n(), n = n())

## Make a cool plot of the data based on species and treatment
ggplot(data_clean, aes(x =species, y = activity)) +
  geom_violin(aes(fill = treatment, col = treatment), alpha = 0.8) +
  geom_sina(aes(fill = treatment, col = treatment), colour = "black", alpha=0.2) +
  labs(x = "Species", y = "Activity (s)", fill = "Treatment", col = "Treatment") +
  ylim(0, 65) +
  geom_errorbar(data = summary_data, aes(x = species, y = x, ymin = x+2*se, ymax = x-2*se), position=position_dodge2(width = 0.01, padding = 0.5), size = 0.5) +
  geom_point(data = summary_data, aes(x = species, y = x, fill = treatment, col = treatment), position = position_dodge(width = 0.9), colour = "black") +
  theme_bw() +
  annotate("text", x = seq(from = 0.75, to = 6.5, by= 0.5), y = 63, label = paste0("n = ", summary_data$n)) +
  scale_fill_manual(breaks = c("CO2", "control"),
                    values=c("#fc9197", "#dde6d5")) +
  scale_color_manual(breaks = c("CO2", "control"),
                  values=c("#fc9197", "#dde6d5")) +
  coord_flip()

ggsave(filename = "species_trt", dpi = 600, device = "png", path = "./output/figs")

  geom_pointrange(data = summary_data, aes(x = species, y = x, ymin = x+se, ymax = x-se), position = position_jitterdodge(dodge.width = 0.85))

  annotate("text", x = species, y = 60, label = summary_data$n)




