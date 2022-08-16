
# Load tidyverse
library(tidyverse); library(ggforce); library(magick); library(png)

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

# Load in some pictures that are provided first, in case you want to pretty it up! Stretch task!
  chromas <- readPNG("./pics/chromis.png")
whitedams <- readPNG("./pics/whitedams.png")
    lemon <- readPNG("./pics/lemon.png")
   humbug <- readPNG("./pics/humbug.png")
  acantho <- readPNG("./pics/acantho.png")
    ambon <- readPNG("./pics/ambon.png")

# Now lets make a plot, boy this took me a while, and I'm still not happy with it, but it's looking cool!
ggplot(data_clean, aes(x =species, y = activity)) +
  geom_violin(aes(fill = treatment, col = treatment), alpha = 0.8) +
  geom_sina(data = . %>% filter(activity >= 25), aes(fill = treatment, col = treatment), colour = "orange", alpha=0.2) + # Here's a COOL feature. You can wrangle data directly from data frame to highly something neat, like activities all above and below 25. The '.' is meant to represent the tibble data frame and just like we would manipulate there we can do the same here.
  geom_sina(data = . %>% filter(activity <= 25), aes(fill = treatment, col = treatment), colour = "blue", alpha=0.2) +
  labs(x = "Species", y = "Activity (s)", fill = "Treatment", col = "Treatment") +
  ylim(0, 100) +  # Make axes large as we may want to add in pictures of the species
  geom_errorbar(data = summary_data, aes(x = species, y = x, ymin = x+2*se, ymax = x-2*se), position=position_dodge2(width = 0.01, padding = 0.5), size = 0.5) +   # We've already created the means and se data so we can feed this new tibble in
  geom_point(data = summary_data, aes(x = species, y = x, fill = treatment, col = treatment), position = position_dodge(width = 0.9), colour = "black") +
  theme_classic() +
  annotate("text", x = seq(from = 0.75, to = 6.5, by= 0.5), y = 63, label = paste0("n = ", summary_data$n)) +
  scale_fill_manual(breaks = c("CO2", "control"),
                    values=c("#fc9197", "#dde6d5")) +
  scale_color_manual(breaks = c("CO2", "control"),
                  values=c("#fc9197", "#dde6d5")) +
  coord_flip() +
  annotation_raster(chromas, xmin = 2.5, xmax = 3.5, ymin = 75, ymax = 100) +
  annotation_raster(whitedams, xmin =5.5, xmax = 6.5, ymin = 75, ymax = 100) +
  annotation_raster(lemon, xmin =4.5, xmax = 5.5, ymin = 75, ymax = 100) +
  annotation_raster(humbug, xmin =3.5, xmax = 4.5, ymin = 75, ymax = 100) +
  annotation_raster(acantho, xmin =0.5, xmax = 1.5, ymin = 75, ymax = 100) +
  annotation_raster(ambon, xmin =1.6, xmax = 2.5, ymin = 75, ymax = 100)



ggsave(filename = "species_trt", dpi = 600, device = "png", path = "./output/figs")

  geom_pointrange(data = summary_data, aes(x = species, y = x, ymin = x+se, ymax = x-se), position = position_jitterdodge(dodge.width = 0.85))

  annotate("text", x = species, y = 60, label = summary_data$n)




