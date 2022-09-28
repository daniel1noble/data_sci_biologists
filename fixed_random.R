set.seed(86) # Set see so that we all get the same simulated results


# We will have 5 studies
stdy  <- 1:5


# We know the variance for each effect
Ves   <- c(0.05, 0.10, 0.02, 0.10, 0.09)


# We'll need this later but these are weights
W     <- 1 / Ves


# We assume they are sampled from a normal distribution  with a mean effect size of 2
es    <- rnorm(length(Ves), 2, sqrt(Ves))


# Data for our fixed effect meta-analysis
dataFE <- data.frame(stdy = stdy-0.2, es, Ves, type = "FE")

esRE <- rnorm(length(Ves), 2, sqrt(Ves + 0.8))


# Data for our random effect meta-analysis
dataRE <-  data.frame(stdy = stdy + 0.2, es = esRE, Ves, type = "RE")

full <- rbind(dataFE, dataRE)



ggplot(full %>% filter(type == "FE"), aes(y = stdy, x = es, colour = type)) +
  geom_point() + geom_errorbar(aes(xmin = es - 2*sqrt(Ves), xmax = es + 2*sqrt(Ves)), width = 0.125) +
  ylim(0, 6) + xlim(-0.1, 4) + geom_vline(xintercept = 2, linetype =2) +
  labs(y = "Study", x = "Effect Size (+95% CI)", colour = "Model Type") + theme_classic()
