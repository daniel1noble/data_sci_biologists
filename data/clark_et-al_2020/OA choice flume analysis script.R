# R script for analyses of predator odour avoidance in choice flumes #

#setwd('')

library(car)
library(MASS)


# import data
flume <- read.csv('OA_flumedat_20190302.csv', header=TRUE)

{
  flume$comment <- as.character(flume$comment)
  flume$animal_id <- as.character(flume$animal_id)
  flume$animal_id <- as.integer(flume$animal_id)
  flume$loc <- as.character(flume$loc)
  flume$species <- as.character(flume$species)
  flume$species <- tolower(flume$species)
  flume$treatment <- relevel(flume$treatment, 'control')
  flume$size <- as.character(flume$size)
  flume$size <- ifelse(flume$size=="", NA, flume$size)
}

str(flume)


with(flume, table(species, loc))









# linear models for humbugs first 

# 2014
humdat <- flume[flume$species=='humbug',]

foo <- humdat[humdat$loc=='LIRS 2014',]
hum21 <- lm(attraction ~ treatment*SL, data=foo)
drop1(hum21, test='F')
hum22 <- lm(attraction ~ treatment + SL, data=foo)
summary(hum22)
drop1(hum22, test='F')
hum23 <- lm(attraction ~ treatment, data=foo)
summary(hum23) #final model



#2016

foo <- humdat[humdat$loc=='LIRS 2016',]
hum25 <- lm(attraction ~ treatment*SL, data=foo)
drop1(hum25, test='F')
hum26 <- lm(attraction ~ treatment + SL, data=foo)
summary(hum26) # final model


# visual checks of model assumptions
# replace model name as needed
# ensure dataframe (e.g., w/ temporary name 'foo') is same one used to fit model
E1 <- resid(hum26)
F1 <- fitted(hum26)

qqp(hum26$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ foo$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = foo$SL, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)











# acanthochromis

acdat <- flume[flume$species=='acantho',]


# 2015 AIMS
foo <- acdat[acdat$loc=='AIMS 2015',]

ac12 <- lm(attraction ~ treatment*SL, data=foo)
summary(ac12) 
drop1(ac12, test='F')
ac12.1 <- lm(attraction ~ treatment + SL, data=foo)
summary(ac12.1)
drop1(ac12.1, test='F')
ac12.2 <- lm(attraction ~ treatment, data=foo)
summary(ac12.2) # FINAL MODEL FOR AIMS 2015 


# 2016 lizard island
foo <- acdat[acdat$loc=='LIRS 2016',]
ac21 <- lm(attraction ~ treatment*SL, data=foo) 
summary(ac21)
drop1(ac21, test='F')
ac21.2 <- lm(attraction ~ treatment + SL, data=foo)
summary(ac21.2) #p = 0.9
drop1(ac21.2, test='F')
ac21.3 <- lm(attraction ~ treatment, data=foo)
summary(ac21.3) # FINAL MODEL FOR ACANTHO 2016



# visual checks of model assumptions
# replace model name as needed
# ensure dataframe (e.g., w/ temporary name 'foo') is same one used to fit model
E1 <- resid(ac21.3)
F1 <- fitted(ac21.3)

qqp(ac21.3$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ foo$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = foo$SL, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)











# white damsels

wdat <- flume[flume$species=='whitedams',]

wd11 <- lm(attraction ~ treatment*SL, data=wdat)
drop1(wd11, test='F')

wd12 <- lm(attraction ~ treatment + SL, data=wdat)
summary(wd12)

hist(wdat$SL) # white damsels actually parse out into two distinct size groups, so we'll try with size class too

wdat$size <- as.factor(wdat$size)
wdat$size <- relevel(wdat$size, ref='small')

wd13 <- lm(attraction ~ treatment + size, data=wdat)
drop1(wd13, test='F') 
# indeed, size class is a stronger effect than standard length here, so we'll present that as our model
# no effect of treatment regardless of how size is treated
summary(wd13)


# visual checks of model assumptions
# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(wd13)
F1 <- fitted(wd13)

qqp(wd13$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ wdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = wdat$size, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)














# ambons

ambdat <- flume[flume$species=='ambon',]

amb11 <- lm(attraction ~ SL*treatment, data=ambdat)
drop1(amb11, test='F')
amb12 <- lm(attraction ~ SL + treatment, data=ambdat)
drop1(amb12, test='F')
amb13 <- lm(attraction ~ treatment, data=ambdat)
summary(amb13) # FINAL MODEL FOR AMBOS (2014)

# visual checks of model assumptions
# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(amb13)
F1 <- fitted(amb13)

qqp(amb13$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ ambdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = ambdat$SL, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)













# lemons

lemdat <- flume[flume$species=='lemon',]

lem11 <- lm(attraction ~ SL*treatment, data=lemdat)
drop1(lem11, test='F')
lem12 <- lm(attraction ~ SL + treatment, data=lemdat)
drop1(lem12, test='F')
lem13 <- lm(attraction ~ treatment, data=lemdat)
drop1(lem13, test='F')
summary(lem13)


# visual checks of model assumptions
# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(lem13)
F1 <- fitted(lem13)

qqp(lem13$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ lemdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = lemdat$SL, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)














# chromis
chdat <- flume[flume$species=='chromis',]

ch11 <- lm(attraction ~ SL*treatment, data=chdat)
drop1(ch11, test='F')
ch12 <- lm(attraction ~ SL + treatment, data=chdat)
drop1(ch12, test='F') # p = 0.056 for effect of SL
ch13 <- lm(attraction ~ treatment, data=chdat)
summary(ch13)


# visual checks of model assumptions
# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(ch13)
F1 <- fitted(ch13)

qqp(ch13$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ chdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = chdat$SL, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)








# end modeling of choice flume data







