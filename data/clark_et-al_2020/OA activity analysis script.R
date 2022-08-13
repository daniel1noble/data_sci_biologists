# R script for analyses of activity levels #

#setwd('')


# import data
oaa <- read.csv('OA_activitydat_20190302.csv', header=TRUE)

library(car)
library(MASS)


# format dataframe
{
oaa$comment <- as.character(oaa$comment)
oaa$animal_id <- as.character(oaa$animal_id)
oaa$animal_id <- as.integer(oaa$animal_id)
oaa$loc <- as.character(oaa$loc)
oaa$species <- as.character(oaa$species)
oaa$species <- tolower(oaa$species)
oaa$treatment <- relevel(oaa$treatment, 'control')
}

#inspect data
str(oaa)
with(oaa, table(species, loc, size))












# linear models looking at fish length (sl = standard length) and treatment
# separate models by year (= loccation) and species


# acanthochromis

acdat <- oaa[oaa$species=='acantho',]
with(acdat, table(loc, treatment))

#2014 data
foo <- acdat[acdat$loc=='LIRS 2014',]
ac100 <- lm(activity ~ treatment*sl, data=foo)
drop1(ac100, test='F')
ac101 <- lm(activity ~ treatment + sl, data=foo)
summary(ac101) # final model for acanthos 2014


#2015 data
foo <- acdat[acdat$loc=='AIMS 2015',]
ac110 <- lm(activity ~ treatment*sl, data=foo)
drop1(ac110, test='F')
summary(ac110) # final model


#2016 data
foo <- acdat[acdat$loc=='LIRS 2016',]
hist(foo$sl)
ac120 <- lm(activity ~ treatment*sl, data=foo)
drop1(ac120, test='F')
ac121 <- lm(activity ~ treatment + sl, data=foo)
summary(ac121)


# visual checks of model assumptions
# replace model name as needed
# ensure dataframe (w/ temporary name 'foo') is same one used to fit model
E1 <- resid(ac121)
F1 <- fitted(ac121)

qqp(ac121$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ foo$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = foo$sl, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)











# humbugs

humdat <- oaa[oaa$species=='humbug',]

# 2014
foo <- humdat[humdat$loc=='LIRS 2014',]
hum101 <- lm(activity ~ treatment*sl, data=foo)
summary(hum101) # final model
drop1(hum101, test='F') 


# 2016
foo <- humdat[humdat$loc=='LIRS 2016',]
hum121 <- lm(activity ~ treatment*sl, data=foo)
summary(hum121)
hum122 <- lm(activity ~ treatment + sl, data=foo)
summary(hum122)
drop1(hum122, test='F')
hum123 <- lm(activity ~ treatment, data=foo)
summary(hum123) # final model


# visual checks of model assumptions
# replace model name as needed
# ensure dataframe (w/ temporary name 'foo') is same one used to fit model
E1 <- resid(hum123)
F1 <- fitted(hum123)

qqp(hum123$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ foo$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = foo$sl, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)














# lemons
lemdat <- oaa[oaa$species=='lemon',]

lem1 <- lm(activity ~ sl*treatment, data=lemdat)
drop1(lem1, test='F')

lem2 <- lm(activity ~ sl + treatment, data=lemdat)
drop1(lem2, test='F')

lem3 <- lm(activity ~ treatment, data=lemdat)
drop1(lem3, test='F')
summary(lem3) # final model

# visual checks of model assumptions
# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(lem3)
F1 <- fitted(lem3)

qqp(lem3$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ lemdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = lemdat$sl, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)



















# white damsels 

wdat <- oaa[oaa$species=='whitedams',]

dams1 <- lm(activity ~ sl*treatment, data=wdat)
drop1(dams1, test='F')

dams2 <- lm(activity ~ treatment + sl, data=wdat)
drop1(dams2, test='F')

summary(dams2) # final model

# visual checks of model assumptions
# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(dams2)
F1 <- fitted(dams2)

qqp(dams2$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ wdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = wdat$sl, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)













# chromis

chrdat <- oaa[oaa$species=='chromis',]

chr1 <- lm(activity ~ sl*treatment, data=chrdat)
drop1(chr1, test='F')

chr2 <- lm(activity ~ sl + treatment, data=chrdat)
drop1(chr2, test='F')

chr3 <- lm(activity ~ treatment, data=chrdat)
drop1(chr3, test='F')
summary(chr3) # final model


# visual checks of model assumptions
# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(chr3)
F1 <- fitted(chr3)

qqp(chr3$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ chrdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = chrdat$sl, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)














# ambons

ambdat <- oaa[oaa$species=='ambon',]
with(ambdat, table(size, loc))

# linear models
amb1 <- lm(activity ~ sl*treatment, data=ambdat)
drop1(amb1, test='F')

amb2 <- lm(activity ~ sl + treatment, data=ambdat)
drop1(amb2, test='F')

amb3 <- lm(activity ~ treatment, data=ambdat)
drop1(amb3, test='F')
summary(amb3)

# nothing significant.

# replace model name as needed
# ensure dataframe is same one used to fit model
E1 <- resid(amb3)
F1 <- fitted(amb3)

qqp(amb3$residuals)

plot(x = F1, y = E1, xlab = 'fitted values', ylab = 'residuals') 
abline(h = 0)

boxplot(E1 ~ ambdat$treatment, main = "treatment", ylab = 'residuals')
abline(h = 0)

plot(x = ambdat$sl, y = E1, xlab = "standard length", ylab = 'residuals') 
abline(h = 0)











# statistics end





















# figure for paper

library(scales)
library(plotrix)


file_nm <- 'OA activity scatter_confint_20190304'

# run one of these, then run from the brackets below
jpeg(paste(file_nm, '.jpg', sep=""), width=14, height=6, units='in', res=600)
pdf(paste(file_nm, '.pdf', sep=""), width=14, height=6)


{
  ymin <- 0
  ymax <- 60
  ylab <- seq(from=ymin, to=ymax, by=10)
  
  layout(matrix(c(1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6), 1, 29))
  
  par(mar=c(1, 0, 0, 1), bty='n', oma=c(7, 7, 1, 0))
  
  # acanthos
  plot(NA, ylim=c(-3, ymax), xlim=c(0.4, 3.6), xaxt = 'n', yaxt = 'n')
  
  rect(-5000, -5000, 10000, 10000, col='grey92')
  abline(v = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2), lwd=3.6, col='white')
  abline(h = c(ylab), lwd=3.6, col='white')
  
  #2014
  foo <- oaa[oaa$species=='acantho' & oaa$loc=='LIRS 2014',]
  
  points(activity ~ jitter(as.numeric(treatment)-0.2, 11), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment)-0.8, 5.5), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  axis(2, at = ylab, labels = ylab, lty=0, tick=NA, cex.axis=2, las=2)
  mtext(expression(paste('Activity duration (', 's min'^'-1', ')')), side=2, line=3.25, cex=2)
  
  #mtext(expression(italic('A. polyacanthus')), cex=1.5, side=3, line=1)
  
  text(0.6, 58, labels = 'A', cex=3.4, font=2)
  
  points(x = 0.8, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 0.8, x1 = 0.8, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.65, x1 = 0.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.65, x1 = 0.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 1.2, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.05, x1 = 1.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.05, x1 = 1.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 1.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  axis(1, at = c(0.8, 1.2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  
  
  #2015
  foo <- oaa[oaa$species=='acantho' & oaa$loc=='AIMS 2015',]
  
  points(activity ~ jitter(as.numeric(treatment)+0.8, 3), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment)+0.2, 3), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  points(x = 1.8, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 1.8, x1 = 1.8, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 1.65, x1 = 1.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 1.65, x1 = 1.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  
  segments(x0 = 2.2, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 2.05, x1 = 2.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 2.05, x1 = 2.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 2.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  axis(1, at = c(1.8, 2.2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  
  
  # 2016 
  foo <- oaa[oaa$species=='acantho' & oaa$loc=='LIRS 2016',]
  
  points(activity ~ jitter(as.numeric(treatment)+1.8, 2), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment)+1.2, 2), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  points(x = 2.8, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 2.8, x1 = 2.8, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 2.65, x1 = 2.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 2.65, x1 = 2.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  #points(x = 3.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=16, col = 'dodgerblue2')
  segments(x0 = 3.2, x1 = 3.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 3.05, x1 = 3.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 3.05, x1 = 3.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 3.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  text(1,-3, labels='2014', cex=2)
  text(2,-3, labels='2015', cex=2)
  text(3,-3, labels='2016', cex=2)
  
  axis(1, at = c(2.8, 3.2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  
  
  
  
  # humbugs
  plot(NA, ylim=c(-3, ymax), xlim=c(0.4, 2.6), xaxt = 'n', yaxt = 'n')
  
  rect(-5000, -5000, 10000, 10000, col='grey92')
  abline(v = c(0.8, 1.2, 1.8, 2.2), lwd=3.6, col='white')
  abline(h = c(ylab), lwd=3.6, col='white')
  
  
  #2014
  
  foo <- oaa[oaa$species=='humbug' & oaa$loc=='LIRS 2014',]
  
  points(activity ~ jitter(as.numeric(treatment)-0.2, 11), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment)-0.8, 5.5), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  #mtext(expression(italic('D. aruanus')), cex=1.5, side=3, line=1)
  
  text(0.6, 58, labels = 'B', cex=3.4, font=2)
  
  points(x = 0.8, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 0.8, x1 = 0.8, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.65, x1 = 0.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.65, x1 = 0.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  #points(x = 1.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=16, col = 'dodgerblue2')
  segments(x0 = 1.2, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.05, x1 = 1.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.05, x1 = 1.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 1.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  
  
  # 2016
  
  foo <- oaa[oaa$species=='humbug' & oaa$loc=='LIRS 2016',]
  
  points(activity ~ jitter(as.numeric(treatment)+0.8, 3), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment)+0.2, 3), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  
  points(x = 1.8, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 1.8, x1 = 1.8, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 1.65, x1 = 1.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 1.65, x1 = 1.95, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  #points(x = 2.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=16, col = 'dodgerblue2')
  segments(x0 = 2.2, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 2.05, x1 = 2.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 2.05, x1 = 2.35, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 2.2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  
  text(1,-3, labels='2015', cex=2)
  text(2,-3, labels='2016', cex=2)
  
  axis(1, at = c(1.8, 2.2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  axis(1, at = c(0.8, 1.2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  
  
  
  # chromis
  plot(NA, ylim=c(-3, ymax), xlim=c(0.4, 2.6), xaxt = 'n', yaxt = 'n')
  
  rect(-5000, -5000, 10000, 10000, col='grey92')
  abline(v = c(1:2), lwd=3.6, col='white')
  abline(h = c(ylab), lwd=3.6, col='white')
  #points(activity ~ jitter(as.numeric(treatment), 1.2), data=oaa[oaa$species=='chromis',], lwd=1, pch=1)
  foo <- oaa[oaa$species=='chromis',]
  
  points(activity ~ jitter(as.numeric(treatment), 10), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment), 5), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  #mtext(expression(italic('C. atripectoralis')), cex=1.5, side=3, line=1)
  axis(1, at = c(1,2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  
  text(1.5,-3, labels='2014', cex=2)
  text(0.6, 58, labels = 'C', cex=3.4, font=2)
  
  points(x = 1, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 1, x1 = 1, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  #points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=16, col = 'dodgerblue2')
  segments(x0 = 2, x1 = 2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  
  # ambon
  plot(NA, ylim=c(-3, ymax), xlim=c(0.4, 2.6), xaxt = 'n', yaxt = 'n')
  
  rect(-5000, -5000, 10000, 10000, col='grey92')
  abline(v = c(1:2), lwd=3.6, col='white')
  abline(h = c(ylab), lwd=3.6, col='white')
  #points(activity ~ jitter(as.numeric(treatment), 1.2), data=oaa[oaa$species=='ambon',], lwd=1, pch=1)
  foo <- oaa[oaa$species=='ambon',]
  
  points(activity ~ jitter(as.numeric(treatment), 10), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment), 5), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  mtext(expression(italic('P. amboinensis')), cex=1.5, side=3, line=1)
  axis(1, at = c(1,2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  text(1.5,-3, labels='2014', cex=2)
  text(0.6, 58, labels = 'D', cex=3.4, font=2)
  
  
  points(x = 1, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 1, x1 = 1, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  #points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=16, col = 'dodgerblue2')
  segments(x0 = 2, x1 = 2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  
  # lemon damsel
  plot(NA, ylim=c(-3, ymax), xlim=c(0.4, 2.6), xaxt = 'n', yaxt = 'n')
  
  rect(-5000, -5000, 10000, 10000, col='grey92')
  abline(v = c(1:2), lwd=3.6, col='white')
  abline(h = c(ylab), lwd=3.6, col='white')
  #points(activity ~ jitter(as.numeric(treatment), 1.2), data=oaa[oaa$species=='lemon',], lwd=1, pch=1)
  foo <- oaa[oaa$species=='lemon',]
  
  points(activity ~ jitter(as.numeric(treatment), 10), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment), 5), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  text(0.6, 58, labels = 'E', cex=3.4, font=2)
  
  # mtext(expression(italic('P. moluccensis')), cex=1.5, side=3, line=1)
  axis(1, at = c(1,2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  text(1.5,-3, labels='2014', cex=2)
  
  
  points(x = 1, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 1, x1 = 1, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  #points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=16, col = 'dodgerblue2')
  segments(x0 = 2, x1 = 2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')  
  points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  
  
  # white damsels
  plot(NA, ylim=c(-3, ymax), xlim=c(0.4, 2.6), xaxt = 'n', yaxt = 'n')
  
  rect(-5000, -5000, 10000, 10000, col='grey92')
  abline(v = c(1:2), lwd=3.6, col='white')
  abline(h = c(ylab), lwd=3.6, col='white')
  #points(activity ~ jitter(as.numeric(treatment), 1.2), data=oaa[oaa$species=='whitedams',], lwd=1, pch=1)
  foo <- oaa[oaa$species=='whitedams',]
  
  points(activity ~ jitter(as.numeric(treatment), 10), data=foo[foo$treatment=='control',], lwd=1, pch=18, col = 'grey40', cex=1.1)
  points(activity ~ jitter(as.numeric(treatment), 5), data=foo[foo$treatment=='CO2',], lwd=1, pch=16, col = 'dodgerblue2')
  
  #mtext(expression(italic('D. perspicillatus')), cex=1.5, side=3, line=1)
  axis(1, at = c(1,2), labels = c('Control', expression(paste('CO'[2]))), lty=0, tick=NA, cex.axis=2, las=2)
  
  text(0.6, 58, labels = 'F', cex=3.4, font=2)
  text(1.5,-3, labels='2016', cex=2)
  
  points(x = 1, y = mean(foo[foo$treatment=='control',]$activity), cex=2.6, pch=16, col = 'grey40')
  segments(x0 = 1, x1 = 1, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  segments(x0 = 0.8, x1 = 1.2, 
           y0 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)), 
           y1 = mean(foo[foo$treatment=='control',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='control',]$activity)),
           lwd = 4, col = 'grey40')
  
  #points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=16, col = 'dodgerblue2')
  segments(x0 = 2, x1 = 2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) - (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  
  segments(x0 = 1.8, x1 = 2.2, 
           y0 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)), 
           y1 = mean(foo[foo$treatment=='CO2',]$activity) + (qnorm(0.975)*std.error(foo[foo$treatment=='CO2',]$activity)),
           lwd = 4, col = 'dodgerblue2')
  points(x = 2, y = mean(foo[foo$treatment=='CO2',]$activity), cex=2.6, pch=21, col = 'dodgerblue2', bg='white', lwd=4)
  
  
  
  
  dev.off()
}



















# visualize models that include effects of standard length
# first up, showing acantho models for 2014 and 2015 (nothing significant in 2016)



# acantho stuff first. excluding 2016 cuz nothing going on there.


foo <- oaa[oaa$species=='acantho',]

summary(ac110) # acanthos 2015
summary(ac101) # acanthos 2014


range(foo[foo$loc=='LIRS 2014',]$sl)

# script for predicting from models # 2014 acanthos based on SL range
SL_seq <- seq(from=50, to=92, by=1) #SL
co2_seq <- rep('CO2', length(SL_seq))
control_seq <- rep('control', length(SL_seq))

treatment_seq <- c(control_seq, co2_seq)
sl_seq <- rep(SL_seq, 2)

length(sl_seq)
length(treatment_seq)

fake_data <- data.frame(treatment = treatment_seq, sl = sl_seq) #merge it all into a dataframe
str(fake_data) #cool, a fake dataset to predict from.
fake_data$treatment <- relevel(fake_data$treatment, ref = 'control')

pdat <- predict(ac101, fake_data, se.fit=TRUE)
activity_fit <- data.frame(treatment = treatment_seq, sl=sl_seq, fit = pdat$fit, se = pdat$se.fit)
activity_fit$lci <- activity_fit$fit - (1.96*activity_fit$se) #calculate lower and upper confidence intervals from the standard errors
activity_fit$uci <- activity_fit$fit + (1.96*activity_fit$se)

ac2014fit <- activity_fit



# script for predicting from models # 2015 acanthos.
SL_seq <- seq(from=8, to=14, by=0.5) #SL
co2_seq <- rep('CO2', length(SL_seq))
control_seq <- rep('control', length(SL_seq))

treatment_seq <- c(control_seq, co2_seq)
sl_seq <- rep(SL_seq, 2)

length(sl_seq)
length(treatment_seq)

fake_data <- data.frame(treatment = treatment_seq, sl = sl_seq) #merge it all into a dataframe
str(fake_data) #cool, a fake dataset to predict from.
fake_data$treatment <- relevel(fake_data$treatment, ref = 'control')

pdat <- predict(ac110, fake_data, se.fit=TRUE)
activity_fit <- data.frame(treatment = treatment_seq, sl=sl_seq, fit = pdat$fit, se = pdat$se.fit)
activity_fit$lci <- activity_fit$fit - (1.96*activity_fit$se) #calculate lower and upper confidence intervals from the standard errors
activity_fit$uci <- activity_fit$fit + (1.96*activity_fit$se)

ac2015fit <- activity_fit






# time to make a plot.

jpeg('acantho activity models_20190304.jpg', width=14, height=7, units='in', res=600)

{
  #define y and x label sequences
  xlab <- seq(from=5, to=95, by=10)
  
  #locations where y-labels will go.
  ylab <- seq(from=0, to=60, by=10)
  
  limX <- c(45, 95)
  limY <- c(-5, 60)
  
  par(mar=c(2,0.5,0.5,0.5)) #plot margins (bottom, left, top, right)
  par(mfrow=c(1,2)) #number of rows and columns in for having murbtiple plots
  par(oma=c(6,8,0,0.5)) #outer margin areas (bottom, left, top, right)
  
  # LIRS 2014
  #empty plot area
  plot(NA, xlim=limX, ylim=limY, main=NULL, xlab=NULL, ylab=NULL, axes=FALSE)
  
  polygon(c(ac2014fit[ac2014fit$treatment=='control',]$sl, 
            rev(ac2014fit[ac2014fit$treatment=='control',]$sl)),
          c(ac2014fit[ac2014fit$treatment=='control',]$uci,
            rev(ac2014fit[ac2014fit$treatment=='control',]$lci)), col = alpha("grey75", 0.5), border = NA)
  
  polygon(c(ac2014fit[ac2014fit$treatment=='CO2',]$sl, 
            rev(ac2014fit[ac2014fit$treatment=='CO2',]$sl)),
          c(ac2014fit[ac2014fit$treatment=='CO2',]$uci,
            rev(ac2014fit[ac2014fit$treatment=='CO2',]$lci)), col = alpha("grey75", 0.5), border = NA)
  
  points(activity ~ sl, data=foo[foo$treatment=='control' & foo$loc=='LIRS 2014',], cex=1.2, pch=18, col = 'grey40')
  
  points(activity ~ sl, data=foo[foo$treatment=='CO2' & foo$loc=='LIRS 2014',], cex=1.1, pch=16, col = 'dodgerblue2')
  
  lines(x = ac2014fit[ac2014fit$treatment=='CO2',]$sl, 
        y = ac2014fit[ac2014fit$treatment=='CO2',]$fit, lwd=3, col = 'dodgerblue2')
  lines(x = ac2014fit[ac2014fit$treatment=='control',]$sl, 
        y = ac2014fit[ac2014fit$treatment=='control',]$fit, lwd=3, col = 'grey40')
  
  axis(1, at = xlab, labels=xlab, cex.axis=1.6)
  axis(2, at = ylab, labels=ylab, cex.axis=1.6, las=2)
  
  mtext(expression(paste('Activity duration (', 's min'^'-1', ')')), side=2, line=3.25, cex=2)
  mtext('Standard length (mm)', side=1, at = 95, line=4, cex=2)
  
  text(x = 85, y = 55, labels = 'LIRS 2014', cex=2)
  
  box(lwd=1.6)
  
  
  
  # 2016
  
  limX <- c(8,14)
  
  plot(NA, xlim=limX, ylim=limY, main=NULL, xlab=NULL, ylab=NULL, axes=FALSE)
  
  polygon(c(ac2015fit[ac2015fit$treatment=='control',]$sl, 
            rev(ac2015fit[ac2015fit$treatment=='control',]$sl)),
          c(ac2015fit[ac2015fit$treatment=='control',]$uci,
            rev(ac2015fit[ac2015fit$treatment=='control',]$lci)), col = alpha("grey75", 0.5), border = NA)
  
  polygon(c(ac2015fit[ac2015fit$treatment=='CO2',]$sl, 
            rev(ac2015fit[ac2015fit$treatment=='CO2',]$sl)),
          c(ac2015fit[ac2015fit$treatment=='CO2',]$uci,
            rev(ac2015fit[ac2015fit$treatment=='CO2',]$lci)), col = alpha("grey75", 0.5), border = NA)
  
  points(activity ~ sl, data=foo[foo$treatment=='control' & foo$loc=='AIMS 2015',], cex=1.2, pch=18, col = 'grey40')
  points(activity ~ sl, data=foo[foo$treatment=='CO2' & foo$loc=='AIMS 2015',], cex=1.1, pch=16, col = 'dodgerblue2')
  
  lines(x = ac2015fit[ac2015fit$treatment=='CO2',]$sl, 
        y = ac2015fit[ac2015fit$treatment=='CO2',]$fit, lwd=3, col = 'dodgerblue2')
  lines(x = ac2015fit[ac2015fit$treatment=='control',]$sl, 
        y = ac2015fit[ac2015fit$treatment=='control',]$fit, lwd=3, col = 'grey40')
  
  axis(1, at = c(8:14), labels=c(8:14), cex.axis=1.6)
  text(x = 13, y = 55, labels = 'AIMS 2015', cex=2)
  
  box(lwd=1.6)
  
  dev.off()
}












# now, plot of humbug 2014 model

foo <- oaa[oaa$species=='humbug',]
foo <- foo[foo$loc=='LIRS 2014',]
summary(hum101)

range(foo$sl)

# script for predicting from 2014 humbugs
SL_seq <- seq(from=17, to=71, by=1) #SL
co2_seq <- rep('CO2', length(SL_seq))
control_seq <- rep('control', length(SL_seq))

treatment_seq <- c(control_seq, co2_seq)
sl_seq <- rep(SL_seq, 2)

length(sl_seq)
length(treatment_seq)

fake_data <- data.frame(treatment = treatment_seq, sl = sl_seq) #merge it all into a dataframe
str(fake_data) #cool, a fake dataset to predict from.
fake_data$treatment <- relevel(fake_data$treatment, ref = 'control')

pdat <- predict(hum101, fake_data, se.fit=TRUE)
activity_fit <- data.frame(treatment = treatment_seq, sl=sl_seq, fit = pdat$fit, se = pdat$se.fit)
activity_fit$lci <- activity_fit$fit - (1.96*activity_fit$se) #calculate lower and upper confidence intervals from the standard errors
activity_fit$uci <- activity_fit$fit + (1.96*activity_fit$se)

hum2014fit <- activity_fit




# humbug plot

jpeg('humbug activity model 2014_20190304.jpg', width=10, height=7, units='in', res=600)

{
  #define y and x label sequences
  xlab <- seq(from=20, to=70, by=10)
  
  #locations where y-labels will go.
  ylab <- seq(from=0, to=60, by=10)
  
  limX <- c(17, 71)
  limY <- c(0, 60)
  
  par(mar=c(2,0.5,0.5,0.5)) #plot margins (bottom, left, top, right)
  par(mfrow=c(1,1)) #number of rows and columns in for having murbtiple plots
  par(oma=c(6,8,0,0.5)) #outer margin areas (bottom, left, top, right)
  
  # LIRS 2014
  #empty plot area
  plot(NA, xlim=limX, ylim=limY, main=NULL, xlab=NULL, ylab=NULL, axes=FALSE)
  
  polygon(c(hum2014fit[hum2014fit$treatment=='control',]$sl, 
            rev(hum2014fit[hum2014fit$treatment=='control',]$sl)),
          c(hum2014fit[hum2014fit$treatment=='control',]$uci,
            rev(hum2014fit[hum2014fit$treatment=='control',]$lci)), col = alpha("grey75", 0.5), border = NA)
  
  polygon(c(hum2014fit[hum2014fit$treatment=='CO2',]$sl, 
            rev(hum2014fit[hum2014fit$treatment=='CO2',]$sl)),
          c(hum2014fit[hum2014fit$treatment=='CO2',]$uci,
            rev(hum2014fit[hum2014fit$treatment=='CO2',]$lci)), col = alpha("grey75", 0.5), border = NA)
  
  points(activity ~ sl, data=foo[foo$treatment=='control' & foo$loc=='LIRS 2014',], cex=1.2, pch=18, col = 'grey40')
  
  points(activity ~ sl, data=foo[foo$treatment=='CO2' & foo$loc=='LIRS 2014',], cex=1.1, pch=16, col = 'dodgerblue2')
  
  lines(x = hum2014fit[hum2014fit$treatment=='CO2',]$sl, 
        y = hum2014fit[hum2014fit$treatment=='CO2',]$fit, lwd=3, col = 'dodgerblue2')
  lines(x = hum2014fit[hum2014fit$treatment=='control',]$sl, 
        y = hum2014fit[hum2014fit$treatment=='control',]$fit, lwd=3, col = 'grey40')
  
  axis(1, at = xlab, labels=xlab, cex.axis=1.6)
  axis(2, at = ylab, labels=ylab, cex.axis=1.6, las=2)
  
  mtext(expression(paste('Activity duration (', 's min'^'-1', ')')), side=2, line=3.25, cex=2)
  mtext('Standard length (mm)', side=1, line=4, cex=2)
  
  box(lwd=1.6)
  
  dev.off()
}







# end

