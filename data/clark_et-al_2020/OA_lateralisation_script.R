# R script for analyses of lateralisation #

rm(list=ls())

#setwd("")

require(lme4)
require(ggplot2)
require(ggpubr)
require(plyr)

# define fixed parameters (can be changed by the user)
 n=10 #number of trials/individual
 P=0.5 #theoretical probability of success (null hypothesis)


#1.Acanthochromis polyacanthus (AIMS 2015 NO OFFSET)
###################################################################################################################

#import data
 aims<-read.csv("AIMS 2015 lat data.csv", header=T)
 head(aims)
 length(aims$treatment)
#subset number of right turns for control fish
 Control_1<-aims[aims$treatment=="Control",c(3,7)]
 colnames(Control_1)=c("ind_1c","X_1c")
 Control_1$ind_1c=as.factor(Control_1$ind_1c)
 X1_1c<-tapply(Control_1$X_1c,Control_1$ind_1c,sum)
 length(X1_1c)
 hist(X1_1c)
 mean(X1_1c)
#subset number of right turns for treatment fish
 CO2_1<-aims[aims$treatment=="CO2",c(3,7)]
 colnames(CO2_1)=c("ind_1t","X_1t")
 CO2_1$ind_1t=as.factor(CO2_1$ind_1t)
 X1_1t<-tapply(CO2_1$X_1t,CO2_1$ind_1t,sum)
 length(X1_1t)
 hist(X1_1t)
 mean(X1_1t)

#Figure
fig1_a<-c(rep("control",120), rep("high CO2",104))
fig1_b<-c(as.vector(X1_1c), as.vector(X1_1t))
fig_1_dat<-data.frame(cbind(fig1_a, fig1_b))
	colnames(fig_1_dat)=c("group","Rturns")
	fig_1_dat$Rturns<-as.numeric(as.character(fig_1_dat$Rturns))   
	range(fig_1_dat$Rturns)
mu_1<-ddply(fig_1_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_1_dat[fig_1_dat$group=="control",]$Rturns)
table(fig_1_dat[fig_1_dat$group=="high CO2",]$Rturns)
table(fig_1_dat$Rturns)
fig_1<-ggplot(fig_1_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=F) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,30), breaks=c(0,10,20,30), labels=c("0","10","20","30"), name="Number of individuals") +
	geom_vline(data=mu_1, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 28, label = "A. polyacanthus") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_1.eps", width = 2.5, height = 2.5)
fig_1 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

###############
###Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_1c=((length(X1_1c)-1)*var(X1_1c)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_1c<-pchisq(chi_sq_1c,df=(length(X1_1c)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_1c=(mean(X1_1c/n)-P)/(.5*.5/(length(X1_1c)*n))^.5			#compute a Z score for the mean
 pZ_1c=pnorm(abs(Z_1c),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_1c<-glmer(X_1c~1+(1|ind_1c),data=Control_1,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_1c<-glm(X_1c~1,data=Control_1,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_1c=summary(ga_1c)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_1c=anova(ga_1c,gb_1c)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_1c<-round(matrix(c(pr_m_1c,pr_v_1c,pZ_1c,sig_v_1c),ncol=4,nrow=1),3); colnames(out_1c)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_1c

###############
###Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_1t=((length(X1_1t)-1)*var(X1_1t)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_1t<-pchisq(chi_sq_1t,df=(length(X1_1t)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_1t=(mean(X1_1t/n)-P)/(.5*.5/(length(X1_1t)*n))^.5			#compute a Z score for the mean
 pZ_1t=pnorm(abs(Z_1t),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_1t<-glmer(X_1t~1+(1|ind_1t),data=CO2_1,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_1t<-glm(X_1t~1,data=CO2_1,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_1t=summary(ga_1t)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_1t=anova(ga_1t,gb_1t)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_1t<-matrix(c(pr_m_1t,pr_v_1t,pZ_1t,sig_v_1t),ncol=4,nrow=1); colnames(out_1t)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_1t

# N.B.for a comparison of wild and RHQ (i.e. aquarium reared) A. polyacanthus, see the end of this script.





#2.Acanthochromis polyacanthus (AIMS 2015 OFFSET)
###################################################################################################################

#import data
 aims<-read.csv("AIMS 2015 lat data.csv", header=T)
#subset number of right turns for control fish
 Control_2<-aims[aims$treatment=="Control",c(3,6)]
 colnames(Control_2)=c("ind_2c","X_2c")
 Control_2$ind_2c=as.factor(Control_2$ind_2c)
 X1_2c<-tapply(Control_2$X_2c,Control_2$ind_2c,sum)
 length(X1_2c)
 hist(X1_2c)
 mean(X1_2c)
 range(X1_2c)
#subset number of right turns for treatment fish
 CO2_2<-aims[aims$treatment=="CO2",c(3,6)]
 colnames(CO2_2)=c("ind_2t","X_2t")
 CO2_2$ind_2t=as.factor(CO2_2$ind_2t)
 X1_2t<-tapply(CO2_2$X_2t,CO2_2$ind_2t,sum)
 length(X1_2t)
 hist(X1_2t)
 mean(X1_2t)
 range(X1_2t)

#Figure
fig2_a<-c(rep("control",120), rep("high CO2",104))
fig2_b<-c(as.vector(X1_2c), as.vector(X1_2t))
fig_2_dat<-data.frame(cbind(fig2_a, fig2_b))
	colnames(fig_2_dat)=c("group","Rturns")
	fig_2_dat$Rturns<-as.numeric(as.character(fig_2_dat$Rturns))   
	range(fig_2_dat$Rturns)
mu_2<-ddply(fig_2_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_2_dat[fig_2_dat$group=="control",]$Rturns)
table(fig_2_dat[fig_2_dat$group=="high CO2",]$Rturns)
fig_2<-ggplot(fig_2_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=T) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,30), breaks=c(0,10,20,30), labels=c("0","10","20","30"), name="Number of individuals") +
	geom_vline(data=mu_2, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) + 
	annotate("text", x = 8, y = 28, label = "A. polyacanthus (offset barrier)") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_2.eps", width = 2.5, height = 2.5)
fig_2 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

###############
###Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_2c=((length(X1_2c)-1)*var(X1_2c)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_2c<-pchisq(chi_sq_2c,df=(length(X1_2c)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_2c=(mean(X1_2c/n)-P)/(.5*.5/(length(X1_2c)*n))^.5			#compute a Z score for the mean
 pZ_2c=pnorm(abs(Z_2c),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_2c<-glmer(X_2c~1+(1|ind_2c),data=Control_2,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_2c<-glm(X_2c~1,data=Control_2,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_2c=summary(ga_2c)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_2c=anova(ga_2c,gb_2c)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_2c<-matrix(c(pr_m_2c,pr_v_2c,pZ_2c,sig_v_2c),ncol=4,nrow=1); colnames(out_2c)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_2c

###############
###Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_2t=((length(X1_2t)-1)*var(X1_2t)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_2t<-pchisq(chi_sq_2t,df=(length(X1_2t)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_2t=(mean(X1_2t/n)-P)/(.5*.5/(length(X1_2t)*n))^.5			#compute a Z score for the mean
 pZ_2t=pnorm(abs(Z_2t),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_2t<-glmer(X_2t~1+(1|ind_2t),data=CO2_2,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_2t<-glm(X_2t~1,data=CO2_2,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_2t=summary(ga_2t)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_2t=anova(ga_2t,gb_2t)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_2t<-matrix(c(pr_m_2t,pr_v_2t,pZ_2t,sig_v_2t),ncol=4,nrow=1); colnames(out_2t)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_2t





#3.Pomacentrus amboinensis (LIRS 2014)
###################################################################################################################

#import data
 lirs<-read.csv("LIRS 2014 lat data.csv", header=T)
#subset for Pomacentrus amboinensis
 lirs_Pa<-lirs[lirs$species=="Ambon",]
#subset number of right turns for control fish
 Control_3<-lirs_Pa[lirs_Pa$treatment=="Control",c(3,15)]
 colnames(Control_3)=c("ind_3c","X_3c")
 Control_3$ind_3c=as.factor(Control_3$ind_3c)
 X1_3c<-tapply(Control_3$X_3c,Control_3$ind_3c,sum)
 length(X1_3c)
 hist(X1_3c)
 mean(X1_3c)
#subset number of right turns for treatment fish
 CO2_3<-lirs_Pa[lirs_Pa$treatment=="CO2",c(3,15)]
 colnames(CO2_3)=c("ind_3t","X_3t")
 CO2_3$ind_3t=as.factor(CO2_3$ind_3t)
 X1_3t<-tapply(CO2_3$X_3t,CO2_3$ind_3t,sum)
 length(X1_3t)
 hist(X1_3t)
 mean(X1_3t)

#Figure
fig3_a<-c(rep("control",21), rep("high CO2",22))
fig3_b<-c(as.vector(X1_3c), as.vector(X1_3t))
fig_3_dat<-data.frame(cbind(fig3_a, fig3_b))
	colnames(fig_3_dat)=c("group","Rturns")
	fig_3_dat$Rturns<-as.numeric(as.character(fig_3_dat$Rturns))   
	range(fig_3_dat$Rturns)
mu_3<-ddply(fig_3_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_3_dat[fig_3_dat$group=="control",]$Rturns)
table(fig_3_dat[fig_3_dat$group=="high CO2",]$Rturns)
fig_3<-ggplot(fig_3_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=T) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,10), breaks=c(0,2,4,6,8,10), labels=c("0","2","4","6","8","10"), name="Number of individuals") +
	geom_vline(data=mu_3, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 9, label = "P. amboinensis") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_3.eps", width = 2.5, height = 2.5)
fig_3 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

###############
###Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_3c=((length(X1_3c)-1)*var(X1_3c)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_3c<-pchisq(chi_sq_3c,df=(length(X1_3c)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_3c=(mean(X1_3c/n)-P)/(.5*.5/(length(X1_3c)*n))^.5			#compute a Z score for the mean
 pZ_3c=pnorm(abs(Z_3c),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_3c<-glmer(X_3c~1+(1|ind_3c),data=Control_3,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_3c<-glm(X_3c~1,data=Control_3,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_3c=summary(ga_3c)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_3c=anova(ga_3c,gb_3c)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_3c<-matrix(c(pr_m_3c,pr_v_3c,pZ_3c,sig_v_3c),ncol=4,nrow=1); colnames(out_3c)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_3c

###############
###Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_3t=((length(X1_3t)-1)*var(X1_3t)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_3t<-pchisq(chi_sq_3t,df=(length(X1_3t)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_3t=(mean(X1_3t/n)-P)/(.5*.5/(length(X1_3t)*n))^.5			#compute a Z score for the mean
 pZ_3t=pnorm(abs(Z_3t),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_3t<-glmer(X_3t~1+(1|ind_3t),data=CO2_3,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_3t<-glm(X_3t~1,data=CO2_3,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_3t=summary(ga_3t)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_3t=anova(ga_3t,gb_3t)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_3t<-matrix(c(pr_m_3t,pr_v_3t,pZ_3t,sig_v_3t),ncol=4,nrow=1); colnames(out_3t)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_3t






#4.Chromis atripectoralis (LIRS 2014)
###################################################################################################################

#import data
 lirs<-read.csv("LIRS 2014 lat data.csv", header=T)
#subset for Chromis
 lirs_Cv<-lirs[lirs$species=="Chromis",]
#subset number of right turns for control fish
 Control_4<-lirs_Cv[lirs_Cv$treatment=="Control",c(3,15)]
 colnames(Control_4)=c("ind_4c","X_4c")
 Control_4$ind_4c=as.factor(Control_4$ind_4c)
 X1_4c<-tapply(Control_4$X_4c,Control_4$ind_4c,sum)
 length(X1_4c)
 hist(X1_4c)
 mean(X1_4c)
#subset number of right turns for treatment fish
 CO2_4<-lirs_Cv[lirs_Cv$treatment=="CO2",c(3,15)]
 colnames(CO2_4)=c("ind_4t","X_4t")
 CO2_4$ind_4t=as.factor(CO2_4$ind_4t)
 X1_4t<-tapply(CO2_4$X_4t,CO2_4$ind_4t,sum)
 length(X1_4t)
 hist(X1_4t)
 mean(X1_4t)

#Figure
fig4_a<-c(rep("control",26), rep("high CO2",17))
fig4_b<-c(as.vector(X1_4c), as.vector(X1_4t))
fig_4_dat<-data.frame(cbind(fig4_a, fig4_b))
	colnames(fig_4_dat)=c("group","Rturns")
	fig_4_dat$Rturns<-as.numeric(as.character(fig_4_dat$Rturns))   
	range(fig_4_dat$Rturns)
mu_4<-ddply(fig_4_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_4_dat[fig_4_dat$group=="control",]$Rturns)
table(fig_4_dat[fig_4_dat$group=="high CO2",]$Rturns)
fig_4<-ggplot(fig_4_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=F) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,8), breaks=c(0,2,4,6,8), labels=c("0","2","4","6","8"), name="Number of individuals") +
	geom_vline(data=mu_4, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 7.5, label = "C. atripectoralis") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_4.eps", width = 2.5, height = 2.5)
fig_4 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

###############
###Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_4c=((length(X1_4c)-1)*var(X1_4c)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_4c<-pchisq(chi_sq_4c,df=(length(X1_4c)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_4c=(mean(X1_4c/n)-P)/(.5*.5/(length(X1_4c)*n))^.5			#compute a Z score for the mean
 pZ_4c=pnorm(abs(Z_4c),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_4c<-glmer(X_4c~1+(1|ind_4c),data=Control_4,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_4c<-glm(X_4c~1,data=Control_4,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_4c=summary(ga_4c)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_4c=anova(ga_4c,gb_4c)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_4c<-matrix(c(pr_m_4c,pr_v_4c,pZ_4c,sig_v_4c),ncol=4,nrow=1); colnames(out_4c)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_4c

###############
###Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_4t=((length(X1_4t)-1)*var(X1_4t)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_4t<-pchisq(chi_sq_4t,df=(length(X1_4t)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_4t=(mean(X1_4t/n)-P)/(.5*.5/(length(X1_4t)*n))^.5			#compute a Z score for the mean
 pZ_4t=pnorm(abs(Z_4t),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_4t<-glmer(X_4t~1+(1|ind_4t),data=CO2_4,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_4t<-glm(X_4t~1,data=CO2_4,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_4t=summary(ga_4t)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_4t=anova(ga_4t,gb_4t)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_4t<-matrix(c(pr_m_4t,pr_v_4t,pZ_4t,sig_v_4t),ncol=4,nrow=1); colnames(out_4t)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_4t






#5.Dascyllus aruanus (LIRS 2014)
###################################################################################################################

#import data
 lirs<-read.csv("LIRS 2014 lat data.csv", header=T)
#subset for Dascyllus aruanus
 lirs_Da<-lirs[lirs$species=="Humbug",]
#subset number of right turns for control fish
 Control_5<-lirs_Da[lirs_Da$treatment=="Control",c(3,15)]
 colnames(Control_5)=c("ind_5c","X_5c")
 Control_5$ind_5c=as.factor(Control_5$ind_5c)
 X1_5c<-tapply(Control_5$X_5c,Control_5$ind_5c,sum)
 length(X1_5c)
 hist(X1_5c)
 mean(X1_5c)
#subset number of right turns for treatment fish
 CO2_5<-lirs_Da[lirs_Da$treatment=="CO2",c(3,15)]
 colnames(CO2_5)=c("ind_5t","X_5t")
 CO2_5$ind_5t=as.factor(CO2_5$ind_5t)
 X1_5t<-tapply(CO2_5$X_5t,CO2_5$ind_5t,sum)
 length(X1_5t)
 hist(X1_5t)
 mean(X1_5t)

#Figure
fig5_a<-c(rep("control",19), rep("high CO2",21))
fig5_b<-c(as.vector(X1_5c), as.vector(X1_5t))
fig_5_dat<-data.frame(cbind(fig5_a, fig5_b))
	colnames(fig_5_dat)=c("group","Rturns")
	fig_5_dat$Rturns<-as.numeric(as.character(fig_5_dat$Rturns))   
	range(fig_5_dat$Rturns)
mu_5<-ddply(fig_5_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_5_dat[fig_5_dat$group=="control",]$Rturns)
table(fig_5_dat[fig_5_dat$group=="high CO2",]$Rturns)
fig_5<-ggplot(fig_5_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=F) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,8), breaks=c(0,2,4,6,8), labels=c("0","2","4","6","8"), name="Number of individuals") +
	geom_vline(data=mu_5, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 7.5, label = "D. aruanus") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_5.eps", width = 2.5, height = 2.5)
fig_5 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

###############
###Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_5c=((length(X1_5c)-1)*var(X1_5c)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_5c<-pchisq(chi_sq_5c,df=(length(X1_5c)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_5c=(mean(X1_5c/n)-P)/(.5*.5/(length(X1_5c)*n))^.5			#compute a Z score for the mean
 pZ_5c=pnorm(abs(Z_5c),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_5c<-glmer(X_5c~1+(1|ind_5c),data=Control_5,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_5c<-glm(X_5c~1,data=Control_5,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_5c=summary(ga_5c)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_5c=anova(ga_5c,gb_5c)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_5c<-matrix(c(pr_m_5c,pr_v_5c,pZ_5c,sig_v_5c),ncol=4,nrow=1); colnames(out_5c)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_5c

###############
###Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_5t=((length(X1_5t)-1)*var(X1_5t)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_5t<-pchisq(chi_sq_5t,df=(length(X1_5t)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_5t=(mean(X1_5t/n)-P)/(.5*.5/(length(X1_5t)*n))^.5			#compute a Z score for the mean
 pZ_5t=pnorm(abs(Z_5t),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_5t<-glmer(X_5t~1+(1|ind_5t),data=CO2_5,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_5t<-glm(X_5t~1,data=CO2_5,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_5t=summary(ga_5t)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_5t=anova(ga_5t,gb_5t)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_5t<-matrix(c(pr_m_5t,pr_v_5t,pZ_5t,sig_v_5t),ncol=4,nrow=1); colnames(out_5t)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_5t






#6.Pomacentrus moluccensis (LIRS 2014)
###################################################################################################################

#import data
 lirs<-read.csv("LIRS 2014 lat data.csv", header=T)
#subset for Pomacentrus moluccensis
 lirs_Pm<-lirs[lirs$species=="Lemon",]
#subset number of right turns for control fish
 Control_6<-lirs_Pm[lirs_Pm$treatment=="Control",c(3,15)]
 colnames(Control_6)=c("ind_6c","X_6c")
 Control_6$ind_6c=as.factor(Control_6$ind_6c)
 X1_6c<-tapply(Control_6$X_6c,Control_6$ind_6c,sum)
 length(X1_6c)
 hist(X1_6c)
 mean(X1_6c)
#subset number of right turns for treatment fish
 CO2_6<-lirs_Pm[lirs_Pm$treatment=="CO2",c(3,15)]
 colnames(CO2_6)=c("ind_6t","X_6t")
 CO2_6$ind_6t=as.factor(CO2_6$ind_6t)
 X1_6t<-tapply(CO2_6$X_6t,CO2_6$ind_6t,sum)
 length(X1_6t)
 hist(X1_6t)
 mean(X1_6t)

#Figure
fig6_a<-c(rep("control",29), rep("high CO2",20))
fig6_b<-c(as.vector(X1_6c), as.vector(X1_6t))
fig_6_dat<-data.frame(cbind(fig6_a, fig6_b))
	colnames(fig_6_dat)=c("group","Rturns")
	fig_6_dat$Rturns<-as.numeric(as.character(fig_6_dat$Rturns))   
	range(fig_6_dat$Rturns)
mu_6<-ddply(fig_6_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_6_dat[fig_6_dat$group=="control",]$Rturns)
table(fig_6_dat[fig_6_dat$group=="high CO2",]$Rturns)
fig_6<-ggplot(fig_6_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=F) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,8), breaks=c(0,2,4,6,8), labels=c("0","2","4","6","10"), name="Number of individuals") +
	geom_vline(data=mu_6, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 7.5, label = "P. moluccensis") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_6.eps", width = 2.5, height = 2.5)
fig_6 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

###############
###Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_6c=((length(X1_6c)-1)*var(X1_6c)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_6c<-pchisq(chi_sq_6c,df=(length(X1_6c)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_6c=(mean(X1_6c/n)-P)/(.5*.5/(length(X1_6c)*n))^.5			#compute a Z score for the mean
 pZ_6c=pnorm(abs(Z_6c),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_6c<-glmer(X_6c~1+(1|ind_6c),data=Control_6,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_6c<-glm(X_6c~1,data=Control_6,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_6c=summary(ga_6c)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_6c=anova(ga_6c,gb_6c)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_6c<-matrix(c(pr_m_6c,pr_v_6c,pZ_6c,sig_v_6c),ncol=4,nrow=1); colnames(out_6c)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_6c

###############
###Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_6t=((length(X1_6t)-1)*var(X1_6t)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_6t<-pchisq(chi_sq_6t,df=(length(X1_6t)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_6t=(mean(X1_6t/n)-P)/(.5*.5/(length(X1_6t)*n))^.5			#compute a Z score for the mean
 pZ_6t=pnorm(abs(Z_6t),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_6t<-glmer(X_6t~1+(1|ind_6t),data=CO2_6,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_6t<-glm(X_6t~1,data=CO2_6,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_6t=summary(ga_6t)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_6t=anova(ga_6t,gb_6t)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_6t<-matrix(c(pr_m_6t,pr_v_6t,pZ_6t,sig_v_6t),ncol=4,nrow=1); colnames(out_6t)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_6t






#7.Pomacentrus amboinensis REDONE (LIRS 2014)
###################################################################################################################

#import data
 lirs_PaRed<-read.csv("LIRS 2014 lat data ambo redone.csv", header=T)
#subset number of right turns for control fish
 Control_7<-lirs_PaRed[lirs_PaRed$treatment=="Control",c(3,13)]
 colnames(Control_7)=c("ind_7c","X_7c")
 Control_7$ind_7c=as.factor(Control_7$ind_7c)
 X1_7c<-tapply(Control_7$X_7c,Control_7$ind_7c,sum)
 length(X1_7c)
 hist(X1_7c)
 mean(X1_7c)
#subset number of right turns for treatment fish
 CO2_7<-lirs_PaRed[lirs_PaRed$treatment=="CO2",c(3,13)]
 colnames(CO2_7)=c("ind_7t","X_7t")
 CO2_7$ind_7t=as.factor(CO2_7$ind_7t)
 X1_7t<-tapply(CO2_7$X_7t,CO2_7$ind_7t,sum)
 length(X1_7t)
 hist(X1_7t)
 mean(X1_7t)

#Figure
fig7_a<-c(rep("control",15), rep("high CO2",15))
fig7_b<-c(as.vector(X1_7c), as.vector(X1_7t))
fig_7_dat<-data.frame(cbind(fig7_a, fig7_b))
	colnames(fig_7_dat)=c("group","Rturns")
	fig_7_dat$Rturns<-as.numeric(as.character(fig_7_dat$Rturns))   
	range(fig_7_dat$Rturns)
mu_7<-ddply(fig_7_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_7_dat[fig_7_dat$group=="control",]$Rturns)
table(fig_7_dat[fig_7_dat$group=="high CO2",]$Rturns)
fig_7<-ggplot(fig_7_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=T) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,6), breaks=c(0,2,4,6), labels=c("0","2","4","6"), name="Number of individuals") +
	geom_vline(data=mu_7, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 5.5, label = "P. amboinensis (retested)") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_7.eps", width = 2.5, height = 2.5)
fig_7 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

###############
###Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_7c=((length(X1_7c)-1)*var(X1_7c)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_7c<-pchisq(chi_sq_7c,df=(length(X1_7c)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_7c=(mean(X1_7c/n)-P)/(.5*.5/(length(X1_7c)*n))^.5			#compute a Z score for the mean
 pZ_7c=pnorm(abs(Z_7c),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_7c<-glmer(X_7c~1+(1|ind_7c),data=Control_7,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_7c<-glm(X_7c~1,data=Control_7,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_7c=summary(ga_7c)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_7c=anova(ga_7c,gb_7c)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_7c<-matrix(c(pr_m_7c,pr_v_7c,pZ_7c,sig_v_7c),ncol=4,nrow=1); colnames(out_7c)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_7c

###############
###Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_7t=((length(X1_7t)-1)*var(X1_7t)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_7t<-pchisq(chi_sq_7t,df=(length(X1_7t)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_7t=(mean(X1_7t/n)-P)/(.5*.5/(length(X1_7t)*n))^.5			#compute a Z score for the mean
 pZ_7t=pnorm(abs(Z_7t),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_7t<-glmer(X_7t~1+(1|ind_7t),data=CO2_7,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_7t<-glm(X_7t~1,data=CO2_7,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_7t=summary(ga_7t)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_7t=anova(ga_7t,gb_7t)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_7t<-matrix(c(pr_m_7t,pr_v_7t,pZ_7t,sig_v_7t),ncol=4,nrow=1); colnames(out_7t)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_7t










#A polyacanthus wild Vs RHQ comparison (AIMS 2015)
###################################################################################################################

##########
#WILD fish
##########

#subset WILD FISH 
 aimsW <- aims[aims$population=="Wild", ]
#subset number of right turns for control fish
 Control_W<-aimsW[aimsW$treatment=="Control",c(3,7)]
 colnames(Control_W)=c("ind_Wc","X_Wc")
 Control_W$ind_1c=as.factor(Control_W$ind_Wc)
 XW_Wc<-tapply(Control_W$X_Wc,Control_W$ind_Wc,sum)
 length(XW_Wc)
 hist(XW_Wc)
 mean(XW_Wc)
#subset number of right turns for treatment fish
 CO2_W<-aimsW[aimsW$treatment=="CO2",c(3,7)]
 colnames(CO2_W)=c("ind_Wt","X_Wt")
 CO2_W$ind_Wt=as.factor(CO2_W$ind_Wt)
 XW_Wt<-tapply(CO2_W$X_Wt,CO2_W$ind_Wt,sum)
 length(XW_Wt)
 hist(XW_Wt)
 mean(XW_Wt)

#Figure
figW_a<-c(rep("control",54), rep("high CO2",42))
figW_b<-c(as.vector(XW_Wc), as.vector(XW_Wt))
fig_W_dat<-data.frame(cbind(figW_a, figW_b))
	colnames(fig_W_dat)=c("group","Rturns")
	fig_W_dat$Rturns<-as.numeric(as.character(fig_W_dat$Rturns))   
	range(fig_W_dat$Rturns)
mu_W<-ddply(fig_W_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_W_dat[fig_W_dat$group=="control",]$Rturns)
table(fig_W_dat[fig_W_dat$group=="high CO2",]$Rturns)
table(fig_W_dat$Rturns)
fig_W<-ggplot(fig_W_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=F) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,15), breaks=c(0,5,10,15), labels=c("0","5","10","15"), name="Number of individuals") +
	geom_vline(data=mu_W, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 14, label = "WILD A. polyacanthus") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_W.eps", width = 2.5, height = 2.5)
fig_W + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

#Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_Wc=((length(XW_Wc)-1)*var(XW_Wc)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_Wc<-pchisq(chi_sq_Wc,df=(length(XW_Wc)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_Wc=(mean(XW_Wc/n)-P)/(.5*.5/(length(XW_Wc)*n))^.5			#compute a Z score for the mean
 pZ_Wc=pnorm(abs(Z_Wc),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_Wc<-glmer(X_Wc~1+(1|ind_Wc),data=Control_W,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_Wc<-glm(X_Wc~1,data=Control_W,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_Wc=summary(ga_Wc)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_Wc=anova(ga_Wc,gb_Wc)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_Wc<-round(matrix(c(pr_m_Wc,pr_v_Wc,pZ_Wc,sig_v_Wc),ncol=4,nrow=1),3); colnames(out_Wc)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_Wc

#Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (CO2 fish)
 chi_sq_Wt=((length(XW_Wt)-1)*var(XW_Wt)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_Wt<-pchisq(chi_sq_Wt,df=(length(XW_Wt)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_Wt=(mean(XW_Wt/n)-P)/(.5*.5/(length(XW_Wt)*n))^.5			#compute a Z score for the mean
 pZ_Wt=pnorm(abs(Z_Wt),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (CO2 fish)
 ga_Wt<-glmer(X_Wt~1+(1|ind_Wt),data=CO2_W,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_Wt<-glm(X_Wt~1,data=CO2_W,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_Wt=summary(ga_Wt)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_Wt=anova(ga_Wt,gb_Wt)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_Wt<-matrix(c(pr_m_Wt,pr_v_Wt,pZ_Wt,sig_v_Wt),ncol=4,nrow=1); colnames(out_Wt)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_Wt


##########
#RHQ fish
##########

#subset WILD FISH 
length(aimsRHQ[,1])
 aimsRHQ <- aims[aims$population=="RHQ", ]
#subset number of right turns for control fish
 Control_RHQ<-aimsRHQ[aimsRHQ$treatment=="Control",c(3,7)]
 colnames(Control_RHQ)=c("ind_RHQc","X_RHQc")
 Control_RHQ$ind_1c=as.factor(Control_RHQ$ind_RHQc)
 XRHQ_RHQc<-tapply(Control_RHQ$X_RHQc,Control_RHQ$ind_RHQc,sum)
 length(XRHQ_RHQc)
 hist(XRHQ_RHQc)
 mean(XRHQ_RHQc)
#subset number of right turns for treatment fish
 CO2_RHQ<-aimsRHQ[aimsRHQ$treatment=="CO2",c(3,7)]
 colnames(CO2_RHQ)=c("ind_RHQt","X_RHQt")
 CO2_RHQ$ind_RHQt=as.factor(CO2_RHQ$ind_RHQt)
 XRHQ_RHQt<-tapply(CO2_RHQ$X_RHQt,CO2_RHQ$ind_RHQt,sum)
 length(XRHQ_RHQt)
 hist(XRHQ_RHQt)
 mean(XRHQ_RHQt)

#Figure
figRHQ_a<-c(rep("control",66), rep("high CO2",62))
figRHQ_b<-c(as.vector(XRHQ_RHQc), as.vector(XRHQ_RHQt))
fig_RHQ_dat<-data.frame(cbind(figRHQ_a, figRHQ_b))
	colnames(fig_RHQ_dat)=c("group","Rturns")
	fig_RHQ_dat$Rturns<-as.numeric(as.character(fig_RHQ_dat$Rturns))   
	range(fig_RHQ_dat$Rturns)
mu_RHQ<-ddply(fig_RHQ_dat, "group", summarise, grp.mean=mean(Rturns))

table(fig_RHQ_dat[fig_RHQ_dat$group=="control",]$Rturns)
table(fig_RHQ_dat[fig_RHQ_dat$group=="high CO2",]$Rturns)
table(fig_RHQ_dat$Rturns)
fig_RHQ<-ggplot(fig_RHQ_dat, aes(x=Rturns, fill=group)) + geom_bar(position = position_dodge2(preserve = "single"), show.legend=F) +
	scale_x_continuous(limits = c(-0.5,10.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10), labels=c("0","1","2","3","4","5","6","7","8","9","10"), name="Right turns") +
	scale_y_continuous(limits = c(0,16), breaks=c(0,5,10,15), labels=c("0","5","10","15"), name="Number of individuals") +
	geom_vline(data=mu_RHQ, aes(xintercept=grp.mean, color=group), linetype="dashed", size=1, show.legend = F) +
	annotate("text", x = 8, y = 14, label = "RHQ A. polyacanthus") +
	scale_color_grey() + scale_fill_grey() + theme_bw()
#postscript("Fig_RHQ.eps", width = 2.5, height = 2.5)
fig_RHQ + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + theme(axis.text.y=element_text(size=10)) + theme(axis.text.x=element_text(size=10), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=10)) + theme(legend.position = c(0.15, 0.9), legend.title=element_blank())
#dev.off()

#Control fish (statistics)
###############

#chi-square test for variance (control fish)
 chi_sq_RHQc=((length(XRHQ_RHQc)-1)*var(XRHQ_RHQc)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_RHQc<-pchisq(chi_sq_RHQc,df=(length(XRHQ_RHQc)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_RHQc=(mean(XRHQ_RHQc/n)-P)/(.5*.5/(length(XRHQ_RHQc)*n))^.5			#compute a Z score for the mean
 pZ_RHQc=pnorm(abs(Z_RHQc),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (control fish)
 ga_RHQc<-glmer(X_RHQc~1+(1|ind_RHQc),data=Control_RHQ,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_RHQc<-glm(X_RHQc~1,data=Control_RHQ,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_RHQc=summary(ga_RHQc)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_RHQc=anova(ga_RHQc,gb_RHQc)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_RHQc<-round(matrix(c(pr_m_RHQc,pr_v_RHQc,pZ_RHQc,sig_v_RHQc),ncol=4,nrow=1),3); colnames(out_RHQc)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_RHQc

#Treatment (CO2) fish (statistics)
###############

#chi-square test for variance (CO2 fish)
 chi_sq_RHQt=((length(XRHQ_RHQt)-1)*var(XRHQ_RHQt)/(n*.5*.5)) 			#compute a chi-square value: numerator is the observed variance, denominator is the expected variance assuming a normal approximation
 sig_v_RHQt<-pchisq(chi_sq_RHQt,df=(length(XRHQ_RHQt)-1),lower.tail=F)	#compute a p-value using the chi-square test with degrees of freedom = N-1 
#z-test for mean (control fish)
 Z_RHQt=(mean(XRHQ_RHQt/n)-P)/(.5*.5/(length(XRHQ_RHQt)*n))^.5			#compute a Z score for the mean
 pZ_RHQt=pnorm(abs(Z_RHQt),lower.tail=F)*2 					#2 tailed z-test of significance for the mean 
#GLMM test for mean and variance (CO2 fish)
 ga_RHQt<-glmer(X_RHQt~1+(1|ind_RHQt),data=CO2_RHQ,family="binomial")	#specify a generalized (binomial) linear random-effects model using the data 'dat', with intercept = grand mean of the data
 gb_RHQt<-glm(X_RHQt~1,data=CO2_RHQ,family="binomial")			#specify a generalized (binomial) linear model using the data 'dat', with intercept = grand mean of the data
 pr_m_RHQt=summary(ga_RHQt)$coefficients[4] 					#extract p-value for model intercept (i.e. test for mean different from 0.5)
 pr_v_RHQt=anova(ga_RHQt,gb_RHQt)[2,8] 						#likelihood ratio test to examine significance of random effects (not very powerful)
#table of results
 out_RHQt<-matrix(c(pr_m_RHQt,pr_v_RHQt,pZ_RHQt,sig_v_RHQt),ncol=4,nrow=1); colnames(out_RHQt)=c("glmer_mean","glmer_var","Z_mean","Chi_var")
 out_RHQt


