require(lme4)
require(stringr) #required to extract stimuli # for both prime and target stimuli strings
require(Deducer) #required for recode function somehow
require(reshape) #required for some convoluted workaround to convert from LONG to WIDE format
#Mixed-effects models for my Weapon Identification Task (WIT) data w/ goal of identifying clusters of 
# (1) Ps showing anti-Black bias vs. (2) not showing vs. (3) showing pro-White bias
#####################################################################################################################################
#Read in replication Sample #1 data
#####################################################################################################################################
#s1.WIT = read.csv("https://raw.githubusercontent.com/eplebel/intra-individual-MS/master/data/wit_agg_detailed_s1.csv")
s1.WIT = read.csv("C:/Users/Etienne/Dropbox/independent replication - Correll (2008, JPSP)/replication #1/data/wit_agg/detailed/wit_agg_detailed_s1.csv")
primeStimNum = substring(s1.WIT$Stim.2,3,3) #create new variables needed for later
targetStimNum = substring(s1.WIT$Stim.3, str_length(s1.WIT$Stim.3)-4, str_length(s1.WIT$Stim.3)-4)
RT.log = log(s1.WIT$RT)
WIT.cond = recode(s1.WIT$Cond, "c(1,2)='1';c(3,4)='2';c(5,6)='3'") #1=control, 2=avoid race, 3=use race (between-subjects conditions)
s1.WIT = cbind(s1.WIT, primeStimNum, targetStimNum, RT.log,WIT.cond) #add new variables to the full data set
s1.WIT.guns = subset(s1.WIT,(Trial>29 & tool=="gun")) #all gun trials, all instruction conditions
s1.WIT.guns.correct = subset(s1.WIT.guns,(Correct=="True")) #correct guns trials only (RT analyses), all conditions
s1.WIT.guns.correct.CONTROL = subset(s1.WIT.guns.correct,(WIT.cond==1)) #correct guns trials, CONTROL CONDITION
s1.WIT.guns.correct.AVOID = subset(s1.WIT.guns.correct,(WIT.cond==2)) #correct guns trials, AVOID RACE CONDITION
s1.WIT.guns.correct.USE = subset(s1.WIT.guns.correct,(WIT.cond==3)) #correct guns trials, USE RACE CONDITION
#####################################################################################################################################
#Read in replication Sample #2 data
#####################################################################################################################################
s2.WIT = read.csv("https://raw.githubusercontent.com/eplebel/intra-individual-MS/master/data/wit_agg_detailed_s2.csv")
primeStimNum = substring(s2.WIT$Stim.2,3,3) #create new variables needed for later
targetStimNum = substring(s2.WIT$Stim.3, str_length(s2.WIT$Stim.3)-4, str_length(s2.WIT$Stim.3)-4)
RT.log = log(s2.WIT$RT)
WIT.cond = recode(s2.WIT$Cond, "c(1,2)='1';c(3,4)='2';c(5,6)='3'") #1=control, 2=avoid race, 3=use race (between-subjects conditions)
s2.WIT = cbind(s2.WIT, primeStimNum, targetStimNum, RT.log, WIT.cond) #add new variables to the full data set
s2.WIT.guns = subset(s2.WIT,(Trial>29 & tool=="gun")) #all gun trials, all instruction conditions
s2.WIT.guns.correct = subset(s2.WIT.guns,(Correct=="True")) #correct guns trials only (RT analyses), all conditions
s2.WIT.guns.correct.CONTROL = subset(s2.WIT.guns.correct,(WIT.cond==1)) #correct guns trials, CONTROL CONDITION
s2.WIT.guns.correct.AVOID = subset(s2.WIT.guns.correct,(WIT.cond==2)) #correct guns trials, AVOID RACE CONDITION
s2.WIT.guns.correct.USE = subset(s2.WIT.guns.correct,(WIT.cond==3)) #correct guns trials, USE RACE CONDITION

######################################################################################################################################
#RT analysis (only correct trials considered following Payne (2001) & Correll (2008), and general convention)
#RT racial bias on GUN trials (i.e., faster to categorize Black-guns compared to White-guns)
#######################
# Replication SAMPLE #1      
#######################
s1.RT.model = lmer(RT.log ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s1.WIT.guns.correct, REML=FALSE) #w/out instruction condition
N=length(coef(s1.RT.model)$Subj[,2])
plot(fitted(s1.RT.model),residuals(s1.RT.model)) #checking model assumptions; things seem OK?
hist(residuals(s1.RT.model)) 
qqnorm(residuals(s1.RT.model))
#adding between-subjects instruction manipulation
s1.RT.model.instrCOND = lmer(RT.log ~ race*WIT.cond + (1+race|Subj) + (1+race|targetStimNum), data=s1.WIT.guns.correct, REML=FALSE)
plot(fitted(s1.RT.model.instrCOND),residuals(s1.RT.model.instrCOND)) #checking model assumptions;
hist(residuals(s1.RT.model.instrCOND)) 
qqnorm(residuals(s1.RT.model.instrCOND))
anova(s1.RT.model,s1.RT.model.instrCOND) #going with simpler model, given adding instruction condition doesn't improve model fit
#plotting histogram of Person-specific race-on-RT slopes
hist(coef(s1.RT.model)$Subj[,2]) 
frequencies(coef(s1.RT.model)$Subj[,2]) #As can be seen, heterogeneity in individual-specific slopes!
#CIs around individual slopes; reshuffling technique inspired by Baron (2010) 
set.seed(4) #make this reproducible
myiter=100 #increase this, though already very slow (i.e., ~3 mins) for 100 iterations
s1.slope.array<-matrix(,N,myiter) #to store iteration specific null slopes
for (i in 1:myiter){
  if (i==i) print(i) #just to show progress 
  s1.WIT.guns.correct$race.reshuffled = recode(rbinom(nrow(s1.WIT.guns.correct),1,.5), "c(0)='white'; else='black'") #reshuffle each Ps trials so that White/Black trials are completely randomized 
  s1.RT.model.reshuffled = lmer(RT.log ~ race.reshuffled + (1+race.reshuffled|Subj) + (1+race.reshuffled|targetStimNum), data=s1.WIT.guns.correct, REML=FALSE) #re-run analysis w/ "race.reshuffled" instead of "race"
  s1.slope.array[,i] = coef(s1.RT.model.reshuffled)$Subj[,2] #store participant-specific null slopes for ith iteration
}
#now calculate Subj-specific 95% CIs for Subj-specific race slopes *AND* determine E%/N%/O% clusters
s1.slope.array.CIs<-data.frame(matrix(,N,5)) #had added these at the end of the slope.array matrix, but wasn't working somehow...
colnames(s1.slope.array.CIs)[1]="Subj" #Subject ID#
colnames(s1.slope.array.CIs)[2]="WIT.cond" #WIT instruction condition (1=control, 2=avoid, 3=use)
colnames(s1.slope.array.CIs)[3]="Slope.Estimate" #having this here will make it easier to plot later
colnames(s1.slope.array.CIs)[4]="CI.LB" #set lower bound variable name
colnames(s1.slope.array.CIs)[5]="CI.UB" #set upper bound variable name
s1.subj.cond.codes<-cast(melt(s1.WIT.guns.correct,id=c("Subj","WIT.cond"),c("Trial")),Subj+WIT.cond~variable,mean) #convoluted workaround to get the Subj IDs and condition #s
s1.slope.array.CIs$Subj=s1.subj.cond.codes$Subj
s1.slope.array.CIs$WIT.cond=s1.subj.cond.codes$WIT.cond
s1.slope.array.CIs$Slope.Estimate=coef(s1.RT.model)$Subj[,2]
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  subject.CIs=quantile(s1.slope.array[i,], c(.025, .975)) #get the 95% C.I. cut-offs
  s1.slope.array.CIs$CI.LB[i] = subject.CIs[1]
  s1.slope.array.CIs$CI.UB[i] = subject.CIs[2]
  if (s1.slope.array.CIs$Slope.Estimate[i]<s1.slope.array.CIs$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s1.slope.array.CIs$Slope.Estimate[i]>s1.slope.array.CIs$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
#Full Sample 1 Clusters
cat(paste("Full Sample 1 Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes):  ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 53% (78)
#N% (Non-signif Slopes): 41% (60)
#O% (Pro-White Slopes):  7% (10)
###########################################################
#CLUSTERS in different instruction conditions
###########################################################
#Sample 1 CONTROL Condition Clusters
s1.slope.array.CIs.CONTROL = subset(s1.slope.array.CIs,s1.slope.array.CIs$WIT.cond==1)
N=nrow(s1.slope.array.CIs.CONTROL)
sig.neg.slope.count=0
sig.pos.slope.count=0
i=1
for (i in 1:N) {
  if (s1.slope.array.CIs.CONTROL$Slope.Estimate[i]<s1.slope.array.CIs.CONTROL$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s1.slope.array.CIs.CONTROL$Slope.Estimate[i]>s1.slope.array.CIs.CONTROL$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 1 CONTROL Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 40% (19)
#N% (Non-signif Slopes): 50% (24)
#O% (Pro-White Slopes ): 10% (5)
#Sample 1 AVOID Condition Clusters
s1.slope.array.CIs.AVOID = subset(s1.slope.array.CIs,s1.slope.array.CIs$WIT.cond==2)
N=nrow(s1.slope.array.CIs.AVOID)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (s1.slope.array.CIs.AVOID$Slope.Estimate[i]<s1.slope.array.CIs.AVOID$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s1.slope.array.CIs.AVOID$Slope.Estimate[i]>s1.slope.array.CIs.AVOID$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 1 AVOID Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 62% (31)
#N% (Non-signif Slopes): 34% (17)
#O% (Pro-White Slopes ): 4% (2)
#Sample 1 USE Condition Clusters
s1.slope.array.CIs.USE = subset(s1.slope.array.CIs,s1.slope.array.CIs$WIT.cond==3)
N=nrow(s1.slope.array.CIs.USE)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (s1.slope.array.CIs.USE$Slope.Estimate[i]<s1.slope.array.CIs.USE$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s1.slope.array.CIs.USE$Slope.Estimate[i]>s1.slope.array.CIs.USE$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 1 USE Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 56% (28)
#N% (Non-signif Slopes): 38% (19)
#O% (Pro-White Slopes ): 6% (3)
###########################################################

#plotting subject-specific racial bias slopes
library(ggplot2) 
#anti-Black bias (faster on Black-gun compared to White-gun trials; thus, positive slope (increase in RT from black to white))
subj12data = subset(s1.WIT.guns.correct,(Subj==12)) #correct gun trials 
coef(lm(RT ~ race, data = subj12data))
p<-qplot(race, RT, data=subj12data, main="Subject 12", xlab="", ylab="RT (correct trials)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = 335.2, slope=146.8, colour="red",size=1)  + ylim(0,5000)
#pro-White bias (slower on Black-gun compared to White-gun trials; thus, negative slope (decrease in RT from black to white))
subj29data = subset(s1.WIT.guns.correct,(Subj==29)) #correct gun trials 
coef(lm(RT ~ race, data = subj29data))
p<-qplot(race, RT, data=subj29data, main="Subject 29", xlab="", ylab="RT (correct trials)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = 300.6, slope=-23.9, colour="red",size=1)  + ylim(0,1000)

#######################
# Replication SAMPLE #2      
#######################
s2.RT.model = lmer(RT.log ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s2.WIT.guns.correct, REML=FALSE) #w/out instruction condition
N=length(coef(s2.RT.model)$Subj[,2])
plot(fitted(s2.RT.model),residuals(s2.RT.model)) #checking model assumptions; 
hist(residuals(s2.RT.model)) 
qqnorm(residuals(s2.RT.model))
#adding between-subjects instruction manipulation
s2.RT.model.instrCOND = lmer(RT.log ~ race*WIT.cond + (1+race|Subj) + (1+race|targetStimNum), data=s2.WIT.guns.correct, REML=FALSE)
plot(fitted(s2.RT.model.instrCOND),residuals(s2.RT.model.instrCOND)) #checking model assumptions; 
hist(residuals(s2.RT.model.instrCOND)) 
qqnorm(residuals(s2.RT.model.instrCOND))
anova(s2.RT.model,s2.RT.model.instrCOND) #going with simpler model, given adding instruction condition doesn't improve model fit
#plotting histogram of Person-specific race-on-RT slopes
hist(coef(s2.RT.model)$Subj[,2]) 
frequencies(coef(s2.RT.model)$Subj[,2]) #As can be seen, heterogeneity in individual-specific slopes!
#CIs around individual slopes; reshuffling technique inspired by Baron (2010) 
set.seed(4) #make this reproducible
myiter=100 #increase this, though already very slow (i.e., ~3 mins) for 100 iterations
s2.slope.array<-matrix(,N,myiter) #to store iteration specific null Slopes
for (i in 1:myiter){
  if (i==i) print(i) #just to show progress 
  s2.WIT.guns.correct$race.reshuffled = recode(rbinom(nrow(s2.WIT.guns.correct),1,.5), "c(0)='white'; else='black'") #reshuffle each Ps trials so that White/Black trials are completely randomized 
  s2.RT.model.reshuffled = lmer(RT.log ~ race.reshuffled + (1+race.reshuffled|Subj) + (1+race.reshuffled|targetStimNum), data=s2.WIT.guns.correct, REML=FALSE) #re-run analysis w/ "race.reshuffled" instead of "race"
  s2.slope.array[,i] = coef(s2.RT.model.reshuffled)$Subj[,2] #store participant-specific null slopes for ith iteration
}
#now calculate Subj-specific 95% CIs for Subj-specific race slopes *AND* determine E%/N%/O% clusters
s2.slope.array.CIs<-data.frame(matrix(,N,5)) #had added these at the end of the slope.array matrix, but wasn't working somehow...
colnames(s2.slope.array.CIs)[1]="Subj" #Subject ID#
colnames(s2.slope.array.CIs)[2]="WIT.cond" #WIT instruction condition (1=control, 2=avoid, 3=use)
colnames(s2.slope.array.CIs)[3]="Slope.Estimate" #having this here will make it easier to plot later
colnames(s2.slope.array.CIs)[4]="CI.LB" #set lower bound variable name
colnames(s2.slope.array.CIs)[5]="CI.UB" #set upper bound variable name
s2.subj.cond.codes<-cast(melt(s2.WIT.guns.correct,id=c("Subj","WIT.cond"),c("Trial")),Subj+WIT.cond~variable,mean) #convoluted workaround to get the Subj IDs and condition #s
s2.slope.array.CIs$Subj=s2.subj.cond.codes$Subj
s2.slope.array.CIs$WIT.cond=s2.subj.cond.codes$WIT.cond
s2.slope.array.CIs$Slope.Estimate=coef(s2.RT.model)$Subj[,2]
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  subject.CIs=quantile(s2.slope.array[i,], c(.025, .975)) #get the 95% C.I. cut-offs
  s2.slope.array.CIs$CI.LB[i] = subject.CIs[1]
  s2.slope.array.CIs$CI.UB[i] = subject.CIs[2]
  if (s2.slope.array.CIs$Slope.Estimate[i]<s2.slope.array.CIs$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s2.slope.array.CIs$Slope.Estimate[i]>s2.slope.array.CIs$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
#Full Sample 2 Clusters
cat(paste("Full Sample 2 Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 62% (92)
#N% (Non-signif Slopes): 32% (48)
#O% (Pro-White Slopes ): 6% (9)
###########################################################
#CLUSTERS in different instruction conditions
###########################################################
#Sample 2 CONTROL Condition Clusters
s2.slope.array.CIs.CONTROL = subset(s2.slope.array.CIs,s2.slope.array.CIs$WIT.cond==1)
N=nrow(s2.slope.array.CIs.CONTROL)
sig.neg.slope.count=0
sig.pos.slope.count=0
i=1
for (i in 1:N) {
  if (s2.slope.array.CIs.CONTROL$Slope.Estimate[i]<s2.slope.array.CIs.CONTROL$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s2.slope.array.CIs.CONTROL$Slope.Estimate[i]>s2.slope.array.CIs.CONTROL$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 2 CONTROL Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 54% (27)
#N% (Non-signif Slopes): 40% (20)
#O% (Pro-White Slopes ): 6% (3)
#Sample 2 AVOID Condition Clusters
s2.slope.array.CIs.AVOID = subset(s2.slope.array.CIs,s2.slope.array.CIs$WIT.cond==2)
N=nrow(s2.slope.array.CIs.AVOID)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (s2.slope.array.CIs.AVOID$Slope.Estimate[i]<s2.slope.array.CIs.AVOID$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s2.slope.array.CIs.AVOID$Slope.Estimate[i]>s2.slope.array.CIs.AVOID$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 2 AVOID Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 62% (31)
#N% (Non-signif Slopes): 34% (17)
#O% (Pro-White Slopes ): 4% (2)
#Sample 2 USE Condition Clusters
s2.slope.array.CIs.USE = subset(s2.slope.array.CIs,s2.slope.array.CIs$WIT.cond==3)
N=nrow(s2.slope.array.CIs.USE)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (s2.slope.array.CIs.USE$Slope.Estimate[i]<s2.slope.array.CIs.USE$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s2.slope.array.CIs.USE$Slope.Estimate[i]>s2.slope.array.CIs.USE$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 2 USE Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 69% (34)
#N% (Non-signif Slopes): 22% (11)
#O% (Pro-White Slopes ): 8% (4)
###########################################################






#NEXT STEPS:
#1-form 95% CIs around individual-level slope estimates (because current "95% CIs" are actually null sampling distributions *not* CIs!!)
#2-figure out issue w/ crossed random effects vs. non-crossed random effects... and also whether to also model prime stimuli as random effects???)
#3-resolve the singular & false convergence warning messages
#4-figure out better way to analyze RT data via LMER *without* log-transforming data??





#new approach to determining & plotting 95% CIs around slopes from Fox (2002)
#this works, but it doesn't allow me to include the full mixed-effects model equations, so values are different!
require(lattice)
require(nlme)
RT.tools.list <- lmList(RT ~ race | Subj, data=s1.WIT.guns.correct)
plot(intervals(RT.tools.list), main="95% CIs")

#playing around with the lmerTest package, as suggested by Fred, 
#for individual-level slope significance but doesn't seem very promising
require(lmerTest)
summary(s1.RT.model)
anova(s1.RT.model,ddf="Kenward-Roger")
st<-step(s1.RT.model)
plot(st)

#browse through specific subjects
#subject=1 #actually case/row #  (not actual Subject ID #)
#hist(slope.array[subject,],breaks=100)
#paste(round(coef(correll.RT.model)$Subj[subject,2],4), ", 95% CI = [",round(slope.array.CIs[subject,4],4), ",", round(slope.array.CIs[subject,5],4), "]",sep="")
