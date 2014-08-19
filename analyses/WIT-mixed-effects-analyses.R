require(lme4)
require(stringr) #required to extract stimuli # for both prime and target stimuli strings
require(Deducer) #required for recode function somehow
require(reshape) #required for some convoluted workaround to convert from LONG to WIDE format
#Mixed-effects models for my Weapon Identification Task (WIT) data w/ goal of identifying clusters of 
# (1) Ps showing anti-Black bias vs. (2) not showing vs. (3) showing pro-White bias
#####################################################################################################################################
#Read in replication Sample #1 data
#####################################################################################################################################
s1.WIT = read.csv("https://raw.githubusercontent.com/eplebel/intra-individual-MS/master/data/wit_agg_detailed_s1.csv")
primeStimNum = substring(s1.WIT$Stim.2,3,3) #create new variables needed for later
targetStimNum = substring(s1.WIT$Stim.3, str_length(s1.WIT$Stim.3)-4, str_length(s1.WIT$Stim.3)-4)
RT.log = log(s1.WIT$RT)
WIT.cond = recode(s1.WIT$Cond, "c(1,2)='1';c(3,4)='2';c(5,6)='3'") #1=control, 2=avoid race, 3=use race (between-subjects conditions)
s1.WIT = cbind(s1.WIT, primeStimNum, targetStimNum, RT.log,WIT.cond) #add new variables to the full data set
s1.WIT.tools = subset(s1.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))) #all tool trials, all instruction conditions
s1.WIT.tools.correct = subset(s1.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")) #correct tool trials only (RT analyses), all conditions
s1.WIT.tools.correct.CONTROL = subset(s1.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(WIT.cond==1)) #correct tool trials, CONTROL CONDITION
s1.WIT.tools.correct.AVOID = subset(s1.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(WIT.cond==2)) #correct tool trials, AVOID RACE CONDITION
s1.WIT.tools.correct.USE = subset(s1.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(WIT.cond==3)) #correct tool trials, USE RACE CONDITION
#####################################################################################################################################
#Read in replication Sample #2 data
#####################################################################################################################################
s2.WIT = read.csv("https://raw.githubusercontent.com/eplebel/intra-individual-MS/master/data/wit_agg_detailed_s2.csv")
primeStimNum = substring(s2.WIT$Stim.2,3,3) #create new variables needed for later
targetStimNum = substring(s2.WIT$Stim.3, str_length(s2.WIT$Stim.3)-4, str_length(s2.WIT$Stim.3)-4)
RT.log = log(s2.WIT$RT)
WIT.cond = recode(s2.WIT$Cond, "c(1,2)='1';c(3,4)='2';c(5,6)='3'") #1=control, 2=avoid race, 3=use race (between-subjects conditions)
s2.WIT = cbind(s2.WIT, primeStimNum, targetStimNum, RT.log, WIT.cond) #add new variables to the full data set
s2.WIT.tools = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))) #all tool trials, all instruction conditions
s2.WIT.tools.correct = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")) #correct tool trials only (RT analyses), all conditions
s2.WIT.tools.correct.CONTROL = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(WIT.cond==1)) #correct tool trials, CONTROL CONDITION
s2.WIT.tools.correct.AVOID = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(WIT.cond==2)) #correct tool trials, AVOID RACE CONDITION
s2.WIT.tools.correct.USE = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(WIT.cond==3)) #correct tool trials, USE RACE CONDITION

######################################################################################################################################
#RT analysis (only correct trials considered following Payne (2001) & Correll (2008), and general convention)
#Only considering RT racial bias on tool trials as starting point (i.e., slower to categorize Black-tools compared to White-tools)
#######################
# Replication SAMPLE #1      
#######################
s1.RT.model = lmer(RT.log ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s1.WIT.tools.correct, REML=FALSE) #w/out instruction condition
N=length(coef(s1.RT.model)$Subj[,2])
plot(fitted(s1.RT.model),residuals(s1.RT.model)) #checking model assumptions; things seem OK?
hist(residuals(s1.RT.model)) 
qqnorm(residuals(s1.RT.model))
#adding between-subjects instruction manipulation
s1.RT.model.instrCOND = lmer(RT.log ~ race*WIT.cond + (1+race|Subj) + (1+race|targetStimNum), data=s1.WIT.tools.correct, REML=FALSE)
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
slope.array<-matrix(,N,myiter) #to store iteration specific null slopes
for (i in 1:myiter){
  if (i==i) print(i) #just to show progress 
  s1.WIT.tools.correct$race.reshuffled = recode(rbinom(nrow(s1.WIT.tools.correct),1,.5), "c(0)='white'; else='black'") #reshuffle each Ps trials so that White/Black trials are completely randomized by shifting down vector (more efficient than invoking rbinom each iteration)
  s1.RT.model.reshuffled = lmer(RT.log ~ race.reshuffled + (1+race.reshuffled|Subj) + (1+race.reshuffled|targetStimNum), data=s1.WIT.tools.correct, REML=FALSE) #re-run analysis w/ "race.reshuffled" instead of "race"
  slope.array[,i] = coef(s1.RT.model.reshuffled)$Subj[,2] #store participant-specific null slopes for ith iteration
}
#now calculate Subj-specific 95% CIs for Subj-specific race slopes *AND* determine E%/N%/O% clusters
slope.array.CIs<-data.frame(matrix(,N,5)) #had added these at the end of the slope.array matrix, but wasn't working somehow...
colnames(slope.array.CIs)[1]="Subj" #Subject ID#
colnames(slope.array.CIs)[2]="WIT.cond" #WIT instruction condition (1=control, 2=avoid, 3=use)
colnames(slope.array.CIs)[3]="Slope.Estimate" #having this here will make it easier to plot later
colnames(slope.array.CIs)[4]="CI.LB" #set lower bound variable name
colnames(slope.array.CIs)[5]="CI.UB" #set upper bound variable name
subj.cond.codes<-cast(melt(s1.WIT.tools.correct,id=c("Subj","WIT.cond"),c("Trial")),Subj+WIT.cond~variable,mean) #convoluted workaround to get the Subj IDs and condition #s
slope.array.CIs$Subj=subj.cond.codes$Subj
slope.array.CIs$WIT.cond=subj.cond.codes$WIT.cond
slope.array.CIs$Slope.Estimate=coef(s1.RT.model)$Subj[,2]
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  subject.CIs=quantile(slope.array[i,], c(.025, .975)) #get the 95% C.I. cut-offs
  slope.array.CIs$CI.LB[i] = subject.CIs[1]
  slope.array.CIs$CI.UB[i] = subject.CIs[2]
  if (slope.array.CIs$Slope.Estimate[i]<slope.array.CIs$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs$Slope.Estimate[i]>slope.array.CIs$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
#Full Sample 1 Clusters
cat(paste("Full Sample 1 Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes):  ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 78% (116)
#N% (Non-signif Slopes): 22% (32)
#O% (Pro-White Slopes):  0% (0)
###########################################################
#CLUSTERS in different instruction conditions
###########################################################
#Sample 1 CONTROL Condition Clusters
slope.array.CIs.CONTROL = subset(slope.array.CIs,slope.array.CIs$WIT.cond==1)
N=nrow(slope.array.CIs.CONTROL)
sig.neg.slope.count=0
sig.pos.slope.count=0
i=1
for (i in 1:N) {
  if (slope.array.CIs.CONTROL$Slope.Estimate[i]<slope.array.CIs.CONTROL$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs.CONTROL$Slope.Estimate[i]>slope.array.CIs.CONTROL$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 1 CONTROL Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 75% (36)
#N% (Non-signif Slopes): 25% (12)
#O% (Pro-White Slopes ): 0% (0)
#Sample 1 AVOID Condition Clusters
slope.array.CIs.AVOID = subset(slope.array.CIs,slope.array.CIs$WIT.cond==2)
N=nrow(slope.array.CIs.AVOID)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (slope.array.CIs.AVOID$Slope.Estimate[i]<slope.array.CIs.AVOID$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs.AVOID$Slope.Estimate[i]>slope.array.CIs.AVOID$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 1 AVOID Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 86% (43)
#N% (Non-signif Slopes): 14% (7)
#O% (Pro-White Slopes ): 0% (0)
#Sample 1 USE Condition Clusters
slope.array.CIs.USE = subset(slope.array.CIs,slope.array.CIs$WIT.cond==3)
N=nrow(slope.array.CIs.USE)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (slope.array.CIs.USE$Slope.Estimate[i]<slope.array.CIs.USE$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs.USE$Slope.Estimate[i]>slope.array.CIs.USE$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 1 USE Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 74% (37)
#N% (Non-signif Slopes): 26% (13)
#O% (Pro-White Slopes ): 0% (0)
###########################################################

#######################
# Replication SAMPLE #2      
#######################
s2.RT.model = lmer(RT.log ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s2.WIT.tools.correct, REML=FALSE) #w/out instruction condition
N=length(coef(s2.RT.model)$Subj[,2])
plot(fitted(s2.RT.model),residuals(s2.RT.model)) #checking model assumptions; 
hist(residuals(s2.RT.model)) 
qqnorm(residuals(s2.RT.model))
#adding between-subjects instruction manipulation
s2.RT.model.instrCOND = lmer(RT.log ~ race*WIT.cond + (1+race|Subj) + (1+race|targetStimNum), data=s2.WIT.tools.correct, REML=FALSE)
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
slope.array<-matrix(,N,myiter) #to store iteration specific null Slopes
for (i in 1:myiter){
  if (i==i) print(i) #just to show progress 
  s2.WIT.tools.correct$race.reshuffled = recode(rbinom(nrow(s2.WIT.tools.correct),1,.5), "c(0)='white'; else='black'") #reshuffle each Ps trials so that White/Black trials are completely randomized by shifting down vector (more efficient than invoking rbinom each iteration)
  s2.RT.model.reshuffled = lmer(RT.log ~ race.reshuffled + (1+race.reshuffled|Subj) + (1+race.reshuffled|targetStimNum), data=s2.WIT.tools.correct, REML=FALSE) #re-run analysis w/ "race.reshuffled" instead of "race"
  slope.array[,i] = coef(s2.RT.model.reshuffled)$Subj[,2] #store participant-specific null slopes for ith iteration
}
#now calculate Subj-specific 95% CIs for Subj-specific race slopes *AND* determine E%/N%/O% clusters
slope.array.CIs<-data.frame(matrix(,N,5)) #had added these at the end of the slope.array matrix, but wasn't working somehow...
colnames(slope.array.CIs)[1]="Subj" #Subject ID#
colnames(slope.array.CIs)[2]="WIT.cond" #WIT instruction condition (1=control, 2=avoid, 3=use)
colnames(slope.array.CIs)[3]="Slope.Estimate" #having this here will make it easier to plot later
colnames(slope.array.CIs)[4]="CI.LB" #set lower bound variable name
colnames(slope.array.CIs)[5]="CI.UB" #set upper bound variable name
subj.cond.codes<-cast(melt(s2.WIT.tools.correct,id=c("Subj","WIT.cond"),c("Trial")),Subj+WIT.cond~variable,mean) #convoluted workaround to get the Subj IDs and condition #s
slope.array.CIs$Subj=subj.cond.codes$Subj
slope.array.CIs$WIT.cond=subj.cond.codes$WIT.cond
slope.array.CIs$Slope.Estimate=coef(s2.RT.model)$Subj[,2]
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  subject.CIs=quantile(slope.array[i,], c(.025, .975)) #get the 95% C.I. cut-offs
  slope.array.CIs$CI.LB[i] = subject.CIs[1]
  slope.array.CIs$CI.UB[i] = subject.CIs[2]
  if (slope.array.CIs$Slope.Estimate[i]<slope.array.CIs$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs$Slope.Estimate[i]>slope.array.CIs$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
#Full Sample 2 Clusters
cat(paste("Full Sample 2 Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 74% (110)
#N% (Non-signif Slopes): 25% (37)
#O% (Pro-White Slopes ): 1% (2)
###########################################################
#CLUSTERS in different instruction conditions
###########################################################
#Sample 2 CONTROL Condition Clusters
slope.array.CIs.CONTROL = subset(slope.array.CIs,slope.array.CIs$WIT.cond==1)
N=nrow(slope.array.CIs.CONTROL)
sig.neg.slope.count=0
sig.pos.slope.count=0
i=1
for (i in 1:N) {
  if (slope.array.CIs.CONTROL$Slope.Estimate[i]<slope.array.CIs.CONTROL$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs.CONTROL$Slope.Estimate[i]>slope.array.CIs.CONTROL$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 2 CONTROL Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 74% (37)
#N% (Non-signif Slopes): 24% (12)
#O% (Pro-White Slopes ): 2% (1)
#Sample 2 AVOID Condition Clusters
slope.array.CIs.AVOID = subset(slope.array.CIs,slope.array.CIs$WIT.cond==2)
N=nrow(slope.array.CIs.AVOID)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (slope.array.CIs.AVOID$Slope.Estimate[i]<slope.array.CIs.AVOID$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs.AVOID$Slope.Estimate[i]>slope.array.CIs.AVOID$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 2 AVOID Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 84% (42)
#N% (Non-signif Slopes): 14% (7)
#O% (Pro-White Slopes ): 2% (1)
#Sample 2 USE Condition Clusters
slope.array.CIs.USE = subset(slope.array.CIs,slope.array.CIs$WIT.cond==3)
N=nrow(slope.array.CIs.USE)
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  if (slope.array.CIs.USE$Slope.Estimate[i]<slope.array.CIs.USE$CI.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (slope.array.CIs.USE$Slope.Estimate[i]>slope.array.CIs.USE$CI.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
cat(paste("Sample 2 USE Condition Clusters: \nE% (Anti-Black Slopes): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Pro-White Slopes ): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 63% (31)
#N% (Non-signif Slopes): 37% (18)
#O% (Pro-White Slopes ): 0% (0)
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
RT.tools.list <- lmList(RT ~ race | Subj, data=s1.WIT.tools.correct)
plot(intervals(RT.tools.list), main="95% CIs")



#browse through specific subjects
subject=1 #actually case/row #  (not actual Subject ID #)
hist(slope.array[subject,],breaks=100)
paste(round(coef(correll.RT.model)$Subj[subject,2],4), ", 95% CI = [",round(slope.array.CIs[subject,4],4), ",", round(slope.array.CIs[subject,5],4), "]",sep="")
