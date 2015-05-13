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
#s1.WIT = read.csv("C:/Users/Etienne/Dropbox/independent replication - Correll (2008, JPSP)/replication #1/data/wit_agg/detailed/wit_agg_detailed_s1.csv")
primeStimNum = substring(s1.WIT$Stim.2,3,3) #create new variables needed for later
targetStimNum = substring(s1.WIT$Stim.3, str_length(s1.WIT$Stim.3)-4, str_length(s1.WIT$Stim.3)-4)
RT.log = log(s1.WIT$RT)
WIT.cond = recode(s1.WIT$Cond, "c(1,2)='1';c(3,4)='2';c(5,6)='3'") #1=control, 2=avoid race, 3=use race (between-subjects conditions)
s1.WIT = cbind(s1.WIT, primeStimNum, targetStimNum, RT.log,WIT.cond) #add new variables to the full data set
s1.WIT = s1.WIT[order(s1.WIT$Subj,s1.WIT$Trial),]
s1.WIT.guns = subset(s1.WIT,(Trial>29 & tool=="gun")) #all gun trials, all instruction conditions
s1.WIT.guns.correct = subset(s1.WIT.guns,(Correct=="True")) #correct guns trials only (RT analyses), all conditions

#####################################################################################################################################
#Read in replication Sample #2 data
#####################################################################################################################################
s2.WIT = read.csv("https://raw.githubusercontent.com/eplebel/intra-individual-MS/master/data/wit_agg_detailed_s2.csv")
primeStimNum = substring(s2.WIT$Stim.2,3,3) #create new variables needed for later
targetStimNum = substring(s2.WIT$Stim.3, str_length(s2.WIT$Stim.3)-4, str_length(s2.WIT$Stim.3)-4)
RT.log = log(s2.WIT$RT)
WIT.cond = recode(s2.WIT$Cond, "c(1,2)='1';c(3,4)='2';c(5,6)='3'") #1=control, 2=avoid race, 3=use race (between-subjects conditions)
s2.WIT = cbind(s2.WIT, primeStimNum, targetStimNum, RT.log, WIT.cond) #add new variables to the full data set
s2.WIT = s2.WIT[order(s2.WIT$Subj,s2.WIT$Trial),]
s2.WIT.guns = subset(s2.WIT,(Trial>29 & tool=="gun")) #all gun trials, all instruction conditions
s2.WIT.guns.correct = subset(s2.WIT.guns,(Correct=="True")) #correct guns trials only (RT analyses), all conditions

######################################################################################################################################
#RT analysis (only correct trials considered following Payne (2001) & Correll (2008), and general convention)
#RT racial bias on GUN trials (i.e., faster to categorize Black-guns compared to White-guns)
#######################
# Replication SAMPLE #1      
#######################
s1.RT.model = lmer(RT.log ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s1.WIT.guns.correct, REML=FALSE) 
N=length(coef(s1.RT.model)$Subj[,2])
plot(fitted(s1.RT.model),residuals(s1.RT.model)) #checking model assumptions; things seem OK?
hist(residuals(s1.RT.model)) 
qqnorm(residuals(s1.RT.model))
#plotting histogram of Person-specific race-on-RT slopes
hist(coef(s1.RT.model)$Subj[,2]) 
frequencies(coef(s1.RT.model)$Subj[,2]) #As can be seen, heterogeneity in individual-specific slopes!
#Calculating stat sign. of individual slopes via reshuffling technique (Baron, 2010) 
set.seed(4) #make this reproducible
myiter=100 #increase this, though already very slow (i.e., ~3 mins) for 100 iterations
s1.slope.array<-matrix(,N,myiter) #to store iteration specific null slopes
for (i in 1:myiter){
  if (i==i) print(i) #just to show progress 
  s1.WIT.guns.correct$race.reshuffled = recode(rbinom(nrow(s1.WIT.guns.correct),1,.5), "c(0)='white'; else='black'") #reshuffle each Ps trials so that White/Black trials are completely randomized 
  s1.RT.model.reshuffled = lmer(RT.log ~ race.reshuffled + (1+race.reshuffled|Subj) + (1+race.reshuffled|targetStimNum), data=s1.WIT.guns.correct, REML=FALSE) #re-run analysis w/ "race.reshuffled" instead of "race"
  s1.slope.array[,i] = coef(s1.RT.model.reshuffled)$Subj[,2] #store participant-specific null slopes for ith iteration
}
#now calculate Subj-specific 95% null cutoffs for Subj-specific race slopes *AND* determine E%/N%/O% clusters
s1.slope.array.CIs<-data.frame(matrix(,N,7)) #had added these at the end of the slope.array matrix, but wasn't working somehow...
colnames(s1.slope.array.CIs)[1]="Subj" #Subject ID#
colnames(s1.slope.array.CIs)[2]="WIT.cond" #WIT instruction condition (1=control, 2=avoid, 3=use)
colnames(s1.slope.array.CIs)[3]="Slope.Estimate" #having this here will make it easier to plot later
colnames(s1.slope.array.CIs)[4]="null.LB" #null lower bound
colnames(s1.slope.array.CIs)[5]="null.UB" #null upper bound
colnames(s1.slope.array.CIs)[6]="CI.LB" #95% C.I. lower bound
colnames(s1.slope.array.CIs)[7]="CI.UB" #95% C.I. upper bound
s1.subj.cond.codes<-cast(melt(s1.WIT.guns.correct,id=c("Subj","WIT.cond"),c("Trial")),Subj+WIT.cond~variable,mean) #convoluted workaround to get the Subj IDs and condition #s
s1.slope.array.CIs$Subj=s1.subj.cond.codes$Subj
s1.slope.array.CIs$WIT.cond=s1.subj.cond.codes$WIT.cond
s1.slope.array.CIs$Slope.Estimate=coef(s1.RT.model)$Subj[,2]
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  subject.null.cutoffs=quantile(s1.slope.array[i,], c(.025, .975)) #get the 95% null cut-offs
  s1.slope.array.CIs$null.LB[i] = subject.null.cutoffs[1]
  s1.slope.array.CIs$null.UB[i] = subject.null.cutoffs[2]
  if (s1.slope.array.CIs$Slope.Estimate[i]<s1.slope.array.CIs$null.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s1.slope.array.CIs$Slope.Estimate[i]>s1.slope.array.CIs$null.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
#Full Sample 1 Clusters
cat(paste("Full Sample 1 Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Anti-White Slopes):  ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 53% (78)
#N% (Non-signif Slopes): 41% (60)
#O% (Anti-White Slopes):  7% (10)

#Calculating 95% C.I.s of individual-slopes via bootstrapping
myiter=100 #increase this to 1000
s1.slope.arrayCIs<-matrix(,N,myiter) #to store iteration specific boostrapped slopes
RT.log.resamp=rep(0,nrow(s1.WIT.guns.correct)) 
corr=as.numeric(s1.WIT.guns$Correct) 
corr=recode(corr, "c(1)='0';c(2)='1'")
s1.WIT.guns = cbind(s1.WIT.guns,corr)
num.correct.array = aggregate(s1.WIT.guns$corr, by = list(subj = s1.WIT.guns$Subj), FUN=mean) #num correct trials for each subj; 
num.corr=rep(0,nrow(s1.WIT.guns.correct)) 
s1.WIT.guns.correct = cbind(s1.WIT.guns.correct,num.corr,RT.log.resamp)
#add number correct as an additional column in the main data frame (repeated but that's fine)
for (i in 1:nrow(s1.WIT.guns.correct)){
  s1.WIT.guns.correct[i,c("num.corr")] = num.correct.array[s1.WIT.guns.correct[i,c("Subj")],c("x")]*100  #slow, but it works!!
}
bootstrap.loop.vector <- array(0, dim=nrow(num.correct.array))
bootstrap.loop.vector[1] = 1
for (s in 2:nrow(num.correct.array)){
  bootstrap.loop.vector[s] = bootstrap.loop.vector[s-1] + num.correct.array[num.correct.array$subj==4,c("x")]*100 #OK now sorted by Subj #s! 
}
#bootstrap.loop.vector = c(1, 101, 201,...) #construct vector required to step through next for-loop (via another loop!)
for (i in 1:myiter){
  for (j in bootstrap.loop.vector) #this will step through bootstrap.loop.vector (e.g., 1, 101, 201)
    { if (j==j) print(j) #just to show progress 
      #resample w/ replacement at subject level (# of correct trials)  #next line not correct..>@#$@#$@#$ too complicated
      s1.WIT.guns.correct$RT.log.resamp[j:j+num.correct[j]-1] = sample(s1.WIT.guns.correct$RT.log[j:j+num.correct[j]-1],num.correct[j],replace=TRUE)   
  
      [1:100]  #ok i think that works! [though still need to test it]
      [101:200]
      [201:300]
    }
  s1.RT.model.bootstrapped = lmer(RT.log.resamp ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s1.WIT.guns.correct, REML=FALSE) #re-run analysis
  s1.slope.arrayCIs[,i] = coef(s1.RT.model.bootstrapped)$Subj[,2] #store participant-specific bootstrapped slopes for ith iteration
}

#plotting subject-specific racial bias slopes
library(ggplot2) 
#anti-Black bias (faster on Black-gun compared to White-gun trials; thus, positive slope (increase in RT from black to white))
subj=100
subjdata = subset(s1.WIT.guns.correct,(Subj==subj)) #correct gun trials 
coefs <- coef(lm(RT ~ race, data = subjdata)) 
p<-qplot(race, RT, data=subjdata, main=paste("Subject ", subj, sep=""), xlab="", ylab="RTs on GUN trials (correct responses)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = coefs[1], slope=coefs[2], colour="red",size=1)  + ylim(0,1300)

#no bias
subj=6
subjdata = subset(s1.WIT.guns.correct,(Subj==subj)) #correct gun trials 
coefs <- coef(lm(RT ~ race, data = subjdata)) 
p<-qplot(race, RT, data=subjdata, main=paste("Subject ", subj, sep=""), xlab="", ylab="RTs on GUN trials (correct responses)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = coefs[1], slope=coefs[2], colour="red",size=1)  + ylim(0,1300)

#anti-White bias (slower on Black-gun compared to White-gun trials; thus, negative slope (decrease in RT from black to white))
subj=29
subjdata = subset(s1.WIT.guns.correct,(Subj==subj)) #correct gun trials 
coefs <- coef(lm(RT ~ race, data = subjdata)) 
p<-qplot(race, RT, data=subjdata, main=paste("Subject ", subj, sep=""), xlab="", ylab="RTs on GUN trials (correct responses)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = coefs[1], slope=coefs[2], colour="red",size=1)  + ylim(0,1300)






#######################
# Replication SAMPLE #2      
#######################
s2.RT.model = lmer(RT.log ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s2.WIT.guns.correct, REML=FALSE) #w/out instruction condition
N=length(coef(s2.RT.model)$Subj[,2])
plot(fitted(s2.RT.model),residuals(s2.RT.model)) #checking model assumptions; 
hist(residuals(s2.RT.model)) 
qqnorm(residuals(s2.RT.model))
#plotting histogram of Person-specific race-on-RT slopes
hist(coef(s2.RT.model)$Subj[,2]) 
frequencies(coef(s2.RT.model)$Subj[,2]) #As can be seen, heterogeneity in individual-specific slopes!
#Calculating stat sign. of individual slopes via reshuffling technique (Baron, 2010) 
set.seed(4) #make this reproducible
myiter=100 #increase this, though already very slow (i.e., ~3 mins) for 100 iterations
s2.slope.array<-matrix(,N,myiter) #to store iteration specific null slopes
for (i in 1:myiter){
  if (i==i) print(i) #just to show progress 
  s2.WIT.guns.correct$race.reshuffled = recode(rbinom(nrow(s2.WIT.guns.correct),1,.5), "c(0)='white'; else='black'") #reshuffle each Ps trials so that White/Black trials are completely randomized 
  s2.RT.model.reshuffled = lmer(RT.log ~ race.reshuffled + (1+race.reshuffled|Subj) + (1+race.reshuffled|targetStimNum), data=s2.WIT.guns.correct, REML=FALSE) #re-run analysis w/ "race.reshuffled" instead of "race"
  s2.slope.array[,i] = coef(s2.RT.model.reshuffled)$Subj[,2] #store participant-specific null slopes for ith iteration
}
#now calculate Subj-specific 95% null cutoffs for Subj-specific race slopes *AND* determine E%/N%/O% clusters
s2.slope.array.CIs<-data.frame(matrix(,N,7)) #had added these at the end of the slope.array matrix, but wasn't working somehow...
colnames(s2.slope.array.CIs)[1]="Subj" #Subject ID#
colnames(s2.slope.array.CIs)[2]="WIT.cond" #WIT instruction condition (1=control, 2=avoid, 3=use)
colnames(s2.slope.array.CIs)[3]="Slope.Estimate" #having this here will make it easier to plot later
colnames(s2.slope.array.CIs)[4]="null.LB" #null lower bound
colnames(s2.slope.array.CIs)[5]="null.UB" #null upper bound
colnames(s2.slope.array.CIs)[6]="CI.LB" #95% C.I. lower bound
colnames(s2.slope.array.CIs)[7]="CI.UB" #95% C.I. upper bound
s2.subj.cond.codes<-cast(melt(s2.WIT.guns.correct,id=c("Subj","WIT.cond"),c("Trial")),Subj+WIT.cond~variable,mean) #convoluted workaround to get the Subj IDs and condition #s
s2.slope.array.CIs$Subj=s2.subj.cond.codes$Subj
s2.slope.array.CIs$WIT.cond=s2.subj.cond.codes$WIT.cond
s2.slope.array.CIs$Slope.Estimate=coef(s2.RT.model)$Subj[,2]
sig.neg.slope.count=0
sig.pos.slope.count=0
for (i in 1:N) {
  subject.null.cutoffs=quantile(s2.slope.array[i,], c(.025, .975)) #get the 95% C.I. cut-offs
  s2.slope.array.CIs$null.LB[i] = subject.null.cutoffs[1]
  s2.slope.array.CIs$null.UB[i] = subject.null.cutoffs[2]
  if (s2.slope.array.CIs$Slope.Estimate[i]<s2.slope.array.CIs$null.LB[i]) sig.neg.slope.count=sig.neg.slope.count+1 #stat sign neg slope
  if (s2.slope.array.CIs$Slope.Estimate[i]>s2.slope.array.CIs$null.UB[i]) sig.pos.slope.count=sig.pos.slope.count+1 #stat sign pos slope
}
#Full Sample 2 Clusters
cat(paste("Full Sample 2 Clusters: \nE% (Anti-Black Slopes): ",round((sig.pos.slope.count/N)*100,2),"% (",sig.pos.slope.count,")\n","N% (Non-signif Slopes): ",round(((N-(sig.neg.slope.count+sig.pos.slope.count))/N)*100,2),"% (",N-(sig.neg.slope.count+sig.pos.slope.count),")\n","O% (Anti-White Slopes ): ",round((sig.neg.slope.count/N)*100,2),"% (",sig.neg.slope.count,")\n",sep=""))
#E% (Anti-Black Slopes): 62% (92)
#N% (Non-signif Slopes): 32% (48)
#O% (Anti-White Slopes ): 6% (9)







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

require(lattice)
s1.rawRT.model = lmer(RT ~ race + (1+race|Subj) + (1+race|targetStimNum), data=s1.WIT.guns.correct, REML=FALSE) #w/out instruction condition
tess <- ranef(s1.rawRT.model, condVar = TRUE)
dotplot(tess) #kind of works, but not sure how CIs are calculated...

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
