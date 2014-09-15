require(lme4)
require(plyr)
require(reshape2)
require(lmerTest)
require(languageR)
require(lattice)

# READ DATA FROM GITHUB ---------------------------------------------------
library(RCurl)

# Copy RAW url from GitHub repository
s1.WIT <- read.csv(textConnection(getURL("https://raw.githubusercontent.com/eplebel/intra-individual-MS/master/data/wit_agg_detailed_s1.csv")))            
s2.WIT <- read.csv(textConnection(getURL("https://raw.githubusercontent.com/eplebel/intra-individual-MS/master/data/wit_agg_detailed_s2.csv")))

# Get the experimental trials, order them by Subject and Timestamp
s1.TOT <- subset(s1.WIT,subset = s1.WIT$Block==5)
s1.TOT <- s1.TOT[ order(s1.TOT['Subj'],s1.TOT['Started.']), ]

# Log transform RT and make False responses NA
s1.TOT$RT.log   <- log(s1.TOT$RT)
s1.TOT$RT.log.c <- s1.TOT$RT.log 
s1.TOT$RT.log.c[s1.TOT$Correct=="False"] <- NA

# Create and re-level some factors
# Tool and Race
s1.TOT$tool <- relevel(s1.TOT$tool, ref="tool")
s1.TOT$race <- relevel(s1.TOT$race, ref="white")
# Factor Subj
s1.TOT$Subj <- factor(s1.TOT$Subj)
# Factor Condition (Between subjects)
s1.TOT$Condition <- factor(s1.TOT$File,labels=c("Avoid Race","Control","Use Race"))
s1.TOT$Condition <- relevel(s1.TOT$Condition, ref="Control")
# Factor assigning unique number to Prime:Target combinations
s1.TOT$PrimeTarget         <- with(s1.TOT, interaction(Stim.2,Stim.3,drop=T))
levels(s1.TOT$PrimeTarget) <- gsub(".bmp","",levels(s1.TOT$PrimeTarget))
# Factor for 0-based 'time' vector: Trialnumber
s1.TOT$TrialNum <- unlist(llply(unique(s1.TOT$Subj),function(s) return(0:(length(which(s==s1.TOT$Subj))-1))))
s1.TOT$TrialTime <- s1.TOT$TrialNum/max(s1.TOT$TrialNum)
# Factor for unique tool:race combinations
s1.TOT$Prime <- with(s1.TOT, interaction(race,tool,drop=TRUE))


# CHECK NESTED STRUCTURE --------------------------------------------------

# Check nesting of Prime:Target within Subjects...
with(s1.TOT, isNested(PrimeTarget,Subj))
# Prime:Target combinations are NOT nested within Subjects
xtabs(~Subj+PrimeTarget,s1.TOT,drop=T,sparse=T)

with(s1.TOT, isNested(Prime,Subj)) # Not nested, 25 duplicates for each Subject
xtabs(~tool+race,s1.TOT,drop=T,sparse=T)

with(s1.TOT, isNested(PrimeTarget,Prime)) # Prime:Target is nested within Prime
xtabs(~PrimeTarget+Prime,s1.TOT,drop=T,sparse=T)

# This will be the random effect for different stimulus combinations
s1.TOT$Stims <- with(s1.TOT, interaction(Prime,PrimeTarget,drop=TRUE))

# PLOT SLOPES -------------------------------------------------------------
xyplot(scale(RT.log.c,scale=F) ~ TrialTime | Subj, data=s1.TOT[which(s1.TOT$Subj%in%c(107,100,sample(unique(s1.TOT$Subj),14))), ], layout=c(4,4), groups=Prime, xlab = "Time (a.u.)", ylab = "Grand mean centered log(RT) of correct trials", main="WIT: 16 Random Participants\nRandom Slope Model?", 
       prepanel = function(x, y, groups) prepanel.lmline(x, y),
       panel=function(x,y,groups,subscripts){
         panel.xyplot(x,y,groups=groups,subscripts = subscripts)
         #panel.superpose(x,y,panel.groups=panel.loess, groups=groups, subscripts = subscripts, lwd=2,alpha=.6)
         panel.superpose(x,y,panel.groups=panel.lmline, groups=groups, subscripts = subscripts, lwd=3)},
       ylim=c(-1.5,1.5),as.table=T,auto.key = list(points = FALSE, lines = TRUE, columns = 2))

# FULL MODEL --------------------------------------------------------------


#rndID <- which(s1.TOT$Subj%in%sample(unique(s1.TOT$Subj),50)) # Random sample to save time
rndID <- 1:nrow(s1.TOT) # Total sample

# s1.M01 is the "empty" model to compare against when using the Multilevel model for change
# PrimeTarget combinations were administered randomly between Subjects and TrialTime indicates the passage of time
s1.M00 <- lmer(RT.log.c ~ 1 + (1|Subj) + (1|Stims), data=s1.TOT[rndID, ], REML=F) 
s1.M01 <- lmer(RT.log.c ~ 1 + TrialTime + (1|Subj) + (1|Stims), data=s1.TOT[rndID, ], REML=F) 
s1.M02 <- lmer(RT.log.c ~ 1 + TrialTime + (TrialTime|Subj) + (1|Stims), data=s1.TOT[rndID, ], REML=F) 
s1.M03 <- lmer(RT.log.c ~ 1 + TrialTime + (TrialTime|Subj) + (TrialTime|Stims), data=s1.TOT[rndID, ], REML=F) 
anova(s1.M00,s1.M01,s1.M02,s1.M03)
# Use M02, refit with REML
s1.M02m <- lmer(RT.log.c ~ TrialTime + (TrialTime|Subj) + (1|Stims), data=s1.TOT[rndID, ])
(s1.M02m.sum <- summary(s1.M02m))

s1.M10 <- lmer(RT.log.c ~ TrialTime * tool * race + (TrialTime|Subj) + (1|Stims), data=s1.TOT[rndID, ], REML=F) 
s1.M11 <- lmer(RT.log.c ~ TrialTime * tool * race + (TrialTime|Subj) + (tool+race|Subj) + (1|Stims), data=s1.TOT[rndID, ], REML=F)
#s1.M12 <- lmer(RT.log.c ~ TrialTime * Prime + (TrialTime|Subj) + (Prime|Subj) + (1|Stims) + (Prime|Stims), data=s1.TOT[rndID, ], REML=F) 
anova(s1.M02,s1.M10,s1.M11)
# Use M11
s1.M11m <- lmer(RT.log.c ~ TrialTime * tool * race + (TrialTime|Subj) + (tool+race|Subj) + (1|Stims), data=s1.TOT[rndID, ])
#,control=lmerControl(optimizer="Nelder_Mead"))
(s1.M11m.sum <- summary(s1.M11m))

interaction.plot(TrialTime, Subj, s1.M11m, xlab="Time", ylab="log(RT)")

s1.M20 <- lmer(RT.log.c ~ TrialTime + Prime + Condition + (TrialTime|Subj) + (Prime|Subj) + (1|Stims), data=s1.TOT[rndID, ], REML=F) 
s1.M21 <- lmer(RT.log.c ~ TrialTime + Prime * Condition + (TrialTime|Subj) + (Prime|Subj) + (1|Stims), data=s1.TOT[rndID, ], REML=F) 
anova(s1.M11,s1.M20,s1.M21)

str(rr1 <- ranef(s1.M11m))
dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subj"]]


which(rr1$Subj==min(rr1$Subj[[1]]))

if(FALSE) )  ## default

plot(resid(s1.M11m, scaled=T))
pr.M1 <- profile(s1.M1, which="beta_")

xyplot(pr.M1, absVal=TRUE)

fitted <- s1.M1@frame

dotplot(ranef(s1.M1, condVar = TRUE))

print(confint(pr.M1))
xyplot(pr.M1)
densityplot(pr.M1)
splom(pr.M1)

min(
fitted <- s1.M1@frame

dotplot(ranef(s1.M1, condVar = TRUE))

print(confint(pr.M1))
xyplot(pr.M1)
densityplot(pr.M1)
splom(pr.M1)

s1.M2 = lmer(RT.log ~ TrialNum + race + WIT.cond.f + (race|Subj) + (race|PrimeTarget), data=s1.WIT.tools.correct[rndID, ], REML=FALSE) 
summary(s1.M2)

pr.M2 <- profile(s1.M2, optimizer="Nelder_Mead", which="beta_")
xyplot(pr.M2)
densityplot(pr.M2)
splom(pr.M2)
print(confint(pr.M2))

