#plotting subject-specific racial bias slopes
library(ggplot2) 

#anti-Black bias
subj26data = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(Subj==26)) #correct tool trials 
coef(lm(RT ~ race, data = subj26data))
p<-qplot(race, RT, data=subj26data, main="Subject 26", xlab="", ylab="RT (correct trials)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = 247.9, slope=-82.6, colour="red",size=1)  + ylim(0,2250)

#no bias
subj45data = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(Subj==45)) #correct tool trials 
coef(lm(RT ~ race, data = subj45data))
p<-qplot(race, RT, data=subj45data, main="Subject 45", xlab="", ylab="RT (correct trials)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = 296.1, slope=31.6, colour="red",size=1) + ylim(0,2250)


#pro-White bias
subj97data = subset(s2.WIT,((Trial>29 & Trial<80)|(Trial>129 & Trial<180))&(Correct=="True")&(Subj==97)) #correct tool trials 
coef(lm(RT ~ race, data = subj97data))
p<-qplot(race, RT, data=subj97data, main="Subject 97", xlab="", ylab="RT (correct trials)", color=race, shape=race, geom=c("jitter"))
p + geom_abline(intercept = 410.3, slope=143, colour="red",size=1)  + ylim(0,2250)

