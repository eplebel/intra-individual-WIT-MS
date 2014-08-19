#power calculations for small 2x2 interaction effect (for p.9 of intra-individual manuscript)
#
#Written by Uri Simonsohn, March 2014
#
#
#In DataColada[17] I propose that 2x2 interaction studies need 2x the sample size
#http://datacolada.org/2014/03/10/17-no-way-interactions
#In a companion ,pdf I show the simple math behind it
#
#
#Simulations are often more persuasive than math, so here it goes.
#I run simulations that compute power for 2 and 4 cell design, the latter testing the interaction
###################################################################################################

#Create function that computes power of Studies 1 and 2, where Study 1  has 2 cells and tests a simple effect
#and Study 2 has 4 cells and tests the interaction

colada17=function(d1,d2,n1,n2,simtot)
{
  #n1: sample size, per cell, study 1
  #n2: sample size, per cell, study 2
  #d1: simple effect M1-M2
  #d2: moderated effect M3-M4, full elimination of effect implies d2=0
  #simtot: how many simulations to run
  
  
  #Here we will store results
  p1=c()    #p-values for Study 1
  p2=c()    #p-values for Study 2
  
  
  for(i in 1:simtot) {
    #draw data 4 samples
    y1=rnorm(n=max(n1,n2),mean=d1)
    y2=rnorm(n=max(n1,n2))
    y3=rnorm(n=max(n1,n2),mean=d2)
    y4=rnorm(n=max(n1,n2))
    
    #GET DATA READY FOR ANOVA  
    y=c(y1,y2,y3,y4)          #the d.v.
    nrep=rep(n2,4)          
    A=rep(c(1,1,0,0),times=nrep) 
    B=rep(c(1,0,1,0),times=nrep)
    
    #STUDY 1
    p1.k=t.test(y1[1:n1],y2[1:n1],var.equal=TRUE)$p.value  #Do a t-test on the first n1 observations
    
    #STUDY 2
    p2.k=anova(lm(y ~ A * B))["A:B", "Pr(>F)"]             #Do anova, keep p-value of the interaction
    
    #Store the results
    p1=c(p1,p1.k)
    p2=c(p2,p2.k)
    
  }
  
  #What share off comparisons are significant
  power1=sum(p1<=.05)/simtot  #Simple test using estimate of variance from 2 cells only
  power2=sum(p2<=.05)/simtot  #Interaction
  
  cat("\nStudy 1 is powered to:",round(power1,2))
  cat("\nStudy 2 is powered to:",round(power2,2))
  
}


#underpowered if run with the same n
colada17(simtot=2000, n1=20,n2=20,d1=1,d2=0)  

#Same power for 2n regardless of n and d
colada17(simtot=2000, n1=20,n2=40,d1=1,d2=0)  
colada17(simtot=2000, n1=50,n2=100,d1=.3,d2=0)
colada17(simtot=2000, n1=150,n2=300,d1=.25,d2=0)

#Need 4n if effect is 70% attenuated
colada17(simtot=2000, n1=25,n2=100,d1=.5, d2=.3*.5)
colada17(simtot=2000, n1=50,n2=200,d1=.5, d2=.3*.5)
colada17(simtot=2000, n1=22,n2=88,d1=.41, d2=.3*.41)





#power calculations for small 2x2 interaction effect (for p.9 of current manuscript)
colada17(simtot=2000, n1=650,n2=1300,d1=.2,d2=0)