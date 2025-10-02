#Teitsworth et al. 


#Set working directory
setwd("*")

# Packages used
library(ggplot2)
library(unmarked)
library(corrplot)
library(AICcmodavg)
library(dplyr)


# Detection data: "0" or "1" in 20 Visits - 4 visits in each of 5 primary periods

# Read in raw data
Data<-read.csv("Data.csv", header=TRUE)
det.data<-FullData[, c(2:21)] #Detection history
surv.cov<-FullData[, c(31:70)] #Detection Covs
site.cov<-FullData[, c(22:30)] #Static Covs
years.cov<-read.csv("YearlySiteCovs.csv", header=TRUE) #Dynamic Covs

##Calling out specific covs to use

#static
sub<-scale(site.cov$Bottom.Substrate)  #scaled Bottom substrate score from habitat assessment
cover<-scale(site.cov$Cover.Score) #scaled Cover score from habitat assessment
TS<-scale(site.cov$Total) #scaled Total Score from habitat assessment
devel<-site.cov$devel # % developed land cover (HAiFLS : Peterson and Pearse 2017)
grass<-site.cov$grass # % herbaceous/pasture land cover (HAiFLS : Peterson and Pearse 2017)
crop<-site.cov$crop # % crop land cover (HAiFLS : Peterson and Pearse 2017)
wet<-site.cov$wetland # % wetland land cover (HAiFLS : Peterson and Pearse 2017)

#Dynamic
TQ<- years.cov$TQ #Yearly TQmean
DI<- years.cov$DI #Yearly Max Drought Index
DI<-as.factor(DI)


#Setting up initial occupancy (psi) covariates
sc<-cbind(sub, cover, TS, devel, grass, crop, wet)
sc<-as.data.frame(sc) #convert to dataframe that unmarked can use  
colnames(sc)<-c("sub","cover", "TS", "devel", "grass", "crop", "wet")

#Setting up colonization/extinction (gamma/epsilon) covariates 
ysc<-data.frame(TQ, DI)

#Setting up detection (p) covariates
discharge<-as.matrix(surv.cov[,1:20]) #mean daily discharge
discharge<-log(discharge)
discharge<-t(discharge) #transpose to 20 rows (total visits) with 176 columns (sites)
discharge<-matrix(discharge, dimnames=list(t(outer(colnames(discharge), rownames(discharge), FUN=paste)), NULL)) #stack all sites into single column (ordered site 1 discharge, site 1 discharge 2, etc.)
bait<-as.matrix(FullData[,51:70]) #bait age, in days
bait<-t(bait) #transpose 
bait<-matrix(bait, dimnames=list(t(outer(colnames(bait), rownames(bait), FUN=paste)), NULL)) #stack into one column

obs<-cbind(discharge, bait) #combine the two detection cov columns
colnames(obs)<-c("discharge", "bait")
obs<- as.data.frame(obs) #convert to dataframe that unmarked can use



##### Define the Unmarked dataframe #####
umf<-unmarkedMultFrame(y=det.data,siteCovs=sc, numPrimary=5, 
                       yearlySiteCovs=ysc,
                       obsCovs=obs)      



####Exploration/hypothesis testing####
# models

s1<- colext(psiformula= ~sub, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s2<- colext(psiformula= ~cover, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s3<- colext(psiformula= ~sub + cover, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s4<- colext(psiformula= ~TS, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s5<- colext(psiformula= ~TS + devel, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s6<- colext(psiformula= ~TS + grass, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s7<- colext(psiformula= ~TS + crop, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s8<- colext(psiformula= ~TS + wet, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s9<- colext(psiformula= ~TS*devel, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s10<- colext(psiformula= ~TS*grass, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s11<- colext(psiformula= ~TS*wet, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s12<- colext(psiformula= ~TS*crop, gammaformula= ~1, epsilonformula= ~1, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s13<- colext(psiformula= ~TS, gammaformula= ~1, epsilonformula= ~TS, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s14<- colext(psiformula= ~TS, gammaformula= ~1, epsilonformula= ~TQ, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s15<- colext(psiformula= ~TS, gammaformula= ~1, epsilonformula= ~DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s16<- colext(psiformula= ~TS, gammaformula= ~1, epsilonformula= ~devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s17<- colext(psiformula= ~TS, gammaformula= ~1, epsilonformula= ~grass, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s18<- colext(psiformula=~TS, gammaformula= ~1, epsilonformula= ~TS + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s19<- colext(psiformula=~TS, gammaformula= ~1, epsilonformula= ~TS + grass, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s20<- colext(psiformula=~TS, gammaformula= ~1, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s21<- colext(psiformula=~TS, gammaformula= ~1, epsilonformula= ~TS + TQ, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s22<- colext(psiformula=~TS, gammaformula= ~1, epsilonformula= ~TS + DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s23<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~TS + DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s24<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s25<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~TS, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s26<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~TS + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s27<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s28<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s29<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s30<- colext(psiformula=~TS, gammaformula= ~TS, epsilonformula= ~devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s31<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~TS + DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s32<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s33<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~TS, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s34<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~TS + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s35<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s36<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s37<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s38<- colext(psiformula=~TS, gammaformula= ~devel, epsilonformula= ~devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s39<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~TS + DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s40<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s41<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~TS, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s42<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~TS + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s43<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s44<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s45<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s46<- colext(psiformula=~TS, gammaformula= ~grass, epsilonformula= ~devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s47<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~TS + DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s48<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s49<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~TS, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s50<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~TS + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s51<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s52<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s53<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s54<- colext(psiformula=~TS, gammaformula= ~DI, epsilonformula= ~devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s55<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~TS + DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s56<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s57<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~TS, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s58<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~TS + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s59<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~DI + devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s60<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~TS + DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s61<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~DI, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))
s62<- colext(psiformula=~TS, gammaformula= ~TQ, epsilonformula= ~devel, pformula= ~discharge, umf, se=TRUE, control=list(maxit=1000))


#### global model
global.model<- colext(psiformula = ~TS + devel + grass + crop + wet,
                      gammaformula =  ~TQ + DI + TS + devel + grass,
                      epsilonformula = ~TQ + DI + TS + devel + grass,
                      pformula = ~discharge + bait,
                      umf, se=TRUE)  
gof<- mb.gof.test(global.model, nsim=999, parallel = FALSE) # Assess goodness of fit
#c-hat=0.23  -> underdispersed (i.e., more covariates than needed given the observed data)




#### Model selection ####
ml2<- list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20,
           s21, s22, s23, s24, s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s36, s37, s38, s39, s40,
           s41, s42, s43, s44, s45, s46, s47, s48, s49, s50, s51, s52, s53, s54, s55, s56, s57, s58, s59, s60,
           s61, s62)
tm<-aictab(ml2, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1) 


#### Results ####

#Top ranked model is s31
#Psi - unsupported positive effect of total score
#gamma - unsupported positive effect of developed
#epsilon - supported negative effect of total score
#        - supported positive effect of drought
#        - supported positive effect of developed
#p - supported positive effect of discharge

#Second ranked model is s22 (delta AICc 2.33)
#same as s31, but with no effect of developed on gamma


# Analyzing model s31 

# getting 95% confidence intervals of each parameter
confint(s31, type= "psi", level=0.95)
confint(s31, type= "col", level=0.95)
confint(s31, type= "ext", level=0.95)
confint(s31, type= "det", level=0.95)

# Bootstrap the standard errors for smoothed trajectory (i.e., the finite sample estimates)
s31 <- nonparboot(s31, B = 100)  # This takes a while!
s31@smoothed.mean.bsse


# Save Data
s31_col <- predict(s31, type="col", se=TRUE)
s31_col1<- s31_col[seq(1,nrow(s31_col),5),] #Only selecting every 5th row (i.e., Year1 estimate) because colonization is constant across years, but site-specific
write.csv(s31_col1,"*/s31_col.csv")

s31_ext <- predict(s31, type="ext", se=TRUE)
write.csv(s31_ext,"*/s31_ext.csv")

s31_p <- predict(s31, type="det", se=TRUE)
write.csv(s31_p,"*/s31_p.csv")

s31_psi <- predict(s31, type="psi", se=TRUE)
write.csv(s31_psi,"*/s31_psi.csv")

# Pull out likelihood of occurrence (finite sample estimate. fixed to 1.0 if true detection)
s31_smoothed <- matrix(s31@smoothed[2,,], 176, 5, byrow=TRUE) 
  # This takes each array slice in smoothed (i.e., data for each site) and converts
  # it to a single matrix where each site is a row and each year is a column
write.csv(s31_smoothed, "*/s31_smoothed.csv")

# Smoothed Predictions of annual site occurrence under the top model (s31) #
s84@smoothed.mean 
sitebyyear<-s84@smoothed



#### Example of plotting the effects under top model #####

# Total Score score on epsilon

# Create a data set for each level of drought "DI"
newdata<- data.frame(TS = seq(from = min(TS), #xaxis variable
                              to = max(TS),
                              length.out=880), DI = "0",devel=mean(devel))
newdata1<- data.frame(TS = seq(from = min(TS), #xaxis variable
                               to = max(TS),
                               length.out=880), DI = "1",devel=mean(devel))
newdata2<- data.frame(TS = seq(from = min(TS), #xaxis variable
                               to = max(TS),
                               length.out=880), DI = "2",devel=mean(devel))
newdata3<- data.frame(TS = seq(from = min(TS), #xaxis variable
                               to = max(TS),
                               length.out=880), DI = "3",devel=mean(devel))

# Make fit predictions based on newdata
pred.tse<-predict(s31, newdata, se.fit=TRUE, type="ext", backTransform=TRUE)
pred.tse1<-predict(s31, newdata1, se.fit=TRUE, type="ext", backTransform=TRUE)
pred.tse2<-predict(s31, newdata2, se.fit=TRUE, type="ext", backTransform=TRUE)
pred.tse3<-predict(s31, newdata3, se.fit=TRUE, type="ext", backTransform=TRUE)

# Plot
ggplot(data=pred.tse, aes(seq(from = min(FullData$Total), #xaxis variable
                              to = max(FullData$Total),
                              length.out=880), y=Predicted)) +
  geom_line(size=1.5, color="gray40") +
  labs(x = "Total Habitat Score", y = "Extinction Probability", col=NULL)+
  ylim(-0.01, 1.0) +
  xlim(35, 100) +
  geom_ribbon(data=pred.tse, aes(x=seq(from = min(FullData$Total), #xaxis variable
                                       to = max(FullData$Total),
                                       length.out=880),ymin=pred.tse$lower, ymax=pred.tse$upper), alpha=0.15) +
  geom_ribbon(data=pred.tse1, aes(x=seq(from = min(FullData$Total), #xaxis variable
                                        to = max(FullData$Total),
                                        length.out=880),ymin=pred.tse1$lower, ymax=pred.tse1$upper), alpha=0.15, fill="cyan") +
  geom_ribbon(data=pred.tse2, aes(x=seq(from = min(FullData$Total), #xaxis variable
                                        to = max(FullData$Total),
                                        length.out=880),ymin=pred.tse2$lower, ymax=pred.tse2$upper), alpha=0.15, fill="magenta") +
  geom_ribbon(data=pred.tse3, aes(x=seq(from = min(FullData$Total), #xaxis variable
                                        to = max(FullData$Total),
                                        length.out=880),ymin=pred.tse3$lower, ymax=pred.tse3$upper), alpha=0.15, fill="orange2") +
  geom_line(data=pred.tse1,size=1.5, color="cyan3")+
  geom_line(data=pred.tse2,size=1.5, color="magenta")+
  geom_line(data=pred.tse3,size=1.5, color="orange")+
  theme_classic(base_size = 12) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        axis.title.x=element_text(margin=margin(t=20)),
        axis.title.y=element_text(margin=margin(r=20)))
