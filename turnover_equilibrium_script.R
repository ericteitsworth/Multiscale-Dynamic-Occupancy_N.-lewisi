#Teitsworth et al. 

### Estimating turnover probability ###

# vignette from CRAN unmarked 
#(https://cran.r-project.org/web/packages/unmarked/vignettes/colext.html)


# Copied example of the unmarked code for model
#turnover <- function(fm) {
#  psi.hat <- plogis(coef(fm, type="psi"))
#  if(length(psi.hat) > 1)
#    stop("this function only works if psi is scalar")
#  T <- getData(fm)@numPrimary
#  tau.hat <- numeric(T-1)
#  gamma.hat <- plogis(coef(fm, type="col"))
#  phi.hat <- 1 - plogis(coef(fm, type="ext"))
#  if(length(gamma.hat) != T-1 | length(phi.hat) != T-1)
#    stop("this function only works if gamma and phi T-1 vectors")
#  for(t in 2:T) {
#    psi.hat[t] <- psi.hat[t-1]*phi.hat[t-1] +
#      (1-psi.hat[t-1])*gamma.hat[t-1]
#    tau.hat[t-1] <- gamma.hat[t-1]*(1-psi.hat[t-1]) / psi.hat[t]
#  }
#  return(tau.hat)
#}


# Our adaptation, since the model has covariates

# Estimating turnover by subpopulation (i.e., management unit)
# (EXAMPLE)

# Eno/Flat
turnover_FE <- function(fm) {
  psi.hat <- predict(fm, type="psi") # site-specific estimated initial occupancies
  psi.hat <- data.frame(psi.hat, FullData$MU) # Add subpopulation (management unit) for subsetting purposes
  psi.hat_FE <- mean(subset(psi.hat, FullData.MU == "FE")$Predicted) #mean of eno/flat sites
  T <- getData(fm)@numPrimary
  tau.hat_FE <- numeric(T-1)
  # Pull out "Predicted" from backtransformed site-specific, year-specific estimated colonization probs
  gamma.hat <- array(predict(fm, type="col")$Predicted, dim = c(5, 176))
  # Pull out Eno/Flat sites
  gamma.hat_FE <- gamma.hat[,c(1:8,11,126,130:132,135,139,155:156)]
  # Global mean colonization each year for eno/flat sites
  mean.gamma.hat_FE <- c(mean(gamma.hat_FE[1,]),mean(gamma.hat_FE[2,]),mean(gamma.hat_FE[3,]),mean(gamma.hat_FE[4,]))
  # Do the same with extinction probabilities
  eps.hat <- array(predict(fm, type="ext")$Predicted, dim = c(5, 176))
  eps.hat_FE <- eps.hat[,c(1:8,11,126,130:132,135,139,155:156)]
  mean.eps.hat_FE <- c(mean(eps.hat_FE[1,]),mean(eps.hat_FE[2,]),mean(eps.hat_FE[3,]),mean(eps.hat_FE[4,]))
  mean.phi.hat_FE <- 1-mean.eps.hat_FE
  
  # Calculate turnover
  for(t in 2:T){
    psi.hat_FE[t] <- psi.hat_FE[t-1]*mean.phi.hat_FE[t-1] + (1-psi.hat_FE[t-1])*mean.gamma.hat_FE[t-1]
    tau.hat_FE[t-1] <- mean.gamma.hat_FE[t-1]*(1-psi.hat_FE[t-1]) / psi.hat_FE[t]
  }
  return(tau.hat_FE)
}  

# Parametric bootstrapping to get 95% confidence intervals of turnover probabilities  
pb_FE <- parboot(s31, statistic=turnover_FE, nsim=500)
turnCI_FE <- cbind(pb_FE@t0,
                   t(apply(pb_FE@t.star, 2, quantile, probs=c(0.025, 0.975))))
colnames(turnCI_FE) <- c("tau", "lower", "upper")
turnCI_FE


# Plotting the eno/flat turnover 
tci_FE <- data.frame(c(1,2,3,4), turnCI_FE)
colnames(tci_FE) <- c("year","tau", "lower", "upper")

ggplot(tci_FE, aes(x = year, y = tau)) +
  ggtitle("Eno/Flat") +
  ylim(0, 1.0) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2)




### Estimating equilibrium occupancy probability (stability) ###

# Eno/Flat
equilibrium_FE <- function(fm) {
  T <- getData(fm)@numPrimary
  psi.eq.hat_FE <- numeric(T-1)
  # Pull out "Predicted" from backtransformed site-specific, year-specific estimated colonization probs
  gamma.hat <- array(predict(fm, type="col")$Predicted, dim = c(5, 176))
  # Pull out Eno/Flat sites
  gamma.hat_FE <- gamma.hat[,c(1:8,11,126,130:132,135,139,155:156)]
  # Global mean colonization each year for eno/flat sites
  mean.gamma.hat_FE <- c(mean(gamma.hat_FE[1,]),mean(gamma.hat_FE[2,]),mean(gamma.hat_FE[3,]),mean(gamma.hat_FE[4,]))
  # Do the same with extinction probabilities
  eps.hat <- array(predict(fm, type="ext")$Predicted, dim = c(5, 176))
  eps.hat_FE <- eps.hat[,c(1:8,11,126,130:132,135,139,155:156)]
  mean.eps.hat_FE <- c(mean(eps.hat_FE[1,]),mean(eps.hat_FE[2,]),mean(eps.hat_FE[3,]),mean(eps.hat_FE[4,]))
  
  # Calculate equilibrium
  for(t in 1:T-1){
    psi.eq.hat_FE[t] <- mean.gamma.hat_FE[t] / (mean.gamma.hat_FE[t] + mean.eps.hat_FE[t])
  }
  return(psi.eq.hat_FE)
}

# Parametric bootstrapping to get 95% confidence intervals of equilibrium probabilities  
eq.pb_FE <- parboot(s31, statistic=equilibrium_FE, nsim=500)
eqCI_FE <- cbind(eq.pb_FE@t0,
                 t(apply(eq.pb_FE@t.star, 2, quantile, probs=c(0.025, 0.975))))
colnames(eqCI_FE) <- c("psi.eq", "lower", "upper")
eqCI_FE


# Plotting the eno/flat equilibrium
eci_FE <- data.frame(c(1.05,2.05,3.05,4.05), eqCI_FE)
colnames(eci_FE) <- c("year","psi.eq", "lower", "upper")

ggplot(eci_FE, aes(x = year, y = psi.eq)) +
  ggtitle("Eno/Flat") +
  ylim(0, 1.0) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2)


# Calculate non-equilibrium (i.e., difference between current occupancy [finite sample estimate] and equilibrium)
# to evaluate if subpopulation is currently above, below, or at the estimate required to maintain stability.

avgeq<- mean(eqCI_FE[,1]) # mean of equilibrium probability estimates across Eno/Flat sites
avg_fs<- mean(sci_FE[,2]) # mean of finite sample estimates across Eno/Flat sites

avg_diff<- avg_fs-avgeq # mean non-equilibrium of Eno/Flat subpopulation
  