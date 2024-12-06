# P Orishaba & Howard Thom August-2024
# Adapted to include correlation between STAMPEDE hazard ratios


# Load the R2OpenBUGS package
library(R2OpenBUGS)
library(readxl)

# Random effects model
model_normal_identity_re <- function()
{
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    
    var[i] <- pow(se[i,2],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i,2] ~ dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- delta[i] # model for linear predictor
    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    
    delta[i] ~ dnorm(md[i], tau) # trial-specific LOR distributions
    md[i] <- d[t[i,2]] - d[t[i,1]] # mean of treat effects distributions (with multi-arm trial correction)
  }
  totresdev <- sum(dev[]) #Total Residual Deviance
  
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  
  # Informative prior on heterogeneity variance - Turner 2019 - LN(−2.59, 1.62^2)
  het.var.prec <- pow(1.62, -2) # code for 1/(1.62*1.62)
  het.var ~ dlnorm(-2.59, het.var.prec) #lognormal distribution
  tau <- pow(het.var, -1)
  
  # ranking on relative scale                 
  for (k in 1:nt) {
    rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    #rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
  }
  
  for(k in 1:nt) {
    bayesian_p[k] <- step(d[k]) # This is 1 if d>=0
  }
}


# Function to generate summary statistics using coda from MCMC simulations
summary_stats <- function(x, n_digits = 2, med = FALSE)
{
  if(med){
    return(paste0(format(median(x), digits = n_digits, nsmall = n_digits), " (", 
                  format(quantile(x, probs = c(0.025)), digits = n_digits, nsmall = n_digits), ", ", 
                  format(quantile(x, probs = c(0.975)), digits = n_digits, nsmall = n_digits), ")"))
  }else{
    return(paste0(format(mean(x), digits = n_digits, nsmall = n_digits)," (",
                  format(quantile(x, probs = c(0.025)), digits = n_digits, nsmall = n_digits),", ",
                  format(quantile(x, probs = c(0.975)), digits = n_digits, nsmall = n_digits),")"))
  }
}

# Function to generate cross relative treatment effects comparison table
# using the d (e.g. log odds, mean differences) from R2OpenBUGS NMA models
cross_effect <- function(bugs_object, t_names, med = FALSE, exp = TRUE)
{
  effect_table <- matrix(NA, nrow = length(t_names), ncol = length(t_names))
  rownames(effect_table) <- colnames(effect_table) <- t_names
  for(i in 2:length(t_names))
  {
    # If the results need to be exponentiated (e.g. for odds ratios or hazard ratios)
    if(exp) {
      # Comparisons with reference
      effect_table[i, 1] <- summary_stats(exp(bugs_object$sims.array[, , 
                                                                     grep("d", rownames(bugs_object$summary))[i - 1]]), med = med)
      effect_table[1, i] <- summary_stats(exp(-bugs_object$sims.array[, , 
                                                                      grep("d", rownames(bugs_object$summary))[i - 1]]), med = med)
      for(j in 2:length(t_names))
      {
        effect_table[i, j] <- summary_stats(exp(
          bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[i - 1]]-
            bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[j - 1]]
        ), med = med)			
      }  
    } else {
      # If results do not need to be exponentiated (e.g. for mean differences)
      effect_table[i, 1] <- summary_stats(bugs_object$sims.array[, , 
                                                                 grep("d", rownames(bugs_object$summary))[i - 1]], med = med)
      effect_table[1, i] <- summary_stats(-bugs_object$sims.array[, , 
                                                                  grep("d", rownames(bugs_object$summary))[i - 1]], med = med)
      for(j in 2:length(t_names))
      {
        effect_table[i, j] <- summary_stats(
          bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[i - 1]]-
            bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[j - 1]]
          , med = med)			
      }
    }
    
  }
  for(i in 1:length(t_names))effect_table[i, i] <- t_names[i]
  
  return(effect_table)
}

draw.rankogram<-function(bugs.object, bugs.data, t.names, cumulative=FALSE)
{
  # Ranking probability indices
  prob.indices<-grep("prob",rownames(bugs.object$summary))
  
  # Probability that each treatment takes each rank plus cumulative
  rank.probs<-rank.cumprobs<-matrix(NA,nrow=bugs.data$nt,ncol=bugs.data$nt)
  # SUCRA	
  sucra<-rep(NA,bugs.data$nt)
  names(sucra)<-rownames(rank.cumprobs)<-rownames(rank.probs)<-t.names
  colnames(rank.cumprobs)<-colnames(rank.probs)<-paste("Rank",c(1:bugs.data$nt))
  plot(c(0,0),col=0,xlim=c(1,bugs.data$nt),ylim=c(0,1),xlab="Rank",ylab="Probability")
  
  for(i in 1:bugs.data$nt)
  {
    # prob[k,i] is probability treatment i is kth best
    rank.probs[i,]<-bugs.object$summary[prob.indices,"mean"][c(0:(bugs.data$nt-1))*bugs.data$nt+i]
    rank.cumprobs[i,]<-cumsum(rank.probs[i,])
    sucra[i]<-(1/(bugs.data$nt-1))*sum(rank.cumprobs[i,1:bugs.data$nt-1])
    if(!cumulative){
      lines(rank.probs[i,],col=i,lty=i,lwd=3)
    }
    if(cumulative){
      lines(rank.cumprobs[i,],col=i,lty=i,lwd=2)
    }
  
}
if(!cumulative){
  legend("topright",legend=t.names,lwd=3,col=c(1:bugs.data$nt),lty=c(1:bugs.data$nt), cex = 0.3)
}
if(cumulative){	
  legend("bottomright",legend=t.names,lwd=2,col=c(1:bugs.data$nt),lty=c(1:bugs.data$nt), cex = 0.3)
}

return("sucra"=sucra)
}

calculate_sucra <- function(bugs_object, bugs_data, t_names)
{
  # Ranking probability indices
  prob.indices <- grep("prob", rownames(bugs_object$summary))
  
  # Probability that each treatment takes each rank plus cumulative
  rank_probs <- rank_cumprobs <- matrix(NA, nrow = bugs_data$nt, ncol = bugs_data$nt)
  # SUCRA	
  sucra <- rep(NA, bugs_data$nt)
  names(sucra) <- rownames(rank_cumprobs) <- rownames(rank_probs) <- t_names
  colnames(rank_cumprobs) <- colnames(rank_probs) <- paste("Rank", c(1:bugs_data$nt))
  for(i in 1:bugs_data$nt)
  {
    # prob[k, i] is probability treatment i is kth best
    rank_probs[i, ] <- bugs_object$summary[prob.indices, "mean"][c(0:(bugs_data$nt - 1))*bugs_data$nt+i]
    rank_cumprobs[i, ] <- cumsum(rank_probs[i, ])
    sucra[i] <- (1/(bugs_data$nt - 1))*sum(rank_cumprobs[i, 1:bugs_data$nt - 1])
  }
  
  return("sucra" = sucra)
}

# Number of MCMC chains and samples
n_chains <- 3
num_sims <- 10000 * n_chains 
burn_in <- 10000 * n_chains	

# Define the bugs data 
ns <- nrow(mCSPC_pain_progression)
t  <- array(c(mCSPC_pain_progression$t1, mCSPC_pain_progression$t2), dim = c(ns, 2)) 
nt <- max(t)
y  <- array(c(rep(0, ns), mCSPC_pain_progression$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), mCSPC_pain_progression$se), dim = c(ns, 2))

bugs_data <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt)


# Create initial values for MCMC simulation 
# initial values according to the number of parameters
inits1 <- list(d=c( NA, rep(0, nt - 1)), het.var= 1)
inits2 <- list(d=c( NA, rep(-1, nt - 1)), het.var = 2)
inits3 <- list(d=c( NA, rep(2, nt - 1)), het.var = 0.5)
bugs_inits <- list(inits1, inits2, inits3)

# Call OpenBUGS

bugs_object_re <- bugs(data = bugs_data, inits = bugs_inits,
                       parameters.to.save = c("d", "totresdev", "rk", "best", "prob", "bayesian_p", "het.var"),
                       model = model_normal_identity_re, clearWD = TRUE, 
                       summary.only = FALSE,
                       n.iter = (num_sims + burn_in), n.burnin = burn_in,
                       n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)


# Look at the raw bugs_object
bugs_object_re$DIC
# There is a difference but it is not huge
# Not that totresdev and DIC can't be compared as data are different
bugs_object_re$summary[, c("mean", "2.5%", "97.5%")]


# Format the odds ratios
cross_hr_os_doublet_pain <- cross_effect(bugs_object = bugs_object_re, t_names = c("ADT", "DAR+ADT", "ENZ+ADT", "APA+ADT", "ABI+ADT"), med = TRUE, exp = TRUE)
write.csv(x = cross_hr_os_doublet_pain, file = "cross_hr_os_doublet_pain")

calculate_sucra(bugs_object = bugs_object_re, bugs_data = bugs_data, t_names = c("ADT", "DAR+ADT", "ENZ+ADT", "APA+ADT", "ABI+ADT"))


# Define your treatment names
treatment_names <- c("ADT", "DAR+ADT", "ENZ+ADT", "APA+ADT", "ABI+ADT")

# Call the function to draw rankograms
sucra_values <- draw.rankogram(bugs_object_re, bugs_data, treatment_names, cumulative = TRUE)

# If you want to print or use the SUCRA values
print(sucra_values)


#mean heterogeneity SD
library(coda)
sd_samples <- sqrt(bugs_object_re$sims.array[, , "het.var"])
mean_sd <- mean(sd_samples)
sd_cri <- quantile(sd_samples, probs = c(0.025, 0.975))
# Print results
mean_sd
sd_cri



#average log hazard
# Extract posterior samples for all d parameters
d_samples <- bugs_object_re$sims.array[, , grep("^d\\[", dimnames(bugs_object_re$sims.array)[[3]])]

# Calculate mean and 95% CrI for each d parameter
mean_d <- apply(d_samples, 2, mean)
cri_d <- apply(d_samples, 2, quantile, probs = c(0.025, 0.975))

# Print results
mean_d
cri_d
