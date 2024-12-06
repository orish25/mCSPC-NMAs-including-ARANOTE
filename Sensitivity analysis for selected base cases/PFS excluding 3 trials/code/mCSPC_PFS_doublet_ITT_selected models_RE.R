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
  
  # Informative prior on heterogeneity variance
  het.var.prec <- pow(1.41, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-4.18, het.var.prec) #lognormal distribution
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

# Normal likelihood, identity link, fixed effects
# With adjustment for Multiarm multistage (MAMS) design
# Assumes only a single MAMS trial with na_mams arms
# Data are
# y_mams
# prec_mams (precision matrix for all arms of MAMS trial)
# t_mams
# na_mams (number of arms in the MAMS trial)


model_mams_adjust_re <- function()
{
  # Non-MAMS designs analysed as usual
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    
    var[i] <- pow(se[i,2],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i,2] ~ dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- delta[i] # model for linear predictor
    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    
    delta[i] ~ dnorm(md[i], tau) # trial-specific LOR distributions
    md[i] <- d[t[i,2]] - d[t[i,1]] # mean of treat effects distributions (with multi-arm trial correction)
  }
  
  
  # MAMS designs
  y_mams[1:na_mams] ~ dmnorm(theta_mams[], prec_mams[ , ])
  for(k in 1:na_mams) {
    theta_mams[k] <- delta_mams[k]
      delta_mams[k] ~ dnorm(md_mams[k], tau)
      md_mams[k] <- d[t_mams[k,2]]-d[t_mams[k,1]]
    dev_mams[k] <-
      (y_mams[k]-theta_mams[k]) * (y_mams[k]-theta_mams[k]) * prec_mams[k, k] #Deviance contribution
  }
  
  # Total residual deviance is sum of non-MAMS and MAMS contributions
  totresdev <- sum(dev[]) + sum(dev_mams[])
  
  # Treatment model same as in unadjusted analysis
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  
  # Informative prior on heterogeneity variance - LN(−3.52, 1.61^2 ) Turner 2019 - composite (mortality or morbidity)
  het.var.prec <- pow(1.61, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-3.52, het.var.prec) #lognormal distribution
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
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(mCSPC_PFS_doublet_ITT_selected_models)
t  <- array(c(mCSPC_PFS_doublet_ITT_selected_models$t1, mCSPC_PFS_doublet_ITT_selected_models$t2), dim = c(ns, 2)) 
nt <- max(t) 
y  <- array(c(rep(0, ns), mCSPC_PFS_doublet_ITT_selected_models$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), mCSPC_PFS_doublet_ITT_selected_models$se), dim = c(ns, 2))

study_names <- gsub("#", "", mCSPC_PFS_doublet_ITT_selected_models$`#ID`)
rownames(t) <- rownames(y) <- rownames(se) <- study_names


# Separate out the STAMPEDE/MAMS studies
mams_indices <- grep("STAMPEDE", rownames(y))
# Data for the MAMS trials
y_mams <- y[mams_indices, 2]
se_mams <- se[mams_indices, ] # Not used as SE may not be correct and need covariance
t_mams <- t[mams_indices, ]
na_mams <- length(y_mams)
# Data for the non-MAMS trials (analysed as usual)
y_non_mams <- y[-mams_indices, ]
se_non_mams <- se[-mams_indices, ]
t_non_mams <- t[-mams_indices, ] 
ns_non_mams <- dim(y_non_mams)[1]

# Covariance matrix for the MAMS trial
# From Section 4.3 of the SAP (FFS/PFS can be found there as well)
var_mams <- matrix(c(0.00810098, 0.00296, 0.00342,
                     0.00296, 0.00538584, 0.00131,
                     0.00342, 0.00131, 0.01905604),
                   nrow = 3)
# Inverse of covariance matrix is precision
prec_mams <- solve(var_mams)


# Bugs data for unadjusted model
bugs_data <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt)

# Bugs data for adjusted model
bugs_data_mams <- list(
  y = y_non_mams,
  se = se_non_mams,
  t = t_non_mams,
  ns = ns_non_mams, 
  y_mams = y_mams,
  prec_mams = prec_mams,
  t_mams = t_mams,
  na_mams = na_mams,
  nt = nt)

# Create initial values for MCMC simulation 
# initial values according to the number of parameters
# These are the same for both adjusted and unadjusted models
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

bugs_object_re_mams <- bugs(data = bugs_data_mams, inits = bugs_inits,
                            parameters.to.save = c("d", "totresdev", "rk", "best", "prob", "bayesian_p", "het.var"),
                            model = model_mams_adjust_re, clearWD = TRUE, 
                            summary.only = FALSE,
                            n.iter = (num_sims + burn_in), n.burnin = burn_in,
                            n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)


# Look at the raw bugs_object
bugs_object_re$DIC
# There is a difference but it is not huge
# Not that totresdev and DIC can't be compared as data are different
bugs_object_re$summary[, c("mean", "2.5%", "97.5%")]
bugs_object_re_mams$summary[, c("mean", "2.5%", "97.5%")]

# Format the odds ratios
cross_hr_doublet_re_pfs_select <- cross_effect(bugs_object = bugs_object_re_mams, t_names = c("DAR+ADT", "ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT"), med = TRUE, exp = TRUE)
write.csv(x = cross_hr_doublet_re_pfs_select, file = "cross_hr_doublet_re_pfs_select")

calculate_sucra(bugs_object = bugs_object_re, bugs_data = bugs_data, t_names = c("DAR+ADT", "ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT"))
calculate_sucra(bugs_object = bugs_object_re_mams, bugs_data = bugs_data_mams, t_names = c("DAR+ADT", "ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT"))


# Define your treatment names
treatment_names <- c("DAR+ADT", "ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT")

# Call the function to draw rankograms
sucra_values <- draw.rankogram(bugs_object_re, bugs_data, treatment_names, cumulative = TRUE)
sucra_values_mams <- draw.rankogram(bugs_object_re_mams, bugs_data_mams, treatment_names, cumulative = TRUE)

# If you want to print or use the SUCRA values
print(sucra_values)
print(sucra_values_mams)


dimnames(bugs_object_re_mams$sims.array)[[3]]  # List all parameter names in the MCMC output



#mean heterogeneity SD
library(coda)
sd_samples <- sqrt(bugs_object_re_mams$sims.array[, , "het.var"])
mean_sd <- mean(sd_samples)
sd_cri <- quantile(sd_samples, probs = c(0.025, 0.975))
print (mean_sd)
print (sd_cri)

#average log hazard
# Extract posterior samples for all d parameters
d_samples <- bugs_object_re_mams$sims.array[, , grep("^d\\[", dimnames(bugs_object_re_mams$sims.array)[[3]])]

# Calculate mean and 95% CrI for each d parameter
mean_d <- apply(d_samples, 2, mean)
cri_d <- apply(d_samples, 2, quantile, probs = c(0.025, 0.975))

# Print results
mean_d
cri_d


# Extract the posterior samples for treatment effects (log hazard ratios)
posterior_samples <- bugs_object_re_mams$sims.array[, , grep("^d\\[", rownames(bugs_object_re_mams$summary))]

# Combine all chains into a single 2D matrix
posterior_samples_matrix <- do.call(rbind, lapply(1:dim(posterior_samples)[2], 
                                                  function(x) posterior_samples[, x, ]))

# Compute the posterior mean for each treatment effect
mean_effects <- apply(posterior_samples_matrix, 2, mean)

# Compute the variance-covariance matrix for the treatment effects
var_cov_matrix <- cov(posterior_samples_matrix)

# Print results
print("Posterior means of treatment effects (log hazard ratios):")
print(mean_effects)

print("Posterior variance-covariance matrix:")
print(var_cov_matrix)


