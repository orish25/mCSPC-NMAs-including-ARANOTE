# P. Orishaba 18/07/2024
# Poisson likelihood, log link, random efffects NMA in R2OpenBUGS

# Load the R2OpenBUGS package
library(R2OpenBUGS)


# Poisson likelihood, log link, random effects

model_poisson_log_re <- function(){ # *** PROGRAM STARTS
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001) # vague priors for all trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      log(lambda[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  #sd ~ dunif(0,5) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  hr <- exp(d[2]) # Hazard ratio for reduced fat dietary intervention
  
  #Informative prior on heterogeneity variance - Turner2019 LN(−1.97, 1.60^2)
  var_d <- pow(1.60, -2) # code for 1/(1.41*1.41)
  var ~ dlnorm(-1.97, var_d) #lognormal distribution
  tau <- pow(var, -1)
  
  # ranking on relative scale                  
  for (k in 1:nt) {
    #rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
  }
} # *** PROGRAM ENDS


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
num_sims <- 100000 * n_chains 
burn_in <- 100000 * n_chains	

# Load the data
#trial_data <- read.csv("Practical1_exercise2_data.csv")
# Store BUGS data in matrices
ns <- nrow(mCSPC_discontinuation_due_to_AE)
nt <- max(mCSPC_discontinuation_due_to_AE$t1)
t  <- array(c(mCSPC_discontinuation_due_to_AE$t1, mCSPC_discontinuation_due_to_AE$t2), dim = c(ns, 2)) 
r  <- array(c(mCSPC_discontinuation_due_to_AE$r1, mCSPC_discontinuation_due_to_AE$r2), dim = c(ns, 2)) 
n  <- array(c(mCSPC_discontinuation_due_to_AE$n1, mCSPC_discontinuation_due_to_AE$n2), dim = c(ns, 2))
E <- array(c(mCSPC_discontinuation_due_to_AE$E1, mCSPC_discontinuation_due_to_AE$E2), dim = c(ns, 2))

# Define the bugs data
bugs_data <- list(
  r = r,
  n = n,
  t = t,
  E = E,
  ns = ns,
  nt = nt,
  na = rep(2, ns))

# Create initial values for MCMC simulation
inits1 <- list(d=c( NA, rep(0, nt - 1)), mu=rep(-0.5, ns), var=1)
inits2 <- list(d=c( NA, rep(-1, nt - 1)), mu=rep(-0.5, ns), var=2)
inits3 <- list(d=c( NA, rep(2, nt - 1)), mu=rep(-3, ns), var=0.5)
bugs_inits <- list(inits1, inits2, inits3)

# Call OpenBUGS
n.thin=10
bugs_object_re <- bugs(data = bugs_data, inits = bugs_inits,
                     parameters.to.save = c("d", "totresdev", "rk", "best", "prob", "var"),
                     model = model_poisson_log_re, clearWD = TRUE, 
                     summary.only = FALSE,
                     n.iter = (num_sims + burn_in), n.burnin = burn_in,
                     n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)

# Look at the raw bugs_object
bugs_object_re$DIC
bugs_object_re$summary

# Format the odds ratios
cross_hr_re_discont_AEs <- cross_effect(bugs_object = bugs_object_re, t_names = c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "ABI+DOC+ADT"), med = TRUE, exp = TRUE)
write.csv(x = cross_hr_re_discont_AEs, file = "cross_hr_re_discont_AEs.csv")

calculate_sucra(bugs_object = bugs_object_re, bugs_data = bugs_data, t_names = c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "ABI+DOC+ADT"))



#mean heterogeneity SD
library(coda)
sd_samples <- sqrt(bugs_object_re$sims.array[, , "var"])
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



