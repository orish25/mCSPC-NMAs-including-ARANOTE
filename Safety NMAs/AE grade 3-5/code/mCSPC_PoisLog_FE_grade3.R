

#install.packages("R2OpenBUGS") # Install the package
library(R2OpenBUGS) # Load the package

#Poisson likelihood with log
# Fixed effects model for multi-arm trials
# ns= number of studies; nt=number of treatments
# t = treatment matrix
# r = matrix number of events
# E = matrix exposure time (= number of patients in arm multiplied by mean follow-up)
# na = vector number of arms 
model_poisson_log_fe <- function() {
  for(i in 1:ns){                      # LOOP THROUGH STUDIES
    mu[i] ~ dnorm(0,.0001)           # vague priors for all trial baselines
    for (k in 1:na[i]) {             # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k])   # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor
      log(lambda[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]]
      #Deviance contribution
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k]))            }
    #  summed residual deviance contribution for this trial
    resdev[i] <- sum(dev[i,1:na[i]])       
  }   
  totresdev <- sum(resdev[])            #Total Residual Deviance
  d[1]<-0       # treatment effect is zero reference treatment
  # vague priors for treatment effects
  for (k in 2:nt){  d[k] ~ dnorm(0,.0001) }
  
  # ranking on relative scale                  #EK added
  for (k in 1:nt) {
    #rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
  }
}                                     # *** PROGRAM ENDS


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
num_sims <- 10000 * n_chains 
burn_in <- 10000 * n_chains	

# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(mCSPC_grade_3_5)
nt <- max(mCSPC_grade_3_5$t1)
t  <- array(c(mCSPC_grade_3_5$t1, mCSPC_grade_3_5$t2), dim = c(ns, 2)) 
r  <- array(c(mCSPC_grade_3_5$r1, mCSPC_grade_3_5$r2), dim = c(ns, 2)) 
n  <- array(c(mCSPC_grade_3_5$n1, mCSPC_grade_3_5$n2), dim = c(ns, 2))
E <- array(c(mCSPC_grade_3_5$E1, mCSPC_grade_3_5$E2), dim = c(ns, 2))


bugs_data <- list(
  r = r,
  n = n,
  t = t,
  E = E,
  ns = ns,
  nt = nt,
  na = rep(2, ns))

bugs_inits <- function(){
  #chain 1
  list(d = c(NA, rep(0, nt - 1)), mu = rep(0.1, ns))
  #chain 2
  list(d = c(NA, rep(-1, nt - 1)), mu = rep(-0.1, ns))
  #chain 3
  list(d = c(NA, rep(2, nt - 1)), mu = rep(-0.5, ns))
  
}

# Call OpenBUGS
n.thin <- 1
bugs_object_fe<-bugs(data = bugs_data, inits = bugs_inits, n.thin=n.thin,
                     parameters.to.save = c("d", "totresdev", "rk", "best", "prob"),
                     model = model_poisson_log_fe, clearWD = TRUE, 
                     summary.only = FALSE,
                     n.iter = (num_sims + burn_in), n.burnin = burn_in,
                     n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)

# Look at the raw bugs_object
bugs_object_fe$DIC
bugs_object_fe$summary

# Format the odds ratios
cross_meandiff_fe_grade3 <- cross_effect(bugs_object = bugs_object_fe, t_names = c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZA+ADT", "ABI+ADT", "APA+ADT", "SNA+ADT", "ABI+DOC+ADT"), med = TRUE, exp = TRUE)
write.csv(x = cross_meandiff_fe_grade3, file = "cross_meandiff_fe_grade3.csv")

calculate_sucra(bugs_object = bugs_object_fe, bugs_data = bugs_data, t_names = c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZA+ADT", "ABI+ADT", "APA+ADT", "SNA+ADT", "ABI+DOC+ADT"))


library(coda)

#average log hazard
# Extract posterior samples for all d parameters
d_samples <- bugs_object_fe$sims.array[, , grep("^d\\[", dimnames(bugs_object_fe$sims.array)[[3]])]

# Calculate mean and 95% CrI for each d parameter
mean_d <- apply(d_samples, 2, mean)
cri_d <- apply(d_samples, 2, quantile, probs = c(0.025, 0.975))

# Print results
mean_d
cri_d