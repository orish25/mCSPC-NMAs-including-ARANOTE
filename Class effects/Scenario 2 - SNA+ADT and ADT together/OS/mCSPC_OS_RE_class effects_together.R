# Based on P Orishaba proportional hazards NMA script for FE OS
# Adapted to include correlation between STAMPEDE hazard ratios
# Howard Thom 30-Dec-2023

# Load the R2OpenBUGS package
library(R2OpenBUGS)
library(readxl)

#Load the data
mCSPC_OS_class_effects_together <- read_excel("mCSPC OS class effects together.xlsx")


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

  
  d[1] <- 0 # ADT remains reference and on its own
  
  
  # Link doublet and triplet treatment effects to their classes
  for (k in 1:n_doublets) {
    d[doublet_indices[k]] ~ dnorm(doublet_mean, doublet_tau)
  }
  for (k in 1:n_triplets) {
    d[triplet_indices[k]] ~ dnorm(triplet_mean, triplet_tau)
  } 
  for (k in 1:n_doc) {
    d[doc_indices[k]] ~ dnorm(doc_mean, doc_tau)
  } 
  for (k in 1:n_adt) {
    d[adt_indices[k]] ~ dnorm(adt_mean, adt_tau)
  }
  
  # Informative prior on heterogeneity variance - Turner 2019 - LN(−4.28, 1.61^2)
  het.var.prec <- pow(1.61, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-4.28, het.var.prec) #lognormal distribution
  tau <- pow(het.var, -1)
  
  doublet_mean ~ dnorm(0, 0.0001)
  triplet_mean ~ dnorm(0, 0.0001)
  doc_mean ~ dnorm(0, 0.0001)
  adt_mean ~ dnorm(0, 0.0001)
  doublet_tau <- pow(doublet_sd, -2) 
  doublet_sd ~ dunif(0, 5) 
  triplet_tau <- pow(triplet_sd, -2) 
  triplet_sd ~ dunif(0, 5) 
  doc_tau <- pow(doc_sd, -2) 
  doc_sd ~ dunif(0, 5) 
  adt_tau <- pow(adt_sd, -2) 
  adt_sd ~ dunif(0, 5)
  
  
  
  # ranking on relative scale                 
  for (k in 1:nt) {
    rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    #rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
  }
}


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
  

  d[1] <- 0 # ADT remains reference and on its own
  
  
  # Link doublet and triplet treatment effects to their classes
  for (k in 1:n_doublets) {
    d[doublet_indices[k]] ~ dnorm(doublet_mean, doublet_tau)
  }
  for (k in 1:n_triplets) {
    d[triplet_indices[k]] ~ dnorm(triplet_mean, triplet_tau)
  } 
  for (k in 1:n_doc) {
    d[doc_indices[k]] ~ dnorm(doc_mean, doc_tau)
  } 
  for (k in 1:n_adt) {
    d[adt_indices[k]] ~ dnorm(adt_mean, adt_tau)
  } 
  
  
  # Informative prior on heterogeneity variance - Turner 2019 - LN(−4.28, 1.61^2)
  het.var.prec <- pow(1.61, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-4.28, het.var.prec) #lognormal distribution
  tau <- pow(het.var, -1)
  
  
  doublet_mean ~ dnorm(0, 0.0001)
  triplet_mean ~ dnorm(0, 0.0001)
  doc_mean ~ dnorm(0, 0.0001)
  adt_mean ~ dnorm(0, 0.0001)
  doublet_tau <- pow(doublet_sd, -2) 
  doublet_sd ~ dunif(0, 5) 
  triplet_tau <- pow(triplet_sd, -2) 
  triplet_sd ~ dunif(0, 5) 
  doc_tau <- pow(doc_sd, -2) 
  doc_sd ~ dunif(0, 5) 
  adt_tau <- pow(adt_sd, -2) 
  adt_sd ~ dunif(0, 5)
  
  
  # ranking on relative scale                
  for (k in 1:nt) {
    rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    #rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
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
num_sims <- 50000 * n_chains 
burn_in <- 50000 * n_chains	

# Define the bugs data 
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(mCSPC_OS_class_effects_together)
t  <- array(c(mCSPC_OS_class_effects_together$t1, mCSPC_OS_class_effects_together$t2), dim = c(ns, 2)) 
nt <- max(t) 
y  <- array(c(rep(0, ns), mCSPC_OS_class_effects_together$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), mCSPC_OS_class_effects_together$se), dim = c(ns, 2))

study_names <- gsub("#", "", mCSPC_OS_class_effects_together$`#ID`)
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
var_mams <- matrix(c(0.01179151, 0.00203, 0.01113,
                     0.00203, 0.00665433, 0.00085,
                     0.01113, 0.00085, 0.0318173),
                   nrow = 3)
# Inverse of covariance matrix is precision
prec_mams <- solve(var_mams)


# Bugs data for unadjusted model
bugs_data <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt,
  doublet_indices = c(2, 5, 6, 7),
  triplet_indices = c(3, 9, 10),
  doc_indices = c(4, 11),
  adt_indices = c(1, 8),
  n_doublets = 4,
  n_triplets = 3,
  n_doc = 2,
  n_adt = 2)

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
  nt = nt,
  doublet_indices = c(2, 5, 6, 7),
  triplet_indices = c(3, 9, 10),
  doc_indices = c(4, 11),
  adt_indices = c(1, 8),
  n_doublets = 4,
  n_triplets = 3,
  n_doc = 2,
  n_adt = 2)

# Create initial values for MCMC simulation 
# initial values according to the number of parameters
# These are the same for both adjusted and unadjusted models
inits1 <- list(d=c( NA, rep(0, nt - 1)), het.var= 1,
               doublet_sd = 0.5,
               triplet_sd = 0.5,
               doc_sd = 0.5,
               adt_sd = 0.5,
               doublet_mean = 0,
               triplet_mean = 0,
               doc_mean = 0,
               adt_mean = 0)
inits2 <- list(d=c( NA, rep(-1, nt - 1)), het.var = 2, 
               doublet_sd = 0.1,
               triplet_sd = 0.1,
               doc_sd = 0.1,
               adt_sd = 0.1,
               doublet_mean = -1,
               triplet_mean = -1,
               doc_mean = -1,
               adt_mean = -1)
inits3 <- list(d=c( NA, rep(2, nt - 1)), het.var = 0.5, 
               doublet_sd = 1,
               triplet_sd = 1,
               doc_sd = 1,
               adt_sd = 1,
               doublet_mean = 2,
               triplet_mean = 2,
               doc_mean = 2,
               adt_mean = 2)
bugs_inits <- list(inits1, inits2, inits3)

# Call OpenBUGS

bugs_object_re <- bugs(data = bugs_data, inits = bugs_inits,
                       parameters.to.save = c("d", "totresdev", "doublet_mean", "triplet_mean", "doc_mean", "adt_mean","doublet_sd", "triplet_sd", "doc_sd", "adt_sd", "rk", "best", "prob"),
                       model = model_normal_identity_re, clearWD = TRUE, 
                       summary.only = FALSE,
                       n.iter = (num_sims + burn_in), n.burnin = burn_in,
                       n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)

bugs_object_re_mams <- bugs(data = bugs_data_mams, inits = bugs_inits,
                            parameters.to.save = c("d", "totresdev", "doublet_mean", "triplet_mean", "doc_mean", "adt_mean", "doublet_sd", "triplet_sd", "doc_sd", "adt_sd", "rk", "best", "prob"),
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
cross_meandiff_re_class_all_merged <- cross_effect(bugs_object = bugs_object_re_mams, t_names = c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "SNA+ADT", "ABI+DOC+ADT", "ENZ+DOC+ADT", "SNA+DOC+ADT"), med = TRUE, exp = TRUE)
write.csv(x = cross_meandiff_re_class_all_merged, file = "cross_meandiff_re_class_all_merged.csv")

calculate_sucra(bugs_object = bugs_object_re, bugs_data = bugs_data, t_names = c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "SNA+ADT", "ABI+DOC+ADT", "ENZ+DOC+ADT", "SNA+DOC+ADT"))
calculate_sucra(bugs_object = bugs_object_re_mams, bugs_data = bugs_data_mams, t_names = c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "SNA+ADT", "ABI+DOC+ADT", "ENZ+DOC+ADT", "SNA+DOC+ADT"))

# Define your treatment names
treatment_names <- c("ADT", "DAR+ADT", "DAR+DOC+ADT", "DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "SNA+ADT", "ABI+DOC+ADT", "ENZ+DOC+ADT", "SNA+DOC+ADT")

# Call the function to draw rankograms
sucra_values <- draw.rankogram(bugs_object_re, bugs_data, treatment_names, cumulative = TRUE)
sucra_values_mams <- draw.rankogram(bugs_object_re_mams, bugs_data_mams, treatment_names, cumulative = TRUE)


print(sucra_values_mams)



# Extract the samples for the desired parameters
samples <- bugs_object_re_mams$sims.array[, , c("triplet_mean", "doublet_mean", 
                                                "doc_mean", "adt_mean", "doublet_sd", 
                                                "triplet_sd", "doc_sd", "adt_sd")]

# Flatten the samples into a matrix with 1000 rows (randomly sampling if needed)
set.seed(123)  # For reproducibility
sample_indices <- sample(1:dim(samples)[1], size = 1000, replace = FALSE)
#samples_flattened <- samples[sample_indices, , drop = TRUE]

samples_flattened <- samples[sample_indices, , ]


# Save the samples to a CSV file
write.csv(samples_flattened, file = "MCMC_samples_OS_RE_1.csv", row.names = FALSE)



# Extract the posterior samples for treatment effects (log hazard ratios)
samples_d <- bugs_object_re_mams$sims.array[, , grep("^d\\[", rownames(bugs_object_re_mams$summary))]

# Flatten the samples into a matrix with 1000 rows (randomly sampling if needed)
set.seed(123)  # For reproducibility
sample_indices_d <- sample(1:dim(samples_d)[1], size = 1000, replace = FALSE)

samples_flattened_d <- samples_d[sample_indices_d, , ]


# Save the samples to a CSV file
write.csv(samples_flattened_d, file = "MCMC_samples_OS_RE_d.csv", row.names = FALSE)

