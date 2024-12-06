library(R2OpenBUGS)
library(readxl)
library(grid)

#mHSPC_PFS_data_PO <- read_excel("mHSPC PFS data PO.xlsx")
#rearranging data so first numbered treatment is first

mCSPC_PFS_doublet_ITT_inconsistency <- transform(mCSPC_PFS_doublet_ITT_inconsistency, y = ifelse(t2 > t1, y, abs(y)))
temp <- mCSPC_PFS_doublet_ITT_inconsistency$t1
mCSPC_PFS_doublet_ITT_inconsistency <- transform(mCSPC_PFS_doublet_ITT_inconsistency, t1 = ifelse(t1 > t2, t2, t1))
mCSPC_PFS_doublet_ITT_inconsistency <- transform(mCSPC_PFS_doublet_ITT_inconsistency, t2 = ifelse(t2 == t1, temp, t2))

# Normal likelihood, identity link, fixed effects
# With adjustment for Multiarm multistage (MAMS) design
# Inconsistency model

model_mams_adjust_fe_ume <- function()
{
  # Non-MAMS designs analysed as usual
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    var[i] <- pow(se[i, 2],2) # calculate variances 
    prec[i] <- 1/var[i] # set precisions
    y[i, 2] ~ dnorm(theta[i],prec[i]) # normal likelihood 
    theta[i] <- d[t[i,1],t[i,2]] # model for linear predictor 
    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
  }
  
  # MAMS designs
  y_mams[1:na_mams] ~ dmnorm(theta_mams[], prec_mams[ , ])
  for(k in 1:na_mams) {
    theta_mams[k] <- d[t_mams[k,1], t_mams[k,2]]
    dev_mams[k] <-
      (y_mams[k]-theta_mams[k]) * (y_mams[k]-theta_mams[k]) * prec_mams[k, k] #Deviance contribution
  }
  
  # Total residual deviance is sum of non-MAMS and MAMS contributions
  totresdev <- sum(dev[]) + sum(dev_mams[])
  
  # Treatment model same as in unadjusted analysis
  for (k in 1:nt) { d[k,k] <- 0 }
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.0001) }
    
    
  }
}

#RE inconsistency model
model_mams_adjust_re_ume <- function()
{
  # Non-MAMS designs analysed as usual
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    
    var[i] <- pow(se[i,2],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i,2] ~ dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- delta[i] # model for linear predictor
    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    
    delta[i] ~ dnorm(md[i], tau) # trial-specific LOR distributions
    md[i] <- d[t[i,1],t[i,2]] # mean of treat effects distributions (with multi-arm trial correction)
  }
  
  # MAMS designs
  y_mams[1:na_mams] ~ dmnorm(theta_mams[], prec_mams[ , ])
  for(k in 1:na_mams) {
    theta_mams[k] <- delta_mams[k]
    delta_mams[k] ~ dnorm(md_mams[k], tau)
    md_mams[k] <- d[t_mams[k,1], t_mams[k,2]]
    dev_mams[k] <-
      (y_mams[k]-theta_mams[k]) * (y_mams[k]-theta_mams[k]) * prec_mams[k, k] #Deviance contribution
  }
  
  # Total residual deviance is sum of non-MAMS and MAMS contributions
  totresdev <- sum(dev[]) + sum(dev_mams[])
  
  # Treatment model same as in unadjusted analysis
  for (k in 1:nt) { d[k,k] <- 0 }
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.0001) }
  }
  
  # Informative prior on heterogeneity variance - LN(−3.52, 1.61^2 ) Turner 2019 - composite (mortality or morbidity)
  het.var.prec <- pow(1.61, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-3.52, het.var.prec) #lognormal distribution
  tau <- pow(het.var, -1)
  
}

# Number of MCMC chains and samples
n_chains <- 3
num_sims <- 50000 * n_chains 
burn_in <- 50000 * n_chains	

# Define the bugs data 
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(mCSPC_PFS_doublet_ITT_inconsistency)
t  <- array(c(mCSPC_PFS_doublet_ITT_inconsistency$t1, mCSPC_PFS_doublet_ITT_inconsistency$t2), dim = c(ns, 2)) 
nt <- max(t) 
y  <- array(c(rep(0, ns), mCSPC_PFS_doublet_ITT_inconsistency$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), mCSPC_PFS_doublet_ITT_inconsistency$se), dim = c(ns, 2))

study_names <- gsub("#", "", mCSPC_PFS_doublet_ITT_inconsistency$X.ID)
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

inits1_ume <- list(
  d = t(structure(.Data = c(
    NA,1,1,1,1,
    1,1,1,1,            NA,            NA,
    1,1,1,1,1,1,
    1,            NA,            NA,            NA,1,
    1,1,1,1,1,            NA,
    NA,            NA,            NA,1,1,
    1,1,1,            NA,            NA,            NA,
    NA,            NA,1,1,1,1,
    NA,            NA,            NA,            NA,            NA,
    NA,1,1,1,            NA,            NA,
    NA,            NA,            NA,            NA,            NA,
    1,1,            NA,            NA,            NA,            NA,
    NA,            NA,            NA,            NA,
    1,            NA,            NA,            NA,            NA,
    NA,            NA,            NA,            NA,            NA),
    .Dim = c(9,9))))

inits2_ume <- list(
  d = t(structure(.Data = c(
    NA,2,2,2,2,
    2,2,2,2,            NA,            NA,
    2,2,2,2,2,
    2,2,            NA,            NA,            NA,2,
    2,2,2,2,2,            NA,
    NA,            NA,            NA,2,2,
    2,2,2,            NA,            NA,            NA,
    NA,            NA,2,2,2,2,
    NA,            NA,            NA,            NA,            NA,
    NA,2,2,2,            NA,            NA,
    NA,            NA,            NA,            NA,            NA,
    2,2,            NA,            NA,            NA,            NA,
    NA,            NA,            NA,            NA, 
    2,            NA,            NA,            NA,            NA,
    NA,            NA,            NA,            NA,           NA),
    .Dim = c(9,9))))


inits3_ume <- list(
  d = t(structure(.Data = c(
    NA,0.5,0.5,0.5,0.5,
    0.5,0.5,0.5,0.5,            NA,            NA,
    0.5,0.5,0.5,0.5,0.5,
    0.5,0.5,            NA,            NA,            NA,0.5,
    0.5,0.5,0.5,0.5,0.5,            NA,
    NA,            NA,            NA,0.5,0.5,
    0.5,0.5,0.5,            NA,            NA,            NA,
    NA,            NA,0.5,0.5,0.5,0.5,
    NA,            NA,            NA,            NA,            NA,
    NA,0.5,0.5,0.5,            NA,            NA,
    NA,            NA,            NA,            NA,            NA,
    0.5,0.5,            NA,            NA,            NA,            NA,
    NA,            NA,            NA,            NA, 
    0.5,            NA,            NA,            NA,            NA,
    NA,            NA,            NA,            NA,            NA),
    .Dim = c(9,9))))
bugs_inits_ume <- list(inits1_ume, inits2_ume, inits3_ume)


bugs_object_fe_mams_ume <- bugs(data = bugs_data_mams, inits = bugs_inits_ume,
                                parameters.to.save = c("d", "totresdev"),
                                model = model_mams_adjust_fe_ume, clearWD = TRUE, 
                                summary.only = FALSE,
                                n.iter = (num_sims + burn_in), n.burnin = burn_in,
                                n.chains = n_chains, bugs.seed = 1, debug = TRUE, save.history = TRUE)

summary_fe_ums <- bugs_object_fe_mams_ume$summary[, c("mean", "2.5%", "97.5%")]
bugs_object_fe_mams_ume$DIC


bugs_object_re_mams <- bugs(data = bugs_data_mams, inits = bugs_inits_ume,
                            parameters.to.save = c("d", "totresdev"),
                            model = model_mams_adjust_re_ume, clearWD = TRUE, 
                            summary.only = FALSE,
                            n.iter = (num_sims + burn_in), n.burnin = burn_in,
                            n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)

summary_re_ums <- bugs_object_re_mams$summary[, c("mean", "2.5%", "97.5%")]
bugs_object_re_mams$DIC

###################
###Node-splitting##
###################

direct_evidence <- as.data.frame(mCSPC_PFS_doublet_ITT_inconsistency[(mCSPC_PFS_doublet_ITT_inconsistency$t1 == 2 & mCSPC_PFS_doublet_ITT_inconsistency$t2 == 4),])
direct_evidence$t1 <- 1
direct_evidence$t2 <- 2
#Treatments 2 and 4 relabelled 1 and 2

# Define the bugs data 
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(direct_evidence)
t  <- array(c(direct_evidence$t1, direct_evidence$t2), dim = c(ns, 2)) 
nt <- max(t) 
y  <- array(c(rep(0, ns), direct_evidence$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), direct_evidence$se), dim = c(ns, 2))

study_names <- gsub("#", "", direct_evidence$`#ID`)
rownames(t) <- rownames(y) <- rownames(se) <- study_names


# Bugs data for unadjusted model
bugs_data_direct <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt)

# Number of MCMC chains and samples
n_chains <- 3
num_sims <- 50000 * n_chains 
burn_in <- 50000 * n_chains	


# Create initial values for MCMC simulation 
# initial values according to the number of parameters
# These are the same for both adjusted and unadjusted models
inits1 <- list(d=c( NA, rep(0, nt - 1)))
inits2 <- list(d=c( NA, rep(-1, nt - 1)))
inits3 <- list(d=c( NA, rep(2, nt - 1)))
bugs_inits <- list(inits1, inits2, inits3)


model_normal_identity_fe <- function()
{
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    var[i] <- pow(se[i, 2],2) # calculate variances 
    prec[i] <- 1/var[i] # set precisions
    y[i, 2] ~ dnorm(theta[i],prec[i]) # normal likelihood 
    theta[i] <- d[t[i,2]]-d[t[i,1]] # model for linear predictor 
    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
  }
  
  totresdev <- sum(dev[]) 
  
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  
  # ranking on relative scale                 
  for (k in 1:nt) {
    rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    #rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
  }
}

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
}

bugs_object_fe_direct <- bugs(data = bugs_data_direct, inits = bugs_inits,
                              parameters.to.save = c("d", "totresdev"), #EK added
                              model = model_normal_identity_fe, clearWD = TRUE, 
                              summary.only = FALSE,
                              n.iter = (num_sims + burn_in), n.burnin = burn_in,
                              n.chains = n_chains, bugs.seed = 1, debug = TRUE, save.history = TRUE)

bugs_object_re_direct <- bugs(data = bugs_data_direct, inits = bugs_inits,
                              parameters.to.save = c("d", "totresdev"), #EK added
                              model = model_normal_identity_re, clearWD = TRUE, 
                              summary.only = FALSE,
                              n.iter = (num_sims + burn_in), n.burnin = burn_in,
                              n.chains = n_chains, bugs.seed = 1, debug = TRUE, save.history = TRUE)

summary_direct_evidence_fe <- bugs_object_fe_direct$summary[, c("mean", "2.5%", "97.5%")]
summary_direct_evidence_re <- bugs_object_re_direct$summary[, c("mean", "2.5%", "97.5%")]


##########################################################################################

indirect_evidence <- as.data.frame(mCSPC_PFS_doublet_ITT_inconsistency[(mCSPC_PFS_doublet_ITT_inconsistency$t2 == 5),])
indirect_evidence$t2 <- 3
indirect_evidence$t1 <- ifelse(indirect_evidence$t1 == 2, 1, 2)
#Treatments 2, 4 and 5 rerabelled 1, 2 and 3

# Define the bugs data 
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(indirect_evidence)
t  <- array(c(indirect_evidence$t1, indirect_evidence$t2), dim = c(ns, 2)) 
nt <- max(t) 
y  <- array(c(rep(0, ns), indirect_evidence$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), indirect_evidence$se), dim = c(ns, 2))

study_names <- gsub("#", "", indirect_evidence$`#ID`)
rownames(t) <- rownames(y) <- rownames(se) <- study_names


# Bugs data for unadjusted model
bugs_data_indirect <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt)

inits1 <- list(d=c( NA, rep(0, nt - 1)))
inits2 <- list(d=c( NA, rep(-1, nt - 1)))
inits3 <- list(d=c( NA, rep(2, nt - 1)))
bugs_inits <- list(inits1, inits2, inits3)

bugs_object_fe_indirect <- bugs(data = bugs_data_indirect, inits = bugs_inits,
                                parameters.to.save = c("d", "totresdev"), #EK added
                                model = model_normal_identity_fe, clearWD = TRUE, 
                                summary.only = FALSE,
                                n.iter = (num_sims + burn_in), n.burnin = burn_in,
                                n.chains = n_chains, bugs.seed = 1, debug = TRUE, save.history = TRUE)

bugs_object_re_indirect <- bugs(data = bugs_data_indirect, inits = bugs_inits,
                                parameters.to.save = c("d", "totresdev"), #EK added
                                model = model_normal_identity_re, clearWD = TRUE, 
                                summary.only = FALSE,
                                n.iter = (num_sims + burn_in), n.burnin = burn_in,
                                n.chains = n_chains, bugs.seed = 1, debug = TRUE, save.history = TRUE)

summary_indirect_evidence_fe <- bugs_object_fe_indirect$summary[, c("mean", "2.5%", "97.5%")]
summary_indirect_evidence_re <- bugs_object_re_indirect$summary[, c("mean", "2.5%", "97.5%")]


diff_fe <- bugs_object_fe_indirect$sims.array[,,1][,1] - bugs_object_fe_direct$sims.array[,,1][,1] 
mean <- round(mean(diff_fe), digits = 2)
lower <- round(quantile(diff_fe, probs = c(0.025)), digits = 2)
upper <- round(quantile(diff_fe, probs = c(0.975)), digits = 2)
mean_diff_fe <- paste0(mean, "(", lower, ", ", upper, ")")
#p_value_fe <- sum(diff_fe > 0)/length(diff_fe)

diff_re <- bugs_object_re_indirect$sims.array[,,1][,1] - bugs_object_re_direct$sims.array[,,1][,1] 
mean <- round(mean(diff_re), digits = 2)
lower <- round(quantile(diff_re, probs = c(0.025)), digits = 2)
upper <- round(quantile(diff_re, probs = c(0.975)), digits = 2)
mean_diff_re <- paste0(mean, "(", lower, ", ", upper, ")")

library(ggplot2)
data <- data.frame(evidence = c("Direct", "Indirect"),
                   index = 1:2,
                   Estimate  = c(summary_direct_evidence_re[1, "mean"], summary_indirect_evidence_re[1, "mean"]), 
                   lower = c(summary_direct_evidence_re[1, "2.5%"], summary_indirect_evidence_re[1, "2.5%"]),
                   upper = c(summary_direct_evidence_re[1, "97.5%"], summary_indirect_evidence_re[1, "97.5%"]))

ggplot(data=data, aes(y=index, x=Estimate, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  coord_fixed() +
  theme_classic() +
  xlim(-1.5, 1) +
  ggtitle("Random effects") +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  annotation_custom(grobTree(textGrob(paste("Mean difference", mean_diff_re), x = 0.3, y= 0.5))) +
  scale_y_continuous(name = "", breaks=1:nrow(data), labels=c("Direct", "Indirect"))


data <- data.frame(evidence = c("Direct", "Indirect"),
                   index = 1:2,
                   Estimate  = c(summary_direct_evidence_fe[1, "mean"], summary_indirect_evidence_fe[1, "mean"]), 
                   lower = c(summary_direct_evidence_fe[1, "2.5%"], summary_indirect_evidence_fe[1, "2.5%"]),
                   upper = c(summary_direct_evidence_fe[1, "97.5%"], summary_indirect_evidence_fe[1, "97.5%"]))

ggplot(data=data, aes(y=index, x=Estimate, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  coord_fixed() +
  theme_classic() +
  xlim(-1.5, 1) +
  ggtitle("Fixed effects") +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  annotation_custom(grobTree(textGrob(paste("Mean difference", mean_diff_fe), x = 0.3, y= 0.5))) +
  scale_y_continuous(name = "", breaks=1:nrow(data), labels=c("Direct", "Indirect"))