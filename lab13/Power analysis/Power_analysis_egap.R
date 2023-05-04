############################################################################
# Power analysis example
#
# from https://egap.org/resource/10-things-to-know-about-statistical-power/
############################################################################

rm(list = ls())

set.seed(123)

possible.ns <- seq(from=100, to=2000, by=40) # The sample sizes we'll be considering
stopifnot(all( (possible.ns %% 2)==0 )) ## require even number of experimental pool
powers <- rep(NA, length(possible.ns)) # Empty object to collect simulation estimates 
alpha <- 0.05 # Standard significance level 
sims <- 500 # Number of simulations to conduct for each N 
#### Outer loop to vary the number of subjects #### 
for (j in 1:length(possible.ns)){ N <- possible.ns[j] # Pick the jth value for N 
Y0 <- rnorm(n=N, mean=60, sd=20) # control potential outcome 
tau <- 5 # Hypothesize treatment effect 
Y1 <- Y0 + tau # treatment potential outcome                                   
significant.experiments <- rep(NA, sims) # Empty object to count significant experiments 

#### Inner loop to conduct experiments "sims" times over for each N #### 
Y0 <- rnorm(n=N, mean=60, sd=20) # control potential outcome 
tau <- 5 # Hypothesize treatment effect 
Y1 <- Y0 + tau # treatment potential outcome 
for (i in 1:sims){ 
  ## Z.sim <- rbinom(n=N, size=1, prob=.5) # Do a random assignment  by coin flip
  Z.sim <- sample(rep(c(0,1),N/2)) ## Do a random assignment ensuring equal sized groups
  Y.sim <- Y1*Z.sim + Y0*(1-Z.sim) # Reveal outcomes according to assignment 
  fit.sim <- lm(Y.sim ~ Z.sim) # Do analysis (Simple regression) 
  p.value <- summary(fit.sim)$coefficients[2,4] # Extract p-values 
  significant.experiments[i] <- (p.value <= alpha) # Determine significance according to p <= 0.05
}
powers[j] <- mean(significant.experiments) # store average success rate (power) for each N 
} 
powers 

hist(powers)
