##########################################
# Design Effects Exercise
# 
# from Cyrus to the masses
##########################################

set.seed(4444)
library(estimatr)

################################
##### Simulated population #####
################################

# N. of clusters
N_clusters <- 100

# Obs for each cluster
cluster_size <- 10

# Population size
N_pop <- N_clusters*cluster_size

# IDs for clusters and obs
cluster_index <- (1:100)%x%rep(1, cluster_size)

unit_index <- rep(1, N_clusters)%x%(1:cluster_size)

# SD of between-cluster REs
between_sd <- 2

# Generate the between-cluster REs
cluster_random_effect <- rnorm(N_clusters, sd=between_sd)

# Attribute to each obs the RE of its cluster
cluster_random_effect_l <- cluster_random_effect%x%rep(1, cluster_size)

# Create a cluster-level covariate, constant for units within the same cluster 
x <- (0.25*cluster_random_effect + rnorm(N_clusters, sd=1))%x%rep(1, cluster_size)

# Potential outcomes that exhibit high intra-cluster correlation
within_sd <- 1
y0 <- 0.25*cluster_random_effect_l + 0.25*x + rnorm(N_pop, sd=within_sd)
y1 <- y0 + .1*x + rnorm(N_pop, sd=within_sd)
dt_cluster <- data.frame( id = 1:N_pop,
                          cluster_index,
                          unit_index,
                          x,
                          y0,
                          y1)



## Simulation setup:

# N. of simulations
nsim <- 10000

# Empty objects to store the estimates for each iteration
b_block_2 <- vhat_block_2 <- 
  b_block_1 <- vhat_block_1 <-
  b_block_nostrat <- vhat_block_nostrat <-
  b_cluster <- vhat_cluster <- 
  b_unit <- vhat_unit <- 
  rep(NA, nsim)


######################################
##### (1) Unit random assignment #####
######################################

for(i in 1:nsim){
  # Create treatment indicator
  dt_cluster$z <- 0
  
  # Randomly assign half obs to treatment
  dt_cluster$z[dt_cluster$id %in% sample(dt_cluster$id, N_pop/2)] <- 1
  
  # Observed outcome for treated and control units
  dt_cluster$y <- dt_cluster$z*dt_cluster$y1 + (1-dt_cluster$z)*dt_cluster$y0
  
  # Compute diff in means
  fit <- lm_robust(y~z, data=dt_cluster)
  
  # Store diff in means estimate
  b_unit[i] <- coef(fit)["z"]
  
  # Store estimated var of the estimator
  vhat_unit[i] <- diag(vcov(fit))["z"]
}


## Unbiasednes and variance:

# True DiM
mean(y1 - y0)

# Expected DiM from the randomization design
mean(b_unit)

# True variance of the DiM estimator
var(b_unit)

# Expected variance of the DiM estimator from the randomization design
mean(vhat_unit)

## --> the Neyman variance is conservative wrt the true value



#########################################
##### (2) Cluster random assignment #####
#########################################


for(i in 1:nsim){
  # Generate treatment indicator
  dt_cluster$z <- 0
  
  # Randomly assign half of the **CLUSTERS** to the treatment
  dt_cluster$z[dt_cluster$cluster_index %in% sample(unique(dt_cluster$cluster_index), N_clusters/2)] <- 1
  
  # Observed outcome for treated and control obs
  dt_cluster$y <- dt_cluster$z*dt_cluster$y1 + (1-dt_cluster$z)*dt_cluster$y0
  
  # Compute difference in means (w/ clustered SEs)
  fit <- lm_robust(y~z,
                   data=dt_cluster,
                   clusters=cluster_index)
  
  # Store DiM estimate
  b_cluster[i] <- coef(fit)["z"]
  
  # Store estimated DiM variance
  vhat_cluster[i] <- diag(vcov(fit))["z"]
}

## Unbiasedness and variance:
mean(y1 - y0)

mean(b_cluster)

var(b_cluster)

mean(vhat_cluster)

var(b_cluster)/var(b_unit)

## --> Cluster randomized design has 3+ higher variance than the simple randomized one
##     That's because the effective sample is smaller



########################################
##### (3) Block cluster assignment #####
########################################

# N. of blocks
N_blocks <- 10

# N. clusters per block
clusters_per_block <- N_clusters/N_blocks

# Assign clusters to their block
dt_cluster <- dt_cluster[order(dt_cluster$x),]
dt_cluster$block <- (1:N_blocks%x%rep(1, clusters_per_block))%x%rep(1, cluster_size)


for(i in 1:nsim){
  
  ### Preliminaries: treatment assignment ###
  
  # Generate treatment indicator
  dt_cluster$z <- 0
  
  # Within each block:
  for(j in unique(dt_cluster$block)){
    # Assign half of the clusters to treatment
    clus_treat <- sample(unique(dt_cluster[dt_cluster$block==j,"cluster_index"]), 
                         clusters_per_block/2)
    dt_cluster[dt_cluster$block==j&dt_cluster$cluster_index%in% clus_treat,"z"] <- 1
  }
  
  # Observed outcome for each obs
  dt_cluster$y <- dt_cluster$z*dt_cluster$y1 + (1-dt_cluster$z)*dt_cluster$y0
  
  ### Estimation 1: ignore stratification ###
  fit <- lm_robust(y~z,
                   data=dt_cluster,
                   clusters = cluster_index)
  b_block_nostrat[i] <- coef(fit)["z"]
  vhat_block_nostrat[i] <- diag(vcov(fit))["z"]
  
  # Preliminaries to estimation 2: construct inverse probability weights for 
  # the block-cluster randomized experiment
  
  # Means of covariate and treatment indicator by cluster (invariant)
  z_aggregated <- aggregate(dt_cluster[,c("x","z","block")],
                            by=list(cluster=dt_cluster$cluster_index),
                            mean)
  z_aggregated <- z_aggregated[order(z_aggregated$x),]
  
  # Means of cluster-specific treatment indicators within block
  p_treatment <- tapply(z_aggregated$z,
                        z_aggregated$block,
                        mean)
  
  # Create IPW score for each block
  dt_cluster$ipw <- NA
  for(k in unique(dt_cluster$block)){
    dt_cluster[dt_cluster$block==k,"ipw"]  <- dt_cluster[dt_cluster$block==k,
                                                         "z"]*(1/p_treatment[as.numeric(names(p_treatment))==k]) +
      (1-dt_cluster[dt_cluster$block==k,
                    "z"])*(1/(1-p_treatment[as.numeric(names(p_treatment))==k]))
  }
  
  ### Estimation 2: Weighted Block FE ###
  # note that the weights in our simulation are uniform, so not actually needed
  fit <- lm_robust(y~z+as.factor(block),
                   data=dt_cluster,
                   clusters=cluster_index,
                   weights=ipw)
  
  b_block_1[i] <- coef(fit)["z"]
  
  vhat_block_1[i] <- diag(vcov(fit))["z"]
  
  ### Estimation 3: Block FE with centered interaction ###
  fit <- lm_lin(y~z,
                covariates=~as.factor(block),
                clusters=cluster_index,
                data=dt_cluster)
  
  b_block_2[i] <- coef(fit)["z"]
  
  vhat_block_2[i] <- diag(vcov(fit))["z"]
}


## Unbiasedness and variance:

mean(y1 - y0)

mean(b_block_nostrat)

mean(b_block_1)

mean(b_block_2)

var(b_block_nostrat)

var(b_block_1)

var(b_block_2)

mean(vhat_block_nostrat)

mean(vhat_block_1)

mean(vhat_block_2)

var(b_block_2)/var(b_cluster)
## --> block-cluster randomization with interaction specification has 0.5*variance than the simple cluster

var(b_block_2)/var(b_unit)
## --> still, block-cluster randomization has ~2*variance than simple randomization 