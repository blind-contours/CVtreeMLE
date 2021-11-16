library(dplyr)
library(tidyr)
library(fastDummies)

#################
### Backwards Histogram Regression Simulation for Mixtures
################

# number observations to create
## create multivariate normal distribution




calc_truth_rev_his_m3 <- function(n_obs,
                                     splits,
                                     mins,
                                     maxs,
                                     w1_betas,
                                     w2_betas,
                                     mix_subspace_betas,
                                     subspace_assoc_strength_betas,
                                     marginal_impact_betas,
                                     eps_sd,
                                     binary) {

  ## create a covariate
  W <- rnorm(n_obs, 0, 0.5)
  W2 <- rnorm(n_obs, 0, 0.5)

  ## probabilities
  b0i <- mix_subspace_betas
  b1i <- w1_betas
  b2i <- w2_betas

  probs_list <- c()

  for (i in seq(W)){

    w <- W[i]
    w2 <- W2[i]

    denominator <- sum(exp(b0i[2] + (b1i[2]*w) + (b2i[2]*w2)),
                       exp(b0i[3] + (b1i[3]*w) + (b2i[3]*w2)),
                       exp(b0i[4] + (b1i[4]*w) + (b2i[4]*w2)),
                       exp(b0i[5] + (b1i[5]*w) + (b2i[5]*w2)),
                       exp(b0i[6] + (b1i[6]*w) + (b2i[6]*w2)),
                       exp(b0i[7] + (b1i[7]*w) + (b2i[7]*w2)),
                       exp(b0i[8] + (b1i[8]*w) + (b2i[8]*w2)))

    a0 <- exp(b0i[1] + (b1i[1] * w) + (b2i[1] * w2)) / (1 + denominator)
    a1 <- exp(b0i[2] + (b1i[2] * w) + (b2i[2] * w2)) / (1 + denominator)
    a2 <- exp(b0i[3] + (b1i[3] * w) + (b2i[3] * w2)) / (1 + denominator)
    a3 <- exp(b0i[4] + (b1i[4] * w) + (b2i[4] * w2)) / (1 + denominator)
    a4 <- exp(b0i[5] + (b1i[5] * w) + (b2i[5] * w2)) / (1 + denominator)
    a5 <- exp(b0i[6] + (b1i[6] * w) + (b2i[6] * w2)) / (1 + denominator)
    a6 <- exp(b0i[7] + (b1i[7] * w) + (b2i[7] * w2)) / (1 + denominator)
    a7 <- 1 - sum(a0,a1,a2,a3,a4,a5,a6)

    probs_a <- c(a0, a1, a2, a3, a4, a5, a6, a7)

    probs_list[[i]] <- probs_a

  }

    probs_df <- as.data.frame(do.call(rbind, probs_list))
  colnames(probs_df) <- paste0("p", seq_len(ncol(probs_df)))

  probs <- as.matrix(probs_df %>% dplyr::select(p1:p8))

  res <- probs_df %>%
    dplyr::mutate(rcat = Hmisc::rMultinom(probs, 1))

  res <- as.data.frame(res)

  mixture_section_indicator <- expand.grid(c(0,1), c(0,1), c(0,1))
  colnames(mixture_section_indicator) <- c("M1", "M2", "M3")

  Ms <- as.data.frame(matrix(data = NA, ncol = 3, nrow = n_obs))
  colnames(Ms) <- c("M1", "M2", "M3")

  for (i in 1:nrow(mixture_section_indicator)) {
    ## iteration through the subspaces
    mix_space <- mixture_section_indicator[i,]

    ## 0 or 1 for for each mixture
    M1_01 <- mix_space[1]
    M2_01 <- mix_space[2]
    M3_01 <- mix_space[3]

    ## set high and low for M1
    if (M1_01 == 0) {
      m1_min <- mins[1]
      m1_max <- splits[1]
    } else{
      m1_min <- splits[1]
      m1_max <- maxs[1]
    }

    ## set high and low for M2
    if (M2_01 == 0) {
      m2_min <- mins[2]
      m2_max <- splits[2]
    } else{
      m2_min <- splits[2]
      m2_max <- maxs[2]
    }

    ## set high and low for M3
    if (M3_01 == 0) {
      m3_min <- mins[3]
      m3_max <- splits[3]
    } else{
      m3_min <- splits[3]
      m3_max <- maxs[3]
    }

    M1_sec <- runif(n_obs, min = m1_min, max = m1_max)
    M2_sec <- runif(n_obs, min = m2_min, max = m2_max)
    M3_sec <- runif(n_obs, min = m3_min, max = m3_max)

    subspace_data <- cbind(M1_sec, M2_sec, M3_sec)

    Ms[res$rcat == paste("p", i, sep = ""),] <- subspace_data[res$rcat == paste("p", i, sep = ""),]
  }

  ## marginals

  m1_marg <- ifelse(Ms$M1 > splits[1], 1,0)
  m2_marg <- ifelse(Ms$M2 > splits[2], 1,0)
  m3_marg <- ifelse(Ms$M3 > splits[3], 1,0)


  if(binary == TRUE) {
    y <-
      plogis(subspace_assoc_strength_betas[1] +
               subspace_assoc_strength_betas[2] * as.numeric(res$rcat == "p2") +
               subspace_assoc_strength_betas[3] * as.numeric(res$rcat == "p3") +
               subspace_assoc_strength_betas[4] * as.numeric(res$rcat == "p4") +
               subspace_assoc_strength_betas[5] * as.numeric(res$rcat == "p5") +
               subspace_assoc_strength_betas[6] * as.numeric(res$rcat == "p6") +
               subspace_assoc_strength_betas[7] * as.numeric(res$rcat == "p7") +
               subspace_assoc_strength_betas[8] * as.numeric(res$rcat == "p8") +
               marginal_impact_betas[1] * m1_marg +
               marginal_impact_betas[2] * m2_marg +
               marginal_impact_betas[3] * m2_marg)
               #W +
               #W2 +
               #rnorm(length(res$rcat), mean = 0, sd = eps_sd))

    y <- ifelse(y > 0.50, 1,0)

  } else{

    ## mix space counterfactual setting
    dummys <- fastDummies::dummy_cols(res$rcat)
    dummys <- dummys[,-1]
    internvention_vars <- which(subspace_assoc_strength_betas != 0 )

    for (i in internvention_vars) {
      dummys[,i] <- 1
    }

    ## marginal counterfactuals
    marginals <- cbind(m1_marg, m2_marg, m3_marg)
    marg_int_vars <- which(marginal_impact_betas != 0 )

    for (i in marg_int_vars) {
      marginals[,i] <- 1
    }

    y1 <-
      ## section for mixture subspaces
      subspace_assoc_strength_betas[1] +
      subspace_assoc_strength_betas[2] * dummys[,2] +
      subspace_assoc_strength_betas[3] * dummys[,3] +
      subspace_assoc_strength_betas[4] * dummys[,4] +
      subspace_assoc_strength_betas[5] * dummys[,5] +
      subspace_assoc_strength_betas[6] * dummys[,6] +
      subspace_assoc_strength_betas[7] * dummys[,7] +
      subspace_assoc_strength_betas[8] * dummys[,8] +
      marginal_impact_betas[1] * marginals[,1] +
      marginal_impact_betas[2] * marginals[,2] +
      marginal_impact_betas[3] * marginals[,3]
      #W +
      #W2 +
      #rnorm(length(res$rcat), mean = 0, sd = eps_sd)

    ## setting to 0 for relevant sections of the mixture space and marginal space
    for (i in internvention_vars) {
      dummys[,i] <- 0
    }

    for (i in marg_int_vars) {
      marginals[,i] <- 0
    }

    y0 <-
      ## section for mixture subspaces
      subspace_assoc_strength_betas[1] +
      subspace_assoc_strength_betas[2] * dummys[,2] +
      subspace_assoc_strength_betas[3] * dummys[,3] +
      subspace_assoc_strength_betas[4] * dummys[,4] +
      subspace_assoc_strength_betas[5] * dummys[,5] +
      subspace_assoc_strength_betas[6] * dummys[,6] +
      subspace_assoc_strength_betas[7] * dummys[,7] +
      subspace_assoc_strength_betas[8] * dummys[,8] +
      marginal_impact_betas[1] * marginals[,1] +
      marginal_impact_betas[2] * marginals[,2] +
      marginal_impact_betas[3] * marginals[,3]
      #W +
      #W2 +
      #rnorm(length(res$rcat), mean = 0, sd = eps_sd)

    counterfactual_y <- y1 - y0

    #counterf_y_scaled <- ((counterfactual_y - min(counterfactual_y)) / (max(counterfactual_y) - min(counterfactual_y)))
    true_ATE <- mean(counterfactual_y)

  }


  return(true_ATE)
}
