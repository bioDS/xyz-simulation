#!/usr/bin/Rscript
require(Matrix)
require(dplyr)

verbose <- TRUE
lambda_min_ratio = 5e-2
if (verbose) cat("Reading data set\n")
Q <- readRDS("./Q1_binary.rds") # possibly not large enough for p > 1000


args <- commandArgs(trailingOnly = TRUE)
#args <- c("1000", "100", "10", "50", "20", "0")

print(args)
n <- args[1] %>% as.numeric
p <- args[2] %>% as.numeric
SNR <- args[3] %>% as.numeric
num_bi <- args[4] %>% as.numeric
num_bij <- args[5] %>% as.numeric
perc_viol <- args[6] %>% as.numeric

#' Perturbation matrix
#'
#' Sample a n \times p perturbation matrix from a frequency vector
#' using the Bernoulli distribution
#'
#' @author Fabian Schmich
#' @param n Number of observations
#' @param p Number of genes
#' @return Perturbation matrix
perturbation_matrix <- function(n, p) {
  colnames(Q) <- rownames(Q) <- NULL
  X <- Q[sample(1:nrow(Q), n), sample(1:ncol(Q), p)] %>% as.matrix
  dup <- which(duplicated(t(X)) == TRUE)
  while(length(dup) > 0) {
    repl <- Q[sample(1:nrow(Q), n), sample(1:ncol(Q), length(dup)), drop = FALSE] %>% as.matrix
    for (i in 1:ncol(repl)) {
      X[,dup[i]] <- repl[,i]
    }
    dup <- which(duplicated(t(X)) == TRUE)
  }
  # X <- matrix(0, nrow = n, ncol = p)
  # for (i in 1:nrow(X)) {
  #   X[i,] <- rbinom(n = p, size = 1, prob = sample(x = freqs, size = 1))
  # }
  return(X)
}

## Simulate
# Perturbation matrix
# o00: the number of times gene_i and gene_j are both false
# o10: number of times gene_i is false, gene_j is true, ....
# o11: both i and j are true.
if (verbose) cat("Make perturbation matrix\n")
X <- perturbation_matrix(n = n, p = p)
obs <- lapply(2:p, function(i) {
  apply(X[,1:(i-1), drop = FALSE], 2, function(col) {
    data.frame(o00 = (!col & !X[,i]) %>% sum,
               o10 = (col & !X[,i]) %>% sum,
               o01 = (!col & X[,i]) %>% sum,
               o11 = (col & X[,i]) %>% sum)
    }) %>% do.call("rbind", .) %>%
    tbl_df %>%
    mutate(gene_j = i, gene_i = 1:(i-1))
}) %>% do.call("rbind", .) %>%
  tbl_df %>%
  arrange(gene_i, gene_j) %>%
  select(gene_i, gene_j, o00, o01, o10, o11) %>%
  group_by(gene_i, gene_j, o00, o01, o10, o11) %>% 
  summarise(omin = min(o00, o10, o01, o11)) %>%
  ungroup

## Interaction terms
# choose num_bij distinct pairs of genes that will be used for interactions.
# choose their strengths at random.
if (verbose) cat("Sampling interaction fx\n")
# s <- seq(from = max(0, median(bij_obs) - floor(num_bij/2)) + 1, to = min(median(bij_obs) + floor(num_bij/2), max(bij_obs)))
bij_ind <- obs %>%
  filter(omin %in% sample(setdiff(obs$omin, 0) %>% unique, size = num_bij)) %>%
  group_by(omin) %>%
  sample_n(size = 1) %>%
  ungroup %>%
  rowwise %>%
  mutate(coef = -100)
#  mutate(coef = rnorm(1, mean = 0, sd = 2)) # double amplitude sd

## Main effects
# both sides of the interaction effects. If we use all of these, it will impose a strong hierarchy (since every main effect will be present for the interaction effects).
# To avoid this, we will remove some main effects, leaving the interactions in place.
if (verbose) cat("Sampling main fx\n")
#bi_ind <- c(bij_ind[["gene_i"]], bij_ind[["gene_j"]]) %>% unique %>% sort

### Add violations
## main effects that violate (what?)
#num_viols <- perc_viol / 100 * nrow(bij_ind)
#bi_ind <- cbind(sample(1:nrow(bij_ind), size = num_viols, replace = FALSE), 
#                sample(1:2, size = num_viols, replace = TRUE)) %>%
#  as.matrix(bij_ind[,1:2])[.] %>%
#  setdiff(bi_ind, .) %>% sort
## Add additional main effects
#bi_main_ind <- sample(setdiff(1:p, bi_ind), size = min(p - length(bi_ind), num_bi), replace = FALSE)
#bi_ind <- data.frame(gene_i = c(bi_ind, bi_main_ind) %>% unique, gene_j = NA) %>%
#  rowwise %>%
#  mutate(coef = rnorm(1, mean = 0, sd = 1))
#  
## remove some main effects at random, to avoid strong hierarchy.
#bi_ind <- sample_n(bi_ind, round(runif(1, 0, nrow(bij_ind))), replace=FALSE)

## Fitness
# as a result of only the main effects, and their interactions.
# if we want to add some explicit synthetic lethals, we need to set up bij_ind so that Y[somewhere] == 0.
# in the interests of maintaining strictly positive fitness values, let's take e^(fitness) later. near-zero fitness can be obtained with a large enough negative coefficient.
if (verbose) cat("Sampling fitness\n")
#Y <- X[,bi_ind[["gene_i"]], drop = FALSE] %*% bi_ind[["coef"]]
for (i in 1:nrow(bij_ind)) {
  Y <- 1 + (X[,bij_ind[i,][["gene_i"]], drop = FALSE] * X[,bij_ind[i,][["gene_j"]], drop = FALSE]) %*% bij_ind[i,][["coef"]]
}

## add noise
# If fitness is the average number of ofspring per generation, then a mean of 0 would suggest the vast majority of things do not survive many generations. Surely this should be 1?
noise <- (rnorm(n = nrow(Y), mean = 0, sd = 1))
Y <- Y + sqrt(var(Y[,1])/(SNR * var(noise))) * noise

# If interaction coefficients are an indication of the number of the generations offspring the combination would kill, then a more negative result would be, at most, dead faster.
#Y[Y<0] <- 0

# make fitness values positive (e.g. representing remaining living fraction)
#TODO: find a better way of doing this, it wreaks havok on estimated coefficients.
#Y <- exp(Y)

## Write out
if (verbose) cat("Saving\n")
saveRDS(list(X=X, Y=Y,
             bij_ind = bij_ind,
#             bi_ind = bi_ind,
             obs = obs),
        file = sprintf("./simulated_lethal_data/n%d_p%d_SNR%d_nbi%d_nbij%d_viol%d_%d_lethal.rds",
                       n, p, SNR, num_bi, num_bij, perc_viol, (runif(1) * 1e5) %>% floor))
