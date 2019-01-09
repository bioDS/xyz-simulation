#!/usr/bin/Rscript
require(Matrix)
require(dplyr)
require(xyz)
#library(Rcpp)
#source('./xyz/R/regression.R')
#source('./xyz/R/search.R')
#source('./xyz/R/xyz.R')
#sourceCpp('./xyz/src/core.cpp')
#require(glinternet)

#setwd("~/Projects/epistasis")

verbose <- TRUE
lambda_min_ratio = 5e-2
if (verbose) cat("Reading data set\n")
Q <- readRDS("./Q1_binary.rds")


#args <- commandArgs(trailingOnly = TRUE)
args <- c("1000", "100", "10", "50", "20", "0", "0.9")

print(args)
n <- args[1] %>% as.numeric
p <- args[2] %>% as.numeric
SNR <- args[3] %>% as.numeric
num_bi <- args[4] %>% as.numeric
num_bij <- args[5] %>% as.numeric
perc_viol <- args[6] %>% as.numeric
regression_alpha <- args[7] %>% as.numeric
# pm_freqs <- readRDS(file = "results/simulation/perturbation_matrix/QIAGEN_targetfreqs.rds") 

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

# Interaction terms
if (verbose) cat("Sampling interaction fx\n")
# s <- seq(from = max(0, median(bij_obs) - floor(num_bij/2)) + 1, to = min(median(bij_obs) + floor(num_bij/2), max(bij_obs)))
bij_ind <- obs %>%
  filter(omin %in% sample(setdiff(obs$omin, 0) %>% unique, size = num_bij)) %>%
  group_by(omin) %>%
  sample_n(size = 1) %>%
  ungroup %>%
  rowwise %>%
  mutate(coef = rnorm(1, mean = 0, sd = 2)) # double amplitude sd

## Main effects
if (verbose) cat("Sampling main fx\n")
bi_ind <- c(bij_ind[["gene_i"]], bij_ind[["gene_j"]]) %>% unique %>% sort
# Add violations
num_viols <- perc_viol / 100 * nrow(bij_ind)
bi_ind <- cbind(sample(1:nrow(bij_ind), size = num_viols, replace = FALSE), 
                sample(1:2, size = num_viols, replace = TRUE)) %>%
  as.matrix(bij_ind[,1:2])[.] %>%
  setdiff(bi_ind, .) %>% sort
# Add additional main effects
bi_main_ind <- sample(setdiff(1:p, bi_ind), size = min(p - length(bi_ind), num_bi), replace = FALSE)
bi_ind <- data.frame(gene_i = c(bi_ind, bi_main_ind) %>% unique, gene_j = NA) %>%
  rowwise %>%
  mutate(coef = rnorm(1, mean = 0, sd = 1))

# Fitness
if (verbose) cat("Sampling fitness\n")
Y <- X[,bi_ind[["gene_i"]], drop = FALSE] %*% bi_ind[["coef"]]
for (i in 1:nrow(bij_ind)) {
  Y <- Y + (X[,bij_ind[i,][["gene_i"]], drop = FALSE] * X[,bij_ind[i,][["gene_j"]], drop = FALSE]) %*% bij_ind[i,][["coef"]]
}
## add noise
noise <- (rnorm(n = nrow(Y), mean = 0, sd = 1))
Y <- Y + sqrt(var(Y[,1])/(SNR * var(noise))) * noise

## Fit model, wip: use xyz
if (verbose) cat("Fitting model\n")
regression_results <- xyz_regression(X, Y %>% as.numeric, standardize=TRUE, standardize_response=TRUE, alpha=regression_alpha)

#q()

# let's not worry about anything past here for the moment.

if (verbose) cat("Collecting stats\n")

## Glinternet version:
#cf <- coef(fit, lambdaType = "lambdaHat") #lambdaIndex = 50)#

## Collect coefficients
#fx_main <- data.frame(gene_i = cf$mainEffects$cont,
#                      effect = cf$mainEffectsCoef$cont %>% lapply(., function(x) x[[1]]) %>% unlist) %>%
#  arrange(gene_i) %>%
#  mutate(type = "main", gene_j = NA, TP = gene_i %in% bi_ind[["gene_i"]]) %>%
#  select(gene_i, gene_j, type, effect, TP) %>%
#  arrange(desc(TP)) %>%
#  tbl_df
#fx_main

# Collect coefficients
fx_main <- data.frame(gene_i = regression_results[[1]][[10]]) %>%
                     # effect = cf$mainEffectsCoef$cont %>% lapply(., function(x) x[[1]]) %>% unlist) %>%
  arrange(gene_i) %>%
  mutate(type = "main", gene_j = NA, TP = gene_i %in% bi_ind[["gene_i"]]) %>%
  select(gene_i, gene_j, type, TP) %>%
  arrange(desc(TP)) %>%
  tbl_df
fx_main

# fx_int <- data.frame(gene_i = cf$interactions$catcat[,1], gene_j = cf$interactions$catcat[,2],
#                      effect = cf$interactionsCoef$catcat %>% lapply(., function(x) x[1,1]) %>% unlist) %>%
#fx_int <- data.frame(gene_i = cf$interactions$contcont[,1], gene_j = cf$interactions$contcont[,2],
#                     effect = cf$interactionsCoef$contcont %>% unlist) %>%
#  arrange(gene_i) %>%
#  left_join(., obs, by = c("gene_i", "gene_j")) %>%
#  mutate(type = "interaction") %>%
#  rowwise %>%
#  left_join(., bij_ind, by = c("gene_i", "gene_j")) %>%
#  ungroup %>%
#  mutate(TP = !is.na(coef)) %>%
#  select(gene_i, gene_j, type, effect, TP) %>%
#  arrange(desc(TP)) %>%
#  tbl_df
#fx_int %>% data.frame

# remove reflexive results
first = unlist(split(regression_results[[3]][[10]], 1:2)[1])
second = unlist(split(regression_results[[3]][[10]], 1:2)[2])

nfirst = first[first != second]
nsecond = second[first != second]

fx_int <- data.frame (gene_i = nfirst, 
                      gene_j = nsecond) %>%
  arrange(gene_i) %>%
  left_join(., obs, by = c("gene_i", "gene_j")) %>%
  mutate(type = "interaction") %>%
  rowwise %>%
  left_join(., bij_ind, by = c("gene_i", "gene_j")) %>%
  ungroup %>%
  mutate(TP = !is.na(coef)) %>%
  select(gene_i, gene_j, type, TP) %>%
  arrange(desc(TP)) %>%
  tbl_df
fx_int %>% data.frame

## Statistical test if b_i and b_ij are sig. > 0
Z <- cbind(X[,fx_main[["gene_i"]]])
for (i in 1:nrow(fx_int)) {
  Z <- cbind(Z, X[,fx_int[i,][["gene_i"]], drop = FALSE] * X[,fx_int[i,][["gene_j"]], drop = FALSE])
}
Z <- as.matrix(Z)
colnames(Z) <- rownames(Z) <- NULL
Ynum <- as.numeric(Y)
fit_red <- lm(Ynum ~ Z)

#TODO: probably shouldn't do this, it seems bad.
#fit_red[[1]][is.na(fit_red[[1]])] <- 0


pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
  filter(!is.na(coef)) %>%
  data.frame(., pval = summary(fit_red)$coef[-1,4]) %>%
  tbl_df
smry <- left_join(rbind(fx_main, fx_int) %>% data.frame(id = 1:nrow(.), .), pvals, by = "id") %>%
  mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  rename(coef.est = coef) %>%
  left_join(., obs, by = c("gene_i", "gene_j")) %>%
   mutate(pval = p.adjust(pval, method = "BH")) %>%
   filter(pval < 0.05)

# smry <- summary(fit_red)$coef %>%
#   data.frame %>% 
#   tbl_df %>%
#   .[-1,c(1,4)] %>%
#   tbl_df %>%
#   rename(coef = Estimate, pval = Pr...t..) %>%
#   mutate(pval = p.adjust(pval, method = "BH")) %>%
#   cbind(rbind(fx_main, fx_int), .) %>%
#   filter(pval < 0.01) %>%
#   left_join(., bij_ind, by = c("gene_i", "gene_j"))


## Write out
if (verbose) cat("Saving\n")
saveRDS(list(fit = regression_results,
             bij = bij_ind,
             bi = bi_ind,
             obs = obs,
             fx_int = fx_int,
             fx_main = fx_main,
             fit_red = fit_red,
             smry = smry),
        file = sprintf("./fits_proper/n%d_p%d_SNR%d_nbi%d_nbij%d_viol%d_%d.rds",
                       n, p, SNR, num_bi, num_bij, perc_viol, (runif(1) * 1e5) %>% floor))

# now that we've printed glinternet restults, let's see how xyz compares
#print(xyz_results)

#TODO: make this whole section look a bit more like R
#first <- unlist(split(unlist(xyz_results[1]), 1:2)[1]) %in% unlist(bij_ind[1])
#second <- unlist(split(unlist(xyz_results[1]), 1:2)[2]) %in% unlist(bij_ind[2])

#tp_indices = list()
#tp_check = first & second
#for (i in c(1:length(tp_check))) {
#	if (tp_check[i] == TRUE) {
#		tp_indices <- append(tp_indices, i)
#	}
#}

#found_results <- 0
#for (i in tp_indices) {
#	for (j in which(unlist(bij_ind[1]) %in% unlist(xyz_results[1])[(i-1)*2 + 1])) {
#			if (unlist(bij_ind[2])[j] == unlist(xyz_results[1])[(i-1)*2 + 2]) {
#			cat(unlist(xyz_results[1])[(i-1)*2 + 1],
#				unlist(xyz_results[1])[(i-1)*2 + 2],
#				unlist(xyz_results[2])[i], "\n")
#			found_results <- found_results + 1
#		}
#	}
#}
#cat("found pairs: ", found_results, "\n")

#Xg <- X
#Xg[Xg == 0] <- -1

#regression_results <- xyz_regression(Xg, Y %>% as.numeric, standardize=TRUE, standardize_response=TRUE, alpha=regression_alpha)

#print(regression_results)

#bij_ind

##for (i in c(1:length(regression_results[[3]][[10]]) / 2))
##for (i in regression_results[[3]][[10]]) {
#	#print(i)
##}

#first <- unlist(split(unlist(regression_results[[3]][[10]]), 1:2)[1]) %in% unlist(bij_ind[1])
#first <- first | unlist(split(unlist(regression_results[[3]][[10]]), 1:2)[1]) %in% unlist(bij_ind[2])
#second <- unlist(split(unlist(regression_results[[3]][[10]]), 1:2)[2]) %in% unlist(bij_ind[1]) | unlist(split(unlist(regression_results[[3]][[10]]), 1:2)[2]) %in% unlist(bij_ind[2])

#tp_indices = list()
#reflexive_indices = list()
#reflexive_results <- 0
#tp_check = first & second
#for (i in c(1:length(tp_check))) {
#	if (regression_results[[3]][[10]][[2*i - 1]] == regression_results[[3]][[10]][[2*i]]) {
#		cat("Ignoring result (", regression_results[[3]][[10]][[2*i - 1]], ",",
#					regression_results[[3]][[10]][[2*i]], ")\n")
#		reflexive_results <- reflexive_results + 1
#		reflexive_indices <- append(reflexive_indices, i)
#	} else if (tp_check[i] == TRUE) {
#		tp_indices <- append(tp_indices, i)
#	}
#}

#found_results <- 0
#for (i in tp_indices) {
#	for (j in which(unlist(bij_ind[1]) %in% unlist(regression_results[[3]][[10]])[(i-1)*2 + 1])) {
#			if (unlist(bij_ind[2])[j] == unlist(regression_results[[3]][[10]])[(i-1)*2 + 2]) {
#			cat(unlist(regression_results[[3]][[10]])[(i-1)*2 + 1],
#				unlist(regression_results[[3]][[10]])[(i-1)*2 + 2],
#				unlist(regression_results[[4]][[10]])[i], "\n")
#			found_results <- found_results + 1
#		}
#	}
#}
##TODO: if there are no results this will still return 1. It appears to be correct for >0 results, however.
#interactions_found <- length(regression_results[[3]][[10]]) / 2 - reflexive_results

#cat("found pairs: ", found_results, " out of ", interactions_found, "\n")

#cat("recall: ", found_results/num_bij, " precision: ", found_results/interactions_found, "\n")

#coef <- regression_results[[2]][[10]]
#for (i in reflexive_indices) {
#    coef <- coef[-i]
#}

#if (size(coef) == 0) {
#    coef = list()
#}


### Write out
#if (verbose) cat("Saving\n")
#saveRDS(list(bij = bij_ind,
#             bi = bi_ind,
#             obs = obs,
#             tp = found_results,
#             count = interactions_found,
#             coef.est = coef.est),
#        file = sprintf("fits_proper/n%d_p%d_SNR%d_nbi%d_nbij%d_viol%d_alpha%f_%d.rds",
#                       n, p, SNR, num_bi, num_bij, perc_viol, regression_alpha, (runif(1) * 1e5) %>% floor))
