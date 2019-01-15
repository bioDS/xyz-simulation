#!/usr/bin/Rscript
require(Matrix)
require(dplyr)
require(xyz)

## Uncomment to use modified xyz
#library(Rcpp)
#source('./xyz/R/regression.R')
#source('./xyz/R/search.R')
#source('./xyz/R/xyz.R')
#sourceCpp('./xyz/src/core.cpp')

verbose <- TRUE

args <- commandArgs(trailingOnly = TRUE)

print(args)
f <- args[1]
regression_alpha <- 0.9

n <- regmatches(x = f, m = regexpr(f, pattern = "(?<=n)\\d+(?=_)", perl = TRUE)) %>% as.numeric
p <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_p)\\d+(?=_)", perl = TRUE)) %>% as.numeric
SNR <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_SNR)\\d+(?=_)", perl = TRUE)) %>% as.numeric
num_bi <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbi)\\d+(?=_)", perl = TRUE)) %>% as.numeric
num_bij <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbij)\\d+(?=_)", perl = TRUE)) %>% as.numeric
perc_viol <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_viol)\\d+(?=_)", perl = TRUE)) %>% as.numeric
ID <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_)\\d+(?=\\.rds)", perl = TRUE)) %>% as.numeric


## not really necessary, but X is neater than data$X
data <- readRDS(paste("simulated_data/", f, sep=''))
X <- data$X
Y <- data$Y
obs <- data$obs
bi_ind <- data$bi_ind
bij_ind <- data$bij_ind

data <- NULL
gc()

## Fit model, wip: use xyz
if (verbose) cat("Fitting model\n")
time <- system.time(regression_results <- xyz_regression(X, Y %>% as.numeric, standardize=TRUE, standardize_response=TRUE, alpha=regression_alpha, L=min(ceiling(sqrt(n)),100)))


if (verbose) cat("Collecting stats\n")

# Collect coefficients
fx_main <- data.frame(gene_i = regression_results[[1]][[10]]) %>%
  arrange(gene_i) %>%
  mutate(type = "main", gene_j = NA, TP = gene_i %in% bi_ind[["gene_i"]]) %>%
  select(gene_i, gene_j, type, TP) %>%
  arrange(desc(TP)) %>%
  tbl_df
fx_main

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


pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
  filter(!is.na(coef)) %>%
  data.frame(., pval = summary(fit_red)$coef[-1,4]) %>%
  tbl_df
  
smry <- left_join(rbind(fx_main, fx_int) %>% data.frame(id = 1:nrow(.), .), pvals, by = "id") %>%
  mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  rename(coef.est = coef) %>%
  left_join(., obs, by = c("gene_i", "gene_j")) 

# Print time taken to actuall run xyz
if (verbose)
    time

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
                       n, p, SNR, num_bi, num_bij, perc_viol, ID))
