#!/usr/bin/env Rscript
require(Matrix)
require(dplyr)
require(glinternet)


args <- commandArgs(trailingOnly = TRUE)

verbose <- TRUE

print(args)
f <- args[1]

# ran from 10:40 am to 06:05 am the following day.
time = 72300
fit = readRDS("./whitedwarf_full_gl_output.rds")
X = readRDS("./whitedwarf_full_X.rds")
Y = readRDS("./whitedwarf_full_Y.rds")

if (verbose) cat("Collecting stats\n")
cf <- coef(fit, lambdaType = "lambdaHat") #lambdaIndex = 50)#

## Collect coefficients
fx_main <- data.frame(gene_i = cf$mainEffects$cont,
                      effect = cf$mainEffectsCoef$cont %>% lapply(., function(x) x[[1]]) %>% unlist) %>%
  arrange(gene_i) %>%
  mutate(type = "main", gene_j = NA) %>%
  select(gene_i, gene_j, type) %>%
  tbl_df
fx_main

fx_int <- data.frame(gene_i = cf$interactions$contcont[,1], gene_j = cf$interactions$contcont[,2],
                     effect = cf$interactionsCoef$contcont %>% unlist) %>%
  arrange(gene_i) %>%
  mutate(type = "interaction") %>%
  rowwise %>%
  ungroup %>%
  select(gene_i, gene_j, type) %>%
  tbl_df

fx_int %>% data.frame

## Fit main effects only for comparison
fit_main_only = lm(Y ~ X)

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
  rename(coef.est = coef)

# Print time taken to actuall run xyz
if (verbose)
    time

## Write out
if (verbose) cat("Saving\n")
saveRDS(list(fit = fit,
             fit_main_only = fit_main_only,
             fx_int = fx_int,
             fx_main = fx_main,
             fit_red = fit_red,
             smry = smry),
        file = sprintf("glinternet_ifx_full_stats.rds"))
