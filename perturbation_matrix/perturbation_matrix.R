#!/usr/bin/Rscript
# Libraries
require(reshape2)
require(Matrix)
require(dplyr)
require(ggplot2)
require(gridExtra)


## General settings
setwd("..")
#source("~/Projects/R/fs_.R")

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



p <- 100
Q <- readRDS("Q1_binary.rds")

## Simulate
# Perturbation matrix
dat <- lapply(c(400, 1000), function(n) {
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
  
  M <- matrix(0, p, p)
  M[upper.tri(M)] <- obs$omin
  M <- M + t(M)
  row.order <- hclust(dist(M))$order
  col.order <- hclust(dist(t(M)))$order
  M <- M[row.order, col.order]
  M[lower.tri(M)] <- diag(M) <- NA
  # M[M == 0] <- NA
  melt(M) %>% 
    data.frame %>%
    tbl_df %>%
    mutate(n = n)
}) %>% do.call("rbind", .) %>%
  tbl_df %>%
  mutate(n = factor(n)) %>% 
  mutate(n = factor(n, levels = levels(n), labels = paste0("n = ", levels(n))))

for (numrows in c(400, 1000)) {
  dat.n <- filter(dat, n == sprintf("n = %d", numrows))
  
  ident_perc <- dat.n %>% filter(!is.na(value)) %>% 
    mutate(Identifiable = (value != 0)) %>%
    group_by(Identifiable) %>% 
    summarise(count = n()) %>% 
    mutate(freq = count / sum(count)) %>%
    filter(Identifiable == TRUE) %>% .[["freq"]]
  
  
  pl.hist <- ggplot(dat.n, aes(x = value)) +
    geom_histogram(aes(y = ..count../sum(..count..)), bins = 25) +
    facet_wrap(~n) +
    xlab("Observations of double knockdown") +
    ylab("Frequency") +
    ggtitle("A") +
#    theme_fs() +
	theme_bw() +
    theme(plot.title = element_text(hjust = -0.1, vjust = -0.2, size = 18))
  pl.hist
  
  pl.heatmap <- ggplot(dat.n, aes(x = Var1, y = Var2)) +
    geom_raster(aes(fill = value)) +
    scale_fill_distiller(name = "Observations\nof double knockdown", palette = "YlOrRd", trans = "log10", direction = 1, na.value = "white") +
    scale_x_continuous(breaks = c(1, seq(20, p, by = 20), p)) +
    scale_y_continuous(breaks = c(1, seq(20, p, by = 20), p)) +
    xlab("Gene") +
    ylab("Gene") +
    ggtitle("B") +
    facet_wrap(~n) +
    annotate(geom = "text", x = Inf, y = -Inf, 
             label = sprintf("%.0f%% non-zero", ident_perc * 100), 
             hjust = 1.2, vjust = -1.5, colour = "black") +
#    theme_fs() +
	theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = -0.1, vjust = -0.2, size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  pl.heatmap
  
  
  pdf(sprintf("perturbation_matrix/pmatrix_n%d.pdf", numrows), width = 11, height = 6.5)
  grid.arrange(pl.hist, pl.heatmap, ncol = 2)
  dev.off()
}










## Create off-target frequencies from Qiagen libraries
# Q <- readRDS("data/QIAGEN_HELA.rds")
# nzs <- lapply( split(1:nrow(Q), 1:10 %>% as.factor), function(rows) {
#   apply(Q[rows,], 1, function(x) sum(x != 0))
# })
# freqs <- nzs %>% unlist %>% sapply(function(x) x / ncol(Q)) %>% as.numeric
# saveRDS(freqs, file = "results/simulation/perturbation_matrix/QIAGEN_targetfreqs.rds")


# Number of observations per pairwise interaction term \beta_{ij}
# freqs <- readRDS(file = "results/simulation/perturbation_matrix/QIAGEN_targetfreqs.rds")
# n <- 500
# p <- 100

# ANS <- lapply(1:10, function(k) {
#   X <- perturbation_matrix(n = n, p = p)
#   num.obs <- lapply(1:p, function(i) {
#     if (i > 1) {
#       apply(X[,1:(i-1), drop = FALSE], 2, function(x) (x * X[,i]) %>% sum)
#     }
#   }) %>% unlist
#     stopifnot(length(num.obs) == choose(p, 2))
#     print(k)
#     return(num.obs)
# })
# # saveRDS(ANS, file = "results/simulation/perturbation_matrix/obs_per_interaction.rds")
# 
# 
# ## Create figure
# freqs <- readRDS(file = "results/simulation/perturbation_matrix/QIAGEN_targetfreqs.rds")
# pl.freqs <- data.frame(Value = freqs * 100) %>%
#   ggplot(aes(x = Value)) +
#   geom_histogram(aes(y = ..density.. * 100), binwidth = 1, fill = "darkgrey") +
#   xlab("Off-targeted genes [%]") +
#   ylab("Observations [%]") +
#   ggtitle("A") +
#   theme_fs() +
#   theme(plot.title = element_text(hjust = -0.2, vjust = -0.2, size = 18),
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank(),
#       panel.grid.minor.y = element_blank())
# 
# ANS <- readRDS(file = "results/simulation/perturbation_matrix/obs_per_interaction.rds")
# pl.bij <- lapply(1:length(ANS), function(i) {
#   data.frame(rep = i, obs = ANS[[i]])
# }) %>% do.call("rbind", .) %>%
#   tbl_df %>%
#   mutate(rep = factor(rep), obs = factor(obs)) %>%
#   group_by(rep, obs) %>%
#   summarise(perc = n() / choose(p, 2) * 100) %>%
#   group_by(obs) %>%
#   summarise(mean = mean(perc), sd = sd(perc)) %>%
#   ggplot(aes(x = obs, y = mean)) +
#   geom_bar(stat = "identity", fill = "darkgrey") + ##99d8c9
#   geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd), width = 0.5) +
#   scale_x_discrete(breaks = seq(0, 305, by = 5), labels = seq(0, 305, by = 5)) +
#   xlab("Double knock-downs [#]") +
#   ylab("Gene pairs [%]") +
#   ggtitle("B") +
#   theme_fs() +
#   theme(plot.title = element_text(hjust = -0.2, vjust = -0.2, size = 18),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_blank())
# 
# pdf("results/simulation/perturbation_matrix/perturbation_matrix.pdf", width = 7, height = 3)
# grid.arrange(pl.freqs, pl.bij, nrow = 1)
# dev.off()
