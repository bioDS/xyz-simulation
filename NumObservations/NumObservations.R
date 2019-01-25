#!/usr/bin/Rscript
require(ggplot2)
require(RColorBrewer)
require(dplyr)
require(tidyr)
require(reshape2)
#source("~/Projects/R/fs_.R")
setwd("..")

args <- commandArgs(trailingOnly = TRUE)

# Number of observations
 ans <- lapply(list.files(path = "./fits_proper/", pattern = "", full.names = TRUE), function(f) {#, spriadditifffntf("n%d_p%d", n, p)), function(f) {
   ans <- readRDS(f)
   n <- regmatches(x = f, m = regexpr(f, pattern = "(?<=n)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   p <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_p)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   SNR <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_SNR)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   nbi <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbi)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   nbij <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbij)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   nlethals <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nlethals)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   perc_viol <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_viol)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   L <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_L)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   id <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_)\\d+(?=\\.rds)", perl = TRUE)) %>% as.numeric
   smry_int <- ans$smry %>% filter(type == "interaction")
   notest <- data.frame(n = n, p = p, SNR = SNR, nbi = nbi, nbij = nbij, L=L, id = id, nlethals = nlethals,
              left_join(ans$bij, smry_int, by = c("gene_i", "gene_j", "o00", "o01", "o10", "o11", "omin")) %>%
                select(observations = omin, TP) %>%
                mutate(type = ifelse(is.na(TP), "FN", "TP")) %>%
                select(-TP) %>%
                rbind(., filter(smry_int, TP == FALSE) %>% select(observations = omin) %>% mutate(type = "FP")),
              test = "no")
   smry_int <- mutate(smry_int, pval = p.adjust(pval, method = "BH")) %>%
     filter(pval < 0.05)
   test <- data.frame(n = n, p = p, SNR = SNR, nbi = nbi, nbij = nbij, L=L, id = id, nlethals = nlethals,
                        left_join(ans$bij, smry_int, by = c("gene_i", "gene_j", "o00", "o01", "o10", "o11", "omin")) %>%
                          select(observations = omin, TP) %>%
                          mutate(type = ifelse(is.na(TP), "FN", "TP")) %>%
                          select(-TP) %>%
                        rbind(., filter(smry_int, TP == FALSE) %>% select(observations = omin) %>% mutate(type = "FP")),
                        test = "yes")
   rbind(test, notest)
 }) %>% do.call("rbind", .) %>%
   tbl_df %>%
   mutate(n = factor(n),
          p = factor(p),
          SNR = factor(SNR),
          nbi = factor(nbi),
          nbij = factor(nbij),
          nlethals = factor(nlethals),
          L = factor(L),
          type = factor(type))
 saveRDS(ans, file = "NumObservations/dat_numobs.rds")


mult <- args[1] %>% as.numeric

for (numrows in c(1000 * mult)) { #1000, 
  if (numrows == 400) {
    rseq <- 0#c(seq(0, 60, by = 20), 100)
  } else if (numrows == 1000) {
    rseq <- rseq <- c(0, 10, 20, 40, 80, Inf)#c(seq(0, 150, by = 50), 250)
  } else {
    rseq <- rseq <- c(0, 10, 20, 40, 80, Inf) #TODO: what should these be?
  }
  for (t in c("yes", "no")) {
  dat_nobs <- readRDS("NumObservations/dat_numobs.rds") %>% filter(nlethals == 0)
    TP=length(unlist((dat_nobs %>% filter(type=="TP"))[1]))
    FP=length(unlist((dat_nobs %>% filter(type=="FP"))[1]))
    FN=length(unlist((dat_nobs %>% filter(type=="FN"))[1]))
    dat_nobs <- dat_nobs %>%
      filter(n == numrows) %>%
      filter(test == t) %>%
      filter(L == round(sqrt(levels(p) %>% as.numeric))) %>%
      filter(nbi == 20*mult, SNR != 1) %>%
      mutate(SNR = factor(SNR, labels = paste0("SNR = ", levels(factor(SNR))))) %>%
      filter(nbij %in% (c(5*mult, 20*mult, 50*mult, 100*mult))) %>%
      group_by(n, p, SNR, nbi, nbij, observations, type) %>%
      summarise(count = n()) %>% 
      spread(type, count, fill = 0) %>%
      mutate(precision = TP / (TP + FP),
             recall = TP / (TP + FN)) %>%
      mutate(F1 = 2 *  (precision * recall) / (precision + recall)) %>%
      filter(!is.na(precision), !is.na(recall), !is.na(F1)) %>%
      gather(measure, value, precision, recall, F1) %>%
      mutate(measure = factor(measure, levels = c("precision", "recall", "F1"), labels = c("Precision", "Recall", "F1"))) %>%
      # filter(measure == "recall") %>%
      mutate(range = cut(observations, rseq)) %>%
      group_by(n, p, SNR, nbij, measure, range) %>%
      summarise(mean = mean(value, na.rm = TRUE), sem = sd(value, na.rm = TRUE) / sqrt(n()))
    
    pl <- ggplot(dat_nobs, aes(x = range,#observations %>% as.character %>% as.numeric,
                               y = mean, 
                               group = nbij,
                               ymax = mean + sem, 
                               ymin = mean - sem)) +
      geom_line(aes(colour = nbij), position = position_dodge(.35)) +
      geom_point(aes(colour = nbij), position = position_dodge(.35), size = 1) +
      geom_errorbar(colour = "darkgrey", width = 0.3, position = position_dodge(.35)) +
      facet_grid(measure~SNR) +
      scale_color_discrete(name = "True interactions") +
      scale_x_discrete(labels = c("(0, 10]", "(10, 20]", "(20, 40]", "(40, 80]", expression(paste("(80, ", infinity, ")")))) +
      ylim(c(0,1)) +
      xlab("Observations of double knockdown") +
      ylab("") +
#      theme_fs() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    pl
    ggsave(pl, file = sprintf("NumObservations/NumObservations_n%d_t%s.pdf", numrows, t), width = 5, height = 7)
  }
}




#
#
#
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


n <- 1000
p <- 100
Q <- readRDS("Q1_binary.rds")

## Simulate
# Perturbation matrix
dat <- lapply(1:10, function(i) {
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
    mutate(n = n, rep = i)
}) %>% do.call("rbind", .) %>%
  tbl_df %>%
  mutate(n = factor(n), rep = factor(rep))



# rseq <- c(seq(0, 150, by = 50), Inf)
rseq <- c(0, 10, 20, 40, 80, Inf)

ans <- filter(dat, !is.na(value)) %>% #, value > 0
  select(-n) %>%
  mutate(range = cut(value, rseq)) %>%
  group_by(range, rep) %>%
  summarise(count = n()) %>%
  group_by(rep) %>%
  mutate(perc = count / sum(count) * 100) %>%
  group_by(range) %>%
  summarise(mean = mean(perc), sem = sd(perc) / sqrt(n())) %>%
  filter(!is.na(range))

pl <- ggplot(ans, aes(x = range)) +
  geom_bar(aes(y = mean), stat = "identity", fill = "lightgrey") +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), stat = "identity", width = 0.3) +
  xlab("Observations of double knockdown") +
  ylab("Gene pairs [%]") +
  scale_x_discrete(labels = c("(0, 10]", "(10, 20]", "(20, 40]", "(40, 80]", expression(paste("(80, ", infinity, ")")))) +
#  theme_fs() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(pl, file = sprintf("NumObservations/NumObservations_n%d_percGenes.pdf", numrows), width = 5, height = 3)
