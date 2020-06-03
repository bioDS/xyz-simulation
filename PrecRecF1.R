#!/usr/bin/Rscript
require(ggplot2)
require(dplyr)
require(RColorBrewer)
#source("~/Projects/R/fs_.R")
#setwd("..")

args <- commandArgs(trailingOnly = TRUE)
target_alpha <- args[1] %>% as.numeric

cat("target alpha", target_alpha, "\n")

# Precision, recall and F1 for interaction terms
ans <- lapply(list.files(path = "./fits_proper/", pattern = "", full.names = TRUE), function(f) {
  ans <- readRDS(f)
  n <- regmatches(x = f, m = regexpr(f, pattern = "(?<=n)\\d+(?=_)", perl = TRUE)) %>% as.numeric
  p <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_p)\\d+(?=_)", perl = TRUE)) %>% as.numeric
  SNR <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_SNR)\\d+(?=_)", perl = TRUE)) %>% as.numeric
  nbi <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbi)\\d+(?=_)", perl = TRUE)) %>% as.numeric
  nbij <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbij)\\d+(?=_)", perl = TRUE)) %>% as.numeric
  perc_viol <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_viol)\\d+(?=_)", perl = TRUE)) %>% as.numeric
  regression_alpha <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_alpha)[\\d.]+(?=_)", perl = TRUE)) %>% as.numeric
  ID <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_)\\d+(?=\\.rds)", perl = TRUE)) %>% as.numeric
  tp <- ans$tp
  count <- ans$count
  notest <- data.frame(n = n, p = p, SNR = SNR, regression_alpha, nbi = nbi, nbij = nbij,
                       precision = tp / count,
                       recall = tp / nbij) %>%
    mutate(F1 = 2 *  (precision * recall) / (precision + recall),
           test = "no")
  test <- data.frame(n = n, p = p, SNR = SNR, regression_alpha, nbi = nbi, nbij = nbij,
                     precision = tp / count,
                     recall = tp / nbij) %>%
    mutate(F1 = 2 *  (precision * recall) / (precision + recall),
           test = "yes")
  rbind(notest, test) %>% data.frame(., id = ID)
}) %>% do.call("rbind", .) %>%
  tbl_df %>%
  mutate(n = factor(n),
         p = factor(p),
         SNR = factor(SNR),
         regression_alpha = factor(regression_alpha),
         nbi = factor(nbi),
         nbij = factor(nbij))
saveRDS(ans, file = "PrecRecF1/dat_precrecf1.rds")

# or we could calculate dat_precrecf1.rds beforehand.

for (numrows in c(1000)) { #400
  for (t in c("yes", "no")) {
    dat_precrecf1 <- readRDS(file = "PrecRecF1/dat_precrecf1.rds") %>%
      filter(n == numrows) %>%
      filter(test == t) %>%
      filter(nbi %in% c("0", "20", "50", "100")) %>%
      filter(nbij %in% c("5", "20", "50", "100")) %>%
      filter(SNR != "1") %>%
      filter(regression_alpha == target_alpha) %>%
      mutate(SNR = factor(SNR, levels = levels(SNR), labels = paste0("SNR = ", levels(SNR)))) %>%
      group_by(n, p, SNR, nbi, nbij, test) #%>% sample_n(10, replace = TRUE)
      #TODO: replace = FALSE when we have a large enough data set.
    
    pl <- group_by(dat_precrecf1, SNR, nbi, nbij) %>%
      summarise(mean = mean(precision, na.rm = TRUE), sd = sd(precision, na.rm = TRUE) / sqrt(n())) %>%
      mutate(Measure = "Precision") %>%
      rbind(., group_by(dat_precrecf1, SNR, nbi, nbij) %>%
              summarise(mean = mean(recall, na.rm = TRUE), sd = sd(recall, na.rm = TRUE) / sqrt(n())) %>%
              mutate(Measure = "Recall")) %>%
      rbind(., group_by(dat_precrecf1, SNR, nbi, nbij) %>%
              summarise(mean = mean(F1, na.rm = TRUE), sd = sd(F1, na.rm = TRUE) / sqrt(n())) %>%
              mutate(Measure = "F1")) %>%
      mutate(Measure = factor(Measure, levels = c("Precision", "Recall", "F1"))) %>%
      ggplot(aes(x = nbij, 
                 y = mean, 
                 group = nbi,
                 ymax = mean + sd, 
                 ymin = mean - sd)) +
      facet_grid(Measure ~ SNR) +
      # geom_bar(aes(fill = nbi), stat = "identity", position = "dodge") +
      geom_line(aes(colour = nbi), position=position_dodge(.35)) +
      geom_point(aes(colour = nbi), position=position_dodge(.35), size = 1) +
      geom_errorbar(colour = "darkgrey", width=.3, position=position_dodge(.35)) +
      scale_color_discrete(name = "True additional\nmain effects") +
      # scale_colour_manual(name = "Additional main effects", values = c("#7fcdbb", "#1d91c0", "#253494")) + #brewer.pal(7, "Set3")[-1]) + #
      ylim(c(0,1)) +
      xlab("True interactions") +
      ylab("") +
#      theme_fs() +
      theme(legend.position = "bottom")
    pl
    ggsave(pl, file = sprintf("PrecRecF1/PrecRecF1_n%d_t%s_a%.2f.pdf", numrows, t, target_alpha), width = 3, height = 4)
  }
}








#
#
#
require(tidyr)
require(gridExtra)
#source("~/Projects/R/fs_.R")
for (numrows in c(400, 1000)) {
  
  dat <- readRDS(file = "PrecRecF1/dat_precrecf1.rds") %>%
    filter(n == numrows) %>%
    mutate(test = factor(test)) %>%
    group_by(n, p, SNR, nbi, nbij, test) %>%
    summarise(precision = mean(precision, na.rm = TRUE), recall = mean(recall, na.rm = TRUE)) %>%
    gather(measure, value, precision:recall) %>%
    ungroup %>%
    spread(test, value) %>%
    filter(SNR == 5) %>%
    filter(nbi %in% c(0, 20, 50, 100)) %>%
    filter(nbij %in% c("5", "20", "50", "100"))
  
  
  pl.prec <- group_by(subset(dat, measure == "precision"), n, p, nbi, nbij, measure) %>%
    summarise(change = yes - no) %>%
    ggplot(aes(x = nbij, y = change)) +
    geom_bar(aes(fill = nbi), stat = "identity", position = "dodge") +
    # facet_wrap(~measure, scale = "free_y", ncol = 1) +
    # scale_fill_manual(name = "Additional\nmain effects", values = c("#7fcdbb", "#1d91c0", "#253494")) +
    scale_fill_discrete(name = "True additional\nmain effects") +
    ylim(c(0, 0.8)) +
    xlab("True interactions") +
    ylab("Precision") +
#    theme_fs() +
    theme(legend.position = "bottom")
  
  pl.rec <- group_by(subset(dat, measure == "recall"), n, p, nbi, nbij, measure) %>%
    summarise(change = yes - no) %>%
    ggplot(aes(x = nbij, y = change)) +
    geom_bar(aes(fill = nbi), stat = "identity", position = "dodge") +
    # facet_wrap(~measure, scale = "free_y", ncol = 1) +
    # scale_fill_manual(name = "Additional\nmain effects", values = c("#7fcdbb", "#1d91c0", "#253494")) +
    scale_fill_discrete(name = "True additional\nmain effects") +
    ylim(c(-0.8, 0)) +
    xlab("True interactions") +
    ylab("Recall") +
#    theme_fs() +
    theme(legend.position = "bottom")
  
  pdf(sprintf("PrecRecF1/test_analysis_n%d.pdf", numrows), width = 3, width = 3)
  grid.arrange(pl.prec, pl.rec, ncol = 1)
  dev.off()
}
