#!/usr/bin/Rscript
require(dplyr)
require(tidyr)
require(ggplot2)
require(RColorBrewer)
#source("~/Projects/R/fs_.R")
#setwd("~/Projects/epistasis/results/simulation")

## Number of observations
 ans <- lapply(list.files(path = "./fits_proper/", pattern = "", full.names = TRUE), function(f) {#, sprintf("n%d_p%d", n, p)), function(f) {
   ans <- readRDS(f)
   n <- regmatches(x = f, m = regexpr(f, pattern = "(?<=n)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   p <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_p)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   SNR <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_SNR)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   nbi <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbi)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   nbij <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbij)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   perc_viol <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_viol)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   regression_alpha <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_alpha)[\\d.]+(?=_)", perl = TRUE)) %>% as.numeric
   id <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_)\\d+(?=\\.rds)", perl = TRUE)) %>% as.numeric
#   smry_int <- ans$smry %>% filter(type == "interaction")
#   results <- ans$results
#   coef <- results[[2]][[10]]
#   coef <- ans$results %>% filter(type == "coefficient")
   coef <- ans$coef
   tp <- ans$tp
   count <- ans$count
   notest <- data.frame(n = n, p = p, SNR = SNR, nbi = nbi, nbij = nbij, id = id, type="TP", coef=coef, #TODO: correct type
                       precision = tp / count,
                       recall = tp / nbij,
#              left_join(ans$bij, smry_int, by = c("gene_i", "gene_j", "o00", "o01", "o10", "o11", "omin")) %>%
#                select(coef = coef, coef_est = coef.est, TP, omin) %>%
#                mutate(type = ifelse(is.na(TP), "FN", "TP")) %>%
#                select(-TP) %>%
#                rbind(., filter(smry_int, TP == FALSE) %>% mutate(coef = NA) %>% select(coef, coef_est = coef.est, omin) %>% mutate(type = "FP")),
              test = "no") %>% tbl_df
#   smry_int <- mutate(smry_int, pval = p.adjust(pval, method = "BH")) %>%
#     filter(pval < 0.05)
   test <- data.frame(n = n, p = p, SNR = SNR, nbi = nbi, nbij = nbij, id = id, type="TP", coef=coef, #TODO: correct type
                       precision = tp / count,
                       recall = tp / nbij,
#                      left_join(ans$bij, smry_int, by = c("gene_i", "gene_j", "o00", "o01", "o10", "o11", "omin")) %>%
#                        select(coef = coef, coef_est = coef.est, TP, omin) %>%
#                        mutate(type = ifelse(is.na(TP), "FN", "TP")) %>%
#                        select(-TP) %>%
#                        rbind(., filter(smry_int, TP == FALSE) %>% mutate(coef = NA) %>% select(coef, coef_est = coef.est, omin) %>% mutate(type = "FP")),
                      test = "yes") %>% tbl_df
   rbind(test, notest)
 }) %>% do.call("rbind", .) %>%
   tbl_df %>%
   mutate(n = factor(n),
          p = factor(p),
          SNR = factor(SNR),
          nbi = factor(nbi),
          nbij = factor(nbij),
          type = factor(type))
 saveRDS(ans, file = "FXstrength/dat_fxstrength.rds")




for (numrows in c(1000)) { #TODO: 400?
  for (t in c("yes", "no")) {
    dat_fxs <- readRDS("FXstrength/dat_fxstrength.rds") %>%
      filter(n == numrows) %>%
      filter(test == t) %>%
      filter(nbi == 0, SNR != 1) %>%
      mutate(SNR = factor(SNR, labels = paste0("SNR = ", levels(factor(SNR))))) %>%
      rowwise %>%
#      mutate(coef = ifelse(is.na(coef), coef_est, coef)) %>%
#      select(-coef_est) %>%
      ungroup %>%
      filter(nbij %in% c(5, 20, 50, 100)) %>%
      group_by(n, p, SNR, nbi, nbij, type, precision, recall) %>%
      mutate(range = cut(coef, c(-Inf, seq(-0.2, 0.2, by = 0.025), Inf))) %>%
      group_by(n, p, SNR, nbij, type, range, precision, recall) %>%
      summarise(count = n()) %>% 
      spread(type, count, fill = 0) %>%
      # filter(FP > 0) %>% # avoid 0 effects for test = yes
#      mutate(precision = TP / (TP + FP),
#             recall = TP / (TP + FN)) %>%
      mutate(F1 = 2 *  (precision * recall) / (precision + recall)) %>%
      filter(!is.na(precision), !is.na(recall), !is.na(F1)) %>%
      gather(measure, value, precision, recall, F1) %>%
      mutate(measure = factor(measure, levels = c("precision", "recall", "F1"), labels = c("Precision", "Recall", "F1")))
    
#    levels(dat_fxs$range) <- c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5)
#    dat_fxs$range <- factor(dat_fxs$range, labels = c("everything"))
    dat_fxs <- dat_fxs %>% group_by(n, p, SNR, nbij, measure, range) %>%
      summarise(mean = mean(value, na.rm = TRUE), sem = sd(value, na.rm = TRUE) / sqrt(n()))
    
    pl <- ggplot(dat_fxs, aes(x = range,#observations %>% as.character %>% as.numeric,
                               y = mean, 
                               group = nbij,
                               ymax = mean + sem, 
                               ymin = mean - sem)) +
      geom_line(aes(colour = nbij), position = position_dodge(.35)) +
      geom_point(aes(colour = nbij), position = position_dodge(.35), size = 1) +
      geom_errorbar(colour = "darkgrey", width = 0.3, position = position_dodge(.35)) +
#      scale_x_discrete(labels = c(expression(paste("(-", infinity, ", -0.3]")), "(-0.3,-0.1]", "(-0.1,0.1]", "(0.1,0.3]", expression(paste("(0.3, ", infinity, ")")))) +
      facet_grid(measure~SNR) +
      scale_color_discrete(name = "True interactions") +
      # ylim(c(0,1)) +
      xlab("True epistasis") +
      ylab("") +
     # theme_fs() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    pl
    ggsave(pl, file = sprintf("FXstrength/FXstrength_PRF_n%d_t%s.pdf", numrows, t), width = 5, height = 7)
    
    
#   pl_wrongdir <- readRDS("FXstrength/dat_fxstrength.rds") %>%
#      filter(n == numrows) %>%
#      filter(test == t) %>%
#      filter(type == "TP") %>%
#      filter(nbi == 0, SNR != 1) %>%
#      filter(nbij %in% c(5, 20, 50, 100)) %>%
#      mutate(SNR = factor(SNR, labels = paste0("SNR = ", levels(factor(SNR))))) %>%
##      mutate(coef = ifelse(is.na(coef), coef_est, coef)) %>%
##      mutate(corrdir = sign(coef) == sign(coef_est)) %>%
#      mutate(range = cut(coef, c(-Inf, seq(-3, 3, by = 2), Inf))) %>%
#      group_by(n, p, SNR, nbi, nbij, range) %>%
##      summarise(wrongdir = 1e-2 + (1 - sum(corrdir) / n()) * 100) %>%
#      ggplot(aes(x = range, y = wrongdir)) +
#      geom_bar(aes(fill = nbij), stat = "identity", position = "dodge") +
#      scale_x_discrete(labels = c(expression(paste("(-", infinity, ", -3]")), "(-3,-1]", "(-1,1]", "(1,3]", expression(paste("(3, ", infinity, ")")))) +
#      facet_grid(.~SNR) +
#      scale_fill_discrete(name = "True interactions") +
#      ylim(c(0,15)) +
#      xlab("True epistasis") +
#      ylab("Incorrect direction [%]") +
##      theme_fs() +
#      theme(legend.position = "bottom",
#            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#    pl_wrongdir
#    ggsave(pl_wrongdir, file = sprintf("FXstrength/FXstrength_direction_n%d_t%s.pdf", numrows, t), width = 4.5, height = 4)
  }
}



