#!/usr/bin/Rscript
require(ggplot2)
require(dplyr)
require(RColorBrewer)
require(gridExtra)
#source("~/Projects/R/fs_.R")
#setwd("~/Projects/epistasis/results/simulation")
setwd("..")

args <- commandArgs(trailingOnly = TRUE)
#args <- c('large', 'no')

SNR_limits <- c(2,5,10)

if (args[1] == 'l') {
    cat("using large data\n")
    large <- TRUE
    graph_numrows <- c(10000)
    graph_nbi <- c("0", "200", "500", "1000")
    graph_nbij <- c("50", "200", "500", "1000")
} else {
    cat("using small data\n")
    graph_numrows <- c(1000)
    graph_nbi <- c("0", "20", "50", "100")
    graph_nbij <- c("5", "20", "50", "100")
    large <- FALSE
}


if (args[2] == 'y') {
# Precision, recall and F1 for interaction terms
 ans <- lapply(list.files(path = "./fits_proper/", pattern = "", full.names = TRUE), function(f) {#, sprintf("n%d_p%d", n, p)), function(f) {
   ans <- readRDS(f)
   n <- regmatches(x = f, m = regexpr(f, pattern = "(?<=n)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   p <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_p)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   SNR <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_SNR)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   nbi <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbi)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   nbij <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbij)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   perc_viol <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_viol)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   L <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_L)\\d+(?=_)", perl = TRUE)) %>% as.numeric
   ID <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_)\\d+(?=\\.rds)", perl = TRUE)) %>% as.numeric
   
   smry_int <- ans$smry %>% filter(type == "interaction")
   notest <- data.frame(n = n, p = p, SNR = SNR, L = L, nbi = nbi, nbij = nbij,
                        TP = sum(smry_int[["TP"]]), FP = nrow(smry_int) - sum(smry_int[["TP"]]),
                        precision = sum(smry_int[["TP"]]) / nrow(smry_int),
                        recall = sum(smry_int[["TP"]]) / nbij) %>%
     mutate(F1 = 2 *  (precision * recall) / (precision + recall),
            test = "no")
   smry_int <- mutate(smry_int, pval = p.adjust(pval, method = "BH")) %>%
     filter(pval < 0.05)
   test <- data.frame(n = n, p = p, SNR = SNR, L = L, nbi = nbi, nbij = nbij,
                      TP = sum(smry_int[["TP"]]), FP = nrow(smry_int) - sum(smry_int[["TP"]]),
                      precision = sum(smry_int[["TP"]]) / nrow(smry_int),
                      recall = sum(smry_int[["TP"]]) / nbij) %>%
     mutate(F1 = 2 *  (precision * recall) / (precision + recall),
            test = "yes")
   rbind(notest, test) %>% data.frame(., id = ID)
 }) %>% do.call("rbind", .) %>%
   tbl_df %>%
   mutate(n = factor(n),
          p = factor(p),
          SNR = factor(SNR),
          L = factor(L),
          nbi = factor(nbi),
          nbij = factor(nbij))
 saveRDS(ans, file = "l_diff/dat_l_diff.rds")
}

for (SNR_limit in SNR_limits) {
    for (numrows in graph_numrows) { #400
      for (t in c("yes", "no")) {
        dat_l_diff <- readRDS(file = "l_diff/dat_l_diff.rds") %>%
          filter(n == numrows) %>%
          filter(test == t) %>%
          filter(nbi %in% graph_nbi) %>%
          filter(nbij %in% graph_nbij) %>%
          filter(SNR == SNR_limit) %>%
          filter(L %in% c(10, 100, 1000)) %>%
          mutate(L = factor(L, levels = levels(L), labels = paste0("L = ", levels(L)))) %>%
          group_by(n, p, L, nbi, nbij, test) #%>% sample_n(10, replace=TRUE) #TODO: improve this?
        
        pl <- group_by(dat_l_diff, L, nbi, nbij) %>%
          summarise(mean = mean(precision, na.rm = TRUE), sd = sd(precision, na.rm = TRUE) / sqrt(n())) %>%
          mutate(Measure = "Precision") %>%
          rbind(., group_by(dat_l_diff, L, nbi, nbij) %>%
                  summarise(mean = mean(recall, na.rm = TRUE), sd = sd(recall, na.rm = TRUE) / sqrt(n())) %>%
                  mutate(Measure = "Recall")) %>%
          rbind(., group_by(dat_l_diff, L, nbi, nbij) %>%
                  summarise(mean = mean(F1, na.rm = TRUE), sd = sd(F1, na.rm = TRUE) / sqrt(n())) %>%
                  mutate(Measure = "F1")) %>%
          mutate(Measure = factor(Measure, levels = c("Precision", "Recall", "F1"))) %>%
          ggplot(aes(x = nbij, 
                     y = mean, 
                     group = nbi,
                     ymax = mean + sd, 
                     ymin = mean - sd)) +
          facet_grid(Measure ~ L) +
          # geom_bar(aes(fill = nbi), stat = "identity", position = "dodge") +
          geom_line(aes(colour = nbi), position=position_dodge(.35)) +
          geom_point(aes(colour = nbi), position=position_dodge(.35), size = 1) +
          geom_errorbar(colour = "darkgrey", width=.3, position=position_dodge(.35)) +
          scale_color_discrete(name = "True additional\nmain effects") +
          # scale_colour_manual(name = "Additional main effects", values = c("#7fcdbb", "#1d91c0", "#253494")) + #brewer.pal(7, "Set3")[-1]) + #
          ylim(c(0,0.25)) +
          xlab("True interactions") +
          ylab("") +
          #theme_fs() +
    #      scale_y_continuous(trans='log10') +
          theme(legend.position = "bottom")
        pl
        ggsave(pl, file = sprintf("l_diff/l_diff_n%d_SNR%d_t%s.pdf", numrows, SNR_limit, t), width = 5, height = 7)
    
      
      pl.prec <- group_by(subset(dat_l_diff), n, p, nbi, L) %>%
        ggplot(aes(x = L, y = TP)) +
        geom_bar(aes(fill = nbij), stat = "identity", position = "dodge") +
#         facet_wrap(~TP, scale = "free_y", ncol = 1) +
#         scale_fill_manual(name = "Additional\nmain effects", values = c("#7fcdbb", "#1d91c0", "#253494")) +
#        scale_fill_discrete(name = "True additional\nmain effects") +
        ylim(c(0, 50)) +
    #    xlab("True interactions") +
        ylab("True Positives") +
        theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
        #theme_fs() +
        theme(legend.position = "none")
      
      pl.rec <- group_by(subset(dat_l_diff), n, p, nbi, nbij) %>%
        ggplot(aes(x = L, y = FP)) +
        geom_bar(aes(fill = nbij), stat = "identity", position = "dodge") +
        scale_y_reverse() +
#         facet_wrap(~FP, scale = "free_y", ncol = 1) +
#         scale_fill_manual(name = "Additional\nmain effects", values = c("#7fcdbb", "#1d91c0", "#253494")) +
        scale_fill_discrete(name = "True interactions") +
        ylim(c(150, 0)) +
        xlab("Number of projections") +
        ylab("False Positives") +
        #theme_fs() +
        theme(legend.position = "bottom")
      
      pdf(sprintf("l_diff/quant_analysis_n%d.pdf", numrows), width = 5, height = 8)
      grid.arrange(pl.prec, pl.rec, ncol = 1)
      dev.off()
      }
   } 
}
