#!/usr/bin/env Rscript
require(ggplot2)
require(dplyr)
require(RColorBrewer)
#source("~/Projects/R/fs_.R")
#setwd("~/Projects/epistasis/results/simulation")
setwd("..")

args <- commandArgs(trailingOnly = TRUE)

if (args[1] == "x") {
	use_xyz = TRUE
	fits_path='./fits_proper/'
} else if (args[1] == "g") {
	use_xyz = FALSE
	fits_path='./fits_glinternet/'
} else {
    cat("[x]yz or [g]linternet required\n")
    q()
}

#if (length(args) >= 3) {
#    append_str = args[3]
#} else {
#    append_str = ''
#}

cat("using lethal data\n")
#graph_numrows <- c(5000)
#graph_nbij <- c("0", "200", "500", "1000")
#graph_nlethals <- c("10", "20", "40", "80", "100")
#large_int <- FALSE
#append_str <- "lethal"
rds_file = sprintf("time_taken/time_taken_xyz%s.rds", use_xyz)

if (args[2] == 'y') {
    # Precision, recall and F1 for interaction terms
     ans <- lapply(list.files(path = fits_path, pattern = "", full.names = TRUE), function(f) {#, sprintf("n%d_p%d", n, p)), function(f) {
       ans <- readRDS(f)
       n <- regmatches(x = f, m = regexpr(f, pattern = "(?<=n)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       p <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_p)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       SNR <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_SNR)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       if (use_xyz)
           L <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_L)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       else
           L <- 0
       nbi <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbi)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       nbij <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nbij)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       nlethals <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_nlethals)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       perc_viol <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_viol)\\d+(?=_)", perl = TRUE)) %>% as.numeric
       ID <- regmatches(x = f, m = regexpr(f, pattern = "(?<=_)\\d+(?=_|\\.rds)", perl = TRUE)) %>% as.numeric
       smry_int <- ans$smry %>% filter(type == "interaction")
       lethal <- ans$smry %>% select(lethal)
       time <- ans$time
       smry_int[["lethal"]][is.na(smry_int[["lethal"]])] <- FALSE
       notest <- data.frame(n = n, p = p, SNR = SNR, L = L, nbi = nbi, nbij = nbij, nlethals = nlethals, time = time[["elapsed"]],
                            precision = sum(smry_int[["lethal"]]) / nrow(smry_int),
                            recall = sum(smry_int[["lethal"]]) / nlethals) %>%
         mutate(F1 = 2 *  (precision * recall) / (precision + recall),
                test = "no")
       smry_int <- mutate(smry_int, pval = p.adjust(pval, method = "BH")) %>%
         filter(pval < 0.05)
       test <- data.frame(n = n, p = p, SNR = SNR, L = L, nbi = nbi, nbij = nbij, nlethals = nlethals, time = time[["elapsed"]],
                          precision = sum(smry_int[["lethal"]]) / nrow(smry_int),
                          recall = sum(smry_int[["lethal"]]) / nlethals) %>%
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
              nbij = factor(nbij),
              nlethals = factor(nlethals),
       	      time = time)
     saveRDS(ans, file = rds_file)
}


dat_time <- readRDS(file = rds_file) %>%
 filter(nbij == as.numeric(as.character(p))/5) %>%
 filter(SNR != "1") %>%
 #filter(lethal == g_lethal) %>%
 filter(p %>% as.character %>% as.numeric == (n %>% as.character %>% as.numeric)/2) %>%
 mutate(frac_lethals = nlethals %>% as.character %>% as.numeric / (p %>% as.character %>% as.numeric)) %>%
 filter(frac_lethals == 0.05) %>%
 mutate(frac_bij = nbij %>% as.character %>% as.numeric / (p %>% as.character %>% as.numeric)) %>%
 filter(frac_bij == 0.2) %>%
 mutate(frac_bi = nbi %>% as.character %>% as.numeric / (p %>% as.character %>% as.numeric)) %>%
 filter(frac_bi == 0.01) %>%
 mutate(frac_lethals = factor(frac_lethals)) %>%
 mutate(pval = p %>% as.character %>% as.numeric)

pl <- group_by(dat_time, pval, frac_lethals) %>%
  summarise(mean = mean(time, na.rm = TRUE), sd = sd(time, na.rm = TRUE) / sqrt(n())) %>%
  ggplot(aes(x = pval, 
             y = mean, 
             group = 1,
             ymax = mean + sd, 
             ymin = mean - sd)) +
  #geom_line(aes(colour = "red")) +
  #geom_errorbar(colour = "darkgrey", width=.3, position=position_dodge(.35)) +
  #geom_point(aes(x=pval,y=mean)) +
  geom_line(aes(colour = frac_lethals), position=position_dodge(.35)) +
  geom_point(aes(colour = frac_lethals), position=position_dodge(.35), size = 1) +
  geom_errorbar(colour = "darkgrey", width=.3, position=position_dodge(.35)) +
  #scale_y_continuous(trans='log2') +
  #scale_x_continuous(trans='log2') +
  xlab("Genes") +
  ylab("Time to discover interactions (seconds)") +
  theme_bw() +
  theme(legend.position = "none")
  #theme_bw(base_size=15)
pl
ggsave(pl, file = sprintf("time_taken/time_taken_xyz%s.pdf", use_xyz), width = 3, height = 3)
