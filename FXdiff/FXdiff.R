#!/usr/bin/env Rscript
require(ggplot2)
require(RColorBrewer)
require(dplyr)
#source("~/Projects/R/fs_.R")
setwd("..")

args <- commandArgs(trailingOnly = TRUE)
mult <- args[1] %>% as.numeric

if (args[2] == "x") {
       use_xyz = TRUE
       fits_path='./fits_proper/'
} else if (args[2] == "g") {
       use_xyz = FALSE
       fits_path='./fits_glinternet/'
} else {
    cat("[x]yz or [g]linternet required\n")
    q()
}
rds_file = sprintf("FXstrength/dat_fxstrength_xyz%s.rds", use_xyz)

for (numrows in c(1000 * mult)) { #400, 1000
  if (numrows == 400) {
    rseq <- 0#c(seq(0, 60, by = 20), 100)
  } else if (numrows == 1000) {
    rseq <- c(0, 10, 20, 40, 80, Inf)#c(seq(0, 150, by = 50), 250)
    #mult<-1
  } else if (numrows == 10000) {
    #mult<-10
    rseq <- mult*c(0, 10, 20, 40, 80, Inf)#c(seq(0, 150, by = 50), 250)
  }
  for (t in c("yes", "no")) {
    dat_dir <- readRDS(file = rds_file)
      if (use_xyz)
        dat_dir <- dat_dir %>% filter(L == round(sqrt(p %>% as.character %>% as.numeric)))
      dat_dir <- dat_dir %>%
      filter(type == "TP") %>%
      filter(n == numrows) %>%
      filter(SNR != 1) %>%
      mutate(SNR = factor(SNR, labels = paste0("SNR = ", levels(SNR)))) %>%
      filter(nbi == 20*mult) %>%
      filter(nbij %in% c(5*mult, 20*mult, 50*mult, 100*mult)) %>%
      filter(test == t) %>%
      mutate(range = cut(observations, rseq)) %>%
      rowwise %>%
      mutate(diff = abs((coef - coef_est) / coef),
             correct_direction = sign(coef) == sign(coef_est)) %>%
      ungroup %>%
      group_by(n, p, SNR, nbi, nbij, range) %>%
      summarise(diff_mean = mean(diff, na.rm = TRUE), diff_sem = sd(diff, na.rm = TRUE) / sqrt(n())) %>%
      data.frame
    
    pl <- ggplot(dat_dir, aes(x = range,#observations %>% as.character %>% as.numeric,
                              y = diff_mean,
                              group = nbij,
                              ymax = diff_mean + diff_sem,
                              ymin = diff_mean - diff_sem)) +
      geom_line(aes(colour = nbij), position = position_dodge(.25)) +
      geom_point(aes(colour = nbij), position = position_dodge(.25), size = 1) +
      geom_errorbar(colour = "darkgrey", width = 0.3, position = position_dodge(.25)) +
      facet_grid(.~SNR) +
      scale_color_discrete(name = "True interactions") +
#      scale_x_discrete(labels = c("(0, 10]", "(10, 20]", "(20, 40]", "(40, 80]", expression(paste("(80, ", infinity, ")")))) +
      xlab("Observations of double knockdown") +
      ylab("Absolute deviation [%]") +
      ylim(c(0, 4)) +
#      theme_fs() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    pl
    ggsave(pl, file = sprintf("FXdiff/FXdiff_n%d_t%s_xyz%s.pdf", numrows, t, use_xyz), width = 5, height = 4)
  }
}



