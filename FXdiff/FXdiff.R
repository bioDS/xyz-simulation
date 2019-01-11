require(ggplot2)
require(RColorBrewer)
source("~/Projects/R/fs_.R")
setwd("~/Projects/epistasis/results/simulation")



for (numrows in c(1000)) { #400, 
  if (numrows == 400) {
    rseq <- 0#c(seq(0, 60, by = 20), 100)
  } else if (numrows == 1000) {
    rseq <- c(0, 10, 20, 40, 80, Inf)#c(seq(0, 150, by = 50), 250)
  }
  for (t in c("yes", "no")) {
    dat_dir <- readRDS("FXstrength/dat_fxstrength.rds") %>%
      filter(type == "TP") %>%
      filter(n == numrows) %>%
      mutate(SNR = factor(SNR, labels = paste0("SNR = ", levels(SNR)))) %>%
      filter(nbi == 20, SNR != "SNR = 1") %>%
      filter(nbij %in% c(5, 20, 50, 100)) %>%
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
      scale_x_discrete(labels = c("(0, 10]", "(10, 20]", "(20, 40]", "(40, 80]", expression(paste("(80, ", infinity, ")")))) +
      xlab("Observations of double knockdown") +
      ylab("Absolute deviation [%]") +
      ylim(c(0, 4)) +
      theme_fs() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    pl
    ggsave(pl, file = sprintf("FXdiff/FXdiff_n%d_t%s.pdf", numrows, t), width = 5, height = 4)
  }
}



