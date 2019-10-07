#!/usr/bin/env Rscript
library(normInfectX)
library(dplyr)
library(glinternet)
library(Matrix)

ifx = readRDS("./IFX_QIAGEN.rds")

cl_output = clean("bartonella", cases = c("Kinome"), libraries = c("Qiagen"))
#cl_output = readRDS("./cl_output.rds")

cl_bartonella = cl_output$bartonella
filtered_cl_bartonella = cl_bartonella %>% select(Catalog_number, eCount_oCells) %>% filter(!is.na(eCount_oCells))

X = ifx[filtered_cl_bartonella$Catalog_number,]
X = X[,colSums(X) != 0]
Y = log2(filtered_cl_bartonella$eCount_oCells / mean(filtered_cl_bartonella$eCount_oCells))

gl_output = glinternet.cv(X, Y, numLevels=rep(1,dim(X)[2]))

saveRDS(gl_output, "gl_output.rds")
