#!/usr/bin/Rscript
library(normInfectX)
library(dplyr)
library(glinternet)

ifx = readRDS("./IFX_QIAGEN.rds")

#cl_output = clean()
cl_output = readRDS("./cl_output.rds")

cl_adeno = cl_output$adeno
filtered_cl_adeno = cl_adeno %>% filter(LIBRARY=="Qiagen") %>% select(Catalog_number, eCount_oCells) %>% filter(!is.na(eCount_oCells))

X = ifx[filtered_cl_adeno$Catalog_number,]
Y = log2(filtered_cl_adeno$eCount_oCells / mean(filtered_cl_adeno$eCount_oCells))

gl_output = glinternet.cv(X, Y, numLevels=rep(1,dim(X)[2]))
