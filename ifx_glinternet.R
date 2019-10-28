#!/usr/bin/env Rscript
library(normInfectX)
library(dplyr)
library(glinternet)
library(Matrix)

ifx = readRDS("./IFX_QIAGEN.rds")

cl_output = clean(cases = c("Kinome"), libraries = c("Qiagen"))
#cl_output = readRDS("./cl_output.rds")

#cl_bartonella = cl_output$bartonella
#filtered_cl_bartonella = cl_bartonella %>% select(Catalog_number, eCount_oCells) %>% filter(!is.na(eCount_oCells))

kinases = c()
catalog_numbers = c()
fitness_measure = c()
for (pathogen in cl_output) {
	kinases = c(kinases, as.character(pathogen$ID))
	pathogen = pathogen %>% filter(!is.na(eCount_oCells))
	catalog_numbers = c(catalog_numbers, pathogen$Catalog_number)
	fitness_measure = c(fitness_measure, pathogen$eCount_oCells)
}
kinases = unique(kinases)
#catalog_numbers = unique(catalog_numbers)

X = ifx[catalog_numbers, kinases]
X[X != 0] <- 1
#X = X[,colSums(X) != 0]
Y = log2(fitness_measure / mean(fitness_measure))

saveRDS(X, "X.rds")
saveRDS(fitness_measure, "fitness_measure.rds")
saveRDS(Y, "Y")
gl_output = glinternet.cv(X, fitness_measure, numLevels=rep(1,dim(X)[2]), verbose=TRUE,  numCores=10)
#
saveRDS(gl_output, "gl_output.rds")
