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
fitness_diff = c()
cl_output$adeno$seed_cells = 700
cl_output$bartonella$seed_cells = 300
cl_output$brucella$seed_cells = 500
cl_output$listeria$seed_cells = 600
cl_output$rhino$seed_cells = 1000
cl_output$salmonella$seed_cells = 550
cl_output$shigella$seed_cells = 600
cl_output$vaccinia$seed_cells = 600
for (pathogen in cl_output) {
	kinases = c(kinases, as.character(pathogen$ID))
	pathogen = pathogen %>% filter(!is.na(eCount_oCells))
	catalog_numbers = c(catalog_numbers, pathogen$Catalog_number)
	fitness_diff = c(fitness_diff, pathogen$eCount_oCells / pathogen$seed_cells)
}
kinases = unique(kinases)
#catalog_numbers = unique(catalog_numbers)

X = ifx[catalog_numbers, kinases]
X[X != 0] <- 1
X = as.matrix(X)
X = X[,!duplicated(t(X))]
#X = X[,colSums(X) != 0]
Y = log2(fitness_diff)

saveRDS(X, "X.rds")
saveRDS(fitness_diff, "fitness_diff.rds")
saveRDS(Y, "Y.rds")
gl_output = glinternet.cv(X, Y, numLevels=rep(1,dim(X)[2]), verbose=TRUE,  numCores=10)
#
saveRDS(gl_output, "gl_output.rds")
