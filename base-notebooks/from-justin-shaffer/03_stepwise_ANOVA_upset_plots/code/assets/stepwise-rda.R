#!/usr/bin/env Rscript

cat(R.version$version.string, "\n")

# load arguments ---------------------------------------------------------
args <- commandArgs(TRUE)
inp.ordination.path  <- args[[1]] # ord table file 
inp.mf.path <- args[[2]] # mf table
output <- args[[3]] # to write to temp file

# load libraries ----------------------------------------------------------
library(vegan)


# load data ---------------------------------------------------------------
if(!file.exists(inp.ordination.path)) {
  errQuit("Input Ordination file does not exist.")
} else {
  Y_16S <- read.csv(inp.ordination.path, sep='\t')
}

if(!file.exists(inp.mf.path)) {
  errQuit("Input Metadata table file does not exist.")
} else {
  X_16S <- read.csv(inp.mf.path, sep='\t')
}

drops <- c("X.SampleID")
Y_16S <- Y_16S[ , !(names(Y_16S) %in% drops)]
X_16S <- X_16S[ , !(names(X_16S) %in% drops)]

# generate effect sizes ---------------------------------------------------_
(dbrda_0 <- rda(Y_16S ~ 1., 
                  X_16S, scale=TRUE)) # Model with intercept only
(dbrda_1 <- rda(Y_16S ~ ., 
                 X_16S, scale=TRUE)) # Model with all explanatory variables
#alias(dbrda_1, names=TRUE)
step.res <- ordiR2step(dbrda_0, dbrda_1,
                       perm.max = 5000,
                       steps = 5000,
                       R2permutations = 5000,
                       permutations = how(nperm = 5000),
                       trace = FALSE,
                       Pin = 0.05,
                       direction = c("both", "forward"),
                       R2scope = FALSE,
                       na.action=na.exclude)
anova_table = step.res$anova

write.table(anova_table,
            file=output,
            quote=FALSE, 
            sep='\t',
            col.names = NA)
