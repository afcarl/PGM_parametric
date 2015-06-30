# Michal Test Case

#################################################
# Read data
#################################################

library(dgRaph)
library(dplyr)
A1BG <- read.delim("tests/examples/A1BG.tab")

#################################################
# Helper functions
#################################################

# Visualize potential
plotPotential <- function(pot){
  require(ggplot2)
  require(reshape2)
  melted <- melt(pot)
  ggplot(melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile()
}

#################################################
# Testing submodel - Beta distributed expression
#################################################

# Expression Modelled as mixture of beta distributions
df_expression <- A1BG %>% 
  mutate(EXPR = read_count / lib_size) %>%
  mutate(H1 = NA) %>% 
  mutate(EXPR = as.integer(cut(EXPR, breaks = seq(0,3e-5,length.out = 101), labels = c(1:100)))) %>%
  select(H1, EXPR)

varDim <- c(2, 100)
facPot <- list(multinomialPotential(dim = c(1,2)),
               betaPotential(dim = c(2, 100), range = c(0,3e-5)) )
facNbs <- c(list(1),
            list(c(1,2)))
dfg_beta <- dfg(varDim, facPot, facNbs)

# Train model
optim <- c('row', 'custombeta')
optimFun <- list(custombeta = betaOptimize(range = c(0, 3e-5)))
train(df_expression, dfg_beta, optim = optim, optimFun = optimFun, threshold = 1e-9)

#################################################
# Testing submodel - Methylation linear regression in expression
#################################################

# Process data
df_methylation <- A1BG %>% 
  mutate(EXPR = read_count / lib_size) %>%
  mutate(EXPR = as.integer(cut(EXPR, breaks = seq(0,3e-5,length.out = 101), labels = c(1:100)))) %>%
  mutate(PR = (pr_1 + pr_2 + pr_3 + pr_4 + pr_5 + pr_6)/6) %>%
  mutate(GB = (gb_1 + gb_2 + gb_3 + gb_4 + gb_5 + gb_6 + gb_7 + gb_8 + gb_9)/9 ) %>%
  mutate(GB = as.integer(cut(GB, breaks = seq(-3,3,length.out = 101), labels = c(1:100)))) %>%
  mutate(PR = as.integer(cut(PR, breaks = seq(0,4,length.out = 101), labels = c(1:100)))) %>%
  select(EXPR, GB, PR)
  
# Build model
varDim <- c(100, 100, 100)
facPot <- c(list(betaPotential(c(1,100))),
            list(linregPotential(dim = c(100, 100))),
            list(linregPotential(dim = c(100, 100))))
facNbs <- c(list(1),
            list(c(1,2)),
            list(c(1,3)))
dfg_methyl <- dfg(varDim, facPot, facNbs)

# Train model
likelihood(df_methylation, dfg_methyl)
optim <- c('custombeta', 'linreg1', 'linreg2')
optimFun <- list(custombeta = betaOptimize(range = c(0, 3e-5)),
                 linreg1 = linregOptimize(range1 = c(0,3e-5), range2 = c(-3,3)),
                 linreg2 = linregOptimize(range1 = c(0,3e-5), range2 = c(0,4)))
train(df_methylation, dfg_methyl, optim = optim, optimFun = optimFun)

# Compare potentials and data

# Beta prior
newFacPot <- potentials(dfg_methyl)
hist(df_methylation$EXPR, breaks = 20)
points(newFacPot[[1]][1,]*150*4)

# Linear reg expr-GB
df_expr_GB = matrix(0, 100, 100)
for(i in 1:nrow(df_methylation)){
  df_expr_GB[df_methylation$EXPR[i], df_methylation$GB[i]] <- df_expr_GB[df_methylation$EXPR[i], df_methylation$GB[i]] + 1
}
library(reshape2)
melt_df_expr_gb <- melt(df_expr_GB)
library(ggplot2)
ggplot(melt_df_expr_gb, aes(x = Var1, y = Var2, fill = as.factor(value))) +
  geom_tile() + 
  scale_fill_manual(values = c('#000000', '#009E73', '#D55E00'))
melt_df_expr_gb_pot <- melt(newFacPot[[2]])
ggplot(melt_df_expr_gb_pot, aes(x = Var1, y = Var2, fill = value)) + geom_tile()

# Linear reg expr-PR
df_expr_PR = matrix(0, 100, 100)
for(i in 1:nrow(df_methylation)){
  df_expr_PR[df_methylation$EXPR[i], df_methylation$PR[i]] <- df_expr_PR[df_methylation$EXPR[i], df_methylation$PR[i]] + 1
}
library(reshape2)
melt_df_expr_pr <- melt(df_expr_PR)
library(ggplot2)
ggplot(melt_df_expr_pr, aes(x = Var1, y = Var2, fill = as.factor(value))) +
  geom_tile() + 
  scale_fill_manual(values = c('#000000', '#009E73', '#D55E00'))
melt_df_expr_pr_pot <- melt(newFacPot[[3]])
ggplot(melt_df_expr_pr_pot, aes(x = Var1, y = Var2, fill = value)) + geom_tile()

#################################################
# Testing submodel - Fixed link normal distributions
#################################################

# Process data
df_probes <- A1BG %>% 
  mutate_each(funs(bin = as.integer(cut(., breaks = seq(-1.5,5,length.out = 101), labels = c(1:100)))), PR = starts_with("pr")) %>%
  mutate(PR = NA) %>%
  select(starts_with("PR", ignore.case = F))

# Build Model
varDim <- rep(100, 7)
facPot <- c(list(fixedNormalPotential()))
facNbs <- c(list(c(7,1)),
            list(c(7,2)),
            list(c(7,3)),
            list(c(7,4)),
            list(c(7,5)),
            list(c(7,6)))
potMap <- rep(1,6)
dfg_probes <- dfg(varDim, facPot, facNbs, potMap = potMap)

# Train
train(df_probes, dfg_probes, optim = rep("noopt", 6))

#################################################
# The whole deal
#################################################

# Process data
df_full <- as.matrix(data_BRCA$RAG1AP1) %>% 
  mutate(EXPR = read_count / lib_size) %>%
  mutate(EXPR = as.integer(cut(EXPR, breaks = seq(0,3e-5,length.out = 101), labels = c(1:100)))) %>%
  mutate(H1 = NA) %>%
  mutate_each(funs(bin = as.integer(cut(., breaks = seq(-1.5,5,length.out = 101), labels = c(1:100)))), PR = starts_with("pr")) %>%
  mutate(PR = NA) %>%
  mutate_each(funs(bin = as.integer(cut(., breaks = seq(-6.5,3.5,length.out = 101), labels = c(1:100)))), GB = starts_with("gb")) %>%
  mutate(GB = NA) %>%
  select(H1, EXPR, starts_with("PR", ignore.case = F), starts_with("GB", ignore.case = F))

# Build Model
varDim <- c(2, rep(100, 19))
facPot <- list(multinomialPotential(dim = c(1,2)), # Prior H1
               betaPotential(c(2,100)), # EXPR | H1
               linregPotential(dim = c(100, 100)), # PR | EXPR
               linregPotential(dim = c(100, 100)), # GB | EXPR
               fixedNormalPotential(), # PR_i | PR
               fixedNormalPotential()) # GB_i | GB
facNbs <- c(list(1), # Prior H1
            list(c(1,2)), # EXPR | H1
            list(c(2,9)), # PR | EXPR
            list(c(2,20)), # GB | EXPR
            lapply(3:8, FUN=function(i){c(9,i)}), # PR_i | PR
            lapply(10:19, FUN=function(i){c(20,i)}) # GB_i | GB
            )

potMap <- c(1, 2, 3, 4, rep(5, 6), rep(6, 10))
dfg_full <- dfg(varDim, facPot, facNbs, potMap, varNames = names(df_full))
plot(dfg_full)

# Train
optimFun <- list(custombeta = betaOptimize(range = c(0, 3e-5)),
                 linreg1 = linregOptimize(range1 = c(0,3e-5), range2 = c(-3,3)),
                 linreg2 = linregOptimize(range1 = c(0,3e-5), range2 = c(0,4)))
train(df_full, dfg_full, optim = c("row", "custombeta", "linreg1", "linreg2", "noopt", "noopt"), optimFun = optimFun, iter.max = 500, threshold = 1e-1)
