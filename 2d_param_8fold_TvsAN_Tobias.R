#################################################
# prepare workspace
#################################################

args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1]) # first ID to process
end <- as.numeric(args[2]) # last consequitve ID to process
load("./data_BRCA.RData") # data to work on, includes folds_an and folds_t

#################################################
# libs
#################################################

library(dgRaph)
library(dplyr)
library(VGAM)
library(bbmle)
library(pROC) # calculates AUC
loglik <- function(alpha, beta){
  -sum(dbetabinom.ab(xi, ni, alpha, beta, log = T))
}

#################################################
# Iterate through jobs
#################################################
for (i in beg:end){
  cat(paste("doing ",i,"\n",sep=""))
  ptm <- proc.time()[3]
  
  # 
  model <- list()
  scores_xval <- vector(length=nrow(data_BRCA[[i]]))
  
  # common data frame with unobserved expression
  df <- data_BRCA[[i]] %>% 
    as.data.frame() %>%
    mutate(EXPR = NA) %>%
    mutate(PR_overall = NA) %>%
    mutate(GB_overall = NA) %>%
    mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_p,max_p,length.out = 101), labels = c(1:100)))), PR = starts_with("pr", ignore.case = F)) %>%
    mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_gb,max_gb,length.out = 101), labels = c(1:100)))), GB = starts_with("gb", ignore.case = F)) %>%
    select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count)
  
  # common data frame with observed expression pointestimate
  df_w_expr <- df %>% 
    mutate(EXPR = data_BRCA[[i]][2]/ data_BRCA[[i]][1] ) %>% # readcount/lib_size
    mutate(EXPR = as.integer(cut(EXPR, breaks = seq(min_e,max_e,length.out = 101), labels = c(1:100))))

  # Number of methylation probes
  n_pr_cpg <- sum(grepl('^PR', names(df)))-1
  n_gb_cpg <- sum(grepl('^GB', names(df)))-1
  
  # upper and lower limits
  min_p <- min(pr)-0.1
  max_p <- max(pr)+0.1
  min_gb <- min(gb)-0.1
  max_gb <- max(gb)+0.1
  min_e <- max(0,min(expr)-0.001*min(expr))
  max_e <- max(expr)+0.001*max(expr)
  
  # Base models
  varDim <- rep(100, 3+n_pr_cpg+n_gb_cpg)
  facPot <- list(linregPotential(dim = c(100, 100)), # PR_overall | EXPR
                 linregPotential(dim = c(100, 100)), # GB_overall | EXPR
                 linregPotential(dim = c(100, 100), range1 = c(min_p, max_p), range2 = c(min_p,max_p), alpha = 1, beta = 0, var = (max_p-min_p)**2/25), # PR_i | PR
                 linregPotential(dim = c(100, 100), range1 = c(min_gb, max_gb), range2 = c(min_gb,max_gb), alpha = 1, beta = 0, var = (max_gb-min_gb)**2/25)) # GB_i | GB
  facNbs <- c(list(c(1,2)), # PR_overall | EXPR
              list(c(1,3)), # GB_overall | EXPR
              lapply(4:(3+n_pr_cpg), FUN=function(i){c(2,i)}), # PR_i | PR_overall
              lapply((1+3+n_pr_cpg):(3+n_pr_cpg+n_gb_cpg), FUN=function(i){c(3,i)}) # GB_i | GB_overall
  )
  potMap <- c(1, 2, rep(3, n_pr_cpg), rep(4, n_gb_cpg))
  optimFun <- list(linreg1 = linregOptimize(range1 = c(min_e,max_e), range2 = c(min_p,max_p)),
                   linreg2 = linregOptimize(range1 = c(min_e,max_e), range2 = c(min_gb,max_gb)))
  dfg_base <- dfg(varDim, facPot, facNbs, potMap, varNames = names(df))
  dfg_normals <- dfg_base
  dfg_tumours <- dfg_base
  
  # dfg with beta prior
  cur_length <- length(varDim)
  varDimPrior <- varDim
  facPotPrior <- facPot
  facPotPrior[[5]] <- matrix(1, 1,100)
  facNbsPrior <- facNbs
  facNbsPrior[[cur_length]] <- 1
  potMapPrior <- potMap
  potMapPrior[cur_length] <- 5
  dfg_prior_base <- dfg(varDimPrior, facPotPrior, facNbsPrior, potMapPrior, varNames = names(df_normals))
  
  # begin x-fold here
  for (fold in 1:8) { # iterate through 8 folds
    # identify training and validation sets
    an_xfold_training <- setdiff((1:82),unlist(folds_an[[fold]]))
    t_xfold_training <- setdiff(82+(1:730),unlist(folds_t[[fold]]))
    df_normals <- df[ an_xfold_training,]
    df_tumours <- df[ t_xfold_training,]
    
    #################################################
    # Normals model
    #################################################
    
    # Setup data
    data_normals_train <- as.data.frame(data_BRCA[[i]][an_xfold_training,])
    
    # Learn betabinomial
    ni <- as.integer(data_normals_train$lib_size)
    xi <- as.integer(data_normals_train$read_count + 1)
    
    m0 <- mle2(minuslogl = loglik, start = list(alpha = 1, beta = 1), method = "L-BFGS-B", lower=c(alpha = 0.0001, beta = 0.0001))
    alpha <- coef(m0)['alpha']
    beta <- coef(m0)['beta']
    
    # Posteriors
    alphaPost <- alpha + data_normals_train$read_count
    betaPost <- beta + data_normals_train$lib_size - data_normals_train$read_count
    
    # Data list
    dataList <- list()
    dataList[[1]] <- lapply(1:nrow(df_normals), FUN=function(i){
      breaks <- seq(min_e, max_e, length.out = 101)
      diff( pbeta(breaks, alphaPost[i], betaPost[i]))
    })
    
    # Train
    dfg_normals <- train(data = df_normals, 
                         dataList = dataList,
                         dfg = dfg_normals, 
                         optim = c("linreg1", "linreg2", "noopt", "noopt"), 
                         optimFun = optimFun, iter.max = 5000, threshold = 1e-5, verbose = T)
    
    # Model with beta prior
    dfg_normals_prior <- dfg_prior_base
    potentials(dfg_normals_prior) <- c(potentials(dfg_normals), list(betaPotential(range=c(min_e,max_e),alphas=alpha, betas=beta)))
    
    #################################################
    # Tumours model
    #################################################
    data_tumours_train <- as.data.frame(data_BRCA[[i]][t_xfold_training,])
    
    # Learn beta binomial
    ni <- as.integer(data_tumours_train$lib_size)
    xi <- as.integer(data_tumours_train$read_count)
    
    m0 <- mle2(minuslogl = loglik, start = list(alpha = 1, beta = 1), method = "L-BFGS-B", lower=c(alpha = 0.0001, beta = 0.0001))
    alpha <- coef(m0)['alpha']
    beta <- coef(m0)['beta']
    
    # Posteriors
    alphaPost <- alpha + data_tumours_train$read_count
    betaPost <- beta + data_tumours_train$lib_size - data_tumours_train$read_count
    
    # Data list
    dataList <- list()
    dataList[[1]] <- lapply(1:nrow(df_tumours), FUN=function(i){
      breaks <- seq(min_e, max_e, length.out = 101)
      diff( pbeta(breaks, alphaPost[i], betaPost[i]))
    })
    
    # Train
    dfg_tumours <- train(data = df_tumours, 
                         dataList = dataList,
                         dfg = dfg_tumours, 
                         optim = c("linreg1", "linreg2", "noopt", "noopt"), 
                         optimFun = optimFun, iter.max = 20, threshold = 1e-5, verbose = T)
    
    # Model with beta prior
    dfg_tumours_prior <- dfg_prior_base
    potentials(dfg_tumours_prior) <- c(potentials(dfg_tumours_prior), list(betaPotential(range=c(min_e,max_e),alphas=alpha, betas=beta)))
    
    #################################################
    # Evaluation and performance
    #################################################
    # training data performance
    df_train <- df_w_expr[c(an_xfold_training,t_xfold_training),]
    df_eval <- df_w_expr[-c(an_xfold_training,t_xfold_training),]

    likelihood1 <- likelihood(dfg = dfg_tumours, data = df_train, log = T)
    likelihood1[which(is.na(likelihood1))] <- vapply(likelihood1[which(is.na(likelihood1))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
    likelihood2 <- likelihood(dfg = dfg_normals, data = df_train, log = T)
    likelihood2[which(is.na(likelihood2))] <- vapply(likelihood2[which(is.na(likelihood2))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
    scores_train <- likelihood1 - likelihood2
    
    likelihood1 <- likelihood(dfg = dfg_tumours, data = df_eval, log = T)
    likelihood1[which(is.na(likelihood1))] <- vapply(likelihood1[which(is.na(likelihood1))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
    likelihood2 <- likelihood(dfg = dfg_normals, data = df_eval, log = T)
    likelihood2[which(is.na(likelihood2))] <- vapply(likelihood2[which(is.na(likelihood2))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
    scores_eval <- likelihood1 - likelihood2
    
    # Metrics
    model$auc_training[fold] <- auc(predictor=scores_train,response=c(rep("N",length(an_xfold_training)),rep("T",length(t_xfold_training))))
    model$auc_evaluation[fold] <- auc(predictor=scores_eval,response=c(rep("N",(82-length(an_xfold_training))),rep("T",(730-length(t_xfold_training)))))
    model$kls[fold] <- kl(dfg_normals,dfg_tumours)
    
    scores_xval[c(folds_an[[fold]],folds_t[[fold]])] <- scores_eval
  }
  # end x-fold here
  model$scores <- scores_xval
  eval(parse(text=paste('save(model, file="./',i,'_model.RData")',sep="")))
  cat(paste("done evaluating for ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
