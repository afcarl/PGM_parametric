#################################################
# prepare workspace
#################################################
args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1]) # first ID to process
end <- as.numeric(args[2]) # last consequitve ID to process
load("./data_BRCA.RData") # data to work on, includes folds_g1 and folds_g2

#################################################
# libs and common functions
#################################################
library(dgRaph)
library(dplyr)
library(VGAM)
library(bbmle)
library(pROC) # calculates AUC
loglik <- function(alpha, beta){-sum(dbetabinom.ab(xi, ni, alpha, beta, log = T))}

# data checks
if (length(folds_g1) != length(folds_g2)) stop("Uneven number of x-folds")
nfolds <- length(folds_g1)
length_g1 <- length(unlist(folds_g1))
length_g2 <- length(unlist(folds_g2))

for (i in beg:end) {
	cat(paste("doing ",i,"\n",sep=""))
	ptm <- proc.time()[3]
	
	model <- list()
	scores_xval <- vector(length=nrow(data_BRCA[[i]]))
	
	# calculate common discretization upper and lower limits
	expr <- as.data.frame(data_BRCA[[i]]) %>% 
		mutate(EXPR = (read_count + 1) / lib_size) %>%
		select(EXPR)
	pr <- as.data.frame(data_BRCA[[i]]) %>% 
		select(starts_with("pr", ignore.case = F))
	n_pr_cpg <- ncol(pr)
	gb <- as.data.frame(data_BRCA[[i]]) %>% 
		select(starts_with("gb", ignore.case = F))
	n_gb_cpg <- ncol(gb)
	
	# upper and lower limits
	min_p <- min(pr)-0.1
	max_p <- max(pr)+0.1
	min_gb <- min(gb)-0.1
	max_gb <- max(gb)+0.1
	min_e <- max(0,min(expr)-0.001*min(expr))
	max_e <- max(expr)+0.001*max(expr)
	
	# common data frame with observed expression point estimate
	df <- data_BRCA[[i]] %>%
		as.data.frame() %>%
		mutate(EXPR = (read_count+1) / lib_size) %>%
		mutate(EXPR = as.integer(cut(EXPR, breaks = seq(min_e,max_e,length.out = 1001), labels = c(1:1000)))) %>%
		mutate(PR_overall = NA) %>%
		mutate(GB_overall = NA) %>%
		mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_p,max_p,length.out = 101), labels = c(1:100)))), PR = starts_with("pr", ignore.case = F)) %>%
		mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_gb,max_gb,length.out = 101), labels = c(1:100)))), GB = starts_with("gb", ignore.case = F)) %>%
		select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count)
	
	# common data frame with unobserved expression
	df_wo_expr <- df %>% 
		mutate(EXPR = NA)
	
	# Base models
	varDim <- c(1000,rep(100, 2+n_pr_cpg+n_gb_cpg))
	facPot <- list(linregPotential(dim = c(1000, 100)), # PR_overall | EXPR
					linregPotential(dim = c(1000, 100)), # GB_overall | EXPR
					linregPotential(dim = c(100, 100), range1 = c(min_p, max_p), range2 = c(min_p,max_p), alpha = 1, beta = 0, var = (max_p-min_p)**2/25), # PR_i | PR
					linregPotential(dim = c(100, 100), range1 = c(min_gb, max_gb), range2 = c(min_gb,max_gb), alpha = 1, beta = 0, var = (max_gb-min_gb)**2/25)) # GB_i | GB
	facNbs <- c(list(c(1,2)), # PR_overall | EXPR
				list(c(1,3)), # GB_overall | EXPR
				lapply(4:(3+n_pr_cpg), FUN=function(i){c(2,i)}), # PR_i | PR_overall
				lapply((1+3+n_pr_cpg):(3+n_pr_cpg+n_gb_cpg), FUN=function(i){c(3,i)})) # GB_i | GB_overall
	potMap <- c(1, 2, rep(3, n_pr_cpg), rep(4, n_gb_cpg))
	optimFun <- list(linreg1 = linregOptimize(range1 = c(min_e,max_e), range2 = c(min_p,max_p)), 
		linreg2 = linregOptimize(range1 = c(min_e,max_e), range2 = c(min_gb,max_gb)), 
		fixedLink1 = fixedlinkOptimize(range1 = c(min_p,max_p), range2 = c(min_p,max_p), alpha = 1, beta = 0), 
		fixedLink2 = fixedlinkOptimize(range1 = c(min_gb,max_gb), range2 = c(min_gb,max_gb), alpha = 1, beta = 0))
	
	dfg_base <- dfg(varDim, facPot, facNbs, potMap, varNames = names(df))
	dfg_g1 <- dfg_base
	dfg_g2 <- dfg_base
	
	
	# dfg with beta prior
	cur_length <- length(varDim)
	varDimPrior <- varDim
	facPotPrior <- facPot
	facPotPrior[[5]] <- matrix(1, 1,1000)
	facNbsPrior <- facNbs
	facNbsPrior[[cur_length]] <- 1
	potMapPrior <- potMap
	potMapPrior[cur_length] <- 5
	dfg_prior_base <- dfg(varDimPrior, facPotPrior, facNbsPrior, potMapPrior, varNames = names(df))
	
	# begin x-fold here
	for (fold in 1:nfolds) { # iterate through nfolds folds
		# identify training and validation sets, assumes g1 and g2 sets are adjacent
		g1_xfold_training <- setdiff((1:length_g1),unlist(folds_g1[[fold]]))
		g2_xfold_training <- setdiff(length_g1+(1:length_g2),unlist(folds_g2[[fold]]))
		
		#################################################
		# G1 model
		#################################################
		#################################################
		# Learn betabinomial
		#################################################
		xini <- as.data.frame(data_BRCA[[i]][g1_xfold_training,1:2])
		ni <- as.integer(xini$lib_size)
		xi <- as.integer(xini$read_count+1)
		m0 <- mle2(minuslogl = loglik, start = list(alpha = 1, beta = 1), method = "L-BFGS-B", lower=c(alpha = 0.0001, beta = 0.0001))
		alpha <- coef(m0)['alpha']
		beta <- coef(m0)['beta']
		# Posteriors
		alphaPost <- alpha + xini$read_count
		betaPost <- beta + xini$lib_size - xini$read_count
		# Data list
		dataList <- list()
		dataList[[1]] <- lapply(1:length(g1_xfold_training), FUN=function(i){
			breaks <- seq(min_e, max_e, length.out = 1001)
			diff( pbeta(breaks, alphaPost[i], betaPost[i]))})
		
		# Train
		dfg_g1 <- train(data = df_wo_expr[g1_xfold_training,], 
						dataList = dataList,
						dfg = dfg_g1, 
						optim = c("linreg1", "linreg2", "fixedLink1", "fixedLink2"), 
						optimFun = optimFun, iter.max = 5000, threshold = 1e-5)
		# Model with beta prior
		dfg_g1_prior <- dfg_prior_base
		potentials(dfg_g1_prior) <- c(potentials(dfg_g1), list(betaPotential(dim=c(1,1000), range=c(min_e,max_e),alphas=alpha, betas=beta)))
		
		#################################################
		# G2 model
		#################################################
		# Learn betabinomial
		xini <- as.data.frame(data_BRCA[[i]][g2_xfold_training,1:2])
		ni <- as.integer(xini$lib_size)
		xi <- as.integer(xini$read_count+1)
		m0 <- mle2(minuslogl = loglik, start = list(alpha = 1, beta = 1), method = "L-BFGS-B", lower=c(alpha = 0.0001, beta = 0.0001))
		alpha <- coef(m0)['alpha']
		beta <- coef(m0)['beta']
		# Posteriors
		alphaPost <- alpha + xini$read_count
		betaPost <- beta + xini$lib_size - xini$read_count
		# Data list
		dataList <- list()
		dataList[[1]] <- lapply(1:length(g2_xfold_training), FUN=function(i){
			breaks <- seq(min_e, max_e, length.out = 1001)
			diff( pbeta(breaks, alphaPost[i], betaPost[i]))})
		
		# Train
		dfg_g2 <- train(data = df_wo_expr[g2_xfold_training,],
						dataList = dataList,
						dfg = dfg_g2,
						optim = c("linreg1", "linreg2", "fixedLink1", "fixedLink2"),
						optimFun = optimFun, iter.max = 5000, threshold = 1e-5)
		# Model with beta prior
		dfg_g2_prior <- dfg_prior_base
		potentials(dfg_g2_prior) <- c(potentials(dfg_g2), list(betaPotential(dim=c(1,1000), range=c(min_e,max_e),alphas=alpha, betas=beta)))
		
		#################################################
		# Evaluation and performance
		#################################################
		# training data scores
		df_train <- df_wo_expr[c(g1_xfold_training,g2_xfold_training),]
		
		likelihood1 <- likelihood(dfg = dfg_g2_prior, data = df_train, log = T)
		likelihood1[which(is.na(likelihood1))] <- vapply(likelihood1[which(is.na(likelihood1))], FUN= function(x) rnorm(n=1,mean=-500,sd=0.01), FUN.VALUE = 500)
		likelihood2 <- likelihood(dfg = dfg_g1_prior, data = df_train, log = T)
		likelihood2[which(is.na(likelihood2))] <- vapply(likelihood2[which(is.na(likelihood2))], FUN= function(x) rnorm(n=1,mean=-500,sd=0.01), FUN.VALUE = 500)
		scores_train <- likelihood1 - likelihood2
		
		# evaluation data scores
		df_eval <- df_wo_expr[-c(g1_xfold_training,g2_xfold_training),]
		
		likelihood1 <- likelihood(dfg = dfg_g2_prior, data = df_eval, log = T)
		likelihood1[which(is.na(likelihood1))] <- vapply(likelihood1[which(is.na(likelihood1))], FUN= function(x) rnorm(n=1,mean=-500,sd=0.01), FUN.VALUE = 500)
		likelihood2 <- likelihood(dfg = dfg_g1_prior, data = df_eval, log = T)
		likelihood2[which(is.na(likelihood2))] <- vapply(likelihood2[which(is.na(likelihood2))], FUN= function(x) rnorm(n=1,mean=-500,sd=0.01), FUN.VALUE = 500)
		scores_eval <- likelihood1 - likelihood2
		
		# Metrics of performance
		model$auc_training[fold] <- auc(predictor=scores_train,response=c(rep("G1",length(g1_xfold_training)),rep("G2",length(g2_xfold_training))))
		model$auc_evaluation[fold] <- auc(predictor=scores_eval,response=c(rep("G1",(length_g1-length(g1_xfold_training))),rep("G2",(length_g2-length(g2_xfold_training)))))
		model$kls[fold] <- kl(dfg_g1_prior,dfg_g2_prior)
		
		scores_xval[c(folds_g1[[fold]],folds_g2[[fold]])] <- scores_eval
	} #  x-validation ends here
	model$scores <- scores_xval
	# save model results
	eval(parse(text=paste('save(model, file="./',i,'_model.RData")',sep="")))
	cat(paste("done evaluating for ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
