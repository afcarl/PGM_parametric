#################################################
# prepare workspace
#################################################
args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1]) # first ID to process
end <- as.numeric(args[2]) # last consequitve ID to process

# data <- args[1]
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

#IDs <- c("RAG1AP1","CPA1","NEK2","RNASEH2A","LOC148145","TMEM63B","TIMM17A","PLK1","RABIF","PTF1A")

# setup variables to hold model and performance metrics across folds
#facPotNormals <- NULL
#facPotTumours <- NULL
#models <- NULL # a list to hold models and performance measures across folds

for (i in beg:end) {#2:length(IDs)) { # iterate through given jobs
	#i <- which(names(data_BRCA) %in% IDs[j]) # find the index to work with
	cat(paste("doing ",i,"\n",sep=""))
	ptm <- proc.time()[3]
	
	#models[[j]] <- list()
	#models[[j]]$gene_name <- IDs[j]
	model <- list()
	#model$gene_name <- IDs[j]
	scores_xval <- vector(length=nrow(data_BRCA[[i]]))
	# calculate common discretization upper and lower limits
	expr <- as.data.frame(data_BRCA[[i]]) %>% 
		mutate(EXPR = read_count / lib_size) %>%
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
	
	# begin x-fold here
	for (fold in 1:8) { # iterate through 8 folds
		# identify training and validation sets
		an_xfold_training <- setdiff((1:82),unlist(folds_an[[fold]]))
		t_xfold_training <- setdiff(82+(1:730),unlist(folds_t[[fold]]))
		
		#################################################
		# Normals model
		#################################################
			data_normals_train <- as.data.frame(data_BRCA[[i]][an_xfold_training,])
			#################################################
			# Learn betabinomial
			#################################################
			ni <- as.integer(data_normals_train$lib_size)
			xi <- as.integer(data_normals_train$read_count+1)
			
			m0 <- mle2(minuslogl = loglik, start = list(alpha = 1, beta = 1), method = "L-BFGS-B", lower=c(alpha = 0.0001, beta = 0.0001))
			alpha <- coef(m0)['alpha']
			beta <- coef(m0)['beta']
			
			# Posteriors
			alphaPost <- alpha + data_normals_train$read_count
			betaPost <- beta + data_normals_train$lib_size - data_normals_train$read_count
			
			#################################################
			# Data-massage
			#################################################
			
			df_normals <- data_normals_train %>% 
			  mutate(EXPR = NA) %>%
			  mutate(PR_overall = NA) %>%
			  mutate(GB_overall = NA) %>%
			  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_p,max_p,length.out = 101), labels = c(1:100)))), PR = starts_with("pr", ignore.case = F)) %>%
			  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_gb,max_gb,length.out = 101), labels = c(1:100)))), GB = starts_with("gb", ignore.case = F)) %>%
			  select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count)
			
			#################################################
			# Building models and training
			#################################################
			
			# Build Model
			varDim <- rep(100, 3+n_pr_cpg+n_gb_cpg)
			if (fold == 1) {
			facPot <- list(linregPotential(dim = c(100, 100)), # PR_overall | EXPR
						linregPotential(dim = c(100, 100)), # GB_overall | EXPR
						linregPotential(dim = c(100, 100), range1 = c(min_p, max_p), range2 = c(min_p,max_p), alpha = 1, beta = 0, var = 0.14**2), # PR_i | PR
						linregPotential(dim = c(100, 100), range1 = c(min_gb, max_gb), range2 = c(min_gb,max_gb), alpha = 1, beta = 0, var = 0.14**2)) # GB_i | GB
			} else facPot <- prev_normals_potentials
			
			facNbs <- c(list(c(1,2)), # PR_overall | EXPR
						list(c(1,3)), # GB_overall | EXPR
						lapply(4:(3+n_pr_cpg), FUN=function(i){c(2,i)}), # PR_i | PR_overall
						lapply((1+3+n_pr_cpg):(3+n_pr_cpg+n_gb_cpg), FUN=function(i){c(3,i)}) # GB_i | GB_overall
			)
			
			potMap <- c(1, 2, rep(3, n_pr_cpg), rep(4, n_gb_cpg))
			dfg_normals <- dfg(varDim, facPot, facNbs, potMap, varNames = names(df_normals))
			
			optimFun <- list(linreg1 = linregOptimize(range1 = c(min_e,max_e), range2 = c(min_p,max_p)),
							 linreg2 = linregOptimize(range1 = c(min_e,max_e), range2 = c(min_gb,max_gb)))
			
			# Data 
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
				  optimFun = optimFun, iter.max = 5000, threshold = 1e-5)
			
			prev_normals_potentials <- potentials(dfg_normals)
			# Hack Expression prior
			cur_length <- length(dfg_normals$facNbs)+1
			facNbs[[cur_length]] <- 1
			potMap[cur_length] <- 5
			facPot <- potentials(dfg_normals)
			facPot[[5]] <- betaPotential(range=c(min_e,max_e),alphas=alpha, betas=beta)
			dfg_normals <- dfg(varDim, facPot, facNbs, potMap, varNames = names(df_normals))
			#facPotNormals[[j]] <- list()
			#facPotNormals[[j]][[fold]] <- potentials(dfg_normals)
			
		#################################################
		# Tumours model
		#################################################
			data_tumours_train <- as.data.frame(data_BRCA[[i]][t_xfold_training,])
			#################################################
			# Learn betabinomial
			#################################################
			ni <- as.integer(data_tumours_train$lib_size)
			xi <- as.integer(data_tumours_train$read_count+1)
			
			m0 <- mle2(minuslogl = loglik, start = list(alpha = 1, beta = 1), method = "L-BFGS-B", lower=c(alpha = 0.0001, beta = 0.0001))
			alpha <- coef(m0)['alpha']
			beta <- coef(m0)['beta']
			
			# Posteriors
			alphaPost <- alpha + data_tumours_train$read_count
			betaPost <- beta + data_tumours_train$lib_size - data_tumours_train$read_count
			
			#################################################
			# Data-massage
			#################################################
			
			df_tumours <- data_tumours_train %>% 
			  mutate(EXPR = NA) %>%
			  mutate(PR_overall = NA) %>%
			  mutate(GB_overall = NA) %>%
			  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_p,max_p,length.out = 101), labels = c(1:100)))), PR = starts_with("pr", ignore.case = F)) %>%
			  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_gb,max_gb,length.out = 101), labels = c(1:100)))), GB = starts_with("gb", ignore.case = F)) %>%
			  select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count)
			
			#################################################
			# Building models and training
			#################################################
			
			# Build Model
			if (fold != 1) facPot <- prev_tumours_potentials
			facNbs[[cur_length]] <- NULL
			potMap <- potMap[-cur_length]
			facPot[[5]] <- NULL
			dfg_tumours <- dfg(varDim, facPot, facNbs, potMap, varNames = names(df_tumours))
			
			# Data 
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
				  optimFun = optimFun, iter.max = 5000, threshold = 1e-5)
			
			prev_tumours_potentials <- potentials(dfg_tumours)
			# Hack Expression prior
			cur_length <- length(dfg_tumours$facNbs)+1
			facNbs[[cur_length]] <- 1
			potMap[cur_length] <- 5
			facPot <- potentials(dfg_tumours)
			facPot[[5]] <- betaPotential(range=c(min_e,max_e),alphas=alpha, betas=beta)
			dfg_tumours <- dfg(varDim, facPot, facNbs, potMap, varNames = names(df_tumours))
			#facPotTumours[[j]] <- list()
			#facPotTumours[[j]][[fold]] <- potentials(dfg_tumours)
			
		#################################################
		# Evaluation and performance
		#################################################
			# training data performance
			data_train <- as.data.frame(data_BRCA[[i]][c(an_xfold_training,t_xfold_training),])
			df_train <- data_train %>% 
				  mutate(EXPR = (read_count+1) / lib_size) %>%
				  mutate(EXPR = as.integer(cut(EXPR, breaks = seq(min_e,max_e,length.out = 101), labels = c(1:100)))) %>%
				  mutate(PR_overall = NA) %>%
				  mutate(GB_overall = NA) %>%
				  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_p,max_p,length.out = 101), labels = c(1:100)))), PR = starts_with("pr", ignore.case = F)) %>%
				  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_gb,max_gb,length.out = 101), labels = c(1:100)))), GB = starts_with("gb", ignore.case = F)) %>%
				  select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count)
			
			likelihood1 <- likelihood(dfg = dfg_tumours, data = df_train, log = T)
			likelihood1[which(is.na(likelihood1))] <- vapply(likelihood1[which(is.na(likelihood1))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
			likelihood2 <- likelihood(dfg = dfg_normals, data = df_train, log = T)
			likelihood2[which(is.na(likelihood2))] <- vapply(likelihood2[which(is.na(likelihood2))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
			scores_train <- likelihood1 - likelihood2
			
			# evaluation data performance
			data_eval <- as.data.frame(data_BRCA[[i]][-c(an_xfold_training,t_xfold_training),])
			df_eval <- data_eval %>% 
				  mutate(EXPR = (read_count+1) / lib_size) %>%
				  mutate(EXPR = as.integer(cut(EXPR, breaks = seq(min_e,max_e,length.out = 101), labels = c(1:100)))) %>%
				  mutate(PR_overall = NA) %>%
				  mutate(GB_overall = NA) %>%
				  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_p,max_p,length.out = 101), labels = c(1:100)))), PR = starts_with("pr", ignore.case = F)) %>%
				  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_gb,max_gb,length.out = 101), labels = c(1:100)))), GB = starts_with("gb", ignore.case = F)) %>%
				  select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count)
			
			#scores_eval <- likelihood(dfg = dfg_tumours, data = df_eval, log = T) - likelihood(dfg = dfg_normals, data = df_eval, log = T)
			likelihood1 <- likelihood(dfg = dfg_tumours, data = df_eval, log = T)
			likelihood1[which(is.na(likelihood1))] <- vapply(likelihood1[which(is.na(likelihood1))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
			likelihood2 <- likelihood(dfg = dfg_normals, data = df_eval, log = T)
			likelihood2[which(is.na(likelihood2))] <- vapply(likelihood2[which(is.na(likelihood2))], FUN= function(x) rnorm(n=1,mean=-500), FUN.VALUE = rnorm(n=1,mean=-500))
			scores_eval <- likelihood1 - likelihood2
			
			# Metrics
			# models[[j]]$auc_training[fold] <- auc(predictor=scores_train,response=c(rep("N",length(an_xfold_training)),rep("T",length(t_xfold_training))))
			# models[[j]]$auc_evaluation[fold] <- auc(predictor=scores_eval,response=c(rep("N",(82-length(an_xfold_training))),rep("T",(730-length(t_xfold_training)))))
			# models[[j]]$kls[fold] <- kl(dfg_normals,dfg_tumours)
			model$auc_training[fold] <- auc(predictor=scores_train,response=c(rep("N",length(an_xfold_training)),rep("T",length(t_xfold_training))))
			model$auc_evaluation[fold] <- auc(predictor=scores_eval,response=c(rep("N",(82-length(an_xfold_training))),rep("T",(730-length(t_xfold_training)))))
			model$kls[fold] <- kl(dfg_normals,dfg_tumours)
			
			scores_xval[c(folds_an[[fold]],folds_t[[fold]])] <- scores_eval
	}
	# end x-fold here
	#models[[j]]$scores <- scores_xval
	model$scores <- scores_xval
	eval(parse(text=paste('save(model, file="./',i,'_model.RData")',sep="")))
	cat(paste("done evaluating for ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
