

i <- 
	###########################################################
	################### prepare data frame ###################
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	ncol = length(IDs_promoter) + length(IDs_body) + 2
	nrow = length(c(G1_train,G2_train))
	temp <- matrix(nrow=nrow,ncol=ncol)
	temp[,1] <- unlist(factors_ls[c(G1_train,G2_train)])*10^6
	temp[,2] <- unlist(counts[workingList_BRCA[i],c(G1_train,G2_train)])
	temp[,3:(2+length(IDs_promoter))] <- t(mmatrix_pc[IDs_promoter,c(G1_train,G2_train)])
	temp[,(3+length(IDs_promoter)):(2+length(IDs_promoter)+length(IDs_body))] <- t(mmatrix_pc[IDs_body,c(G1_train,G2_train)])
	colnames(temp) <- c("lib_size","read_count",paste("pr_",1:length(IDs_promoter),sep=""),paste("gb_",1:length(IDs_body),sep=""))
	rownames(temp) <- c(G1_train,G2_train)
	temp <- as.data.frame(temp)
	
	
	min_p <- min(mmatrix_pc[IDs_promoter,c(G1,G2)])-0.1
	max_p <- max(mmatrix_pc[IDs_promoter,c(G1,G2)])+0.1
	min_gb <- min(mmatrix_pc[IDs_body,c(G1,G2)])-0.1
	max_gb <- max(mmatrix_pc[IDs_body,c(G1,G2)])+0.1
	rpms <- counts[workingList_BRCA[i],c(G1,G2)]/(factors_ls[c(G1,G2)]*10^6)
	min_e <- max(0,min(rpms)-0.01*min(rpms))
	max_e <- max(rpms)+0.01*max(rpms)
	
	# Process data
	df_full <- temp %>% 
	  mutate(EXPR = read_count / lib_size) %>%
	  mutate(EXPR = as.integer(cut(EXPR, breaks = seq(min_e,max_e,length.out = 101), labels = c(1:100)))) %>%
	  mutate(H1 = NA) %>%
	  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_p,max_p,length.out = 101), labels = c(1:100)))), PR = starts_with("pr")) %>%
	  mutate(PR = NA) %>%
	  mutate_each(funs(bin = as.integer(cut(., breaks = seq(min_gb,max_gb,length.out = 101), labels = c(1:100)))), GB = starts_with("gb")) %>%
	  mutate(GB = NA) %>%
	  select(H1, EXPR, starts_with("PR", ignore.case = F), starts_with("GB", ignore.case = F))
	
	
	# Build Model
	nVar <- 3 + length(IDs_body) + length(IDs_promoter)
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
					 linreG1_train = linregOptimize(range1 = c(0,3e-5), range2 = c(-3,3)),
					 linreG2_train = linregOptimize(range1 = c(0,3e-5), range2 = c(0,4)))
	train(df_full, dfg_full, optim = c("row", "custombeta", "linreG1_train", "linreG2_train", "noopt", "noopt"), optimFun = optimFun, iter.max = 500, threshold = 1e-1)
