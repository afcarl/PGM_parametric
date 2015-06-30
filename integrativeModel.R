# the following objects are expected by this script: workingList_BRCA, mmatrix_pc, counts, G1, G2, factors_ls, TSS1500Ind, TSS200Ind, UTR5Ind, EXON1Ind, GENEBODYInd, UTR3Ind

args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1]) # first ID to process
end <- as.numeric(args[2]) # last consequitve ID to process

# data <- args[1]
load("./essentials_allBRCA.RData") # data to work on
IDs_length <- nchar(G1[1])+2
G1 <- G1[c(1,2,5,6,7,8,10,14,18,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,37,38,39,41,42,45,46,47,48,49,51,53,54,56,57,58,59,61,62,64,65,67,68,70,71,75,76,77,78,80,81)]
G2 <- G2[c(1,2,3,5,7,9,10,12,13,14,15,16,19,20,21,22,23,24,25,26,27,30,31,32,34,35,36,37,38,39,41,42,44,46,47,48,52,53,54,55,58,59,61,62,64,65,66,67,69,72,74,76,77,79,81,82,83,87,88,89,90,91,92,93,94,95,97,98,99,101,102,103,104,106,108,112,113,115,116,117,119,120,121,122,125,126,128,129,130,131,132,133,135,136,137,139,141,142,143,145,146,147,148,149,152,155,156,157,158,159,160,164,165,166,168,169,170,172,173,174,176,178,180,181,182,183,184,187,189,190,191,192,193,194,195,197,198,201,202,203,204,205,208,210,212,214,215,217,219,220,221,223,224,225,228,229,232,237,238,239,240,241,244,245,246,247,248,250,251,252,253,254,255,256,259,260,263,264,265,266,268,269,270,272,277,278,279,280,282,283,284,287,289,290,291,292,295,296,297,298,299,300,301,302,303,305,308,309,310,311,313,314,316,317,318,319,320,321,322,324,326,330,332,333,334,335,336,337,339,340,342,345,349,351,352,353,354,355,356,357,358,360,361,362,363,364,365,367,368,369,370,371,372,373,376,377,379,380,382,383,387,390,391,393,397,399,400,402,403,404,405,406,407,409,411,412,413,414,415,417,418,420,422,423,424,426,427,428,429,431,432,433,435,436,437,438,439,442,445,447,448,450,451,452,453,454,455,456,459,461,463,467,468,471,472,473,474,475,477,479,481,482,484,485,486,488,489,490,491,493,494,495,497,500,502,503,505,507,509,510,511,512,513,514,517,518,519,520,521,523,524,525,527,528,529,530,531,532,534,535,536,537,539,540,541,545,546,551,552,553,554,555,556,557,559,561,567,569,570,575,579,580,582,583,584,585,586,589,591,594,596,598,600,602,603,604,605,606,609,610,611,612,613,614,616,618,619,620,621,622,623,625,628,630,631,632,633,634,635,636,638,639,640,641,642,643,645,646,647,648,649,650,652,653,654,658,659,660,661,662,663,664,666,667,669,671,672,673,674,676,677,679,680,681,682,684,685,686,688,689,690,691,696,698,702,703,705,706,707,708,709,711,712,713,714,715,716,718,720,721,722,723,724,726,728,729,730)]

# configure here
res_pr <- 25 # number of bins for promoter methylation
res_gb <- 25 # number of bins for gene body methylation
res_expr <- 25 # number of bins for gene expression
nruns <- 100 # number of permutations to perform for null distribution calculation


for (i in beg:end){
	cat(paste("doing ",i,"\n",sep=""))
	ptm <- proc.time()[3]
	
	###########################################################
	############### identify constitutive CpGs ################
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	ncol = length(IDs_promoter) + length(IDs_body) + 1
	
	###########################################################
	############### define the binning scheme  ################
	all_labels_pr <- as.character(seq(1,res_pr,1))
	all_labels_gb <- as.character(seq(1,res_gb,1))
	all_labels_expr <- as.character(seq(1,res_expr,1))
	promoter_G2 <- matrix(ncol=res_pr,nrow=length(G2))
	promoter_G1 <- matrix(ncol=res_pr,nrow=length(G1))
	body_G2 <- matrix(ncol=res_gb,nrow=length(G2))
	body_G1 <- matrix(ncol=res_gb,nrow=length(G1))
	expr_G2 <- matrix(ncol=res_expr,nrow=length(G2))
	expr_G1 <- matrix(ncol=res_expr,nrow=length(G1))
	pseudo_counts_pr <- matrix(ncol=res_pr,nrow=res_expr,data=rep(1,res_expr*res_pr))
	pseudo_counts_gb <- matrix(ncol=res_gb,nrow=res_expr,data=rep(1,res_expr*res_gb))
	
	# gene body
	x <- mmatrix_pc[IDs_body,c(G2,G1)]
	density <- density(x,bw=0.14,from=-7,to=7,n=2801,na.rm=TRUE)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_gb-1
	for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksBODY <- sort(c(-7.01,breaks,7.01))
	
	# promoter
	x <- mmatrix_pc[IDs_promoter,c(G2,G1)]
	density <- density(x,bw=0.14,from=-7,to=7,n=2801,na.rm=TRUE)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_pr-1
	for (j in 1:noBreaks) {breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
	
	# expression
	temp <- counts[workingList_BRCA[i],c(G2,G1)]
	tempAN <- matrix(ncol=2)
	colnames(tempAN) <- c("cpm","density")
	for (j in 1:length(temp)) {
		lambda <- as.numeric(temp[j])
		if (lambda > 0) X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),0)),round(lambda+(4*lambda*lambda^(-1/2)))) else X <- 0:1
		current <- factors_ls[c(G2,G1)[j]]
		tempAN <- rbind(tempAN,cbind(X/current,dpois(X,lambda=lambda)))
	}
	tempAN <- as.data.frame(tempAN[-1,],)
	tempAN <- tempAN[order(tempAN$cpm),]
	tempAN[,3] <- cumsum(tempAN[,2])
	if (max(tempAN[,3]) != 0) tempAN[,3] <- tempAN[,3]/max(tempAN[,3])
	breaks <- NULL
	noBreaks <- res_expr-1
	for (j in 1:noBreaks) { breaks <- c(breaks, tempAN[which(tempAN[,3] >= j*(1/(1+noBreaks))),1][1])}
	max <- max(temp)
	if (max > 0) max_boundary <- 20+max+(4*max*max^(-1/2)) else max_boundary <- 20
	breaks <- unique(breaks)
	if (length(breaks) == 1) breaks <- seq(max(tempAN[,1]),max_boundary-max(tempAN[,1]),(max_boundary-max(tempAN[,1]))/24) else if (length(breaks) != 24){
		breaks <- breaks[-which(breaks %in% 0)]
		breaks <- c(breaks,seq(max(tempAN[,1])+max(breaks),max_boundary-max(tempAN[,1]),(max_boundary-max(tempAN[,1])-max(breaks))/(24-length(breaks))))
	}
	breaksEXPRESSION <- c(0,breaks,max_boundary)
	########################################################
	#############  model specification  here ###############
	########################################################
	####################### variables map ######################
	varDim <- c(res_expr, res_pr, res_gb)
	########################################################
	############## initial factor potentials ###############
	
	########################################################
	##################### factor graph #####################
	
	###########################################################################
	############################## G1 model ###################################
	## full G1 model developed from here, to obtain likelihoods of G1 ########
	
	# generate FacData for full set of G1
	tempS_G1 <- matrix(ncol=ncol,nrow=length(G1))
	for (current_sample in 1:length(G1)) {
		# expression
		read_count <- trunc(as.numeric(counts[workingList_BRCA[i],G1[current_sample]]))
		lambdas <- breaksEXPRESSION * factors_ls[G1[current_sample]]
		frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
		for (freq in 1:res_expr) {
			frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count,stop.on.error=FALSE)[1]
		}
		frequencies_expr <- unlist(frequencies_expr)
		if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
		frequencies_expr <- frequencies_expr + epsilon_e_G1
		frequencies_expr <- frequencies_expr/sum(frequencies_expr)
		
		# gene body
		cpg_list_gb <- NULL
		for (cpg in 1:length(geneBodyVars)) {
			miu <- mmatrix_pc[IDs_body[cpg],G1[current_sample]]
			if (!is.na(miu)) {
				frequencies_gb <- rep(0,res_gb)
				for (freq in 1:res_gb) {
					frequencies_gb[freq] <- integrate(integrand_m,lower=breaksBODY[freq],upper=breaksBODY[freq+1],mean=miu)$value
				}
				frequencies_gb <- unlist(frequencies_gb) + epsilon_gb_G1
				frequencies_gb <- frequencies_gb/sum(frequencies_gb)
				cpg_list_gb[[cpg]] <- frequencies_gb
			} else cpg_list_gb[[cpg]] <- rep(1/res_gb,res_gb)
		}
		
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix_pc[IDs_promoter[cpg],G1[current_sample]]
			if (!is.na(miu)) {
				frequencies_pr <- rep(0,res_pr)
				for (freq in 1:res_pr) {
					frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
				}
				frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_G1
				frequencies_pr <- frequencies_pr/sum(frequencies_pr)
				cpg_list_pr[[cpg]] <- frequencies_pr
			} else cpg_list_pr[[cpg]] <- rep(1/res_pr,res_pr)
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		tempS_formated[1,1] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="",collapse="")
		cur_ncol <- 1
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="",collapse="")
			cur_ncol <- cur_ncol + 1
		}
		for (element in 1:length(cpg_list_gb)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_gb,']((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep="",collapse="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G1[current_sample,] <- tempS_formated
		
		#start precomputing correct initialization of parameters
		promoter_G1[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
		body_G1[current_sample,] <- apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean))
		expr_G1[current_sample,] <- frequencies_expr
	}
	
	# build and the full model with G1 samples
	# precompute correct initialization of parameters for G1-only model
	prior_pr <- apply(promoter_G1,2,mean)
	prior_gb <- apply(body_G1,2,mean)
	if (smooth_1d > 0) {
		prior_expr <- kernsm(apply(expr_G1,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_G1,2,mean)
	
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_gb,collapse=",")
	geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(body_G1,expr_G1,smooth_h=smooth_gb)
	expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(promoter_G1,expr_G1,smooth_h=smooth_pr)
	expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
	
	potentials <- file(paste("./",i,"/G1_model/all/factorPotentials.txt",sep=""),"w")
	cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G1
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
	rownames(tempFac) <- G1
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G1_model/all/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full G1 model with G1 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/G1_model/all/G1_FacData.tab',sep=""))
	G1_G1model_mlogliks <- as.numeric(substring(string[-1],IDs_length))
	###########################################################################
	############################## G2 model ###################################
	### full G2 model developed from here, to obtain likelihoods of G2 #######

	# generate FacData for full set of G2
	tempS_G2 <- matrix(ncol=ncol,nrow=length(G2))
	for (current_sample in 1:length(G2)) {
		# expression
		read_count <- trunc(as.numeric(counts[workingList_BRCA[i],G2[current_sample]]))
		lambdas <- breaksEXPRESSION * factors_ls[G2[current_sample]]
		frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
		for (freq in 1:res_expr) {
			frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count,stop.on.error=FALSE)[1]
		}
		frequencies_expr <- unlist(frequencies_expr)
		if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
		frequencies_expr <- frequencies_expr + epsilon_e_G2
		frequencies_expr <- frequencies_expr/sum(frequencies_expr)
		
		# gene body
		cpg_list_gb <- NULL
		for (cpg in 1:length(geneBodyVars)) {
			miu <- mmatrix_pc[IDs_body[cpg],G2[current_sample]]
			if (!is.na(miu)) {
				frequencies_gb <- rep(0,res_gb)
				for (freq in 1:res_gb) {
					frequencies_gb[freq] <- integrate(integrand_m,lower=breaksBODY[freq],upper=breaksBODY[freq+1],mean=miu)$value
				}
				frequencies_gb <- unlist(frequencies_gb) + epsilon_gb_G2
				frequencies_gb <- frequencies_gb/sum(frequencies_gb)
				cpg_list_gb[[cpg]] <- frequencies_gb
			} else cpg_list_gb[[cpg]] <- rep(1/res_gb,res_gb)
		}
		
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix_pc[IDs_promoter[cpg],G2[current_sample]]
			if (!is.na(miu)) {
				frequencies_pr <- rep(0,res_pr)
				for (freq in 1:res_pr) {
					frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
				}
				frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_G2
				frequencies_pr <- frequencies_pr/sum(frequencies_pr)
				cpg_list_pr[[cpg]] <- frequencies_pr
			} else cpg_list_pr[[cpg]] <- rep(1/res_pr,res_pr)
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		tempS_formated[1,1] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="",collapse="")
		cur_ncol <- 1
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="",collapse="")
			cur_ncol <- cur_ncol + 1
		}
		for (element in 1:length(cpg_list_gb)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_gb,']((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep="",collapse="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G2[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_G2[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
		body_G2[current_sample,] <- apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean))
		expr_G2[current_sample,] <- frequencies_expr
	}
	
	# build and the full model with G2 samples
	# precompute correct initialization of parameters for G2-only model
	prior_pr <- apply(promoter_G2,2,mean)
	prior_gb <- apply(body_G2,2,mean)
	if (smooth_1d > 0) {
		prior_expr <- kernsm(apply(expr_G2,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_G2,2,mean)
	
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_gb,collapse=",")
	geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(body_G2,expr_G2,smooth_h=smooth_gb)
	expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(promoter_G2,expr_G2,smooth_h=smooth_pr)
	expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
	
	potentials <- file(paste("./",i,"/G2_model/all/factorPotentials.txt",sep=""),"w")
	cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G2
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
	rownames(tempFac) <- G2
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G2_model/all/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with G2 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/G2_model/all/G2_FacData.tab',sep=""))
	G2_G2model_mlogliks <- as.numeric(substring(string[-1],IDs_length))
	##########################################################################
	
	###########################################################################
	######################## All data model ###################################
	## Full model developed from here, to obtain likelihoods of G2 and G1 ####
	# precompute correct initialization of parameters for joint G1-G2 model
	promoter_all <- rbind(promoter_G2,promoter_G1)
	body_all <- rbind(body_G2,body_G1)
	expr_all <- rbind(expr_G2,expr_G1)
	
	prior_pr <- apply(promoter_all,2,mean)
	prior_gb <- apply(body_all,2,mean)
	if (smooth_1d > 0) {
		prior_expr <- kernsm(apply(expr_all,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_all,2,mean)
	
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_gb,collapse=",")
	geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(body_all,expr_all,smooth_h=smooth_gb)
	expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
	
	result <- tensor_product(promoter_all,expr_all,smooth_h=smooth_pr)
	expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
	
	potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
	cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
	close(potentials)
	
	tempFac <- rbind(tempS_G2,tempS_G1)
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
	rownames(tempFac) <- c(G2,G1)
	eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with T and AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/full_model/ -l -n - ./',i,'/full_model/full_VarData.tab ./',i,'/full_model/full_FacData.tab',sep=""))
	allData_jointModel_mlogliks <- as.numeric(substring(string[-1],IDs_length))
	###########################################################################################
	
	###########################################################################
	######################## D calculation ###################################
	D <- 2*(sum(allData_jointModel_mlogliks) - (sum(G1_G1model_mlogliks)+sum(G2_G2model_mlogliks)))
	###########################################################################
	################# P val calculation using null distr. #####################
	Ds <- vector(length=nruns,mode="numeric")
	for (run in 1:nruns) {
		cur <- sample(x=1:(length(G2)+length(G1)),size=length(G2),replace=FALSE)
		# G2
		tempFac_G2 <- tempFac[cur,]
		rownames(tempFac_G2) <- G2
		eval(parse(text = paste('write.table(', paste('tempFac_G2,file = "./',i,'/null/G2_model/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- apply(promoter_all[cur,],2,mean)
		prior_gb <- apply(body_all[cur,],2,mean)
		if (smooth_1d > 0) {
			prior_expr <- kernsm(apply(expr_all[cur,],2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_all[cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_all[cur,],expr_all[cur,],smooth_h=smooth_gb)
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_all[cur,],expr_all[cur,],smooth_h=smooth_pr)
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/null/G2_model/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G2_model/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/null/G2_model/G2_FacData.tab',sep=""))
		G2_G2model_mlogliks <- as.numeric(substring(string[-1],IDs_length))
		
		# G1
		tempFac_G1 <- tempFac[-cur,]
		rownames(tempFac_G1) <- G1
		eval(parse(text = paste('write.table(', paste('tempFac_G1,file = "./',i,'/null/G1_model/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- apply(promoter_all[-cur,],2,mean)
		prior_gb <- apply(body_all[-cur,],2,mean)
		if (smooth_1d > 0) {
			prior_expr <- kernsm(apply(expr_all[-cur,],2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_all[-cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_all[-cur,],expr_all[-cur,],smooth_h=smooth_gb)
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_all[-cur,],expr_all[-cur,],smooth_h=smooth_pr)
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/null/G1_model/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G1_model/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/null/G1_model/G1_FacData.tab',sep=""))
		G1_G1model_mlogliks <- as.numeric(substring(string[-1],IDs_length))
		
		Ds[run] <- 2*(sum(allData_jointModel_mlogliks) - (sum(G1_G1model_mlogliks)+sum(G2_G2model_mlogliks)))
	}
	if (sd(Ds) != 0 & D > 0.1) zscore <- (D - mean(Ds)) / sd(Ds) else zscore <- -6
	pval_zscore <- pnorm(zscore,lower.tail=FALSE)
	###########################################################################################
	
	eval(parse(text=paste('write.table(x=t(c(pval_zscore,D,mean(Ds),sd(Ds),zscore)), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	eval(parse(text=paste('write.table(x=t(breaksEXPRESSION), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	eval(parse(text=paste('write.table(x=t(breaksBODY), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	eval(parse(text=paste('write.table(x=t(breaksPROMOTER), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	system(intern=TRUE,command=paste('tar cf ',i,'.tar ',i,sep=""))
	cat(paste("done ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
