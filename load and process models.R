readModel <- function(x) {eval(parse(text=paste('load("',x,'_model.RData",.GlobalEnv)',sep="")))}
loadModels <- function(x) {if (file.exists(paste(x,"_model.RData",sep=""))) readModel(x) else model <- list()}
models <- list()
for (i in 1:17728) {
	loadModels(i)
	models[[i]] <- model
}
# tumour vs normal
xfold_kls <- matrix(nrow=17728,ncol=8)
for (i in 1:17728) xfold_kls[i,] <- models[[i]]$kls
xfold_kls[is.na(xfold_kls)] <- 0
xfold_kls[xfold_kls == Inf] <- rnorm(n=length(xfold_kls[xfold_kls == Inf]),mean=125) 
workingList_BRCA <- names(data_BRCA)
top10 <- apply(-xfold_kls,2,order)
top10 <- top10[1:10,]
top10_kls <- -apply(-xfold_kls,2,sort)
top10_kls <- top10_kls[1:10,]
top10_names <- as.data.frame(top10)
for (i in 1:14) top10_names[,i] <- workingList_BRCA[top10[,i]]
ranks <- apply(-xfold_kls,2,rank)
mean_ranks <- apply(ranks,1,mean)
top10_meanRank <- (order(mean_ranks))[1:10]
workingList_BRCA[top10_meanRank]
cbind(ranks[top10_meanRank,], workingList_BRCA[top10_meanRank], apply(xfold_kls[top10_meanRank,],2,mean))

# progressing
xfold_kls <- matrix(nrow=17728,ncol=14)
for (i in 1:17728) xfold_kls[i,] <- models[[i]]$kls
xfold_kls[is.na(xfold_kls)] <- 0
xfold_kls[xfold_kls == Inf] <- rnorm(n=length(xfold_kls[xfold_kls == Inf]),mean=125) 
workingList_BRCA <- names(data_BRCA_progressing)
top10 <- apply(-xfold_kls,2,order)
top10 <- top10[1:10,]
top10_kls <- -apply(-xfold_kls,2,sort)
top10_kls <- top10_kls[1:10,]
top10_names <- as.data.frame(top10)
for (i in 1:14) top10_names[,i] <- workingList_BRCA[top10[,i]]
ranks <- apply(-xfold_kls,2,rank)
mean_ranks <- apply(ranks,1,mean)
top10_meanRank <- (order(mean_ranks))[1:10]
workingList_BRCA[top10_meanRank]
cbind(ranks[top10_meanRank,], workingList_BRCA[top10_meanRank], apply(xfold_kls[top10_meanRank,],2,mean))


# performance
library(pROC)
top100 <- apply(-xfold_kls,2,order)
top100 <- top100[1:100,]
xfold_scores <- matrix(nrow=length(models[[1]]$scores), ncol=100)
for (i in 1:100) {
	scores_temp <- vector(length=length(models[[1]]$scores), mode="numeric")
	for (fold in 1:length(folds_g1)) {
		scores_temp[c(folds_g1[[fold]], folds_g2[[fold]])] <- models[[top100[i,fold]]]$scores[c(folds_g1[[fold]], folds_g2[[fold]])]
	}
	xfold_scores[,i] <- scores_temp
}
aucs_single <- vector(length=10,mode="numeric")
for (i in 1:100) {
	aucs_single[i] <- auc(predictor=xfold_scores[,i],response=c(rep("G1",length(unlist(folds_g1))),rep("G2",length(unlist(folds_g2)))))
}

aucs_running_combinations <- vector(length=10,mode="numeric")
for (i in 2:100) {
	aucs_running_combinations[i] <- auc(predictor=apply(xfold_scores[,1:i],1,sum),response=c(rep("G1",length(unlist(folds_g1))),rep("G2",length(unlist(folds_g2)))))
}

# performance using training AUC as selection criteria
xfold_aucs <- matrix(nrow=17728,ncol=length(folds_g1))
for (i in 1:17728) xfold_aucs[i,] <- models[[i]]$auc_training
top <- apply(-xfold_aucs,2,order)
top <- top[1:20,]
xfold_scores <- matrix(nrow=length(models[[1]]$scores), ncol=20)
for (i in 1:20) {
	scores_temp <- vector(length=length(models[[1]]$scores), mode="numeric")
	for (fold in 1:length(folds_g1)) {
		scores_temp[c(folds_g1[[fold]], folds_g2[[fold]])] <- models[[top[i,fold]]]$scores[c(folds_g1[[fold]], folds_g2[[fold]])]
	}
	xfold_scores[,i] <- scores_temp
}
aucs_single <- vector(length=10,mode="numeric")
for (i in 1:20) {
	aucs_single[i] <- auc(predictor=xfold_scores[,i],response=c(rep("G1",length(unlist(folds_g1))),rep("G2",length(unlist(folds_g2)))))
}

aucs_running_combinations <- vector(length=10,mode="numeric")
for (i in 2:20) {
	aucs_running_combinations[i] <- auc(predictor=apply(xfold_scores[,1:i],1,sum),response=c(rep("G1",length(unlist(folds_g1))),rep("G2",length(unlist(folds_g2)))))
}

top <- apply(-xfold_aucs,2,order)
top_names <- as.data.frame(top)
for (i in 1:ncol(top_names)) top_names[,i] <- workingList_BRCA[top[,i]]

# plotting
library(ggplot2)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

plotter <- data.frame(score=xfold_scores[,1], sampleType=c(rep("AN's",length(unlist(folds_g1))),rep("T's",length(unlist(folds_g2)))), top=rep("top-1 single, AUC=0.9827",812))
plotter <- rbind(plotter, data.frame(score=apply(xfold_scores[,1:2],1,sum), sampleType=c(rep("AN's",length(unlist(folds_g1))),rep("T's",length(unlist(folds_g2)))), top=rep("top-2 combined, AUC=0.9903",812)))
plotter <- rbind(plotter, data.frame(score=apply(xfold_scores[,1:3],1,sum), sampleType=c(rep("AN's",length(unlist(folds_g1))),rep("T's",length(unlist(folds_g2)))), top=rep("top-3 combined, AUC=0.9938",812)))
pdf(file="top-3 models combined.pdf",width=11.7,height=4)
ggplot(plotter, aes(x=score, colour=sampleType)) + geom_density() + geom_rug(alpha=0.15) + facet_wrap(~top, scales="free_y") + theme_bw() + scale_colour_manual(values=cb_palette) + theme(legend.position="bottom")
dev.off()

plotter <- rbind(plotter, data.frame(score=apply(xfold_scores[,1:4],1,sum), sampleType=c(rep("AN's",length(unlist(folds_g1))),rep("T's",length(unlist(folds_g2)))), top=rep("top-4 combined, AUC=0.9943",812)))
plotter <- rbind(plotter, data.frame(score=apply(xfold_scores[,1:5],1,sum), sampleType=c(rep("AN's",length(unlist(folds_g1))),rep("T's",length(unlist(folds_g2)))), top=rep("top-5 combined, AUC=0.9914",812)))
plotter <- rbind(plotter, data.frame(score=apply(xfold_scores[,1:6],1,sum), sampleType=c(rep("AN's",length(unlist(folds_g1))),rep("T's",length(unlist(folds_g2)))), top=rep("top-6 combined, AUC=0.9923",812)))
ggplot(plotter, aes(x=score, colour=sampleType)) + geom_density() + geom_rug(alpha=0.15) + facet_wrap(~top, scales="free_y") + theme_bw() + scale_colour_manual(values=cb_palette) + theme(legend.position="bottom")