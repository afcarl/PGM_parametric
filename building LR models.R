library(pROC)
library(dplyr)
#####################################################################
# build logistic regression models using PGM-found genes at each fold

nfolds <- length(folds_g1)
length_g1 <- length(unlist(folds_g1))
length_g2 <- length(unlist(folds_g2))
xfold_LR_probabilities <- matrix(nrow=71,ncol=11, data=0)
for (j in 1:11) {
	for (fold in 1:nfolds) { # iterate through nfolds folds
		g1_xfold_training <- setdiff((1:length_g1),unlist(folds_g1[[fold]]))
		g2_xfold_training <- setdiff(length_g1+(1:length_g2),unlist(folds_g2[[fold]]))
		i <- top[j,fold]
		df <- data_BRCA_progressing[[i]] %>%
			as.data.frame() %>%
			mutate(EXPR = (read_count+1) / lib_size) %>%
			mutate(GB_overall = select(., starts_with("gb")) %>% rowMeans(., na.rm=T)) %>%
			mutate(PR_overall = select(., starts_with("pr")) %>% rowMeans(., na.rm=T)) %>%
			select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count) %>%
			mutate(class = factor(c(rep("G1",14), rep("G2",57))))
			
		model <- glm(as.formula("class ~ EXPR + GB_overall + PR_overall"), data=df[c(g1_xfold_training,g2_xfold_training),], family="binomial")
		xfold_LR_probabilities[unlist(c(folds_g1[fold],folds_g2[fold])), j] <- predict(model,df[-c(g1_xfold_training,g2_xfold_training),],type="response")
	}
}

aucs_single_lr <- vector()
aucs_running_combinations_lr <- vector()
for (j in 1:11) {
	aucs_single_lr[j] <- auc(predictor=xfold_LR_probabilities[,j],response=c(rep("G1",14),rep("G2",57)))
	if (j > 1) aucs_running_combinations_lr[j] <- auc(predictor=apply(xfold_LR_probabilities[,1:j],1,prod), response=c(rep("G1",14),rep("G2",57)))
}

###############################################################
# build models using top genes according to x-fold LR procedure
