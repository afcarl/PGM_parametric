folds_an <- NULL
folds_an[[1]] <- sample(1:82,size=11)
folds_an[[2]] <- sample(setdiff(1:82,unlist(folds_an)),size=11)
folds_an[[3]] <- sample(setdiff(1:82,unlist(folds_an)),size=10)
folds_an[[4]] <- sample(setdiff(1:82,unlist(folds_an)),size=10)
folds_an[[5]] <- sample(setdiff(1:82,unlist(folds_an)),size=10)
folds_an[[6]] <- sample(setdiff(1:82,unlist(folds_an)),size=10)
folds_an[[7]] <- sample(setdiff(1:82,unlist(folds_an)),size=10)
folds_an[[8]] <- sample(setdiff(1:82,unlist(folds_an)),size=10)
# sort(unlist(folds_an))

folds_t <- NULL
folds_t[[1]] <- sample(82+(1:730),size=92)
folds_t[[2]] <- sample(setdiff((82+(1:730)),unlist(folds_t)),size=92)
folds_t[[3]] <- sample(setdiff((82+(1:730)),unlist(folds_t)),size=91)
folds_t[[4]] <- sample(setdiff((82+(1:730)),unlist(folds_t)),size=91)
folds_t[[5]] <- sample(setdiff((82+(1:730)),unlist(folds_t)),size=91)
folds_t[[6]] <- sample(setdiff((82+(1:730)),unlist(folds_t)),size=91)
folds_t[[7]] <- sample(setdiff((82+(1:730)),unlist(folds_t)),size=91)
folds_t[[8]] <- sample(setdiff((82+(1:730)),unlist(folds_t)),size=91)
# sort(unlist(folds_t))



# progressing
folds_g1 <- NULL
for (i in 1:14) {folds_g1[[i]] <- i}

folds_g2 <- NULL
folds_g2[[1]] <- sample(14+(1:57),size=5)
folds_g2[[2]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[3]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[4]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[5]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[6]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[7]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[8]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[9]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[10]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[11]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[12]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[13]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
folds_g2[[14]] <- sample(setdiff(14+(1:57),unlist(folds_g2)),size=4)
