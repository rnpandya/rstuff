#this only runs in Linux, run on erg00lx
setwd("~/AML_drug_regression")
rm(list=ls())

library(parallelRandomForest)
library(parallel)
load("train_test_samples.Rdata")
source("code/plot_functions.r")

rftest <- function(expr, aucs, min_samples=60) {
	aucs[na.omit(match(colnames(expr), rownames(aucs))),] -> aucs_subset
	apply(aucs_subset, 2, function(col) {
		sum(!is.na(col))
	}) -> samples
	aucs_subset <- aucs_subset[, samples > min_samples]
	drug_models <- list()
	Ytest_mat <- aucs_subset[!is.element(rownames(aucs_subset),train_samples),]
	Yhat_mat <- Ytest_mat; Yhat_mat[] <- NA;
	for(drug in colnames(aucs_subset)) {
		y <- na.omit(aucs_subset[,drug])
		x <- expr[,names(y)]
		x <- subset(x, apply(x,1,max) > 1)
		x[] <- log2(x+1)
		x <- t(x)
		train_samples_subset <- train_samples[train_samples %in% names(y)]
		y.train <- y[train_samples_subset]
		x.train <- x[train_samples_subset,]
		y.test <- y[!is.element(names(y),train_samples_subset)];
		x.test <- x[!is.element(names(y),train_samples_subset),];
		randomForest(x=x.train, y=y.train, ntree=1000, keep.forest=T, nthreads=20) -> rf
		drug_models[[drug]] <- rf
		predict(rf, newdata=x.test) -> yhat
		#plot(yhat, y.test)
		Yhat_mat[names(yhat),drug] <- yhat
		print(drug)
	}
	list(drug_models = drug_models, Yhat_mat = Yhat_mat, Ytest_mat = Ytest_mat)
}
readRDS("aucs_filtered.Rds") -> aucs
readRDS("fpkm_UQ.Rds") -> expr

# variables for testing stuff works
#min_samples = 60
#drug = "17-AAG"

fpkm_UQ_res <- rftest(expr, aucs)
save(fpkm_UQ_res, file="randomForest_fpkm_UQ_res.Rdata")


####################
setwd("~/AML_drug_regression")

rm(list=ls())
source("code/plot_functions.r")
library(parallelRandomForest)
load("randomForest_fpkm_UQ_res.Rdata")

Ytest_all <- as.vector(fpkm_UQ_res$Ytest_mat)
Yhat_all <- as.vector(fpkm_UQ_res$Yhat_mat)

Yhat_all[!is.na(Yhat_all)] -> Yhat_all
Ytest_all[!is.na(Ytest_all)] -> Ytest_all

plot_all(Ytest_all, Yhat_all, "randomForest_overall.pdf", "Random Forest overall")
plot_byDrug(fpkm_UQ_res$Yhat_mat, fpkm_UQ_res$Ytest_mat, "randomForest_bydrug.pdf")
nrmse_mean(fpkm_UQ_res$Yhat_mat,fpkm_UQ_res$Ytest_mat)