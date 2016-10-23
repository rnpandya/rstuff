nrmse <- function(Yhat, Yobs) {
	top = sum((Yobs - Yhat)^2)
	mu = mean(Yobs)
	bottom = sum((Yobs - mu)^2)
	return(top/bottom)
}

plot_all <- function(Ytest_all, Yhat_all, pdf_name, title) {
	quantile(Ytest_all,0.75) -> q3
	quantile(Ytest_all,0.25) -> q1
	color <- rep(1,length(Ytest_all))
	color[Yhat_all > q3 & Ytest_all > q3] <- 2 #TP
	color[Yhat_all > q3 & Ytest_all < q1] <- 3 #FP
	color[Yhat_all < q1 & Ytest_all < q1] <- 4 #TN
	color[Yhat_all < q1 & Ytest_all > q3] <- 5 #FN

	accuracy <- (sum(color==2) + sum(color==4))/sum(color %in% 2:5)
	#print(accuracy)

	#ydiff2 <- (Yhat_all - Ytest_all)^2
	#MSE <- sum(ydiff2)/length(Yhat_all)
	#NRMSD <- sqrt(MSE) / (max(Ytest_all) - min(Ytest_all))
	NRMSE <- nrmse(Yhat_all, Ytest_all)
	#CVRMSE <- sqrt(MSE) / (mean(Ytest_all))
	#print(NRMSD)

	pdf(pdf_name)
	#plot(main="fpkm UQ kbmtl combined predictions",  sub = paste0("NRMSD: ",signif(NRMSD,4)), Ytest_all, Yhat_all, pch=16, cex=0.8, col="blue", ylim=c(0,max(Yhat_all, Ytest_all)),xlim=c(0,max(Yhat_all, Ytest_all)))
	#plot(main=title,  sub = paste0("NRMSE: ",signif(NRMSD,4)," , accuracy: ",signif(accuracy,4)), Ytest_all, Yhat_all, pch=16, cex=0.6, col=color, ylim=c(0,max(Yhat_all, Ytest_all)),xlim=c(0,max(Yhat_all, Ytest_all)))
	plot(main=title,  sub = paste0("NRMSE: ",signif(NRMSE,4)), Ytest_all, Yhat_all, pch=16, cex=0.6, col=rgb(0,0,1,0.7), ylim=c(0,max(Yhat_all, Ytest_all)),xlim=c(0,max(Yhat_all, Ytest_all)))
	abline(a=0,b=1, col="red",lty=3)
	dev.off()
}

plot_byDrug <- function(Yhat, Ytest, pdf_name) {

	#if(class(Yhat) == "matrix") {
	pdf(pdf_name)
	par(mfrow=c(2,2))
	for(drug in colnames(Ytest)) {
		yhat_drug <- Yhat[,drug]
		ytest_drug <- Ytest[,drug]
		yhat_drug <- yhat_drug[!is.na(ytest_drug)]
		ytest_drug <- ytest_drug[!is.na(ytest_drug)]
		#ydiff2 <- (yhat_drug - ytest_drug)^2
		#MSE <- sum(ydiff2)/length(yhat_drug)
		#NRMSE <- sqrt(MSE) / (max(ytest_drug) - min(ytest_drug))
		NRMSE <- nrmse(yhat_drug, ytest_drug)
		#CVRMSE <- sqrt(MSE) / (mean(ytest_drug))
		print(paste0(drug," ", NRMSE))
		plot(ytest_drug, yhat_drug, xlim = c(0,max(ytest_drug, yhat_drug)*1.1), ylim = c(0,max(ytest_drug, yhat_drug)*1.1), pch=16, cex=0.8, col="blue", main=paste0("Drug: ",drug), sub = paste0("NRMSE: ",signif(NRMSE,4)))
		#abline(lm(yhat_drug ~ ytest_drug + 0), col="red")
		abline(a=0,b=1, col="red",lty=3)
	}
	dev.off()
	#}
}

nrmse_mean <- function(Yhat, Ytest) {
	nrmse_vector <- vector()
	for(drug in colnames(Ytest)) {
		yhat_drug <- Yhat[,drug]
		ytest_drug <- Ytest[,drug]
		yhat_drug <- yhat_drug[!is.na(ytest_drug)]
		ytest_drug <- ytest_drug[!is.na(ytest_drug)]
		NRMSE <- nrmse(yhat_drug, ytest_drug)
		if(!is.na(NRMSE) & !identical(NRMSE, Inf)) {
			nrmse_vector <- c(NRMSE, nrmse_vector)
		}
	}
	print(mean(nrmse_vector))
	print(sd(nrmse_vector))
}