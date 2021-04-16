library(CountClust)
library(readr)

set.seed(5)

counts_3D <- read_csv('/home/elisa/files/countsORs_3D.zip')

print("df loaded")
genes <- counts_3D$X1
counts_3D <- counts_3D[,-c(1)] #remove the column with gene names
rownames(counts_3D) <- genes
#spatial samples are interpreted as words of the models, genes as documents
log_counts <- log((10 ^ counts_3D) )

min_topics <- 2
max_topics <- 12

step_size <- 1
topics_range <- c(min_topics:max_topics)
#BIC_raw <- c()
BIC_log <- c()
likelihood <- c()
min_BIC_log <- 1000000000000 #arbitrary value, to initiate the Bayesian Information Criterion and get returned the model with the lowest BIC (not necessary when we save the output of each model inferred with k topics) 

for (k in topics_range){

print(paste0("going to infer model with logcounts data", k, " topics"))

name_log <- paste0("./lda_log_",k,"_forhsbm.rda")

lda_log <- FitGoM(log_counts, K = k, tol = 1, num_trials = 3, control=list(tmax=180), path_rda = name_log)  #the fitted model is directly saved in the directory whit pathname = name_log


metrics_log <- compGoM(log_counts,lda_log[["fit"]]) #computes BIC and log_likelihood metrics


BIC_log <- c(BIC_log,metrics_log[["BIC"]])
print(metrics_log[["BIC"]])
likelihood <- c(likelihood, metrics_log[["loglik"]])
	#if (lda_log[["BIC"]] < min_BIC_log){
	#	min_BIC_log <- lda_log[["BIC"]]
	#	print("This is the best log model")
	#	lda_best_log <- lda_log
	#}

#	if (lda_raw[["BIC"]] < min_BIC_raw){
#		min_BIC_raw <- lda_raw[["BIC"]]
#		print("This is the best raw model")
#		lda_best_raw <- lda_raw
#	}

}





#write.table(lda_best_log[["fit"]][["omega"]],'./Countclust_topic-dist.csv', sep = "\t")
 
#write.table(lda_best_log[["fit"]][["theta"]],'./Countclust_topics_word-dist.csv', sep = "\t")



#}

statistics <- as.data.frame(cbind(BIC_log, likelihood ))
write.csv(statistics,'./Countclust_metrics_log.csv')


