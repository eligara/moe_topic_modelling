library(CountClust)
library(readr)

#set.seed(5)

counts_3D<- read_csv('/home/elisa/files/countsORs_3D.zip')

print("df loaded")
genes <- counts_3D$X1
counts_3D <- counts_3D[,-c(1)]
rownames(counts_3D) <- genes
#counts_3D <- t(counts_3D)
#raw_counts <- (10 ^ counts_3D) - 1
log_counts <- log((10 ^ counts_3D))

#min_topics <- 7
#max_topics <- 12

#topics_raw <- 5

#print(paste0("topics to infer", 5))

step_size <- 1
topics_range <- c(24,95)#(min_topics:max_topics)
#BIC_raw <- c()
BIC_log <- c()
likelihood <- c()
#min_BIC_raw <- 1000000000000
min_BIC_log <- 1000000000000

for (k in topics_range){

print(paste0("going to infer model with logcounts data", k, " topics"))

name_log <- paste0("./lda_log_",k,"_forhsbm_1trial.rda")
#name_raw <- paste0("./lda_raw_",k,"_forhsbm.rda")

lda_log <- FitGoM(log_counts,K = k, tol = 1, num_trials = 1,control=list(tmax=180), path_rda = name_log)
#lda_raw <- FitGoM(raw_counts,K = k, tol = 1, num_trials = 3,control=list(tmax=180), path_rda = name_raw)


#metrics_raw <- compGoM(raw_counts,lda_raw[["fit"]])
metrics_log <- compGoM(log_counts,lda_log[["fit"]])

#print(" absolute counts fitted")
#BIC_raw <- c(BIC_raw,metrics_raw[["BIC"]])
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

#lda_8 <- FitGoM(raw_counts,K = 8, tol = 5, num_trials = 3, path_rda = "./lda_8_raw.rda",control=list(tmax=180))
#lda_5_raw <- FitGoM(raw_counts,K = 5, tol = 5, num_trials = 3, path_rda = "./lda_5_raw.rda",control=list(tmax=180))


#write.table(lda_best_raw[["fit"]][["omega"]],'./Countclust_topic-dist.csv', sep = "\t")
 
#write.table(lda_best_log[["fit"]][["theta"]],'./Countclust_topics_word-dist.csv', sep = "\t")



#}

statistics <- as.data.frame(cbind(BIC_log, likelihood ))
write.csv(statistics,'./Countclust_metrics_log.csv')
#write.csv(BIC_raw,'./Countclust_metrics_raw.csv')

