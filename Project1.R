# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

file <- "TCGA_breast_cancer_LumA_vs_Basal_PAM50.tsv"
first10 <- c('NAT1','BIRC5','BAG1','BCL2','BLVRA','CCNB1','CCNE1','CDC6','CDC20','CDH3')
nfold=5
  
header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

LumA <- data[data$sample_id %in% first10,header2=='Luminal A']
Basal <- data[data$sample_id %in% first10,header2=='Basal-like']


# define function cross_valid so we can rerun the cross validataion with various parameters
#cross_validation <- function (nfold, alg="centroid") {

  LumA_groups <- split(sample(colnames(LumA)), 1+(seq_along(colnames(LumA)) %% nfold))
  Basal_groups <- split(sample(colnames(Basal)), 1+(seq_along(colnames(Basal)) %% nfold))
  
  result <- array()
  
  for (test_group in 1:nfold) {
    
    testA <- LumA[,colnames(LumA) %in% unlist(LumA_groups[test_group])]
    testB <- Basal[,colnames(Basal) %in% unlist(Basal_groups[test_group])]
    
    trainingA <- LumA[,!(colnames(LumA) %in% unlist(LumA_groups[test_group]))]
    trainingB <- Basal[,!(colnames(Basal) %in% unlist(Basal_groups[test_group]))]
    
    centroidA <- rowMeans(trainingA)
    centroidB <- rowMeans(trainingB)
    
    misclassifiedA <- sum(sapply(testA, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))>0 }))
    misclassifiedB <- sum(sapply(testB, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))<0 }))
    
    result[test_group] <- (misclassifiedA+misclassifiedB)/(ncol(testA)+ncol(testB))
 }
 
 c(mean(result), sd(result))
#}

#x<-data.frame(three=cross_validation(3), five=cross_validation(5), ten=cross_validation(10))
#rownames(x) <- c('mean','sd')
#x

