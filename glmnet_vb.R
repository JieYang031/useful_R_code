#data process
setwd("/scratch/bcb/jyang32/UniD_pipeline/IDHmut/Data_process/")
load("./GBM_LGG_unid_M.RData")
anno<-read.csv("./IDHmut_anno.csv", header =T ,row.names = 1)

names <- colnames(GBM_LGG)
names2 <- gsub("\\.","-", names)
names3 <- substr(names2, 1,12)
colnames(GBM_LGG)<- names3
inter <- intersect(rownames(anno), names3)
length(inter)
#644
anno2 <- cbind(anno,anno)
anno3 <- anno2[inter,]
anno4 <- anno3[sort(rownames(anno3)),]

GBM_LGG2 <- GBM_LGG[,sort(colnames(GBM_LGG))]
table(colnames(GBM_LGG2) == rownames(anno4))
#TRUE 644

#remove missing value, left 637 sampels
save(X, Y, file = "./GBM_LGG_IDHmut_clean.Rdata")

####variable selection####
var_sel <- function(alpha){
  
  list.of.packages <- c("caret", "mlr", "glmnet", "xgboost")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, require, character.only = TRUE)
  
  load("/scratch/bcb/jyang32/UniD_pipeline/IDHmut/Data_process/GBM_LGG_IDHmut_clean.Rdata")
  
  
  set.seed(2016)
  trainID = createDataPartition(Y, p = 0.6, list = F)
  testID = createDataPartition(Y[-trainID], p = 0.5, list = F)
  
  trainX = X[trainID, ]
  testX  = X[-trainID, ][testID, ]
  validX = X[-trainID, ][-testID, ]
  trainY = Y[trainID]
  testY  = Y[-trainID][testID]
  validY = Y[-trainID][-testID]
  
  
  alpha = alpha;
  K = 5
  foldInds <- createFolds(trainY, k = K, list=TRUE, returnTrain=FALSE)
  
  ## use K fold cross validation
  selectModel = lapply(1:K, function(k) {
    foldID = foldInds[[k]]
    Xtrain = trainX[-foldID, ]
    Ytrain = trainY[-foldID]
    Xtest  = trainX[foldID,  ]
    Ytest  = trainY[foldID]
    model = glmnet(Xtrain, Ytrain, alpha = alpha, nlambda = 200, family = "binomial")
    Ypred = predict(model, newx = as.matrix(Xtest), type = "class", s = model$lambda)
    Acc = apply(Ypred, 2, function(x) mean(x == Ytest))
    
    ## compute probes select rate 
    maxAcc = max(Acc)
    maxLambdaID = which(Acc == maxAcc)
    maxLambdaCoef = as.matrix(coef(model, s = model$lambda[maxLambdaID]))
    selectCoefRate = apply(maxLambdaCoef, 1, function(x) sum(x != 0)) / length(maxLambdaID)
    
    list(k = k, maxLambda = model$lambda[maxLambdaID], 
         selectCoefRate = selectCoefRate, 
         maxLambdaCoef = maxLambdaCoef,
         maxAcc = max(Acc),
         medianAcc = median(Acc))
  })
  
  
  probeRate = apply(Reduce(cbind, lapply(selectModel, function(x) x$selectCoefRate)), 1 , sum) / K
  medianAcc = median(sapply(selectModel, function(x) x$maxAcc))
  maxAcc    = max(sapply(selectModel, function(x) x$maxAcc))
  minAcc    = min(sapply(selectModel, function(x) x$maxAcc))
  medianLambda = median(unlist(lapply(selectModel, function(x) x$maxLambda)))
  maxLambda = max(unlist(lapply(selectModel, function(x) x$maxLambda)))
  minLambda = min(unlist(lapply(selectModel, function(x) x$maxLambda)))
  save(list = c("selectModel", "probeRate", "medianAcc", "maxAcc", "minAcc",
                "medianLambda", "maxLambda", "minLambda"), 
       file = sprintf("elastic_alpha_%s.Rdata", alpha))
  
}

####feature selection####
setwd("/scratch/bcb/jyang32/UniD_pipeline/IDHmut/Variable_selection/result")
load("/scratch/bcb/jyang32/UniD_pipeline/IDHmut/Data_process/GBM_LGG_IDH_clean.Rdata")

Rfiles = grep("\\belastic_alpha_\\w.*.Rdata", list.files(), value = T)

probeRateList = lapply(c(1:length(Rfiles)), function(x) {
  print(x)
  load(Rfiles[x])
  probeRate
})

probeRateMat = Reduce(cbind, probeRateList)[-1,]
probeRateSum = apply(probeRateMat, 1, mean)
sum(probeRateSum !=0)
#In total probes: 1513

top20 = head(sort(probeRateSum, decreasing = T), 20)
top50 = head(sort(probeRateSum, decreasing = T), 50)
top100 = head(sort(probeRateSum, decreasing = T), 100)
top200 = head(sort(probeRateSum, decreasing = T), 200)
top500 = head(sort(probeRateSum, decreasing = T), 500)
top1000 = head(sort(probeRateSum, decreasing = T), 1000)
top1500 = head(sort(probeRateSum, decreasing = T), 1500)

top20Probes = names(top20)
top50Probes = names(top50)
top100Probes = names(top100)
top200Probes = names(top200)
top500Probes = names(top500)
top1000Probes = names(top1000)
top1500Probes = names(top1500)

allSelectedProbes = names(probeRateSum[probeRateSum > 0])

save(list = c("top20Probes" ,"top50Probes", "top100Probes", "top200Probes", "top500Probes",  "top1000Probes", "top1500Probes", "allSelectedProbes"), 
     file  = "IDH_GBM_LGG_features.Rdata")

X20 = X[,top20Probes]
X50 = X[,top50Probes]
X100 = X[,top100Probes]
X200 = X[,top200Probes]
X500 = X[,top500Probes]
X1000 = X[,top1000Probes]
X1500 = X[,top1500Probes]

save(list ="X20" ,"X50", "X100", "X200", "X500","X1000","X1500", "Y",
     file = "IDH_GBM_LGG_Select.Rdata")

#####model_selection
load("/scratch/bcb/jyang32/UniD_pipeline/IDHmut/feature/IDH_GBM_LGG_Select.Rdata")
plot_parameter <- function(input){
  list.of.packages <- c("caret", "mlr", "glmnet", "xgboost", "reshape2", "ggplot2")
  lapply(list.of.packages, require, character.only = TRUE)
  
  
  
  set.seed(2016)
  trainID = createDataPartition(Y, p = 0.6, list = F)
  testID = createDataPartition(Y[-trainID], p = 0.5, list = F)
  
  X = input
  trainX = X[trainID, ]
  testX  = X[-trainID, ][testID, ]
  validX = X[-trainID, ][-testID, ]
  trainY = Y[trainID]
  testY  = Y[-trainID][testID]
  validY = Y[-trainID][-testID]
  
  
  K = 5
  foldInds <- createFolds(trainY, k = K, list=TRUE, returnTrain=FALSE)
  
  alphaSeq = seq(0, 1, by = 0.1)
  lambdaSeq = seq(5, 0, by = -0.05)
  ## use K fold cross validation
  selectParam = lapply(alphaSeq, function(alpha) {
    selectLambda = lapply(1:K, function(k) {
      foldID = foldInds[[k]]
      Xtrain = trainX[-foldID, ]
      Ytrain = trainY[-foldID]
      Xtest  = trainX[foldID,  ]
      Ytest  = trainY[foldID]
      model = glmnet(Xtrain, Ytrain, alpha = alpha, lambda = lambdaSeq, family = "binomial")
      Ypred = predict(model, newx = as.matrix(Xtest), type = "class", s = model$lambda)
      Acc = apply(Ypred, 2, function(x) mean(x == Ytest))
      Acc
    })
    LambdaMat = as.data.frame(Reduce(cbind, selectLambda))
  })
  
  
  selectParamMat = Reduce(rbind, selectParam)
  selectParamDF = data.frame(mean = rowMeans(selectParamMat), 
                             sd = apply(selectParamMat, 1, sd),
                             lambda = rep(lambdaSeq, length(alphaSeq)),
                             alpha = rep(alphaSeq, each = length(lambdaSeq)))
  selectParamDF$alpha = as.factor(selectParamDF$alpha)
  
  limits <- aes(ymax = mean + sd, ymin = mean - sd)
  ggplot(selectParamDF, aes(x = lambda, y = mean, color = alpha)) + geom_line() + 
    geom_point() + geom_pointrange(limits) + 
    theme_bw() + xlab("Lambda") + ylab("CV Acc") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position="right", 
          plot.title = element_text(size = 20,  hjust = 0),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    ggtitle(sprintf("%d folds Cross Validation Accuracy versus Lambda & Alpha (nprobes = %d)", K, ncol(X)))
  
  ggsave(file = sprintf("Elastic_CV_nprobes_%d.pdf", ncol(X)), width = 11, height = 8)
  
}

model_sel <- function(input, alpha , lambda){
  set.seed(2016)
  trainID = createDataPartition(Y, p = 0.6, list = F)
  testID = createDataPartition(Y[-trainID], p = 0.5, list = F)
  
  alpha = alpha
  lambda = lambda
  
  X = input
  trainX = X[trainID, ]
  testX  = X[-trainID, ][testID, ]
  validX = X[-trainID, ][-testID, ]
  trainY = Y[trainID]
  testY  = Y[-trainID][testID]
  validY = Y[-trainID][-testID]
  
  model = glmnet(trainX, trainY, alpha = alpha, family = "binomial")
  
  #draw ROC curve
  library(ROCR)
  trainprob <- predict(model, as.matrix(trainX),type = "response", s = lambda)
  pred <- prediction(trainprob, trainY)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  pdf(paste0("ROC_curve_IDH(GBM_LGG)_",ncol(X),"probe.pdf"))
  plot(perf, col = rainbow(7))
  abline(0,1)
  dev.off()
  
  #calcualte AUC
  library(caTools)
  auc_value<- colAUC(trainprob, trainY, plotROC = T)
  message(paste0("AUC ", auc_value))
  
  YpredTrain = predict(model, newx = as.matrix(trainX), type = "class", s = lambda)
  trainAcc = mean(YpredTrain == trainY)
  message(paste0("trainAcc ", trainAcc))
  
  
  YpredTest = predict(model, newx = as.matrix(testX), type = "class", s = lambda)
  testAcc = mean(YpredTest == testY)
  message(paste0("testAcc ", testAcc))
  
  
  YpredValid = predict(model, newx = as.matrix(validX), type = "class", s = lambda)
  validAcc = mean(YpredValid == validY)
  message(paste0("validAcc ", validAcc))
  
  
  # expore about the difference
  probeNames = colnames(X)
  
  save(list = c("model", "probeNames", "trainAcc","testAcc", "validAcc", "auc_value"), 
       file = sprintf("IDH_probes_%d_alpha_%g_lambda_%g.Rdata", ncol(X), alpha, lambda))
}


# > model_sel(alpha = 0, lambda = 0.1, input=X20)
# AUC 1
# trainAcc 0.997389033942559
# testAcc 1
# validAcc 1
# > model_sel(alpha = 0, lambda = 0.1, input=X50)
# AUC 1
# trainAcc 1
# testAcc 1
# validAcc 1
# > model_sel(alpha = 0, lambda = 0.1, input=X100)
# AUC 1
# trainAcc 1
# testAcc 1
# validAcc 1
# > model_sel(alpha = 0, lambda = 0.1, input=X200)
# AUC 1
# trainAcc 1
# testAcc 1
# validAcc 1
# > model_sel(alpha = 0, lambda = 0.1, input=X500)
# AUC 1
# trainAcc 0.997389033942559
# testAcc 1
# validAcc 1
# > model_sel(alpha = 0, lambda = 0.1, input=X1000)
# AUC 1
# trainAcc 0.997389033942559
# testAcc 1
# validAcc 0.992125984251969
# > model_sel(alpha = 0, lambda = 0.1, input=X1500)
# AUC 1
# trainAcc 1
# testAcc 1
# validAcc 0.992125984251969

set.seed(2016)
trainID = createDataPartition(Y, p = 0.6, list = F)
testID = createDataPartition(Y[-trainID], p = 0.5, list = F)

input = X100
alpha = 0
lambda = 0.1

X = input
trainX = X[trainID, ]
testX  = X[-trainID, ][testID, ]
validX = X[-trainID, ][-testID, ]
trainY = Y[trainID]
testY  = Y[-trainID][testID]
validY = Y[-trainID][-testID]

model = glmnet(trainX, trainY, alpha = alpha, family = "binomial")
probeNames <- colnames(X)

save(model, probeNames, file= "UniD_IDHmut_codel_probe100_alpha_0_lambda_0.1.Rdata")
