library(readxl)
library(dplyr)
library(readr)
library(caret)
library(pROC)
library(e1071)
options(scipen = 999)
#save(list=ls(),file="genentech_v2.rda")
options(readr.default_locale=readr::locale(tz="US/Eastern"))
#-----------------------
# TRIAGE NOTES
#-----------------------
#rm(list=ls(pattern="^triage"))
triage.cases <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CASES_INTERESTED_NOTES_FINAL_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))
triage.controls.1 <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CONTROL_1_INTERESTED_NOTES_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))
triage.controls.3 <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CONTROL_3_INTERESTED_NOTES_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))
#intersect(colnames(tmp.cases),colnames(tmp.controls))
triage.controls.1[setdiff(names(triage.controls.3), names(triage.controls.1))] <- 0
triage.controls.3[setdiff(names(triage.controls.1), names(triage.controls.3))] <- 0
triage.controls <- rbind(triage.controls.1, triage.controls.3)
triage.cases[setdiff(names(triage.controls), names(triage.cases))] <- 0
triage.controls[setdiff(names(triage.cases), names(triage.controls))] <- 0
triage <- rbind(triage.cases, triage.controls ) 
triage.svd <- svd(triage)
sum(diag(triage.svd$d)[-1,])
sum(diag(triage.svd$d)[-1,][1:250,])
triage.approx <- list()
triage.approx[[1]] <- triage.svd$u[,1:31] %*% diag(triage.svd$d)[1:31,1:31] %*% t((triage.svd$v)[1:31,1:31])
triage.approx[[2]] <- triage.svd$u[,1:62] %*% diag(triage.svd$d)[1:62,1:62] %*% t((triage.svd$v)[1:62,1:62])
triage.approx[[3]] <- triage.svd$u[,1:124] %*% diag(triage.svd$d)[1:124,1:124] %*% t((triage.svd$v)[1:124,1:124])
triage.approx[[4]] <- triage.svd$u[,1:332] %*% diag(triage.svd$d)[1:332,1:332] %*% t((triage.svd$v)[1:332,1:332])
triage.approx[[5]] <- triage.svd$u %*% diag(triage.svd$d) %*% t((triage.svd$v))
names(triage.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(triage.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

triage.trainSet <- list()
triage.testSet <- list()
for (i in 1:length(triage.approx)) {
  triage.approx[[i]] <- cbind(triage.approx[[i]], label)%>%
    select(label,everything())
  rownames(triage.approx[[i]]) <- NULL
  index <- caret::createDataPartition(triage.approx[[i]]$label, p = .8, list = F)
  triage.trainSet[[i]] <- triage.approx[[i]][index,]
  triage.testSet[[i]] <- triage.approx[[i]][-index,]
  rm(index)
}
names(triage.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
names(triage.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
paramGrid<-trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        summaryFunction = twoClassSummary,                      # Evaluate performance
                        classProbs = T,                                         # Estimate class probabilities
                        allowParallel = T,
                        search = "random")
ml.train.triage <- list()
startTime <- list()
finishTime <- list()
classifiers <- c("glm", "rf", "svmRadial", "xgbDART")
# Model development
for(i in 1:length(triage.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(triage.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(triage.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.triage[[paste0(names(triage.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                       data = triage.trainSet[[i]],
                                                                                       method = classifiers[[c]],
                                                                                       preProcess = c("center","scale"),
                                                                                       metric = "ROC",
                                                                                       trControl = paramGrid,
                                                                                       tuneLength = 3)
    finishTime[[paste0(names(triage.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(triage.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
fetchResults<-function(x,y){
  z<-as.data.frame(x)
  colnames(z)<-names(y)
  return(z)
}
incrementStart<-function(x){4*x-3}
incrementEnd<-function(x){4*x}
triage.results <- as.data.frame(list())
for (dat in 1:length(triage.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        triage.testSet[[dat]]$label,
        predict(object = ml.train.triage[[c]], triage.testSet[[dat]], type = "prob"),
        predict(object = ml.train.triage[[c]], triage.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    triage.results <- rbind(triage.results,
                            tmp%>%
                              mutate(
                                "Classifier" = names(ml.train.triage[c]),
                                "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                              )%>%
                              select(
                                c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                  "Kappa","Sensitivity","Specificity",
                                  "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                              )
    )
    rm(tmp, cm, predictions)
  }
}


#-----------------------
# PROVIDER NOTES
#-----------------------
#rm(list=ls(pattern="^provider"))
setwd("~/Genentech/NLP")
provider.cases <- read_csv("DEID_STROKE_CASES_INTERESTED_PROVIDER_NOTES_FINAL_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))
provider.controls.1 <- read_csv("DEID_STROKE_CONTROL_1_INTERESTED_PROVIDER_NOTES_FINAL_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))
provider.controls.3 <- read_csv("DEID_STROKE_CONTROL_3_INTERESTED_PROVIDER_NOTES_FINAL_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))
provider.controls.1[setdiff(names(provider.controls.3), names(provider.controls.1))] <- 0
provider.controls.3[setdiff(names(provider.controls.1), names(provider.controls.3))] <- 0
provider.controls <- rbind(provider.controls.1, provider.controls.3)
provider.cases[setdiff(names(provider.controls), names(provider.cases))] <- 0
provider.controls[setdiff(names(provider.cases), names(provider.controls))] <- 0
provider <- rbind(provider.cases, provider.controls ) 
provider.svd <- svd(provider)
sum(diag(provider.svd$d)[-1,])
sum(diag(provider.svd$d)[-1,][1:250,])
provider.approx <- list()
provider.approx[[1]] <- provider.svd$u[,1:31] %*% diag(provider.svd$d)[1:31,1:31] %*% t((provider.svd$v)[1:31,1:31])
provider.approx[[2]] <- provider.svd$u[,1:62] %*% diag(provider.svd$d)[1:62,1:62] %*% t((provider.svd$v)[1:62,1:62])
provider.approx[[3]] <- provider.svd$u[,1:124] %*% diag(provider.svd$d)[1:124,1:124] %*% t((provider.svd$v)[1:124,1:124])
provider.approx[[4]] <- provider.svd$u[,1:291] %*% diag(provider.svd$d)[1:291,1:291] %*% t((provider.svd$v)[1:291,1:291])
provider.approx[[5]] <- provider.svd$u %*% diag(provider.svd$d) %*% t((provider.svd$v))
names(provider.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

provider.trainSet <- list()
provider.testSet <- list()
for (i in 1:length(provider.approx)) {
  provider.approx[[i]] <- cbind(provider.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.approx[[i]]) <- NULL
  index <- caret::createDataPartition(provider.approx[[i]]$label, p = .8, list = F)
  provider.trainSet[[i]] <- provider.approx[[i]][index,]
  provider.testSet[[i]] <- provider.approx[[i]][-index,]
  rm(index)
}
names(provider.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
names(provider.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
ml.train.provider <- list()
# Model development
for(i in 1:length(provider.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(provider.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(provider.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.provider[[paste0(names(provider.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                           data = provider.trainSet[[i]],
                                                                                           method = classifiers[[c]],
                                                                                           preProcess = c("center","scale"),
                                                                                           metric = "ROC",
                                                                                           trControl = paramGrid,
                                                                                           tuneLength = 3)
    finishTime[[paste0(names(provider.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(provider.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
provider.results <- as.data.frame(list())
rm(c,i,dat)
for (dat in 1:length(provider.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.testSet[[dat]]$label,
        predict(object = ml.train.provider[[c]], provider.testSet[[dat]], type = "prob"),
        predict(object = ml.train.provider[[c]], provider.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.provider[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.provider[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.provider[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,        # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                                # Calculate AUPR
    )))
    provider.results <- rbind(provider.results,
                              tmp%>%
                                mutate(
                                  "Classifier" = names(ml.train.provider[c]),
                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                )%>%
                                select(
                                  c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                    "Kappa","Sensitivity","Specificity",
                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                )
    )
    rm(tmp, cm, predictions)
  }
}

#-----------------------
# Triage: Binary
#-----------------------
triage.binary <- triage %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
triage.binary.svd <- svd(triage.binary)
triage.binary.approx <- list()
triage.binary.approx[[1]] <- triage.binary.svd$u[,1:31] %*% diag(triage.binary.svd$d)[1:31,1:31] %*% t((triage.binary.svd$v)[1:31,1:31])
triage.binary.approx[[2]] <- triage.binary.svd$u[,1:62] %*% diag(triage.binary.svd$d)[1:62,1:62] %*% t((triage.binary.svd$v)[1:62,1:62])
triage.binary.approx[[3]] <- triage.binary.svd$u[,1:124] %*% diag(triage.binary.svd$d)[1:124,1:124] %*% t((triage.binary.svd$v)[1:124,1:124])
triage.binary.approx[[4]] <- triage.binary.svd$u[,1:332] %*% diag(triage.binary.svd$d)[1:332,1:332] %*% t((triage.binary.svd$v)[1:332,1:332])
triage.binary.approx[[5]] <- triage.binary.svd$u %*% diag(triage.binary.svd$d) %*% t((triage.binary.svd$v))
names(triage.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(triage.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

triage.binary.trainSet <- list()
triage.binary.testSet <- list()
for (i in 1:length(triage.binary.approx)) {
  triage.binary.approx[[i]] <- cbind(triage.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(triage.binary.approx[[i]]) <- NULL
  index <- caret::createDataPartition(triage.binary.approx[[i]]$label, p = .8, list = F)
  triage.binary.trainSet[[i]] <- triage.binary.approx[[i]][index,]
  triage.binary.testSet[[i]] <- triage.binary.approx[[i]][-index,]
  rm(index)
}
names(triage.binary.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
names(triage.binary.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
ml.train.triage.binary <- list()
# Model development
for(i in 1:length(triage.binary.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(triage.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(triage.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.triage.binary[[paste0(names(triage.binary.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                     data = triage.binary.trainSet[[i]],
                                                                                                     method = classifiers[[c]],
                                                                                                     preProcess = c("center","scale"),
                                                                                                     metric = "ROC",
                                                                                                     trControl = paramGrid,
                                                                                                     tuneLength = 3)
    finishTime[[paste0(names(triage.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(triage.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
incrementStart<-function(x){4*x-3}
incrementEnd<-function(x){4*x}
triage.binary.results <- as.data.frame(list())
for (dat in 1:length(triage.binary.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        triage.binary.testSet[[dat]]$label,
        predict(object = ml.train.triage.binary[[c]], triage.binary.testSet[[dat]], type = "prob"),
        predict(object = ml.train.triage.binary[[c]], triage.binary.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    triage.binary.results <- rbind(triage.binary.results,
                                   tmp%>%
                                     mutate(
                                       "Classifier" = names(ml.train.triage.binary[c]),
                                       "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                     )%>%
                                     select(
                                       c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                         "Kappa","Sensitivity","Specificity",
                                         "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                     )
    )
    rm(tmp, cm, predictions)
  }
}
#-----------------------
# Provider: Binary
#-----------------------
provider.binary <- provider %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
provider.binary.svd <- svd(provider.binary)
provider.binary.approx <- list()
provider.binary.approx[[1]] <- provider.binary.svd$u[,1:31] %*% diag(provider.binary.svd$d)[1:31,1:31] %*% t((provider.binary.svd$v)[1:31,1:31])
provider.binary.approx[[2]] <- provider.binary.svd$u[,1:62] %*% diag(provider.binary.svd$d)[1:62,1:62] %*% t((provider.binary.svd$v)[1:62,1:62])
provider.binary.approx[[3]] <- provider.binary.svd$u[,1:124] %*% diag(provider.binary.svd$d)[1:124,1:124] %*% t((provider.binary.svd$v)[1:124,1:124])
provider.binary.approx[[4]] <- provider.binary.svd$u[,1:291] %*% diag(provider.binary.svd$d)[1:291,1:291] %*% t((provider.binary.svd$v)[1:291,1:291])
provider.binary.approx[[5]] <- provider.binary.svd$u %*% diag(provider.binary.svd$d) %*% t((provider.binary.svd$v))
names(provider.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

provider.binary.trainSet <- list()
provider.binary.testSet <- list()
for (i in 1:length(provider.binary.approx)) {
  provider.binary.approx[[i]] <- cbind(provider.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.binary.approx[[i]]) <- NULL
  index <- caret::createDataPartition(provider.binary.approx[[i]]$label, p = .8, list = F)
  provider.binary.trainSet[[i]] <- provider.binary.approx[[i]][index,]
  provider.binary.testSet[[i]] <- provider.binary.approx[[i]][-index,]
  rm(index)
}
names(provider.binary.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
names(provider.binary.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
ml.train.provider.binary <- list()
# Model development
for(i in 1:length(provider.binary.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(provider.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(provider.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.provider.binary[[paste0(names(provider.binary.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                         data = provider.binary.trainSet[[i]],
                                                                                                         method = classifiers[[c]],
                                                                                                         preProcess = c("center","scale"),
                                                                                                         metric = "ROC",
                                                                                                         trControl = paramGrid,
                                                                                                         tuneLength = 3)
    finishTime[[paste0(names(provider.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(provider.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
provider.binary.results <- as.data.frame(list())
rm(c,i,dat)
for (dat in 1:length(provider.binary.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.binary.testSet[[dat]]$label,
        predict(object = ml.train.provider.binary[[c]], provider.binary.testSet[[dat]], type = "prob"),
        predict(object = ml.train.provider.binary[[c]], provider.binary.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.provider.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.provider.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.provider.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,        # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                                # Calculate AUPR
    )))
    provider.binary.results <- rbind(provider.binary.results,
                                     tmp%>%
                                       mutate(
                                         "Classifier" = names(ml.train.provider.binary[c]),
                                         "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                       )%>%
                                       select(
                                         c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                           "Kappa","Sensitivity","Specificity",
                                           "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                       )
    )
    rm(tmp, cm, predictions)
  }
}

################################################################################
##########################      NEGATED    #####################################
################################################################################
#-----------------------
# TRIAGE NOTES Negated
#-----------------------
#rm(list=ls(pattern="^negated"))
negated.triage.cases <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CASES_INTERESTED_NOTES_FINAL_FORMATTED_NEGATED.csv")%>%
  select(-c('pt_id','enc_dt'))
negated.triage.controls.1 <- read_csv("DEID_STROKE_CONTROL_1_INTERESTED_NOTES_FORMATTED_NEGATED.csv")%>%
  select(-c('pt_id','enc_dt'))
negated.triage.controls.3 <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CONTROL_3_INTERESTED_NOTES_FORMATTED_NEGATED.csv")%>%
  select(-c('pt_id','enc_dt'))
negated.triage.controls.1[setdiff(names(negated.triage.controls.3), names(negated.triage.controls.1))] <- 0
negated.triage.controls.3[setdiff(names(negated.triage.controls.1), names(negated.triage.controls.3))] <- 0
negated.triage.controls <- rbind(negated.triage.controls.1, negated.triage.controls.3)
negated.triage.cases[setdiff(names(negated.triage.controls), names(negated.triage.cases))] <- 0
negated.triage.controls[setdiff(names(negated.triage.cases), names(negated.triage.controls))] <- 0
negated.triage <- rbind(negated.triage.cases, negated.triage.controls ) 
negated.triage.svd <- svd(negated.triage)
negated.triage.approx <- list()
negated.triage.approx[[1]] <- negated.triage.svd$u[,1:31] %*% diag(negated.triage.svd$d)[1:31,1:31] %*% t((negated.triage.svd$v)[1:31,1:31])
negated.triage.approx[[2]] <- negated.triage.svd$u[,1:62] %*% diag(negated.triage.svd$d)[1:62,1:62] %*% t((negated.triage.svd$v)[1:62,1:62])
negated.triage.approx[[3]] <- negated.triage.svd$u[,1:124] %*% diag(negated.triage.svd$d)[1:124,1:124] %*% t((negated.triage.svd$v)[1:124,1:124])
negated.triage.approx[[4]] <- negated.triage.svd$u[,1:332] %*% diag(negated.triage.svd$d)[1:332,1:332] %*% t((negated.triage.svd$v)[1:332,1:332])
negated.triage.approx[[5]] <- negated.triage.svd$u %*% diag(negated.triage.svd$d) %*% t((negated.triage.svd$v))
names(negated.triage.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(negated.triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(negated.triage.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

negated.triage.trainSet <- list()
negated.triage.testSet <- list()
for (i in 1:length(negated.triage.approx)) {
  negated.triage.approx[[i]] <- cbind(negated.triage.approx[[i]], label)%>%
    select(label,everything())
  rownames(negated.triage.approx[[i]]) <- NULL
  index <- caret::createDataPartition(negated.triage.approx[[i]]$label, p = .8, list = F)
  negated.triage.trainSet[[i]] <- negated.triage.approx[[i]][index,]
  negated.triage.testSet[[i]] <- negated.triage.approx[[i]][-index,]
  rm(index)
}
names(negated.triage.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
names(negated.triage.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
ml.train.negated.triage <- list()
# Model development
for(i in 1:length(negated.triage.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(negated.triage.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(negated.triage.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.negated.triage[[paste0(names(negated.triage.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                       data = negated.triage.trainSet[[i]],
                                                                                                       method = classifiers[[c]],
                                                                                                       preProcess = c("center","scale"),
                                                                                                       metric = "ROC",
                                                                                                       trControl = paramGrid,
                                                                                                       tuneLength = 3)
    finishTime[[paste0(names(negated.triage.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(negated.triage.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
incrementStart<-function(x){4*x-3}
incrementEnd<-function(x){4*x}
negated.triage.results <- as.data.frame(list())
for (dat in 1:length(negated.triage.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        negated.triage.testSet[[dat]]$label,
        predict(object = ml.train.negated.triage[[c]], negated.triage.testSet[[dat]], type = "prob"),
        predict(object = ml.train.negated.triage[[c]], negated.triage.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.triage[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.triage[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.triage[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    negated.triage.results <- rbind(negated.triage.results,
                                    tmp%>%
                                      mutate(
                                        "Classifier" = names(ml.train.negated.triage[c]),
                                        "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                      )%>%
                                      select(
                                        c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                          "Kappa","Sensitivity","Specificity",
                                          "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                      )
    )
    rm(tmp, cm, predictions)
  }
}
rm(c,i,dat)

#-----------------------
# Triage: Binary Negated
#-----------------------
negated.triage.binary <- negated.triage %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
negated.triage.binary.svd <- svd(negated.triage.binary)
negated.triage.binary.approx <- list()
negated.triage.binary.approx[[1]] <- negated.triage.binary.svd$u[,1:31] %*% diag(negated.triage.binary.svd$d)[1:31,1:31] %*% t((negated.triage.binary.svd$v)[1:31,1:31])
negated.triage.binary.approx[[2]] <- negated.triage.binary.svd$u[,1:62] %*% diag(negated.triage.binary.svd$d)[1:62,1:62] %*% t((negated.triage.binary.svd$v)[1:62,1:62])
negated.triage.binary.approx[[3]] <- negated.triage.binary.svd$u[,1:124] %*% diag(negated.triage.binary.svd$d)[1:124,1:124] %*% t((negated.triage.binary.svd$v)[1:124,1:124])
negated.triage.binary.approx[[4]] <- negated.triage.binary.svd$u[,1:332] %*% diag(negated.triage.binary.svd$d)[1:332,1:332] %*% t((negated.triage.binary.svd$v)[1:332,1:332])
negated.triage.binary.approx[[5]] <- negated.triage.binary.svd$u %*% diag(negated.triage.binary.svd$d) %*% t((negated.triage.binary.svd$v))
names(negated.triage.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(negated.triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(negated.triage.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

negated.triage.binary.trainSet <- list()
negated.triage.binary.testSet <- list()
for (i in 1:length(negated.triage.binary.approx)) {
  negated.triage.binary.approx[[i]] <- cbind(negated.triage.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(negated.triage.binary.approx[[i]]) <- NULL
  index <- caret::createDataPartition(negated.triage.binary.approx[[i]]$label, p = .8, list = F)
  negated.triage.binary.trainSet[[i]] <- negated.triage.binary.approx[[i]][index,]
  negated.triage.binary.testSet[[i]] <- negated.triage.binary.approx[[i]][-index,]
  rm(index)
}
names(negated.triage.binary.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
names(negated.triage.binary.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
ml.train.negated.triage.binary <- list()
# Model development
for(i in 1:length(negated.triage.binary.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(negated.triage.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(negated.triage.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.negated.triage.binary[[paste0(names(negated.triage.binary.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                     data = negated.triage.binary.trainSet[[i]],
                                                                                                                     method = classifiers[[c]],
                                                                                                                     preProcess = c("center","scale"),
                                                                                                                     metric = "ROC",
                                                                                                                     trControl = paramGrid,
                                                                                                                     tuneLength = 3)
    finishTime[[paste0(names(negated.triage.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(negated.triage.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
incrementStart<-function(x){4*x-3}
incrementEnd<-function(x){4*x}
negated.triage.binary.results <- as.data.frame(list())
for (dat in 1:length(negated.triage.binary.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        negated.triage.binary.testSet[[dat]]$label,
        predict(object = ml.train.negated.triage.binary[[c]], negated.triage.binary.testSet[[dat]], type = "prob"),
        predict(object = ml.train.negated.triage.binary[[c]], negated.triage.binary.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.triage.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.triage.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.triage.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    negated.triage.binary.results <- rbind(negated.triage.binary.results,
                                           tmp%>%
                                             mutate(
                                               "Classifier" = names(ml.train.negated.triage.binary[c]),
                                               "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                             )%>%
                                             select(
                                               c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                 "Kappa","Sensitivity","Specificity",
                                                 "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                             )
    )
    rm(tmp, cm, predictions)
  }
}
rm(c,i,dat)
#-----------------------
# PROVIDER NOTES Negated
#-----------------------
#rm(list=ls(pattern="^provider"))
options(readr.default_locale=readr::locale(tz="US/Eastern"))
negated.provider.cases <- read_csv("DEID_STROKE_CASES_INTERESTED_PROVIDER_NOTES_FINAL_FORMATTED_NEGATED.csv")%>%
  select(-c('pt_id','enc_dt'))
negated.provider.controls.1 <- read_csv("DEID_STROKE_CONTROL_1_INTERESTED_PROVIDER_NOTES_FINAL_FORMATTED_NEGATED.csv")%>%
  select(-c('pt_id','enc_dt'))
negated.provider.controls.3 <- read_csv("DEID_STROKE_CONTROL_3_INTERESTED_PROVIDER_NOTES_FINAL_FORMATTED_NEGATED.csv")%>%
  select(-c('pt_id','enc_dt'))
negated.provider.controls.1[setdiff(names(negated.provider.controls.3), names(negated.provider.controls.1))] <- 0
negated.provider.controls.3[setdiff(names(negated.provider.controls.1), names(negated.provider.controls.3))] <- 0
negated.provider.controls <- rbind(negated.provider.controls.1, negated.provider.controls.3)
negated.provider.cases[setdiff(names(negated.provider.controls), names(negated.provider.cases))] <- 0
negated.provider.controls[setdiff(names(negated.provider.cases), names(negated.provider.controls))] <- 0
negated.provider <- rbind(negated.provider.cases, negated.provider.controls ) 
negated.provider.svd <- svd(negated.provider)
negated.provider.approx <- list()
negated.provider.approx[[1]] <- negated.provider.svd$u[,1:31] %*% diag(negated.provider.svd$d)[1:31,1:31] %*% t((negated.provider.svd$v)[1:31,1:31])
negated.provider.approx[[2]] <- negated.provider.svd$u[,1:62] %*% diag(negated.provider.svd$d)[1:62,1:62] %*% t((negated.provider.svd$v)[1:62,1:62])
negated.provider.approx[[3]] <- negated.provider.svd$u[,1:124] %*% diag(negated.provider.svd$d)[1:124,1:124] %*% t((negated.provider.svd$v)[1:124,1:124])
negated.provider.approx[[4]] <- negated.provider.svd$u[,1:291] %*% diag(negated.provider.svd$d)[1:291,1:291] %*% t((negated.provider.svd$v)[1:291,1:291])
negated.provider.approx[[5]] <- negated.provider.svd$u %*% diag(negated.provider.svd$d) %*% t((negated.provider.svd$v))
names(negated.provider.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(negated.provider.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(negated.provider.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

negated.provider.trainSet <- list()
negated.provider.testSet <- list()
for (i in 1:length(negated.provider.approx)) {
  negated.provider.approx[[i]] <- cbind(negated.provider.approx[[i]], label)%>%
    select(label,everything())
  rownames(negated.provider.approx[[i]]) <- NULL
  index <- caret::createDataPartition(negated.provider.approx[[i]]$label, p = .8, list = F)
  negated.provider.trainSet[[i]] <- negated.provider.approx[[i]][index,]
  negated.provider.testSet[[i]] <- negated.provider.approx[[i]][-index,]
  rm(index)
}
names(negated.provider.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
names(negated.provider.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
ml.train.negated.provider <- list()
# Model development
for(i in 1:length(negated.provider.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(negated.provider.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(negated.provider.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.negated.provider[[paste0(names(negated.provider.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                           data = negated.provider.trainSet[[i]],
                                                                                                           method = classifiers[[c]],
                                                                                                           preProcess = c("center","scale"),
                                                                                                           metric = "ROC",
                                                                                                           trControl = paramGrid,
                                                                                                           tuneLength = 3)
    finishTime[[paste0(names(negated.provider.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(negated.provider.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
negated.provider.results <- as.data.frame(list())
rm(c,i,dat)
for (dat in 1:length(negated.provider.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        negated.provider.testSet[[dat]]$label,
        predict(object = ml.train.negated.provider[[c]], negated.provider.testSet[[dat]], type = "prob"),
        predict(object = ml.train.negated.provider[[c]], negated.provider.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.provider[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.provider[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.provider[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,        # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                                # Calculate AUPR
    )))
    negated.provider.results <- rbind(negated.provider.results,
                                      tmp%>%
                                        mutate(
                                          "Classifier" = names(ml.train.negated.provider[c]),
                                          "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                        )%>%
                                        select(
                                          c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                            "Kappa","Sensitivity","Specificity",
                                            "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                        )
    )
    rm(tmp, cm, predictions)
  }
}
rm(c,i,dat)
#-----------------------
# Provider: Binary Negated
#-----------------------
negated.provider.binary <- negated.provider %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
negated.provider.binary.svd <- svd(negated.provider.binary)
negated.provider.binary.approx <- list()
negated.provider.binary.approx[[1]] <- negated.provider.binary.svd$u[,1:31] %*% diag(negated.provider.binary.svd$d)[1:31,1:31] %*% t((negated.provider.binary.svd$v)[1:31,1:31])
negated.provider.binary.approx[[2]] <- negated.provider.binary.svd$u[,1:62] %*% diag(negated.provider.binary.svd$d)[1:62,1:62] %*% t((negated.provider.binary.svd$v)[1:62,1:62])
negated.provider.binary.approx[[3]] <- negated.provider.binary.svd$u[,1:124] %*% diag(negated.provider.binary.svd$d)[1:124,1:124] %*% t((negated.provider.binary.svd$v)[1:124,1:124])
negated.provider.binary.approx[[4]] <- negated.provider.binary.svd$u[,1:291] %*% diag(negated.provider.binary.svd$d)[1:291,1:291] %*% t((negated.provider.binary.svd$v)[1:291,1:291])
negated.provider.binary.approx[[5]] <- negated.provider.binary.svd$u %*% diag(negated.provider.binary.svd$d) %*% t((negated.provider.binary.svd$v))
names(negated.provider.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(negated.provider.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(negated.provider.controls), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)

negated.provider.binary.trainSet <- list()
negated.provider.binary.testSet <- list()
for (i in 1:length(negated.provider.binary.approx)) {
  negated.provider.binary.approx[[i]] <- cbind(negated.provider.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(negated.provider.binary.approx[[i]]) <- NULL
  index <- caret::createDataPartition(negated.provider.binary.approx[[i]]$label, p = .8, list = F)
  negated.provider.binary.trainSet[[i]] <- negated.provider.binary.approx[[i]][index,]
  negated.provider.binary.testSet[[i]] <- negated.provider.binary.approx[[i]][-index,]
  rm(index)
}
names(negated.provider.binary.trainSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
names(negated.provider.binary.testSet) <- c('approx.31', 'approx.62', 'approx.124', 'approx.291', 'approx.583')
ml.train.negated.provider.binary <- list()
# Model development
for(i in 1:length(negated.provider.binary.trainSet)) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(negated.provider.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    startTime[[paste0(names(negated.provider.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    ml.train.negated.provider.binary[[paste0(names(negated.provider.binary.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                         data = negated.provider.binary.trainSet[[i]],
                                                                                                                         method = classifiers[[c]],
                                                                                                                         preProcess = c("center","scale"),
                                                                                                                         metric = "ROC",
                                                                                                                         trControl = paramGrid,
                                                                                                                         tuneLength = 3)
    finishTime[[paste0(names(negated.provider.binary.trainSet[i]),".",classifiers[[c]])]] <- Sys.time()
    print(paste("finished:",names(negated.provider.binary.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
negated.provider.binary.results <- as.data.frame(list())
rm(c,i,dat)
for (dat in 1:length(negated.provider.binary.testSet)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        negated.provider.binary.testSet[[dat]]$label,
        predict(object = ml.train.negated.provider.binary[[c]], negated.provider.binary.testSet[[dat]], type = "prob"),
        predict(object = ml.train.negated.provider.binary[[c]], negated.provider.binary.testSet[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.provider.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.provider.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.provider.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,        # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                                # Calculate AUPR
    )))
    negated.provider.binary.results <- rbind(negated.provider.binary.results,
                                             tmp%>%
                                               mutate(
                                                 "Classifier" = names(ml.train.negated.provider.binary[c]),
                                                 "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                               )%>%
                                               select(
                                                 c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                   "Kappa","Sensitivity","Specificity",
                                                   "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                               )
    )
    rm(tmp, cm, predictions)
  }
}
rm(c,i,dat)

################################################################################
################################################################################
##########################      ONLY Controls 3    #############################
################################################################################
triage.controls.3[setdiff(names(triage.cases), names(triage.controls.3))] <- 0
triage.cases[setdiff(names(triage.controls.3), names(triage.cases))] <- 0
triage.w.controls3.only <- rbind(
  triage.cases %>%
    mutate(label = 'X1'),
  triage.controls.3 %>%
    mutate(label = 'X0')
)
index <- caret::createDataPartition(triage.w.controls3.only$label, p = .8, list = F)
trainSet.w.controls3.only <- triage.w.controls3.only[index,]
testSet.w.controls3.only <- triage.w.controls3.only[-index,]
testSet.w.controls3.only$label <- factor(testSet.w.controls3.only$label)
paramGrid<-trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        summaryFunction = twoClassSummary,                      # Evaluate performance
                        classProbs = T,                                         # Estimate class probabilities
                        allowParallel = T,
                        search = "random")
ml.train.triage.w.controls3.only <- list()
for (c in 1:length(classifiers)) {
  ml.train.triage.w.controls3.only[[classifiers[[c]]]] <- train(label~.,
                                                                data = trainSet.w.controls3.only,
                                                                method = classifiers[[c]],
                                                                preProcess = c("center","scale"),
                                                                metric = "ROC",
                                                                trControl = paramGrid,
                                                                tuneLength = 3)
}
triage.w.controls3.only.results <- as.data.frame(list())
for (c in 1:length(ml.train.triage.w.controls3.only)) {
  predictions <- setNames(
    data.frame(
      testSet.w.controls3.only$label,
      predict(object = ml.train.triage.w.controls3.only[[c]], testSet.w.controls3.only, type = "prob"),
      predict(object = ml.train.triage.w.controls3.only[[c]], testSet.w.controls3.only, type = "raw")
    ),
    c("obs","X0","X1","pred")
  )
  cm <- confusionMatrix(
    reference = predictions$obs,
    data = predictions$pred,
    mode = "everything",
    positive = "X1"
  )
  tmp <- as.data.frame(t(rbind(
    fetchResults(cm$byClass, ml.train.triage.w.controls3.only[c]),                                                              # Fetch Recall,Specificity,Precision
    fetchResults(cm$overall, ml.train.triage.w.controls3.only[c]),                                                              # Fetch Accuracy,95%CI
    fetchResults(as.data.frame(cm$table)$Freq, ml.train.triage.w.controls3.only[c]),                                             # Fetch TP,FP,FN,TN
    roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,        # Calculate AUROC
    prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                                # Calculate AUPR
  )))
  triage.w.controls3.only.results <- rbind(triage.w.controls3.only.results,
                                           tmp%>%
                                             mutate(
                                               "Classifier" = names(ml.train.triage.w.controls3.only[c]),
                                               "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                             )%>%
                                             select(
                                               c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                 "Kappa","Sensitivity","Specificity",
                                                 "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                             )
  )
  rm(tmp, cm, predictions)
}
#-----------------------
## Simulations with 60,70,80 splits
#-----------------------
#-----------------------
# TRIAGE
#-----------------------
triage.controls.3[setdiff(names(triage.cases), names(triage.controls.3))] <- 0
triage.cases[setdiff(names(triage.controls.3), names(triage.cases))] <- 0
triage.controls3.svd <- svd(rbind(
  triage.cases,
  triage.controls.3))
triage.approx.controls3 <- list()
triage.approx.controls3[[1]] <- triage.controls3.svd$u[,1:31] %*% diag(triage.controls3.svd$d)[1:31,1:31] %*% t((triage.controls3.svd$v)[1:31,1:31])
triage.approx.controls3[[2]] <- triage.controls3.svd$u[,1:62] %*% diag(triage.controls3.svd$d)[1:62,1:62] %*% t((triage.controls3.svd$v)[1:62,1:62])
triage.approx.controls3[[3]] <- triage.controls3.svd$u[,1:124] %*% diag(triage.controls3.svd$d)[1:124,1:124] %*% t((triage.controls3.svd$v)[1:124,1:124])
triage.approx.controls3[[4]] <- triage.controls3.svd$u[,1:306] %*% diag(triage.controls3.svd$d)[1:306,1:306] %*% t((triage.controls3.svd$v)[1:306,1:306])
triage.approx.controls3[[5]] <- triage.controls3.svd$u %*% diag(triage.controls3.svd$d) %*% t((triage.controls3.svd$v))
names(triage.approx.controls3) <- c('approx.31', 'approx.62', 'approx.124', 'approx.306', 'approx.611')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(triage.controls.3), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
for (i in 1:length(triage.approx.controls3)) {
  triage.approx.controls3[[i]] <- cbind(triage.approx.controls3[[i]], label)%>%
    select(label,everything())
  rownames(triage.approx.controls3[[i]]) <- NULL
}
p <- seq(.6, .8, by = .1)
triage.trainSet.controls3.only <- list()
triage.testSet.controls3.only <- list()
for (dat in 1:length(triage.approx.controls3)) {
  for (i in 1:length(p)) {
    index <- caret::createDataPartition(triage.approx.controls3[[dat]]$label, p = p[i], list = F)
    triage.trainSet.controls3.only[[paste0(names(triage.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i])*100)]] <- triage.approx.controls3[[dat]][index,]
    triage.testSet.controls3.only[[paste0(names(triage.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i]*100))]] <- triage.approx.controls3[[dat]][-index,]
    rm(index)
  }
}
rm(label, triage.approx.controls3, triage.cases, triage.controls.3, triage.controls3.svd, dat, i, p, triage.cases, triage.controls.3,
   triage.controls3.svd)
classifiers <- c("glm", "rf", "xgbDART")
ml.train.triage.w.controls3.simulations <- list()
# Model development
for(dat in 1:length(triage.trainSet.controls3.only)) {
  for(c in 1:length(classifiers)) {
    print(paste("started: ",names(triage.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    ml.train.triage.w.controls3.simulations[[paste0(names(triage.trainSet.controls3.only[dat]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                                data = triage.trainSet.controls3.only[[dat]],
                                                                                                                                method = classifiers[[c]],
                                                                                                                                preProcess = c("center","scale"),
                                                                                                                                metric = "ROC",
                                                                                                                                trControl = paramGrid,
                                                                                                                                tuneLength = 3)
    print(paste("finished: ",names(triage.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
rm(triage.trainSet.controls3.only)
# Model validation
fetchResults<-function(x,y){
  z<-as.data.frame(x)
  colnames(z)<-names(y)
  return(z)
}
incrementStart<-function(x){3*x-2}
incrementEnd<-function(x){3*x}
results.triage.w.controls3.simulations <- as.data.frame(list())
for (dat in 1:length(triage.testSet.controls3.only)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        triage.testSet.controls3.only[[dat]]$label,
        predict(object = ml.train.triage.w.controls3.simulations[[c]], triage.testSet.controls3.only[[dat]], type = "prob"),
        predict(object = ml.train.triage.w.controls3.simulations[[c]], triage.testSet.controls3.only[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.w.controls3.simulations[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.w.controls3.simulations[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.w.controls3.simulations[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    results.triage.w.controls3.simulations <- rbind(results.triage.w.controls3.simulations,
                                                    tmp%>%
                                                      mutate(
                                                        "Classifier" = names(ml.train.triage.w.controls3.simulations[c]),
                                                        "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                      )%>%
                                                      select(
                                                        c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                          "Kappa","Sensitivity","Specificity",
                                                          "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                      )
    )
    rm(tmp, cm, predictions)
  }
}
#-----------------------
# TRIAGE: Binary
#-----------------------
triage.controls.3[setdiff(names(triage.cases), names(triage.controls.3))] <- 0
triage.cases[setdiff(names(triage.controls.3), names(triage.cases))] <- 0
triage <- rbind(
  triage.cases%>%
    mutate(label = 'Case'),
  triage.controls.3%>%
    mutate(label = 'Control3'))%>%
  select(label, everything())
triage.binary <- triage[,-1]%>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
triage.binary <- cbind(triage[,1],triage.binary)%>%
  select(label, everything())
colnames(triage.binary[,1]) <- 'label'
triage.controls3.svd <- svd(triage.binary)
triage.binary.approx.controls3 <- list()
triage.binary.approx.controls3[[1]] <- triage.controls3.svd$u[,1:31] %*% diag(triage.controls3.svd$d)[1:31,1:31] %*% t((triage.controls3.svd$v)[1:31,1:31])
triage.binary.approx.controls3[[2]] <- triage.controls3.svd$u[,1:62] %*% diag(triage.controls3.svd$d)[1:62,1:62] %*% t((triage.controls3.svd$v)[1:62,1:62])
triage.binary.approx.controls3[[3]] <- triage.controls3.svd$u[,1:124] %*% diag(triage.controls3.svd$d)[1:124,1:124] %*% t((triage.controls3.svd$v)[1:124,1:124])
triage.binary.approx.controls3[[4]] <- triage.controls3.svd$u[,1:306] %*% diag(triage.controls3.svd$d)[1:306,1:306] %*% t((triage.controls3.svd$v)[1:306,1:306])
triage.binary.approx.controls3[[5]] <- triage.controls3.svd$u %*% diag(triage.controls3.svd$d) %*% t((triage.controls3.svd$v))
names(triage.binary.approx.controls3) <- c('approx.31', 'approx.62', 'approx.124', 'approx.306', 'approx.611')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(triage.controls.3), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
for (i in 1:length(triage.binary.approx.controls3)) {
  triage.binary.approx.controls3[[i]] <- cbind(triage.binary.approx.controls3[[i]], label)%>%
    select(label,everything())
  rownames(triage.binary.approx.controls3[[i]]) <- NULL
}
p <- seq(.6, .8, by = .1)
triage.binary.trainSet.controls3.only <- list()
triage.binary.testSet.controls3.only <- list()
for (dat in 1:length(triage.binary.approx.controls3)) {
  for (i in 1:length(p)) {
    index <- caret::createDataPartition(triage.binary.approx.controls3[[dat]]$label, p = p[i], list = F)
    triage.binary.trainSet.controls3.only[[paste0(names(triage.binary.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i])*100)]] <- triage.binary.approx.controls3[[dat]][index,]
    triage.binary.testSet.controls3.only[[paste0(names(triage.binary.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i]*100))]] <- triage.binary.approx.controls3[[dat]][-index,]
    rm(index)
  }
}
rm(label, triage.binary.approx.controls3, triage.cases, triage.controls.3, triage.binary.controls3.svd, dat, i, p, 
   triage.binary.cases, triage.binary.controls.3, triage.controls3.svd, triage.binary)
ml.train.triage.binary.w.controls3.simulations <- list()
# Model development
for(dat in 1:length(triage.binary.trainSet.controls3.only)) {
  for(c in 1:length(classifiers)) {
    print(paste("started: ",names(triage.binary.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    ml.train.triage.binary.w.controls3.simulations[[paste0(names(triage.binary.trainSet.controls3.only[dat]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                                              data = triage.binary.trainSet.controls3.only[[dat]],
                                                                                                                                              method = classifiers[[c]],
                                                                                                                                              preProcess = c("center","scale"),
                                                                                                                                              metric = "ROC",
                                                                                                                                              trControl = paramGrid,
                                                                                                                                              tuneLength = 3)
    print(paste("finished: ",names(triage.binary.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
results.triage.binary.w.controls3.simulations <- as.data.frame(list())
for (dat in 1:length(triage.binary.testSet.controls3.only)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        triage.binary.testSet.controls3.only[[dat]]$label,
        predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], triage.binary.testSet.controls3.only[[dat]], type = "prob"),
        predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], triage.binary.testSet.controls3.only[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.binary.w.controls3.simulations[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.binary.w.controls3.simulations[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.binary.w.controls3.simulations[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    results.triage.binary.w.controls3.simulations <- rbind(results.triage.binary.w.controls3.simulations,
                                                           tmp%>%
                                                             mutate(
                                                               "Classifier" = names(ml.train.triage.binary.w.controls3.simulations[c]),
                                                               "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                             )%>%
                                                             select(
                                                               c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                                 "Kappa","Sensitivity","Specificity",
                                                                 "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                             )
    )
    rm(tmp, cm, predictions)
  }
}
#-----------------------
# TRIAGE wo Negation
#-----------------------
negated.triage.controls.3[setdiff(names(negated.triage.cases), names(negated.triage.controls.3))] <- 0
negated.triage.cases[setdiff(names(negated.triage.controls.3), names(negated.triage.cases))] <- 0
triage.controls3.svd <- svd(rbind(
  negated.triage.cases,
  negated.triage.controls.3))
triage.approx.controls3 <- list()
triage.approx.controls3[[1]] <- triage.controls3.svd$u[,1:31] %*% diag(triage.controls3.svd$d)[1:31,1:31] %*% t((triage.controls3.svd$v)[1:31,1:31])
triage.approx.controls3[[2]] <- triage.controls3.svd$u[,1:62] %*% diag(triage.controls3.svd$d)[1:62,1:62] %*% t((triage.controls3.svd$v)[1:62,1:62])
triage.approx.controls3[[3]] <- triage.controls3.svd$u[,1:124] %*% diag(triage.controls3.svd$d)[1:124,1:124] %*% t((triage.controls3.svd$v)[1:124,1:124])
triage.approx.controls3[[4]] <- triage.controls3.svd$u[,1:306] %*% diag(triage.controls3.svd$d)[1:306,1:306] %*% t((triage.controls3.svd$v)[1:306,1:306])
triage.approx.controls3[[5]] <- triage.controls3.svd$u %*% diag(triage.controls3.svd$d) %*% t((triage.controls3.svd$v))
names(triage.approx.controls3) <- c('approx.31', 'approx.62', 'approx.124', 'approx.306', 'approx.611')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(negated.triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(negated.triage.controls.3), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
for (i in 1:length(triage.approx.controls3)) {
  triage.approx.controls3[[i]] <- cbind(triage.approx.controls3[[i]], label)%>%
    select(label,everything())
  rownames(triage.approx.controls3[[i]]) <- NULL
}
p <- seq(.6, .8, by = .1)
negated.triage.trainSet.controls3.only <- list()
negated.triage.testSet.controls3.only <- list()
for (dat in 1:length(triage.approx.controls3)) {
  for (i in 1:length(p)) {
    index <- caret::createDataPartition(triage.approx.controls3[[dat]]$label, p = p[i], list = F)
    negated.triage.trainSet.controls3.only[[paste0(names(triage.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i])*100)]] <- triage.approx.controls3[[dat]][index,]
    negated.triage.testSet.controls3.only[[paste0(names(triage.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i]*100))]] <- triage.approx.controls3[[dat]][-index,]
    rm(index)
  }
}
rm(label, triage.approx.controls3, negated.triage.cases, negated.triage.controls.3, triage.controls3.svd, dat, i, p, triage.cases, 
   triage.controls3.svd)
classifiers <- c("glm", "rf", "xgbDART")
ml.train.neg.triage.w.controls3.simulations <- list()
# Model development
for(dat in 1:length(negated.triage.trainSet.controls3.only)) {
  for(c in 1:length(classifiers)) {
    print(paste("started: ",names(negated.triage.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    ml.train.neg.triage.w.controls3.simulations[[paste0(names(negated.triage.trainSet.controls3.only[dat]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                                            data = negated.triage.trainSet.controls3.only[[dat]],
                                                                                                                                            method = classifiers[[c]],
                                                                                                                                            preProcess = c("center","scale"),
                                                                                                                                            metric = "ROC",
                                                                                                                                            trControl = paramGrid,
                                                                                                                                            tuneLength = 3)
    print(paste("finished: ",names(negated.triage.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
results.neg.triage.w.controls3.simulations <- as.data.frame(list())
for (dat in 1:length(negated.triage.testSet.controls3.only)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        negated.triage.testSet.controls3.only[[dat]]$label,
        predict(object = ml.train.neg.triage.w.controls3.simulations[[c]], negated.triage.testSet.controls3.only[[dat]], type = "prob"),
        predict(object = ml.train.neg.triage.w.controls3.simulations[[c]], negated.triage.testSet.controls3.only[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.neg.triage.w.controls3.simulations[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.neg.triage.w.controls3.simulations[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.neg.triage.w.controls3.simulations[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    results.neg.triage.w.controls3.simulations <- rbind(results.neg.triage.w.controls3.simulations,
                                                        tmp%>%
                                                          mutate(
                                                            "Classifier" = names(ml.train.neg.triage.w.controls3.simulations[c]),
                                                            "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                          )%>%
                                                          select(
                                                            c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                              "Kappa","Sensitivity","Specificity",
                                                              "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                          )
    )
    rm(tmp, cm, predictions)
  }
}
#-----------------------
# TRIAGE: Binary wo negation
#-----------------------
negated.triage.controls.3[setdiff(names(negated.triage.cases), names(negated.triage.controls.3))] <- 0
negated.triage.cases[setdiff(names(negated.triage.controls.3), names(negated.triage.cases))] <- 0
triage.binary <- rbind(
  negated.triage.cases,
  negated.triage.controls.3)%>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
triage.controls3.svd <- svd(triage.binary)
triage.binary.approx.controls3 <- list()
triage.binary.approx.controls3[[1]] <- triage.controls3.svd$u[,1:31] %*% diag(triage.controls3.svd$d)[1:31,1:31] %*% t((triage.controls3.svd$v)[1:31,1:31])
triage.binary.approx.controls3[[2]] <- triage.controls3.svd$u[,1:62] %*% diag(triage.controls3.svd$d)[1:62,1:62] %*% t((triage.controls3.svd$v)[1:62,1:62])
triage.binary.approx.controls3[[3]] <- triage.controls3.svd$u[,1:124] %*% diag(triage.controls3.svd$d)[1:124,1:124] %*% t((triage.controls3.svd$v)[1:124,1:124])
triage.binary.approx.controls3[[4]] <- triage.controls3.svd$u[,1:306] %*% diag(triage.controls3.svd$d)[1:306,1:306] %*% t((triage.controls3.svd$v)[1:306,1:306])
triage.binary.approx.controls3[[5]] <- triage.controls3.svd$u %*% diag(triage.controls3.svd$d) %*% t((triage.controls3.svd$v))
names(triage.binary.approx.controls3) <- c('approx.31', 'approx.62', 'approx.124', 'approx.306', 'approx.611')
label <-data.frame(rbind(
  t(data.frame(lapply(1:nrow(negated.triage.cases), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(negated.triage.controls.3), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
for (i in 1:length(triage.binary.approx.controls3)) {
  triage.binary.approx.controls3[[i]] <- cbind(triage.binary.approx.controls3[[i]], label)%>%
    select(label,everything())
  rownames(triage.binary.approx.controls3[[i]]) <- NULL
}
p <- seq(.6, .8, by = .1)
negated.triage.binary.trainSet.controls3.only <- list()
negated.triage.binary.testSet.controls3.only <- list()
for (dat in 1:length(triage.binary.approx.controls3)) {
  for (i in 1:length(p)) {
    index <- caret::createDataPartition(triage.binary.approx.controls3[[dat]]$label, p = p[i], list = F)
    negated.triage.binary.trainSet.controls3.only[[paste0(names(triage.binary.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i])*100)]] <- triage.binary.approx.controls3[[dat]][index,]
    negated.triage.binary.testSet.controls3.only[[paste0(names(triage.binary.approx.controls3[dat]),".split.",p[i]*100,':',(1-p[i]*100))]] <- triage.binary.approx.controls3[[dat]][-index,]
    rm(index)
  }
}
rm(label, triage.binary.approx.controls3, triage.cases, triage.controls.3, triage.binary.controls3.svd, dat, i, p, 
   triage.binary.cases, triage.binary.controls.3, triage.controls3.svd, triage.binary)
ml.train.neg.triage.binary.w.controls3.simulations <- list()
# Model development
for(dat in 1:length(negated.triage.binary.trainSet.controls3.only)) {
  for(c in 1:length(classifiers)) {
    print(paste("started: ",names(negated.triage.binary.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    ml.train.neg.triage.binary.w.controls3.simulations[[paste0(names(negated.triage.binary.trainSet.controls3.only[dat]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                                                          data = negated.triage.binary.trainSet.controls3.only[[dat]],
                                                                                                                                                          method = classifiers[[c]],
                                                                                                                                                          preProcess = c("center","scale"),
                                                                                                                                                          metric = "ROC",
                                                                                                                                                          trControl = paramGrid,
                                                                                                                                                          tuneLength = 3)
    print(paste("finished: ",names(negated.triage.binary.trainSet.controls3.only[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
# Model validation
results.neg.triage.binary.w.controls3.simulations <- as.data.frame(list())
for (dat in 1:length(negated.triage.binary.testSet.controls3.only)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        negated.triage.binary.testSet.controls3.only[[dat]]$label,
        predict(object = ml.train.neg.triage.binary.w.controls3.simulations[[c]], negated.triage.binary.testSet.controls3.only[[dat]], type = "prob"),
        predict(object = ml.train.neg.triage.binary.w.controls3.simulations[[c]], negated.triage.binary.testSet.controls3.only[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.neg.triage.binary.w.controls3.simulations[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.neg.triage.binary.w.controls3.simulations[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.neg.triage.binary.w.controls3.simulations[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    results.neg.triage.binary.w.controls3.simulations <- rbind(results.neg.triage.binary.w.controls3.simulations,
                                                               tmp%>%
                                                                 mutate(
                                                                   "Classifier" = names(ml.train.neg.triage.binary.w.controls3.simulations[c]),
                                                                   "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                                 )%>%
                                                                 select(
                                                                   c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                                     "Kappa","Sensitivity","Specificity",
                                                                     "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                                 )
    )
    rm(tmp, cm, predictions)
  }
}





#-----------------------
# TRIAGE: No SVD
#-----------------------
triage.controls.3[setdiff(names(triage.cases), names(triage.controls.3))] <- 0
triage.controls.3$label <- 'X0'
triage.cases[setdiff(names(triage.controls.3), names(triage.cases))] <- 0
triage.cases$label <- 'X1'
triage <- rbind(triage.cases, triage.controls.3)
triage$label <- as.factor(triage$label)
p <- seq(.6, .8, by = .1)
triage.trainSet.controls3.nosvd <- list()
triage.testSet.controls3.nosvd <- list()
for (i in 1:length(p)) {
  index <- caret::createDataPartition(triage$label, p = p[i], list = F)
  triage.trainSet.controls3.nosvd[[paste0("split.",p[i]*100,':',(1-p[i])*100)]] <- triage[index,]
  triage.testSet.controls3.nosvd[[paste0("split.",p[i]*100,':',(1-p[i])*100)]]  <- triage[-index,]
  rm(index)
}
rm(label, triage.approx.controls3, triage.cases, triage.controls.3, triage.controls3.svd, dat, i, p, triage.cases, triage.controls.3,
   triage.controls3.svd)
classifiers <- c("glm", "rf", "xgbDART")
ml.train.triage.w.controls3.nosvd <- list()
# Model development
for(dat in 1:length(triage.trainSet.controls3.nosvd)) {
  for(c in 1:length(classifiers)) {
    print(paste("started: ",names(triage.trainSet.controls3.nosvd[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    ml.train.triage.w.controls3.nosvd[[paste0(names(triage.trainSet.controls3.nosvd[dat]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                           data = triage.trainSet.controls3.nosvd[[dat]],
                                                                                                                           method = classifiers[[c]],
                                                                                                                           preProcess = c("center","scale"),
                                                                                                                           metric = "ROC",
                                                                                                                           trControl = paramGrid,
                                                                                                                           tuneLength = 3)
    print(paste("finished: ",names(triage.trainSet.controls3.nosvd[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
#-----------------------
# TRIAGE Binary: NOSVD
#-----------------------
triage <- triage %>%
  select(label, everything())
triage.binary <- triage[,-1]%>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
triage.binary <- cbind(triage[,1], triage.binary)
p <- seq(.6, .8, by = .1)
triage.binary.trainSet.controls3.nosvd <- list()
triage.binary.testSet.controls3.nosvd <- list()
for (i in 1:length(p)) {
  index <- caret::createDataPartition(triage.binary$label, p = p[i], list = F)
  triage.binary.trainSet.controls3.nosvd[[paste0("split.",p[i]*100,':',(1-p[i])*100)]] <- triage.binary[index,]
  triage.binary.testSet.controls3.nosvd[[paste0("split.",p[i]*100,':',(1-p[i])*100)]]  <- triage.binary[-index,]
  rm(index)
}
ml.train.triage.binary.w.controls3.nosvd <- list()
# Model development
for(dat in 1:length(triage.binary.trainSet.controls3.nosvd)) {
  for(c in 1:length(classifiers)) {
    print(paste("started: ",names(triage.binary.trainSet.controls3.nosvd[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    ml.train.triage.binary.w.controls3.nosvd[[paste0(names(triage.binary.trainSet.controls3.nosvd[dat]),".",classifiers[[c]])]] <- train(label~.,
                                                                                                                                         data = triage.binary.trainSet.controls3.nosvd[[dat]],
                                                                                                                                         method = classifiers[[c]],
                                                                                                                                         preProcess = c("center","scale"),
                                                                                                                                         metric = "ROC",
                                                                                                                                         trControl = paramGrid,
                                                                                                                                         tuneLength = 3)
    print(paste("finished: ",names(triage.binary.trainSet.controls3.nosvd[dat]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
rm(triage.trainSet.controls3.only)
#-----------------------
# Introduce Noise 
#-----------------------
# Model validation of Triage Binary 80:20 split for ALL cuis
random.index <- as.double(sample(seq_len(ncol(tmp)), 10, replace = F))
noise.10.triage.binary.cui.611.split.80.20 <- cbind(
  data.frame(sdcMicro::addNoise(triage.binary.testSet.controls3.only[[15]]%>%
                                  select(-label), 
                                #variables = random.index,
                                noise = 10, 
                                method = 'additive')$xm)%>%
    rename_all(.funs = funs(sub("X","",.))),
  data.frame(triage.binary.testSet.controls3.only[[15]]$label)%>%
    rename_at(ncol(.),~'label')
)%>%
  select(label, everything())
results.noise.10.triage.binary.cui.611.split.80.20 <- as.data.frame(list())
for (c in 43:45) {
  predictions <- setNames(
    data.frame(
      noise.10.triage.binary.cui.611.split.80.20$label,
      predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], noise.10.triage.binary.cui.611.split.80.20, type = "prob"),
      predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], noise.10.triage.binary.cui.611.split.80.20, type = "raw")
    ),
    c("obs","X0","X1","pred")
  )
  cm <- confusionMatrix(
    reference = predictions$obs,
    data = predictions$pred,
    mode = "everything",
    positive = "X1"
  )
  tmp <- as.data.frame(t(rbind(
    fetchResults(cm$byClass, ml.train.triage.binary.w.controls3.simulations[c]),                                                              # Fetch Recall,Specificity,Precision
    fetchResults(cm$overall, ml.train.triage.binary.w.controls3.simulations[c]),                                                              # Fetch Accuracy,95%CI
    fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.binary.w.controls3.simulations[c]),                                             # Fetch TP,FP,FN,TN
    roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
    prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
  )))
  results.noise.10.triage.binary.cui.611.split.80.20 <- rbind(results.noise.10.triage.binary.cui.611.split.80.20,
                                                              tmp%>%
                                                                mutate(
                                                                  "Classifier" = names(ml.train.triage.binary.w.controls3.simulations[c]),
                                                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                                )%>%
                                                                select(
                                                                  c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                                    "Kappa","Sensitivity","Specificity",
                                                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                                )
  )
  rm(tmp, cm, predictions)
}
results.noise.10.triage.binary.cui.611.split.80.20$noise <- 10
# Loops
noise <- function(dat, noiseper){
  return(
    cbind(
      data.frame(sdcMicro::addNoise(triage.binary.testSet.controls3.only[[dat]]%>%
                                      select(-label), 
                                    #variables = random.index,
                                    noise = noiseper, 
                                    method = 'additive')$xm)%>%
        rename_all(.funs = funs(sub("X","",.))),
      data.frame(triage.binary.testSet.controls3.only[[dat]]$label)%>%
        rename_at(ncol(.),~'label')
    )%>%
      select(label, everything())
  )
}



#-----------------------
# Extract Feature importance
#-----------------------
cui_to_english <- read_excel("cui_to_english.xlsx")
feaImp.binary.nosvd <- list()
for (c in 1:length(ml.train.triage.binary.w.controls3.nosvd)) {
  tmp <- varImp(ml.train.triage.binary.w.controls3.nosvd[[c]])$importance
  tmp$cui <- rownames(tmp)
  tmp$classifier <- names(ml.train.triage.binary.w.controls3.nosvd[c])
  #feaImp.binary.nosvd[[c]] <- tmp
  feaImp.binary.nosvd <- rbind(feaImp.binary.nosvd, tmp)
  #rm(tmp)
}
x <- reduce(feaImp.binary.nosvd, full_join, by = 'cui')
openxlsx::write.xlsx(feaImp.binary.nosvd,"results1.xlsx",row.names=F,quote=F)
feaImp.binary.nosvd <- read_xlsx("feaImpScore_triage.binary.nosvd.xlsx", sheet = 'List')
feaImp.binary.nosvd$split <- as.factor(feaImp.binary.nosvd$split)
feaImp.binary.nosvd$classifier <- as.factor(feaImp.binary.nosvd$classifier)
tmp <- reshape2::melt(feaImp.binary.nosvd)
feaImp.binary.nosvd%>%
  arrange(desc(score))%>%
  group_by(split, classifier)%>%
  top_n(10, score)%>%
  arrange(desc(score),split, classifier)%>%
  ggplot(aes(y = cui, x = classifier, fill = score))+
  geom_tile()+
  facet_wrap(split~.)+
  scale_fill_continuous(trans = 'reverse')

tmp <- triage[,-1] %>%
  mutate_all(funs(ifelse(.==0,0,1)))
tmp <- cbind(triage[,1], tmp)
tmp <- tmp%>%
  group_by(label)%>%
  summarise_all(sum)
tmp <- t(tmp)
colnames(tmp) <- tmp[1,]
tmp <- data.frame(tmp[-1,])
tmp$cui <- rownames(tmp)
tmp$X0 <- as.numeric(tmp$X0)
tmp$X1 <- as.numeric(tmp$X1)
tmp <- left_join(tmp, cui_to_english[,c('cui','text')], by='cui')
tmp <- tmp%>%
  mutate(wtd.avg = (733*X0 + 1054*X1)/1787)
# density plots
tmp2 <- melt(tmp%>%select(cui,X0,X1), id.vars = c('cui'))
ggplot(tmp2, aes(x=value))+
  geom_histogram(aes(y=..density..),
                 binwidth = .5, colour = "blue", fill = "white")+
  geom_density(alpha=.2, fill="#FF6655")+
  facet_wrap(variable~.)

ggplot(tmp2, aes(x=value, fill=variable))+
  geom_density()+
  facet_wrap(variable~.)
# distribution plots
ggplot(tmp2, aes(x=value, fill=variable))+
  geom_histogram(binwidth = 25, alpha =.75, position = "identity")

tmp$scaled.X0 <- scale(tmp$X0, center = FALSE, scale = max(tmp$X0, na.rm = TRUE)/100)
tmp$scaled.X1 <- scale(tmp$X1, center = FALSE, scale = max(tmp$X1, na.rm = TRUE)/100)
tmp1 <- tmp%>%select('cui','scaled.X0','scaled.X1')
tmp1 <- melt(tmp1, id.vars = c('cui'))
melt(tmp)%>%
  ggplot(aes(fill=variable, y=value, x=cui)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(size = 8, angle = 90),
  )
# heatmap with Controls-3
tiff("cui_controls3_heatmap.tiff", units = "in", width = 20, height = 20, res = 300)
gplots::heatmap.2(as.matrix(t(triage.binary[,-1])),
                  dendrogram = c("row"), 
                  Colv = "NA",
                  reorderfun = function(d, w) reorder(d, w),
                  trace = "none", 
                  main = "",
                  key.title = "", 
                  keysize = .5,
                  key.xlab = "CUI value",
                  key.par = list(mar=c(5,.5,1,.5)),
                  density.info = "none",
                  margins = c(2,10.5),
                  xlab = "Patients", 
                  labCol = "",
                  cexRow = .5, 
                  ColSideColors = c(rep("gray",1054), rep("blue",733))
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Cases","Controls3"), # category labels
       col = c("gray","blue"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()
# heatmap with Controls-3 and Controls-1 
triage.cases <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CASES_INTERESTED_NOTES_FINAL_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))%>%
  mutate(label = 'CASES')
triage.controls.1 <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CONTROL_1_INTERESTED_NOTES_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))%>%
  mutate(label = 'CONTROLS1')
triage.controls.3 <- read_csv("C:/Users/vavula1/Desktop/Genentech/v1/DEID_STROKE_CONTROL_3_INTERESTED_NOTES_FORMATTED.csv")%>%
  select(-c('pt_id','enc_dt'))%>%
  mutate(label = 'CONTROLS3')
triage.controls.1[setdiff(names(triage.controls.3), names(triage.controls.1))] <- 0
triage.controls.3[setdiff(names(triage.controls.1), names(triage.controls.3))] <- 0
triage.controls <- rbind(triage.controls.1, triage.controls.3)
triage.cases[setdiff(names(triage.controls), names(triage.cases))] <- 0
triage.controls[setdiff(names(triage.cases), names(triage.controls))] <- 0
triage <- rbind(triage.cases, triage.controls)%>%
  select(label, everything())
triage.binary <- cbind(triage[,1],
                       triage[,-1] %>% 
                         mutate_all(funs(
                           ifelse( .>= 1, 1, 
                                   ifelse(. <= -1, -1, 0)))
                         ))
rm(triage.cases, triage.controls.1, triage.controls.3, triage.controls, triage)
tiff("cui_controls_1_3_heatmap1.tiff", units = "in", width = 20, height = 20, res = 300)
gplots::heatmap.2(as.matrix(t(triage.binary[,-1])),
                  dendrogram = c("row"), 
                  Colv = "NA",
                  trace = "none", 
                  main = "",
                  key.title = "", 
                  keysize = .5,
                  key.xlab = "CUI value",
                  key.par = list(mar=c(5,.5,1,.5)),
                  density.info = "none",
                  margins = c(2,10.5),
                  xlab = "Patients", 
                  labCol = "",
                  cexRow = .5, 
                  ColSideColors = c(rep("gray",1054) ,rep("steelblue",1814), rep("blue",733))
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Cases", "Controls1", "Controls3"), # category labels
       col = c("gray","steelblue", "blue"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()
################################################################################
##########################      ONLY Controls 1    #############################
################################################################################
triage.controls.1[setdiff(names(triage.cases), names(triage.controls.1))] <- 0
triage.cases[setdiff(names(triage.controls.1), names(triage.cases))] <- 0
triage.w.controls1.only <- rbind(
  triage.cases %>%
    mutate(label = 'X1'),
  triage.controls.1 %>%
    mutate(label = 'X0')
)
triage.w.controls1.only$label <- factor(triage.w.controls1.only$label)
index <- caret::createDataPartition(triage.w.controls1.only$label, p = .8, list = F)
trainSet.w.controls1.only <- triage.w.controls1.only[index,]
testSet.w.controls1.only <- triage.w.controls1.only[-index,]
ml.train.triage.w.controls1.only <- list()
for (c in 1:length(classifiers)) {
  ml.train.triage.w.controls1.only[[classifiers[[c]]]] <- train(label~.,
                                                                data = trainSet.w.controls1.only,
                                                                method = classifiers[[c]],
                                                                preProcess = c("center","scale"),
                                                                metric = "ROC",
                                                                trControl = paramGrid,
                                                                tuneLength = 3)
}
triage.w.controls1.only.results <- as.data.frame(list())
for (c in 1:length(ml.train.triage.w.controls1.only)) {
  predictions <- setNames(
    data.frame(
      testSet.w.controls1.only$label,
      predict(object = ml.train.triage.w.controls1.only[[c]], testSet.w.controls1.only, type = "prob"),
      predict(object = ml.train.triage.w.controls1.only[[c]], testSet.w.controls1.only, type = "raw")
    ),
    c("obs","X0","X1","pred")
  )
  cm <- confusionMatrix(
    reference = predictions$obs,
    data = predictions$pred,
    mode = "everything",
    positive = "X1"
  )
  tmp <- as.data.frame(t(rbind(
    fetchResults(cm$byClass, ml.train.triage.w.controls1.only[c]),                                                              # Fetch Recall,Specificity,Precision
    fetchResults(cm$overall, ml.train.triage.w.controls1.only[c]),                                                              # Fetch Accuracy,95%CI
    fetchResults(as.data.frame(cm$table)$Freq, ml.train.triage.w.controls1.only[c]),                                             # Fetch TP,FP,FN,TN
    roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,        # Calculate AUROC
    prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                                # Calculate AUPR
  )))
  triage.w.controls1.only.results <- rbind(triage.w.controls1.only.results,
                                           tmp%>%
                                             mutate(
                                               "Classifier" = names(ml.train.triage.w.controls1.only[c]),
                                               "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                             )%>%
                                             select(
                                               c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                 "Kappa","Sensitivity","Specificity",
                                                 "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                             )
  )
  rm(tmp, cm, predictions)
}

openxlsx::write.xlsx(triage.w.controls3.only.results,"results.xlsx",row.names=F,quote=F)

#-----------------------
# testing
#-----------------------
tmp.cases[setdiff(names(tmp.controls), names(tmp.cases))] <- 0
tmp.controls[setdiff(names(tmp.cases), names(tmp.controls))] <- 0
tmp<-rbind(tmp.cases,tmp.controls)
A <- svd(tmp)
A.test <- A$u%*%diag(A$d)%*%t(A$v)
sum(diag(A$d)[-1,])
sum(diag(A$d)[-1,][1:250,])
A.approx <- A$u[,1:200] %*% diag(A$d)[1:200,1:200] %*% t((A$v)[1:200,1:200])
A.approx <- A$u[,1:100] %*% diag(A$d)[1:100,1:100] %*% t((A$v)[1:100,1:100])
#1
tidytext::tidy(tmp, matrix = 'gamma')

#2
A <- svd(tmp)
xxd <- A$v %*% sqrt(diag(A$d))
x1 <- xxd[,1]
y1 <- xxd[,2]
plot(x1, y1)
text(x1+20, y1, labels = tmp1$column)
plot(A$d^2/sum(A$d^2), xlim = c(0, 400), type = "b", pch = 8, xlab = "principal components", 
     ylab = "variance explained")

################################################################################
##########################  validation: 03/30/21  ##############################
################################################################################
library(readr)
library(dplyr)
library(caret)
library(pROC)
library(openxlsx)
options(readr.default_locale=readr::locale(tz="US/Eastern"))
cases.preCOVID <- read_csv("stroke_aaa_ext_val_cases_nlp_output_pre_covid_final_FORMATTED.csv")%>%
  select(-c('_c0'))
cases.postCOVID <- read_csv("stroke_aaa_ext_val_cases_nlp_output_post_covid_final_FORMATTED.csv")%>%
  select(-c('_c0'))
controls.preCOVID <- read_csv("stroke_aaa_ext_val_control_nlp_output_pre_covid_final_FORMATTED.csv")%>%
  select(-c('_c0'))
controls.postCOVID <- read_csv("stroke_aaa_ext_val_control_nlp_output_post_covid_final_FORMATTED.csv")%>%
  select(-c('_c0'))
#-----------------------
# TRIAGE NOTES
#-----------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.preCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.preCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.preCOVID[setdiff(names(controls.preCOVID), names(cases.preCOVID))] <- 0
controls.preCOVID[setdiff(names(cases.preCOVID), names(controls.preCOVID))] <- 0
preCOVID <- rbind(cases.preCOVID, controls.preCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(preCOVID))
tmp <- data.frame(matrix( , nrow = nrow(preCOVID), 
                          ncol = length(setdiff(colnames(triage), colnames(preCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
preCOVID <- cbind(preCOVID, tmp)
rm(tmp)
preCOVID.svd <- svd(preCOVID)
preCOVID.approx <- list()
preCOVID.approx[[1]] <- data.frame(preCOVID.svd$u[,1:31] %*% diag(preCOVID.svd$d)[1:31,1:31] %*% t((preCOVID.svd$v)[1:31,1:31]))
preCOVID.approx[[2]] <- data.frame(preCOVID.svd$u[,1:62] %*% diag(preCOVID.svd$d)[1:62,1:62] %*% t((preCOVID.svd$v)[1:62,1:62]))
preCOVID.approx[[3]] <- data.frame(preCOVID.svd$u[,1:124] %*% diag(preCOVID.svd$d)[1:124,1:124] %*% t((preCOVID.svd$v)[1:124,1:124]))
preCOVID.approx[[4]] <- data.frame(preCOVID.svd$u[,1:332] %*% diag(preCOVID.svd$d)[1:332,1:332] %*% t((preCOVID.svd$v)[1:332,1:332]))
preCOVID.approx[[5]] <- data.frame(preCOVID.svd$u %*% diag(preCOVID.svd$d) %*% t((preCOVID.svd$v)))
names(preCOVID.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(preCOVID.approx)){
  names(preCOVID.approx[[i]]) <- substring(names(preCOVID.approx[[i]]),2)
  preCOVID.approx[[i]] <- cbind(preCOVID.approx[[i]], label)%>%
    select(label,everything())
  rownames(preCOVID.approx[[i]]) <- NULL
}
preCOVID.results <- as.data.frame(list())
incrementStart<-function(x){4*x-3}
incrementEnd<-function(x){4*x}
for (dat in 1:length(preCOVID.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.approx[[dat]]$label,
        predict(object = ml.train.triage[[c]], preCOVID.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage[[c]], preCOVID.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.results <- rbind(preCOVID.results,
                              tmp%>%
                                mutate(
                                  "Classifier" = names(ml.train.triage[c]),
                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                )%>%
                                select(
                                  c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                    "Kappa","Sensitivity","Specificity",
                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                )
    )
    rm(tmp, cm, predictions)
  }
}
#write.xlsx(preCOVID.results,"preCOVIDresults.xlsx",row.names = F, quote = F)
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.postCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.postCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.postCOVID[setdiff(names(controls.postCOVID), names(cases.postCOVID))] <- 0
controls.postCOVID[setdiff(names(cases.postCOVID), names(controls.postCOVID))] <- 0
postCOVID <- rbind(cases.postCOVID, controls.postCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(postCOVID))
tmp <- data.frame(matrix( , nrow = nrow(postCOVID), 
                          ncol = length(setdiff(colnames(triage), colnames(postCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
postCOVID <- cbind(postCOVID, tmp)
rm(tmp)
postCOVID.svd <- svd(postCOVID)
postCOVID.approx <- list()
postCOVID.approx[[1]] <- data.frame(postCOVID.svd$u[,1:31] %*% diag(postCOVID.svd$d)[1:31,1:31] %*% t((postCOVID.svd$v)[1:31,1:31]))
postCOVID.approx[[2]] <- data.frame(postCOVID.svd$u[,1:62] %*% diag(postCOVID.svd$d)[1:62,1:62] %*% t((postCOVID.svd$v)[1:62,1:62]))
postCOVID.approx[[3]] <- data.frame(postCOVID.svd$u[,1:124] %*% diag(postCOVID.svd$d)[1:124,1:124] %*% t((postCOVID.svd$v)[1:124,1:124]))
postCOVID.approx[[4]] <- data.frame(postCOVID.svd$u[,1:332] %*% diag(postCOVID.svd$d)[1:332,1:332] %*% t((postCOVID.svd$v)[1:332,1:332]))
postCOVID.approx[[5]] <- data.frame(postCOVID.svd$u %*% diag(postCOVID.svd$d) %*% t((postCOVID.svd$v)))
names(postCOVID.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(postCOVID.approx)){
  names(postCOVID.approx[[i]]) <- substring(names(postCOVID.approx[[i]]),2)
  postCOVID.approx[[i]] <- cbind(postCOVID.approx[[i]], label)%>%
    select(label,everything())
  rownames(postCOVID.approx[[i]]) <- NULL
}
postCOVID.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.approx[[dat]]$label,
        predict(object = ml.train.triage[[c]], postCOVID.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage[[c]], postCOVID.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.triage[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.results <- rbind(postCOVID.results,
                               tmp%>%
                                 mutate(
                                   "Classifier" = names(ml.train.triage[c]),
                                   "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                 )%>%
                                 select(
                                   c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                     "Kappa","Sensitivity","Specificity",
                                     "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                 )
    )
    rm(tmp, cm, predictions)
  }
}
#write.xlsx(postCOVID.results,"postCOVIDresults.xlsx",row.names = F, quote = F)
#-----------------------
# TRIAGE NOTES: Binary
#-----------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.preCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.preCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.preCOVID[setdiff(names(controls.preCOVID), names(cases.preCOVID))] <- 0
controls.preCOVID[setdiff(names(cases.preCOVID), names(controls.preCOVID))] <- 0
preCOVID <- rbind(cases.preCOVID, controls.preCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(preCOVID))
tmp <- data.frame(matrix( , nrow = nrow(preCOVID), 
                          ncol = length(setdiff(colnames(triage), colnames(preCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
preCOVID <- cbind(preCOVID, tmp)
rm(tmp)
preCOVID.binary <- preCOVID %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
preCOVID.binary.svd <- svd(preCOVID.binary)
preCOVID.binary.approx <- list()
preCOVID.binary.approx[[1]] <- data.frame(preCOVID.binary.svd$u[,1:31] %*% diag(preCOVID.binary.svd$d)[1:31,1:31] %*% t((preCOVID.binary.svd$v)[1:31,1:31]))
preCOVID.binary.approx[[2]] <- data.frame(preCOVID.binary.svd$u[,1:62] %*% diag(preCOVID.binary.svd$d)[1:62,1:62] %*% t((preCOVID.binary.svd$v)[1:62,1:62]))
preCOVID.binary.approx[[3]] <- data.frame(preCOVID.binary.svd$u[,1:124] %*% diag(preCOVID.binary.svd$d)[1:124,1:124] %*% t((preCOVID.binary.svd$v)[1:124,1:124]))
preCOVID.binary.approx[[4]] <- data.frame(preCOVID.binary.svd$u[,1:332] %*% diag(preCOVID.binary.svd$d)[1:332,1:332] %*% t((preCOVID.binary.svd$v)[1:332,1:332]))
preCOVID.binary.approx[[5]] <- data.frame(preCOVID.binary.svd$u %*% diag(preCOVID.binary.svd$d) %*% t((preCOVID.binary.svd$v)))
names(preCOVID.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(preCOVID.binary.approx)){
  names(preCOVID.binary.approx[[i]]) <- substring(names(preCOVID.binary.approx[[i]]),2)
  preCOVID.binary.approx[[i]] <- cbind(preCOVID.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(preCOVID.binary.approx[[i]]) <- NULL
}
preCOVID.binary.results <- as.data.frame(list())
for (dat in 1:length(preCOVID.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.binary.approx[[dat]]$label,
        predict(object = ml.train.triage.binary[[c]], preCOVID.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage.binary[[c]], preCOVID.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.binary.results <- rbind(preCOVID.binary.results,
                                     tmp%>%
                                       mutate(
                                         "Classifier" = names(ml.train.triage.binary[c]),
                                         "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                       )%>%
                                       select(
                                         c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                           "Kappa","Sensitivity","Specificity",
                                           "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                       )
    )
    rm(tmp, cm, predictions)
  }
}
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.postCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.postCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.postCOVID[setdiff(names(controls.postCOVID), names(cases.postCOVID))] <- 0
controls.postCOVID[setdiff(names(cases.postCOVID), names(controls.postCOVID))] <- 0
postCOVID <- rbind(cases.postCOVID, controls.postCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(postCOVID))
tmp <- data.frame(matrix( , nrow = nrow(postCOVID), 
                          ncol = length(setdiff(colnames(triage), colnames(postCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
postCOVID <- cbind(postCOVID, tmp)
rm(tmp)
postCOVID.binary <- postCOVID %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
postCOVID.binary.svd <- svd(postCOVID.binary)
postCOVID.binary.approx <- list()
postCOVID.binary.approx[[1]] <- data.frame(postCOVID.binary.svd$u[,1:31] %*% diag(postCOVID.binary.svd$d)[1:31,1:31] %*% t((postCOVID.binary.svd$v)[1:31,1:31]))
postCOVID.binary.approx[[2]] <- data.frame(postCOVID.binary.svd$u[,1:62] %*% diag(postCOVID.binary.svd$d)[1:62,1:62] %*% t((postCOVID.binary.svd$v)[1:62,1:62]))
postCOVID.binary.approx[[3]] <- data.frame(postCOVID.binary.svd$u[,1:124] %*% diag(postCOVID.binary.svd$d)[1:124,1:124] %*% t((postCOVID.binary.svd$v)[1:124,1:124]))
postCOVID.binary.approx[[4]] <- data.frame(postCOVID.binary.svd$u[,1:332] %*% diag(postCOVID.binary.svd$d)[1:332,1:332] %*% t((postCOVID.binary.svd$v)[1:332,1:332]))
postCOVID.binary.approx[[5]] <- data.frame(postCOVID.binary.svd$u %*% diag(postCOVID.binary.svd$d) %*% t((postCOVID.binary.svd$v)))
names(postCOVID.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(postCOVID.binary.approx)){
  names(postCOVID.binary.approx[[i]]) <- substring(names(postCOVID.binary.approx[[i]]),2)
  postCOVID.binary.approx[[i]] <- cbind(postCOVID.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(postCOVID.binary.approx[[i]]) <- NULL
}
postCOVID.binary.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.binary.approx[[dat]]$label,
        predict(object = ml.train.triage.binary[[c]], postCOVID.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage.binary[[c]], postCOVID.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.binary[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.triage.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.binary.results <- rbind(postCOVID.binary.results,
                                      tmp%>%
                                        mutate(
                                          "Classifier" = names(ml.train.triage.binary[c]),
                                          "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                        )%>%
                                        select(
                                          c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                            "Kappa","Sensitivity","Specificity",
                                            "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                        )
    )
    rm(tmp, cm, predictions)
  }
}
#-----------------------
# Splits: 60, 70, 80
#-----------------------
incrementStart <- function(x){9*x-8}
incrementEnd <- function(x){9*x}
#----TRIAGE: PRECOVID
preCOVID.splits.triage.results <- as.data.frame(list())
for (dat in 1:length(preCOVID.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.approx[[dat]]$label,
        predict(object = ml.train.triage.w.controls3.simulations[[c]], preCOVID.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage.w.controls3.simulations[[c]], preCOVID.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.splits.triage.results <- rbind(preCOVID.splits.triage.results,
                              tmp%>%
                                mutate(
                                  "Classifier" = names(ml.train.triage.w.controls3.simulations[c]),
                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                )%>%
                                select(
                                  c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                    "Kappa","Sensitivity","Specificity",
                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                )
    )
    rm(tmp, cm, predictions)
  }
}
preCOVID.splits.triage.results$category <- 'preCOVID.splits.triage'
#----TRIAGE: POSTCOVID
postCOVID.splits.triage.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.approx[[dat]]$label,
        predict(object = ml.train.triage.w.controls3.simulations[[c]], postCOVID.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage.w.controls3.simulations[[c]], postCOVID.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.splits.triage.results <- rbind(postCOVID.splits.triage.results,
                                            tmp%>%
                                              mutate(
                                                "Classifier" = names(ml.train.triage.w.controls3.simulations[c]),
                                                "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                              )%>%
                                              select(
                                                c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                  "Kappa","Sensitivity","Specificity",
                                                  "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                              )
    )
    rm(tmp, cm, predictions)
  }
}
postCOVID.splits.triage.results$category <- 'postCOVID.splits.triage'
#----TRIAGE Binary: PRECOVID
preCOVID.splits.triage.binary.results <- as.data.frame(list())
for (dat in 1:length(preCOVID.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.binary.approx[[dat]]$label,
        predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], preCOVID.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], preCOVID.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.binary.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.binary.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage.binary.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.splits.triage.binary.results <- rbind(preCOVID.splits.triage.binary.results,
                                             tmp%>%
                                               mutate(
                                                 "Classifier" = names(ml.train.triage.binary.w.controls3.simulations[c]),
                                                 "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                               )%>%
                                               select(
                                                 c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                   "Kappa","Sensitivity","Specificity",
                                                   "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                               )
    )
    rm(tmp, cm, predictions)
  }
}
preCOVID.splits.triage.binary.results$category <- 'preCOVID.splits.triage.binary'
#----TRIAGE Binary: POSTCOVID
postCOVID.splits.triage.binary.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.binary.approx[[dat]]$label,
        predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], postCOVID.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.triage.binary.w.controls3.simulations[[c]], postCOVID.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage.binary.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage.binary.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq, ml.train.triage.binary.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.splits.triage.binary.results <- rbind(postCOVID.splits.triage.binary.results,
                                                   tmp%>%
                                                     mutate(
                                                       "Classifier" = names(ml.train.triage.binary.w.controls3.simulations[c]),
                                                       "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                     )%>%
                                                     select(
                                                       c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                         "Kappa","Sensitivity","Specificity",
                                                         "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                     )
    )
    rm(tmp, cm, predictions)
  }
}
postCOVID.splits.triage.binary.results$category <- 'postCOVID.splits.triage.binary'
#-----------------------
cases.preCOVID.NEGATED <- read_csv("stroke_aaa_ext_val_cases_nlp_output_pre_covid_final_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
cases.postCOVID.NEGATED <- read_csv("stroke_aaa_ext_val_cases_nlp_output_post_covid_final_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
controls.preCOVID.NEGATED <- read_csv("stroke_aaa_ext_val_control_nlp_output_pre_covid_final_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
controls.postCOVID.NEGATED <- read_csv("stroke_aaa_ext_val_control_nlp_output_post_covid_final_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
#-----------------------
# TRIAGE NOTES NEGATED
#-----------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.preCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.preCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.preCOVID.NEGATED[setdiff(names(controls.preCOVID.NEGATED), names(cases.preCOVID.NEGATED))] <- 0
controls.preCOVID.NEGATED[setdiff(names(cases.preCOVID.NEGATED), names(controls.preCOVID.NEGATED))] <- 0
preCOVID.NEGATED <- rbind(cases.preCOVID.NEGATED, controls.preCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(preCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(preCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(triage), colnames(preCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
preCOVID.NEGATED <- cbind(preCOVID.NEGATED, tmp)
rm(tmp)
preCOVID.NEGATED.svd <- svd(preCOVID.NEGATED)
preCOVID.NEGATED.approx <- list()
preCOVID.NEGATED.approx[[1]] <- data.frame(preCOVID.NEGATED.svd$u[,1:31] %*% diag(preCOVID.NEGATED.svd$d)[1:31,1:31] %*% t((preCOVID.NEGATED.svd$v)[1:31,1:31]))
preCOVID.NEGATED.approx[[2]] <- data.frame(preCOVID.NEGATED.svd$u[,1:62] %*% diag(preCOVID.NEGATED.svd$d)[1:62,1:62] %*% t((preCOVID.NEGATED.svd$v)[1:62,1:62]))
preCOVID.NEGATED.approx[[3]] <- data.frame(preCOVID.NEGATED.svd$u[,1:124] %*% diag(preCOVID.NEGATED.svd$d)[1:124,1:124] %*% t((preCOVID.NEGATED.svd$v)[1:124,1:124]))
preCOVID.NEGATED.approx[[4]] <- data.frame(preCOVID.NEGATED.svd$u[,1:332] %*% diag(preCOVID.NEGATED.svd$d)[1:332,1:332] %*% t((preCOVID.NEGATED.svd$v)[1:332,1:332]))
preCOVID.NEGATED.approx[[5]] <- data.frame(preCOVID.NEGATED.svd$u %*% diag(preCOVID.NEGATED.svd$d) %*% t((preCOVID.NEGATED.svd$v)))
names(preCOVID.NEGATED.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(preCOVID.NEGATED.approx)){
  names(preCOVID.NEGATED.approx[[i]]) <- substring(names(preCOVID.NEGATED.approx[[i]]),2)
  preCOVID.NEGATED.approx[[i]] <- cbind(preCOVID.NEGATED.approx[[i]], label)%>%
    select(label,everything())
  rownames(preCOVID.NEGATED.approx[[i]]) <- NULL
}
preCOVID.NEGATED.results <- as.data.frame(list())
for (dat in 1:length(preCOVID.NEGATED.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.NEGATED.approx[[dat]]$label,
        predict(object = ml.train.negated.triage[[c]], preCOVID.NEGATED.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.triage[[c]], preCOVID.NEGATED.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.triage[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.triage[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.triage[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.NEGATED.results <- rbind(preCOVID.NEGATED.results,
                                      tmp%>%
                                        mutate(
                                          "Classifier" = names(ml.train.negated.triage[c]),
                                          "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                        )%>%
                                        select(
                                          c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                            "Kappa","Sensitivity","Specificity",
                                            "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                        )
    )
    rm(tmp, cm, predictions)
  }
}
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.postCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.postCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.postCOVID.NEGATED[setdiff(names(controls.postCOVID.NEGATED), names(cases.postCOVID.NEGATED))] <- 0
controls.postCOVID.NEGATED[setdiff(names(cases.postCOVID.NEGATED), names(controls.postCOVID.NEGATED))] <- 0
postCOVID.NEGATED <- rbind(cases.postCOVID.NEGATED, controls.postCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(postCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(postCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(triage), colnames(postCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
postCOVID.NEGATED <- cbind(postCOVID.NEGATED, tmp)
rm(tmp)
postCOVID.NEGATED.svd <- svd(postCOVID.NEGATED)
postCOVID.NEGATED.approx <- list()
postCOVID.NEGATED.approx[[1]] <- data.frame(postCOVID.NEGATED.svd$u[,1:31] %*% diag(postCOVID.NEGATED.svd$d)[1:31,1:31] %*% t((postCOVID.NEGATED.svd$v)[1:31,1:31]))
postCOVID.NEGATED.approx[[2]] <- data.frame(postCOVID.NEGATED.svd$u[,1:62] %*% diag(postCOVID.NEGATED.svd$d)[1:62,1:62] %*% t((postCOVID.NEGATED.svd$v)[1:62,1:62]))
postCOVID.NEGATED.approx[[3]] <- data.frame(postCOVID.NEGATED.svd$u[,1:124] %*% diag(postCOVID.NEGATED.svd$d)[1:124,1:124] %*% t((postCOVID.NEGATED.svd$v)[1:124,1:124]))
postCOVID.NEGATED.approx[[4]] <- data.frame(postCOVID.NEGATED.svd$u[,1:332] %*% diag(postCOVID.NEGATED.svd$d)[1:332,1:332] %*% t((postCOVID.NEGATED.svd$v)[1:332,1:332]))
postCOVID.NEGATED.approx[[5]] <- data.frame(postCOVID.NEGATED.svd$u %*% diag(postCOVID.NEGATED.svd$d) %*% t((postCOVID.NEGATED.svd$v)))
names(postCOVID.NEGATED.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(postCOVID.NEGATED.approx)){
  names(postCOVID.NEGATED.approx[[i]]) <- substring(names(postCOVID.NEGATED.approx[[i]]),2)
  postCOVID.NEGATED.approx[[i]] <- cbind(postCOVID.NEGATED.approx[[i]], label)%>%
    select(label,everything())
  rownames(postCOVID.NEGATED.approx[[i]]) <- NULL
}
postCOVID.NEGATED.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.NEGATED.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.NEGATED.approx[[dat]]$label,
        predict(object = ml.train.negated.triage[[c]], postCOVID.NEGATED.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.triage[[c]], postCOVID.NEGATED.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.triage[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.negated.triage[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.triage[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.NEGATED.results <- rbind(postCOVID.NEGATED.results,
                                       tmp%>%
                                         mutate(
                                           "Classifier" = names(ml.train.negated.triage[c]),
                                           "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                         )%>%
                                         select(
                                           c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                             "Kappa","Sensitivity","Specificity",
                                             "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                         )
    )
    rm(tmp, cm, predictions)
  }
}

#------------------------------
# TRIAGE NOTES NEGATED: Binary
#-----------------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.preCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.preCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.preCOVID.NEGATED[setdiff(names(controls.preCOVID.NEGATED), names(cases.preCOVID.NEGATED))] <- 0
controls.preCOVID.NEGATED[setdiff(names(cases.preCOVID.NEGATED), names(controls.preCOVID.NEGATED))] <- 0
preCOVID.NEGATED <- rbind(cases.preCOVID.NEGATED, controls.preCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(preCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(preCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(triage), colnames(preCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
preCOVID.NEGATED <- cbind(preCOVID.NEGATED, tmp)
rm(tmp)
preCOVID.NEGATED.binary <- preCOVID.NEGATED %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
preCOVID.NEGATED.binary.svd <- svd(preCOVID.NEGATED.binary)
preCOVID.NEGATED.binary.approx <- list()
preCOVID.NEGATED.binary.approx[[1]] <- data.frame(preCOVID.NEGATED.binary.svd$u[,1:31] %*% diag(preCOVID.NEGATED.binary.svd$d)[1:31,1:31] %*% t((preCOVID.NEGATED.binary.svd$v)[1:31,1:31]))
preCOVID.NEGATED.binary.approx[[2]] <- data.frame(preCOVID.NEGATED.binary.svd$u[,1:62] %*% diag(preCOVID.NEGATED.binary.svd$d)[1:62,1:62] %*% t((preCOVID.NEGATED.binary.svd$v)[1:62,1:62]))
preCOVID.NEGATED.binary.approx[[3]] <- data.frame(preCOVID.NEGATED.binary.svd$u[,1:124] %*% diag(preCOVID.NEGATED.binary.svd$d)[1:124,1:124] %*% t((preCOVID.NEGATED.binary.svd$v)[1:124,1:124]))
preCOVID.NEGATED.binary.approx[[4]] <- data.frame(preCOVID.NEGATED.binary.svd$u[,1:332] %*% diag(preCOVID.NEGATED.binary.svd$d)[1:332,1:332] %*% t((preCOVID.NEGATED.binary.svd$v)[1:332,1:332]))
preCOVID.NEGATED.binary.approx[[5]] <- data.frame(preCOVID.NEGATED.binary.svd$u %*% diag(preCOVID.NEGATED.binary.svd$d) %*% t((preCOVID.NEGATED.binary.svd$v)))
names(preCOVID.NEGATED.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(preCOVID.NEGATED.binary.approx)){
  names(preCOVID.NEGATED.binary.approx[[i]]) <- substring(names(preCOVID.NEGATED.binary.approx[[i]]),2)
  preCOVID.NEGATED.binary.approx[[i]] <- cbind(preCOVID.NEGATED.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(preCOVID.NEGATED.binary.approx[[i]]) <- NULL
}
preCOVID.NEGATED.binary.results <- as.data.frame(list())
for (dat in 1:length(preCOVID.NEGATED.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.NEGATED.binary.approx[[dat]]$label,
        predict(object = ml.train.negated.triage.binary[[c]], preCOVID.NEGATED.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.triage.binary[[c]], preCOVID.NEGATED.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.triage.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.triage.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.triage.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.NEGATED.binary.results <- rbind(preCOVID.NEGATED.binary.results,
                                             tmp%>%
                                               mutate(
                                                 "Classifier" = names(ml.train.negated.triage.binary[c]),
                                                 "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                               )%>%
                                               select(
                                                 c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                   "Kappa","Sensitivity","Specificity",
                                                   "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                               )
    )
    rm(tmp, cm, predictions)
  }
}
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(cases.postCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(controls.postCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
cases.postCOVID.NEGATED[setdiff(names(controls.postCOVID.NEGATED), names(cases.postCOVID.NEGATED))] <- 0
controls.postCOVID.NEGATED[setdiff(names(cases.postCOVID.NEGATED), names(controls.postCOVID.NEGATED))] <- 0
postCOVID.NEGATED <- rbind(cases.postCOVID.NEGATED, controls.postCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(postCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(postCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(triage), colnames(postCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
postCOVID.NEGATED <- cbind(postCOVID.NEGATED, tmp)
rm(tmp)
postCOVID.NEGATED.binary <- postCOVID.NEGATED%>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
postCOVID.NEGATED.binary.svd <- svd(postCOVID.NEGATED.binary)
postCOVID.NEGATED.binary.approx <- list()
postCOVID.NEGATED.binary.approx[[1]] <- data.frame(postCOVID.NEGATED.binary.svd$u[,1:31] %*% diag(postCOVID.NEGATED.binary.svd$d)[1:31,1:31] %*% t((postCOVID.NEGATED.binary.svd$v)[1:31,1:31]))
postCOVID.NEGATED.binary.approx[[2]] <- data.frame(postCOVID.NEGATED.binary.svd$u[,1:62] %*% diag(postCOVID.NEGATED.binary.svd$d)[1:62,1:62] %*% t((postCOVID.NEGATED.binary.svd$v)[1:62,1:62]))
postCOVID.NEGATED.binary.approx[[3]] <- data.frame(postCOVID.NEGATED.binary.svd$u[,1:124] %*% diag(postCOVID.NEGATED.binary.svd$d)[1:124,1:124] %*% t((postCOVID.NEGATED.binary.svd$v)[1:124,1:124]))
postCOVID.NEGATED.binary.approx[[4]] <- data.frame(postCOVID.NEGATED.binary.svd$u[,1:332] %*% diag(postCOVID.NEGATED.binary.svd$d)[1:332,1:332] %*% t((postCOVID.NEGATED.binary.svd$v)[1:332,1:332]))
postCOVID.NEGATED.binary.approx[[5]] <- data.frame(postCOVID.NEGATED.binary.svd$u %*% diag(postCOVID.NEGATED.binary.svd$d) %*% t((postCOVID.NEGATED.binary.svd$v)))
names(postCOVID.NEGATED.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(postCOVID.NEGATED.binary.approx)){
  names(postCOVID.NEGATED.binary.approx[[i]]) <- substring(names(postCOVID.NEGATED.binary.approx[[i]]),2)
  postCOVID.NEGATED.binary.approx[[i]] <- cbind(postCOVID.NEGATED.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(postCOVID.NEGATED.binary.approx[[i]]) <- NULL
}
postCOVID.NEGATED.binary.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.NEGATED.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.NEGATED.binary.approx[[dat]]$label,
        predict(object = ml.train.negated.triage.binary[[c]], postCOVID.NEGATED.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.triage.binary[[c]], postCOVID.NEGATED.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.triage.binary[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.negated.triage.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.triage.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.NEGATED.binary.results <- rbind(postCOVID.NEGATED.binary.results,
                                              tmp%>%
                                                mutate(
                                                  "Classifier" = names(ml.train.negated.triage.binary[c]),
                                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                )%>%
                                                select(
                                                  c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                    "Kappa","Sensitivity","Specificity",
                                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                )
    )
    rm(tmp, cm, predictions)
  }
}


#-----------------------
# Splits: 60, 70, 80
#-----------------------
incrementStart <- function(x){9*x-8}
incrementEnd <- function(x){9*x}
#----TRIAGE NEGATED: PRECOVID 
preCOVID.neg.splits.triage.results <- as.data.frame(list())
for (dat in 1:length(preCOVID.NEGATED.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.NEGATED.approx[[dat]]$label,
        predict(object = ml.train.neg.triage.w.controls3.simulations[[c]], preCOVID.NEGATED.approx[[dat]], type = "prob"),
        predict(object = ml.train.neg.triage.w.controls3.simulations[[c]], preCOVID.NEGATED.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.neg.triage.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.neg.triage.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq, ml.train.neg.triage.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.neg.splits.triage.results <- rbind(preCOVID.neg.splits.triage.results,
                                            tmp%>%
                                              mutate(
                                                "Classifier" = names(ml.train.neg.triage.w.controls3.simulations[c]),
                                                "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                              )%>%
                                              select(
                                                c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                  "Kappa","Sensitivity","Specificity",
                                                  "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                              )
    )
    rm(tmp, cm, predictions)
  }
}
preCOVID.neg.splits.triage.results$category <- 'preCOVID.neg.splits.triage'
#----TRIAGE NEGATED: POSTCOVID 
postCOVID.neg.splits.triage.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.NEGATED.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.NEGATED.approx[[dat]]$label,
        predict(object = ml.train.neg.triage.w.controls3.simulations[[c]], postCOVID.NEGATED.approx[[dat]], type = "prob"),
        predict(object = ml.train.neg.triage.w.controls3.simulations[[c]], postCOVID.NEGATED.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.neg.triage.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.neg.triage.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq, ml.train.neg.triage.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.neg.splits.triage.results <- rbind(postCOVID.neg.splits.triage.results,
                                                tmp%>%
                                                  mutate(
                                                    "Classifier" = names(ml.train.neg.triage.w.controls3.simulations[c]),
                                                    "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                  )%>%
                                                  select(
                                                    c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                      "Kappa","Sensitivity","Specificity",
                                                      "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                  )
    )
    rm(tmp, cm, predictions)
  }
}
postCOVID.neg.splits.triage.results$category <- 'postCOVID.neg.splits.triage'
#----TRIAGE NEGATED Binary: PRECOVID 
preCOVID.neg.binary.splits.triage.results <- as.data.frame(list())
for (dat in 1:length(preCOVID.NEGATED.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        preCOVID.NEGATED.binary.approx[[dat]]$label,
        predict(object = ml.train.neg.triage.binary.w.controls3.simulations[[c]], preCOVID.NEGATED.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.neg.triage.binary.w.controls3.simulations[[c]], preCOVID.NEGATED.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.neg.triage.binary.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.neg.triage.binary.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq, ml.train.neg.triage.binary.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    preCOVID.neg.binary.splits.triage.results <- rbind(preCOVID.neg.binary.splits.triage.results,
                                                 tmp%>%
                                                   mutate(
                                                     "Classifier" = names(ml.train.neg.triage.binary.w.controls3.simulations[c]),
                                                     "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                   )%>%
                                                   select(
                                                     c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                       "Kappa","Sensitivity","Specificity",
                                                       "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                   )
    )
    rm(tmp, cm, predictions)
  }
}
preCOVID.neg.binary.splits.triage.results$category <- 'preCOVID.neg.binary.splits.triage'
#----TRIAGE NEGATED Binary: POSTCOVID 
postCOVID.neg.binary.splits.triage.results <- as.data.frame(list())
for (dat in 1:length(postCOVID.NEGATED.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        postCOVID.NEGATED.binary.approx[[dat]]$label,
        predict(object = ml.train.neg.triage.binary.w.controls3.simulations[[c]], postCOVID.NEGATED.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.neg.triage.binary.w.controls3.simulations[[c]], postCOVID.NEGATED.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.neg.triage.binary.w.controls3.simulations[c]),                                      # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.neg.triage.binary.w.controls3.simulations[c]),                                      # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq, ml.train.neg.triage.binary.w.controls3.simulations[c]),                     # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    postCOVID.neg.binary.splits.triage.results <- rbind(postCOVID.neg.binary.splits.triage.results,
                                                       tmp%>%
                                                         mutate(
                                                           "Classifier" = names(ml.train.neg.triage.binary.w.controls3.simulations[c]),
                                                           "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                         )%>%
                                                         select(
                                                           c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                             "Kappa","Sensitivity","Specificity",
                                                             "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                         )
    )
    rm(tmp, cm, predictions)
  }
}
postCOVID.neg.binary.splits.triage.results$category <- 'postCOVID.neg.binary.splits.triage'

#-----------------------
# Export results to xlsx
#-----------------------
preCOVID.results$desc <- 'triage'
preCOVID.binary.results$desc <- 'triage.binary'
preCOVID.NEGATED.results$desc <- 'triage.negated' 
preCOVID.NEGATED.binary.results$desc <- 'triage.negated.binary'
write.xlsx(rbind(
  preCOVID.results, preCOVID.binary.results, preCOVID.NEGATED.results, preCOVID.NEGATED.binary.results
),
"preCOVIDresults.xlsx",row.names = F, quote = F)
postCOVID.results$desc <- 'triage'
postCOVID.binary.results$desc <- 'triage.binary'
postCOVID.NEGATED.results$desc <- 'triage.negated' 
postCOVID.NEGATED.binary.results$desc <- 'triage.negated.binary'
write.xlsx(rbind(
  postCOVID.results, postCOVID.binary.results, postCOVID.NEGATED.results, postCOVID.NEGATED.binary.results
),
"postCOVIDresults.xlsx",row.names = F, quote = F)

# write splits
write.xlsx(
  rbind(
    preCOVID.splits.triage.results,
    postCOVID.splits.triage.results,
    preCOVID.splits.triage.binary.results,
    postCOVID.splits.triage.binary.results,
    preCOVID.neg.splits.triage.results,
    postCOVID.neg.splits.triage.results,
    preCOVID.neg.binary.splits.triage.results,
    postCOVID.neg.binary.splits.triage.results
  ),
  "results_covid_splits.xlsx", row.names = F, quote = F
)


# summarize findings
summary.results <- readxl::read_excel('v3Results_1.1.xlsx')
cols.factors <- c("COVID","desc","What Controls?","No.of cuis","Model","Splits" )
summary.results[cols.factors] <- lapply(summary.results[cols.factors], factor)
summary.results.long <- reshape2::melt(summary.results %>% 
                                         select("COVID","desc","What Controls?","No.of cuis","Model","Splits", "AUROC","Accuracy","Sensitivity","Specificity")
                                       #id.vars = 
                                       )
tiff("post_covid_results_summary.tiff", units = "in", width = 15, height = 8, res = 300)
ggplot(data = subset(summary.results.long, COVID == 'post'), aes(x = value, y = desc))+
  geom_point(aes(shape = Model, color = Model), size = 2, alpha = .5)+
  facet_wrap(variable~`No.of cuis`)+
  scale_colour_manual(values = c('green'))+
  theme(
    axis.title = element_blank()
  )+
  scale_color_manual(values = c('dark green', 'steelblue', 'red', 'blue'))
dev.off()

################################################################################
##########################  validation: 04/19/21  ##############################
################################################################################
library(readr)
library(dplyr)
library(caret)
library(pROC)
library(openxlsx)
options(readr.default_locale=readr::locale(tz="US/Eastern"))
provider.cases.preCOVID <- read_csv("stroke_cases_pre_covid_provider_notes_202104161838_FORMATTED.csv")%>%
  select(-c('_c0'))
provider.cases.postCOVID <- read_csv("stroke_cases_post_covid_provider_notes_202104161837_FORMATTED.csv")%>%
  select(-c('_c0'))
provider.controls.preCOVID <- read_csv("stroke_control_pre_covid_provider_notes_202104161838_FORMATTED.csv")%>%
  select(-c('_c0'))
provider.controls.postCOVID <- read_csv("stroke_control_post_covid_provider_notes_202104161838_FORMATTED.csv")%>%
  select(-c('_c0'))
#-----------------------
# PROVIDER NOTES
#-----------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.preCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.preCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.preCOVID[setdiff(names(provider.controls.preCOVID), names(provider.cases.preCOVID))] <- 0
provider.controls.preCOVID[setdiff(names(provider.cases.preCOVID), names(provider.controls.preCOVID))] <- 0
provider.preCOVID <- rbind(provider.cases.preCOVID, provider.controls.preCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.preCOVID))
tmp <- data.frame(matrix( , nrow = nrow(provider.preCOVID), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.preCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.preCOVID <- cbind(provider.preCOVID, tmp)
rm(tmp)
provider.preCOVID.svd <- svd(provider.preCOVID)
provider.preCOVID.approx <- list()
provider.preCOVID.approx[[1]] <- data.frame(provider.preCOVID.svd$u[,1:31] %*% diag(provider.preCOVID.svd$d)[1:31,1:31] %*% t((provider.preCOVID.svd$v)[1:31,1:31]))
provider.preCOVID.approx[[2]] <- data.frame(provider.preCOVID.svd$u[,1:62] %*% diag(provider.preCOVID.svd$d)[1:62,1:62] %*% t((provider.preCOVID.svd$v)[1:62,1:62]))
provider.preCOVID.approx[[3]] <- data.frame(provider.preCOVID.svd$u[,1:124] %*% diag(provider.preCOVID.svd$d)[1:124,1:124] %*% t((provider.preCOVID.svd$v)[1:124,1:124]))
provider.preCOVID.approx[[4]] <- data.frame(provider.preCOVID.svd$u[,1:332] %*% diag(provider.preCOVID.svd$d)[1:332,1:332] %*% t((provider.preCOVID.svd$v)[1:332,1:332]))
provider.preCOVID.approx[[5]] <- data.frame(provider.preCOVID.svd$u %*% diag(provider.preCOVID.svd$d) %*% t((provider.preCOVID.svd$v)))
names(provider.preCOVID.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.preCOVID.approx)){
  names(provider.preCOVID.approx[[i]]) <- substring(names(provider.preCOVID.approx[[i]]),2)
  provider.preCOVID.approx[[i]] <- cbind(provider.preCOVID.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.preCOVID.approx[[i]]) <- NULL
}
incrementStart <- function(x){4*x-3}
incrementEnd <- function(x){4*x}
provider.preCOVID.results <- as.data.frame(list())
for (dat in 1:length(provider.preCOVID.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.preCOVID.approx[[dat]]$label,
        predict(object = ml.train.provider[[c]], provider.preCOVID.approx[[dat]], type = "prob"),
        predict(object = ml.train.provider[[c]], provider.preCOVID.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.provider[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.provider[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.provider[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.preCOVID.results <- rbind(provider.preCOVID.results,
                              tmp%>%
                                mutate(
                                  "Classifier" = names(ml.train.provider[c]),
                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                )%>%
                                select(
                                  c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                    "Kappa","Sensitivity","Specificity",
                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                )
    )
    rm(tmp, cm, predictions)
  }
}
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.postCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.postCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.postCOVID[setdiff(names(provider.controls.postCOVID), names(provider.cases.postCOVID))] <- 0
provider.controls.postCOVID[setdiff(names(provider.cases.postCOVID), names(provider.controls.postCOVID))] <- 0
provider.postCOVID <- rbind(provider.cases.postCOVID, provider.controls.postCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.postCOVID))
tmp <- data.frame(matrix( , nrow = nrow(provider.postCOVID), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.postCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.postCOVID <- cbind(provider.postCOVID, tmp)
rm(tmp)
provider.postCOVID.svd <- svd(provider.postCOVID)
provider.postCOVID.approx <- list()
provider.postCOVID.approx[[1]] <- data.frame(provider.postCOVID.svd$u[,1:31] %*% diag(provider.postCOVID.svd$d)[1:31,1:31] %*% t((provider.postCOVID.svd$v)[1:31,1:31]))
provider.postCOVID.approx[[2]] <- data.frame(provider.postCOVID.svd$u[,1:62] %*% diag(provider.postCOVID.svd$d)[1:62,1:62] %*% t((provider.postCOVID.svd$v)[1:62,1:62]))
provider.postCOVID.approx[[3]] <- data.frame(provider.postCOVID.svd$u[,1:124] %*% diag(provider.postCOVID.svd$d)[1:124,1:124] %*% t((provider.postCOVID.svd$v)[1:124,1:124]))
provider.postCOVID.approx[[4]] <- data.frame(provider.postCOVID.svd$u[,1:332] %*% diag(provider.postCOVID.svd$d)[1:332,1:332] %*% t((provider.postCOVID.svd$v)[1:332,1:332]))
provider.postCOVID.approx[[5]] <- data.frame(provider.postCOVID.svd$u %*% diag(provider.postCOVID.svd$d) %*% t((provider.postCOVID.svd$v)))
names(provider.postCOVID.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.postCOVID.approx)){
  names(provider.postCOVID.approx[[i]]) <- substring(names(provider.postCOVID.approx[[i]]),2)
  provider.postCOVID.approx[[i]] <- cbind(provider.postCOVID.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.postCOVID.approx[[i]]) <- NULL
}
provider.postCOVID.results <- as.data.frame(list())
for (dat in 1:length(provider.postCOVID.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.postCOVID.approx[[dat]]$label,
        predict(object = ml.train.provider[[c]], provider.postCOVID.approx[[dat]], type = "prob"),
        predict(object = ml.train.provider[[c]], provider.postCOVID.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.provider[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.provider[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.provider[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.postCOVID.results <- rbind(provider.postCOVID.results,
                               tmp%>%
                                 mutate(
                                   "Classifier" = names(ml.train.provider[c]),
                                   "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                 )%>%
                                 select(
                                   c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                     "Kappa","Sensitivity","Specificity",
                                     "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                 )
    )
    rm(tmp, cm, predictions)
  }
}

#-----------------------
# PROVIDER NOTES: Binary
#-----------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.preCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.preCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.preCOVID[setdiff(names(provider.controls.preCOVID), names(provider.cases.preCOVID))] <- 0
provider.controls.preCOVID[setdiff(names(provider.cases.preCOVID), names(provider.controls.preCOVID))] <- 0
provider.preCOVID <- rbind(provider.cases.preCOVID, provider.controls.preCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.preCOVID))
tmp <- data.frame(matrix( , nrow = nrow(provider.preCOVID), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.preCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.preCOVID <- cbind(provider.preCOVID, tmp)
rm(tmp)
provider.preCOVID.binary <- provider.preCOVID %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
provider.preCOVID.binary.svd <- svd(provider.preCOVID.binary)
provider.preCOVID.binary.approx <- list()
provider.preCOVID.binary.approx[[1]] <- data.frame(provider.preCOVID.binary.svd$u[,1:31] %*% diag(provider.preCOVID.binary.svd$d)[1:31,1:31] %*% t((provider.preCOVID.binary.svd$v)[1:31,1:31]))
provider.preCOVID.binary.approx[[2]] <- data.frame(provider.preCOVID.binary.svd$u[,1:62] %*% diag(provider.preCOVID.binary.svd$d)[1:62,1:62] %*% t((provider.preCOVID.binary.svd$v)[1:62,1:62]))
provider.preCOVID.binary.approx[[3]] <- data.frame(provider.preCOVID.binary.svd$u[,1:124] %*% diag(provider.preCOVID.binary.svd$d)[1:124,1:124] %*% t((provider.preCOVID.binary.svd$v)[1:124,1:124]))
provider.preCOVID.binary.approx[[4]] <- data.frame(provider.preCOVID.binary.svd$u[,1:332] %*% diag(provider.preCOVID.binary.svd$d)[1:332,1:332] %*% t((provider.preCOVID.binary.svd$v)[1:332,1:332]))
provider.preCOVID.binary.approx[[5]] <- data.frame(provider.preCOVID.binary.svd$u %*% diag(provider.preCOVID.binary.svd$d) %*% t((provider.preCOVID.binary.svd$v)))
names(provider.preCOVID.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.preCOVID.binary.approx)){
  names(provider.preCOVID.binary.approx[[i]]) <- substring(names(provider.preCOVID.binary.approx[[i]]),2)
  provider.preCOVID.binary.approx[[i]] <- cbind(provider.preCOVID.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.preCOVID.binary.approx[[i]]) <- NULL
}
provider.preCOVID.binary.results <- as.data.frame(list())
for (dat in 1:length(provider.preCOVID.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.preCOVID.binary.approx[[dat]]$label,
        predict(object = ml.train.provider.binary[[c]], provider.preCOVID.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.provider.binary[[c]], provider.preCOVID.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.provider.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.provider.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.provider.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.preCOVID.binary.results <- rbind(provider.preCOVID.binary.results,
                                     tmp%>%
                                       mutate(
                                         "Classifier" = names(ml.train.provider.binary[c]),
                                         "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                       )%>%
                                       select(
                                         c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                           "Kappa","Sensitivity","Specificity",
                                           "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                       )
    )
    rm(tmp, cm, predictions)
  }
}
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.postCOVID), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.postCOVID), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.postCOVID[setdiff(names(provider.controls.postCOVID), names(provider.cases.postCOVID))] <- 0
provider.controls.postCOVID[setdiff(names(provider.cases.postCOVID), names(provider.controls.postCOVID))] <- 0
provider.postCOVID <- rbind(provider.cases.postCOVID, provider.controls.postCOVID)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.postCOVID))
tmp <- data.frame(matrix( , nrow = nrow(provider.postCOVID), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.postCOVID)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.postCOVID <- cbind(provider.postCOVID, tmp)
rm(tmp)
provider.postCOVID.binary <- provider.postCOVID %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
provider.postCOVID.binary.svd <- svd(provider.postCOVID.binary)
provider.postCOVID.binary.approx <- list()
provider.postCOVID.binary.approx[[1]] <- data.frame(provider.postCOVID.binary.svd$u[,1:31] %*% diag(provider.postCOVID.binary.svd$d)[1:31,1:31] %*% t((provider.postCOVID.binary.svd$v)[1:31,1:31]))
provider.postCOVID.binary.approx[[2]] <- data.frame(provider.postCOVID.binary.svd$u[,1:62] %*% diag(provider.postCOVID.binary.svd$d)[1:62,1:62] %*% t((provider.postCOVID.binary.svd$v)[1:62,1:62]))
provider.postCOVID.binary.approx[[3]] <- data.frame(provider.postCOVID.binary.svd$u[,1:124] %*% diag(provider.postCOVID.binary.svd$d)[1:124,1:124] %*% t((provider.postCOVID.binary.svd$v)[1:124,1:124]))
provider.postCOVID.binary.approx[[4]] <- data.frame(provider.postCOVID.binary.svd$u[,1:332] %*% diag(provider.postCOVID.binary.svd$d)[1:332,1:332] %*% t((provider.postCOVID.binary.svd$v)[1:332,1:332]))
provider.postCOVID.binary.approx[[5]] <- data.frame(provider.postCOVID.binary.svd$u %*% diag(provider.postCOVID.binary.svd$d) %*% t((provider.postCOVID.binary.svd$v)))
names(provider.postCOVID.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.postCOVID.binary.approx)){
  names(provider.postCOVID.binary.approx[[i]]) <- substring(names(provider.postCOVID.binary.approx[[i]]),2)
  provider.postCOVID.binary.approx[[i]] <- cbind(provider.postCOVID.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.postCOVID.binary.approx[[i]]) <- NULL
}
provider.postCOVID.binary.results <- as.data.frame(list())
for (dat in 1:length(provider.postCOVID.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.postCOVID.binary.approx[[dat]]$label,
        predict(object = ml.train.provider.binary[[c]], provider.postCOVID.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.provider.binary[[c]], provider.postCOVID.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.provider.binary[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.provider.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.provider.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.postCOVID.binary.results <- rbind(provider.postCOVID.binary.results,
                                      tmp%>%
                                        mutate(
                                          "Classifier" = names(ml.train.provider.binary[c]),
                                          "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                        )%>%
                                        select(
                                          c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                            "Kappa","Sensitivity","Specificity",
                                            "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                        )
    )
    rm(tmp, cm, predictions)
  }
}

#-----------------------
provider.cases.preCOVID.NEGATED <- read_csv("stroke_cases_pre_covid_provider_notes_202104161838_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
provider.cases.postCOVID.NEGATED <- read_csv("stroke_cases_post_covid_provider_notes_202104161837_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
provider.controls.preCOVID.NEGATED <- read_csv("stroke_control_pre_covid_provider_notes_202104161838_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
provider.controls.postCOVID.NEGATED <- read_csv("stroke_control_post_covid_provider_notes_202104161838_FORMATTED_NEGATED.csv")%>%
  select(-c('_c0'))
#-----------------------
# PROVIDER NOTES NEGATED
#-----------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.preCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.preCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.preCOVID.NEGATED[setdiff(names(provider.controls.preCOVID.NEGATED), names(provider.cases.preCOVID.NEGATED))] <- 0
provider.controls.preCOVID.NEGATED[setdiff(names(provider.cases.preCOVID.NEGATED), names(provider.controls.preCOVID.NEGATED))] <- 0
provider.preCOVID.NEGATED <- rbind(provider.cases.preCOVID.NEGATED, provider.controls.preCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.preCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(provider.preCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.preCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.preCOVID.NEGATED <- cbind(provider.preCOVID.NEGATED, tmp)
rm(tmp)
provider.preCOVID.NEGATED.svd <- svd(provider.preCOVID.NEGATED)
provider.preCOVID.NEGATED.approx <- list()
provider.preCOVID.NEGATED.approx[[1]] <- data.frame(provider.preCOVID.NEGATED.svd$u[,1:31] %*% diag(provider.preCOVID.NEGATED.svd$d)[1:31,1:31] %*% t((provider.preCOVID.NEGATED.svd$v)[1:31,1:31]))
provider.preCOVID.NEGATED.approx[[2]] <- data.frame(provider.preCOVID.NEGATED.svd$u[,1:62] %*% diag(provider.preCOVID.NEGATED.svd$d)[1:62,1:62] %*% t((provider.preCOVID.NEGATED.svd$v)[1:62,1:62]))
provider.preCOVID.NEGATED.approx[[3]] <- data.frame(provider.preCOVID.NEGATED.svd$u[,1:124] %*% diag(provider.preCOVID.NEGATED.svd$d)[1:124,1:124] %*% t((provider.preCOVID.NEGATED.svd$v)[1:124,1:124]))
provider.preCOVID.NEGATED.approx[[4]] <- data.frame(provider.preCOVID.NEGATED.svd$u[,1:332] %*% diag(provider.preCOVID.NEGATED.svd$d)[1:332,1:332] %*% t((provider.preCOVID.NEGATED.svd$v)[1:332,1:332]))
provider.preCOVID.NEGATED.approx[[5]] <- data.frame(provider.preCOVID.NEGATED.svd$u %*% diag(provider.preCOVID.NEGATED.svd$d) %*% t((provider.preCOVID.NEGATED.svd$v)))
names(provider.preCOVID.NEGATED.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.preCOVID.NEGATED.approx)){
  names(provider.preCOVID.NEGATED.approx[[i]]) <- substring(names(provider.preCOVID.NEGATED.approx[[i]]),2)
  provider.preCOVID.NEGATED.approx[[i]] <- cbind(provider.preCOVID.NEGATED.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.preCOVID.NEGATED.approx[[i]]) <- NULL
}
provider.preCOVID.NEGATED.results <- as.data.frame(list())
for (dat in 1:length(provider.preCOVID.NEGATED.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.preCOVID.NEGATED.approx[[dat]]$label,
        predict(object = ml.train.negated.provider[[c]], provider.preCOVID.NEGATED.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.provider[[c]], provider.preCOVID.NEGATED.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.provider[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.provider[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.provider[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.preCOVID.NEGATED.results <- rbind(provider.preCOVID.NEGATED.results,
                                      tmp%>%
                                        mutate(
                                          "Classifier" = names(ml.train.negated.provider[c]),
                                          "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                        )%>%
                                        select(
                                          c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                            "Kappa","Sensitivity","Specificity",
                                            "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                        )
    )
    rm(tmp, cm, predictions)
  }
}
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.postCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.postCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.postCOVID.NEGATED[setdiff(names(provider.controls.postCOVID.NEGATED), names(provider.cases.postCOVID.NEGATED))] <- 0
provider.controls.postCOVID.NEGATED[setdiff(names(provider.cases.postCOVID.NEGATED), names(provider.controls.postCOVID.NEGATED))] <- 0
provider.postCOVID.NEGATED <- rbind(provider.cases.postCOVID.NEGATED, provider.controls.postCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.postCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(provider.postCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.postCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.postCOVID.NEGATED <- cbind(provider.postCOVID.NEGATED, tmp)
rm(tmp)
provider.postCOVID.NEGATED.svd <- svd(provider.postCOVID.NEGATED)
provider.postCOVID.NEGATED.approx <- list()
provider.postCOVID.NEGATED.approx[[1]] <- data.frame(provider.postCOVID.NEGATED.svd$u[,1:31] %*% diag(provider.postCOVID.NEGATED.svd$d)[1:31,1:31] %*% t((provider.postCOVID.NEGATED.svd$v)[1:31,1:31]))
provider.postCOVID.NEGATED.approx[[2]] <- data.frame(provider.postCOVID.NEGATED.svd$u[,1:62] %*% diag(provider.postCOVID.NEGATED.svd$d)[1:62,1:62] %*% t((provider.postCOVID.NEGATED.svd$v)[1:62,1:62]))
provider.postCOVID.NEGATED.approx[[3]] <- data.frame(provider.postCOVID.NEGATED.svd$u[,1:124] %*% diag(provider.postCOVID.NEGATED.svd$d)[1:124,1:124] %*% t((provider.postCOVID.NEGATED.svd$v)[1:124,1:124]))
provider.postCOVID.NEGATED.approx[[4]] <- data.frame(provider.postCOVID.NEGATED.svd$u[,1:332] %*% diag(provider.postCOVID.NEGATED.svd$d)[1:332,1:332] %*% t((provider.postCOVID.NEGATED.svd$v)[1:332,1:332]))
provider.postCOVID.NEGATED.approx[[5]] <- data.frame(provider.postCOVID.NEGATED.svd$u %*% diag(provider.postCOVID.NEGATED.svd$d) %*% t((provider.postCOVID.NEGATED.svd$v)))
names(provider.postCOVID.NEGATED.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.postCOVID.NEGATED.approx)){
  names(provider.postCOVID.NEGATED.approx[[i]]) <- substring(names(provider.postCOVID.NEGATED.approx[[i]]),2)
  provider.postCOVID.NEGATED.approx[[i]] <- cbind(provider.postCOVID.NEGATED.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.postCOVID.NEGATED.approx[[i]]) <- NULL
}
provider.postCOVID.NEGATED.results <- as.data.frame(list())
for (dat in 1:length(provider.postCOVID.NEGATED.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.postCOVID.NEGATED.approx[[dat]]$label,
        predict(object = ml.train.negated.provider[[c]], provider.postCOVID.NEGATED.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.provider[[c]], provider.postCOVID.NEGATED.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.provider[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.negated.provider[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.provider[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.postCOVID.NEGATED.results <- rbind(provider.postCOVID.NEGATED.results,
                                       tmp%>%
                                         mutate(
                                           "Classifier" = names(ml.train.negated.provider[c]),
                                           "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                         )%>%
                                         select(
                                           c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                             "Kappa","Sensitivity","Specificity",
                                             "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                         )
    )
    rm(tmp, cm, predictions)
  }
}

#------------------------------
# PROVIDER NOTES NEGATED: Binary
#-----------------------------
# pre COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.preCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.preCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.preCOVID.NEGATED[setdiff(names(provider.controls.preCOVID.NEGATED), names(provider.cases.preCOVID.NEGATED))] <- 0
provider.controls.preCOVID.NEGATED[setdiff(names(provider.cases.preCOVID.NEGATED), names(provider.controls.preCOVID.NEGATED))] <- 0
provider.preCOVID.NEGATED <- rbind(provider.cases.preCOVID.NEGATED, provider.controls.preCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.preCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(provider.preCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.preCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.preCOVID.NEGATED <- cbind(provider.preCOVID.NEGATED, tmp)
rm(tmp)
provider.preCOVID.NEGATED.binary <- provider.preCOVID.NEGATED %>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
provider.preCOVID.NEGATED.binary.svd <- svd(provider.preCOVID.NEGATED.binary)
provider.preCOVID.NEGATED.binary.approx <- list()
provider.preCOVID.NEGATED.binary.approx[[1]] <- data.frame(provider.preCOVID.NEGATED.binary.svd$u[,1:31] %*% diag(provider.preCOVID.NEGATED.binary.svd$d)[1:31,1:31] %*% t((provider.preCOVID.NEGATED.binary.svd$v)[1:31,1:31]))
provider.preCOVID.NEGATED.binary.approx[[2]] <- data.frame(provider.preCOVID.NEGATED.binary.svd$u[,1:62] %*% diag(provider.preCOVID.NEGATED.binary.svd$d)[1:62,1:62] %*% t((provider.preCOVID.NEGATED.binary.svd$v)[1:62,1:62]))
provider.preCOVID.NEGATED.binary.approx[[3]] <- data.frame(provider.preCOVID.NEGATED.binary.svd$u[,1:124] %*% diag(provider.preCOVID.NEGATED.binary.svd$d)[1:124,1:124] %*% t((provider.preCOVID.NEGATED.binary.svd$v)[1:124,1:124]))
provider.preCOVID.NEGATED.binary.approx[[4]] <- data.frame(provider.preCOVID.NEGATED.binary.svd$u[,1:332] %*% diag(provider.preCOVID.NEGATED.binary.svd$d)[1:332,1:332] %*% t((provider.preCOVID.NEGATED.binary.svd$v)[1:332,1:332]))
provider.preCOVID.NEGATED.binary.approx[[5]] <- data.frame(provider.preCOVID.NEGATED.binary.svd$u %*% diag(provider.preCOVID.NEGATED.binary.svd$d) %*% t((provider.preCOVID.NEGATED.binary.svd$v)))
names(provider.preCOVID.NEGATED.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.preCOVID.NEGATED.binary.approx)){
  names(provider.preCOVID.NEGATED.binary.approx[[i]]) <- substring(names(provider.preCOVID.NEGATED.binary.approx[[i]]),2)
  provider.preCOVID.NEGATED.binary.approx[[i]] <- cbind(provider.preCOVID.NEGATED.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.preCOVID.NEGATED.binary.approx[[i]]) <- NULL
}
provider.preCOVID.NEGATED.binary.results <- as.data.frame(list())
for (dat in 1:length(provider.preCOVID.NEGATED.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.preCOVID.NEGATED.binary.approx[[dat]]$label,
        predict(object = ml.train.negated.provider.binary[[c]], provider.preCOVID.NEGATED.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.provider.binary[[c]], provider.preCOVID.NEGATED.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.provider.binary[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.negated.provider.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.provider.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.preCOVID.NEGATED.binary.results <- rbind(provider.preCOVID.NEGATED.binary.results,
                                             tmp%>%
                                               mutate(
                                                 "Classifier" = names(ml.train.negated.provider.binary[c]),
                                                 "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                               )%>%
                                               select(
                                                 c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                   "Kappa","Sensitivity","Specificity",
                                                   "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                               )
    )
    rm(tmp, cm, predictions)
  }
}
# post COVID
label <- data.frame(rbind(
  t(data.frame(lapply(1:nrow(provider.cases.postCOVID.NEGATED), function(x) 'X1'))),
  t(data.frame(lapply(1:nrow(provider.controls.postCOVID.NEGATED), function(x) 'X0')))
))
colnames(label) <- 'label'
label$label <- as.factor(label$label)
provider.cases.postCOVID.NEGATED[setdiff(names(provider.controls.postCOVID.NEGATED), names(provider.cases.postCOVID.NEGATED))] <- 0
provider.controls.postCOVID.NEGATED[setdiff(names(provider.cases.postCOVID.NEGATED), names(provider.controls.postCOVID.NEGATED))] <- 0
provider.postCOVID.NEGATED <- rbind(provider.cases.postCOVID.NEGATED, provider.controls.postCOVID.NEGATED)
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(provider.postCOVID.NEGATED))
tmp <- data.frame(matrix( , nrow = nrow(provider.postCOVID.NEGATED), 
                          ncol = length(setdiff(colnames(provider), colnames(provider.postCOVID.NEGATED)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
provider.postCOVID.NEGATED <- cbind(provider.postCOVID.NEGATED, tmp)
rm(tmp)
provider.postCOVID.NEGATED.binary <- provider.postCOVID.NEGATED%>% 
  mutate_all(funs(
    ifelse( .>= 1, 1, 
            ifelse(. <= -1, -1, 0)))
  )
provider.postCOVID.NEGATED.binary.svd <- svd(provider.postCOVID.NEGATED.binary)
provider.postCOVID.NEGATED.binary.approx <- list()
provider.postCOVID.NEGATED.binary.approx[[1]] <- data.frame(provider.postCOVID.NEGATED.binary.svd$u[,1:31] %*% diag(provider.postCOVID.NEGATED.binary.svd$d)[1:31,1:31] %*% t((provider.postCOVID.NEGATED.binary.svd$v)[1:31,1:31]))
provider.postCOVID.NEGATED.binary.approx[[2]] <- data.frame(provider.postCOVID.NEGATED.binary.svd$u[,1:62] %*% diag(provider.postCOVID.NEGATED.binary.svd$d)[1:62,1:62] %*% t((provider.postCOVID.NEGATED.binary.svd$v)[1:62,1:62]))
provider.postCOVID.NEGATED.binary.approx[[3]] <- data.frame(provider.postCOVID.NEGATED.binary.svd$u[,1:124] %*% diag(provider.postCOVID.NEGATED.binary.svd$d)[1:124,1:124] %*% t((provider.postCOVID.NEGATED.binary.svd$v)[1:124,1:124]))
provider.postCOVID.NEGATED.binary.approx[[4]] <- data.frame(provider.postCOVID.NEGATED.binary.svd$u[,1:332] %*% diag(provider.postCOVID.NEGATED.binary.svd$d)[1:332,1:332] %*% t((provider.postCOVID.NEGATED.binary.svd$v)[1:332,1:332]))
provider.postCOVID.NEGATED.binary.approx[[5]] <- data.frame(provider.postCOVID.NEGATED.binary.svd$u %*% diag(provider.postCOVID.NEGATED.binary.svd$d) %*% t((provider.postCOVID.NEGATED.binary.svd$v)))
names(provider.postCOVID.NEGATED.binary.approx) <- c('approx.31', 'approx.62', 'approx.124', 'approx.332', 'approx.665')
for (i in 1:length(provider.postCOVID.NEGATED.binary.approx)){
  names(provider.postCOVID.NEGATED.binary.approx[[i]]) <- substring(names(provider.postCOVID.NEGATED.binary.approx[[i]]),2)
  provider.postCOVID.NEGATED.binary.approx[[i]] <- cbind(provider.postCOVID.NEGATED.binary.approx[[i]], label)%>%
    select(label,everything())
  rownames(provider.postCOVID.NEGATED.binary.approx[[i]]) <- NULL
}
provider.postCOVID.NEGATED.binary.results <- as.data.frame(list())
for (dat in 1:length(provider.postCOVID.NEGATED.binary.approx)){
  for (c in incrementStart(dat):incrementEnd(dat)) {
    predictions <- setNames(
      data.frame(
        provider.postCOVID.NEGATED.binary.approx[[dat]]$label,
        predict(object = ml.train.negated.provider.binary[[c]], provider.postCOVID.NEGATED.binary.approx[[dat]], type = "prob"),
        predict(object = ml.train.negated.provider.binary[[c]], provider.postCOVID.NEGATED.binary.approx[[dat]], type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.negated.provider.binary[c]),                                                              # Fetch Recall,Specificity,postcision
      fetchResults(cm$overall, ml.train.negated.provider.binary[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.negated.provider.binary[c]),                                             # Fetch TP,FP,FN,TN
      roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
      prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
    )))
    provider.postCOVID.NEGATED.binary.results <- rbind(provider.postCOVID.NEGATED.binary.results,
                                              tmp%>%
                                                mutate(
                                                  "Classifier" = names(ml.train.negated.provider.binary[c]),
                                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                                )%>%
                                                select(
                                                  c("Classifier",AUROC = "23",AUPR = "24","Accuracy","95%CI",NIR = "AccuracyNull",
                                                    "Kappa","Sensitivity","Specificity",
                                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                                )
    )
    rm(tmp, cm, predictions)
  }
}


#-----------------------
# Export results to xlsx
#-----------------------
provider.preCOVID.results$desc <- 'provider.preCOVID'
provider.postCOVID.results$desc <- 'provider.postCOVID'
provider.preCOVID.binary.results$desc <- 'provider.preCOVID.binary'
provider.postCOVID.binary.results$desc <- 'provider.postCOVID.binary'
provider.preCOVID.NEGATED.results$desc <- 'provider.preCOVID.negated'
provider.postCOVID.NEGATED.results$desc <- 'provider.postCOVID.negated'
provider.preCOVID.NEGATED.binary.results$desc <- 'provider.preCOVID.negated.binary'
provider.postCOVID.NEGATED.binary.results$desc <- 'provider.postCOVID.negated.binary'
write.xlsx(
  rbind(
    provider.preCOVID.results,
    provider.postCOVID.results,
    provider.preCOVID.binary.results,
    provider.postCOVID.binary.results,
    provider.preCOVID.NEGATED.results,
    provider.postCOVID.NEGATED.results,
    provider.preCOVID.NEGATED.binary.results,
    provider.postCOVID.NEGATED.binary.results
  ),
  "results_covid_provider.xlsx", row.names = F, quote = F
)
#-----------------
# cui-to-text
#-----------------
`cui-to-text` <- read.csv("~/Genentech/v3/_select_cui_text_count_from_select_ROW_NUMBER_OVER_PARTITION_BY__202104281643.csv")
tmp <- left_join(
  read_excel("cui-to-english.xlsx", sheet = 'triage_cases_pre_covid')%>%
    select('cui'),
  `cui-to-text`,
  by = 'cui'
)
write.xlsx(tmp
           , 'tmp.xlsx')

#-----------
# Plotting
#----------
#----
# summary plots
#----
summary.results <- readxl::read_excel('v3Results_04202021.xlsx')
cols.factors <- c("Category","COVID","desc","What Controls?","No.of cuis","Model","Splits" )
summary.results[cols.factors] <- lapply(summary.results[cols.factors], factor)
summary.results.long <- reshape2::melt((summary.results)%>% 
                                         select("Category","COVID","desc","What Controls?","No.of cuis","Model","Splits",
                                                "AUROC","Accuracy","Sensitivity","Specificity")
)
tiff("provider_pre_covid_results_summary.tiff", units = "in", width = 15, height = 8, res = 300)
ggplot(data = subset(summary.results.long, COVID == 'pre' &  Category == 'provider'), aes(x = value, y = desc))+
  #geom_point(aes(shape = Model, color = Model), size = 2, alpha = .5)+
  geom_boxplot() +
  facet_wrap(variable~`No.of cuis`)+
  scale_colour_manual(values = c('green'))+
  theme(
    axis.title = element_blank()
  )+
  scale_color_manual(values = c('dark green', 'steelblue', 'red', 'blue'))
dev.off()
#----
# feature importance 
#----
# Model development
tmp.ml.train.triage <- list()
for(i in 5) {
  for(c in 1:length(classifiers)) {
    print(paste("started:",names(triage.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
    tmp.ml.train.triage[[paste0(names(triage.trainSet[i]),".",classifiers[[c]])]] <- train(label~.,
                                                                                       data = triage.trainSet[[i]],
                                                                                       method = classifiers[[c]],
                                                                                       preProcess = c("center","scale"),
                                                                                       metric = "ROC",
                                                                                       trControl = paramGrid,
                                                                                       tuneLength = 3)
    print(paste("finished:",names(triage.trainSet[i]),".",classifiers[[c]]," at ",Sys.time(), sep = ""))
  }
}
cui_to_english <- read.csv("~/Genentech/v3/_select_cui_text_count_from_select_ROW_NUMBER_OVER_PARTITION_BY__202104281643.csv")
get_feaimp <- function(model, modelNo, dat, modelName, catName){
  tmp <- data.frame(varImp(model[[modelNo]])$importance[,1])
  tmp$rowno <- rownames(tmp)
  tmp1 <- data.frame(colnames(dat))
  tmp1$rowno <- rownames(tmp1)
  tmp <- left_join(tmp, tmp1,
                   by = 'rowno')
  colnames(tmp) <- c('score', 'rowno', 'cui')
  tmp$model <- modelName
  tmp$category <- catName
  return(tmp)
}
feaImp <- rbind(
  get_feaimp(ml.train.triage, 17, triage, 'glm', 'triage'),
  get_feaimp(ml.train.triage, 18, triage, 'rf', 'triage'),
  get_feaimp(ml.train.triage, 19, triage, 'svm', 'triage'),
  get_feaimp(ml.train.triage, 20, triage, 'xgb', 'triage'),
  get_feaimp(ml.train.provider, 17, provider, 'glm', 'provider'),
  get_feaimp(ml.train.provider, 18, provider, 'rf', 'provider'),
  get_feaimp(ml.train.provider, 19, provider, 'svm', 'provider'),
  get_feaimp(ml.train.provider, 20, provider, 'xgb', 'provider'),
  
  get_feaimp(ml.train.triage.binary, 17, triage, 'glm', 'triage.binary'),
  get_feaimp(ml.train.triage.binary, 18, triage, 'rf', 'triage.binary'),
  get_feaimp(ml.train.triage.binary, 19, triage, 'svm', 'triage.binary'),
  get_feaimp(ml.train.triage.binary, 20, triage, 'xgb', 'triage.binary'),
  get_feaimp(ml.train.provider.binary, 17, provider, 'glm', 'provider.binary'),
  get_feaimp(ml.train.provider.binary, 18, provider, 'rf', 'provider.binary'),
  get_feaimp(ml.train.provider.binary, 19, provider, 'svm', 'provider.binary'),
  get_feaimp(ml.train.provider.binary, 20, provider, 'xgb', 'provider.binary'),
  
  get_feaimp(ml.train.negated.triage, 17, triage, 'glm', 'triage.negated'),
  get_feaimp(ml.train.negated.triage, 18, triage, 'rf', 'triage.negated'),
  get_feaimp(ml.train.negated.triage, 19, triage, 'svm', 'triage.negated'),
  get_feaimp(ml.train.negated.triage, 20, triage, 'xgb', 'triage.negated'),
  get_feaimp(ml.train.negated.provider, 17, provider, 'glm', 'provider.negated'),
  get_feaimp(ml.train.negated.provider, 18, provider, 'rf', 'provider.negated'),
  get_feaimp(ml.train.negated.provider, 19, provider, 'svm', 'provider.negated'),
  get_feaimp(ml.train.negated.provider, 20, provider, 'xgb', 'provider.negated'),
  
  get_feaimp(ml.train.negated.triage.binary, 17, triage, 'glm', 'triage.binary.negated'),
  get_feaimp(ml.train.negated.triage.binary, 18, triage, 'rf', 'triage.binary.negated'),
  get_feaimp(ml.train.negated.triage.binary, 19, triage, 'svm', 'triage.binary.negated'),
  get_feaimp(ml.train.negated.triage.binary, 20, triage, 'xgb', 'triage.binary.negated'),
  get_feaimp(ml.train.negated.provider.binary, 17, provider, 'glm', 'provider.binary.negated'),
  get_feaimp(ml.train.negated.provider.binary, 18, provider, 'rf', 'provider.binary.negated'),
  get_feaimp(ml.train.negated.provider.binary, 19, provider, 'svm', 'provider.binary.negated'),
  get_feaimp(ml.train.negated.provider.binary, 20, provider, 'xgb', 'provider.binary.negated')
)
feaImp <- left_join(feaImp,
                 cui_to_english,
                 by = 'cui')
feaImp.avg <- feaImp %>%
  group_by(cui, category) %>%
  summarise(avg = mean(score)) %>%
  left_join(.,
            cui_to_english,
            by = 'cui') %>%
  ungroup()
write.xlsx(feaImp.avg, 'feaimp_avg.xlsx')
feaImp.avg.wide <- reshape2::dcast(feaImp.avg %>% select(c('category', 'avg', 'text')),
                                   text ~ category, value.var = "avg", fun.aggregate = mean)
is.na(feaImp.avg.wide) <- 0
feaImp.avg.wide[is.na(feaImp.avg.wide)] <- 0
rownames(feaImp.avg.wide) <- feaImp.avg.wide$text
tiff("feaImp.tiff", units = "in", width = 15, height = 30, res = 300)
gplots::heatmap.2(as.matrix(t(feaImp.avg.wide[,-1])),
                  dendrogram = c("row"), 
                  #Colv = "NA",
                  trace = "none", 
                  main = "",
                  key.title = "", 
                  keysize = .5,
                  key.xlab = "CUI value",
                  key.par = list(mar=c(5,.5,1,.5)),
                  density.info = "none",
                  margins = c(2,10.5), # b, l, t, r
                  #xlab = "text", 
                  #labCol = "",
                  cexRow = .5,
                  #ColSideColors = c(rep("gray",1054) ,rep("steelblue",1814), rep("blue",733))
)
dev.off()
gplots::heatmap.2(as.matrix(feaImp.avg.wide[,-1]),
                  dendrogram = c("row"), srtCol = 30
)
gplots::heatmap.2(as.matrix(feaImp.avg.wide[,-1]),
                  dendrogram = c("row"), srtCol = 0,
                  cexCol  = 1, cexRow =.75,
                  key.title = ""
                  ,  margins = c(5, 10.5)
                  , adjCol = c(.5, 1)
                  , offsetCol = -.75 , offsetRow = 0
                  , key.par = list(mar = c(7, 2, 2, 4.5))
)
#----
# ROC curve plots
#----
      # comments::::
      # works in R but very clumsy
roc.plot <- function(model, dat, i_dat, c){
  predictions <- setNames(
    data.frame(
      dat[[i_dat]]$label,
      predict(object = model[[c]], dat[[i_dat]], type = "prob")
    ),
    c("obs","X0","X1")
  )
  return(roc(predictor = predictions$X1, response = predictions$obs, levels = rev(levels(predictions$obs))))
  rm(predictions)
}
ggroc(
  list(
    triage.pre.665 = roc.plot(ml.train.triage, preCOVID.approx, 5, 19)
  , triage.post.665 = roc.plot(ml.train.triage, postCOVID.approx, 5, 19)
  , triage.negated.pre.665 = roc.plot(ml.train.negated.triage, preCOVID.NEGATED.approx, 5, 19)
  , triage.negated.post.665 = roc.plot(ml.train.negated.triage, postCOVID.NEGATED.approx, 5, 19)
  , provider.negated.binary.pre.583 = roc.plot(ml.train.negated.provider.binary, provider.preCOVID.NEGATED.binary.approx, 5, 19)
  , provider.negated.binary.post.583 = roc.plot(ml.train.negated.provider.binary, provider.postCOVID.NEGATED.binary.approx, 5, 19)
       )
  , legacy.axes = TRUE
  , alpha = 0.5
#  , aes=c("linetype", "color")
  , size = 1.25
) +
  theme_minimal() + 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = .5)
        ) +
  ggtitle("") + xlab("FPR") + ylab("TPR") +
#+ scale_colour_manual(values = c("red", "blue", "black"))
   geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype="solid") 
  # + scale_colour_brewer(palette="YlOrRd")

# plot(roc(testdata[["outcome"]], testdata[["fi"]]), print.auc=TRUE, col="black", lty=1, lwd=2, legacy.axes = TRUE, print.auc.y=.35, grid=TRUE)

################################################################################
##########################  validation: 07/14/21  ##############################
################################################################################
misdiagnosed.TRIAGE.cases <- read_csv("misdiagnosed_TRIAGE_FORMATTED.csv") %>% select(-c(noteid))
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(triage), colnames(misdiagnosed.TRIAGE.cases))
tmp <- data.frame(matrix( , nrow = nrow(misdiagnosed.TRIAGE.cases), 
                          ncol = length(setdiff(colnames(triage), colnames(misdiagnosed.TRIAGE.cases)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
misdiagnosed.TRIAGE.cases <- cbind(misdiagnosed.TRIAGE.cases, tmp)
rm(tmp)
#misdiagnosed.TRIAGE.cases.svd <- svd(misdiagnosed.TRIAGE.cases)
#misdiagnosed.TRIAGE.cases.approx <- list()
#misdiagnosed.TRIAGE.cases.approx[[5]] <- data.frame(misdiagnosed.TRIAGE.cases.svd$u %*% diag(misdiagnosed.TRIAGE.cases.svd$d) %*% t((misdiagnosed.TRIAGE.cases.svd$v)))
#names(misdiagnosed.TRIAGE.cases.approx) <- c('approx.665')
colnames(misdiagnosed.TRIAGE.cases) <- seq(1, length(misdiagnosed.TRIAGE.cases))
misdiagnosed.TRIAGE.cases$label <- 'X1'
misdiagnosed.TRIAGE.cases$label <- factor(misdiagnosed.TRIAGE.cases$label, levels = c('X0', 'X1'))
misdiagnosed.TRIAGE.cases.results <- as.data.frame(list())
for (c in 17:20) {
    predictions <- setNames(
      data.frame(
        misdiagnosed.TRIAGE.cases$label,
        predict(object = ml.train.triage[[c]], misdiagnosed.TRIAGE.cases, type = "prob"),
        predict(object = ml.train.triage[[c]], misdiagnosed.TRIAGE.cases, type = "raw")
      ),
      c("obs","X0","X1","pred")
    )
    cm <- confusionMatrix(
      reference = predictions$obs,
      data = predictions$pred,
      mode = "everything",
      positive = "X1"
    )
    tmp <- as.data.frame(t(rbind(
      fetchResults(cm$byClass, ml.train.triage[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.triage[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.triage[c])                                             # Fetch TP,FP,FN,TN                                             
    )))
    misdiagnosed.TRIAGE.cases.results <- rbind(misdiagnosed.TRIAGE.cases.results,
                              tmp%>%
                                mutate(
                                  "Classifier" = names(ml.train.triage[c]),
                                  "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                )%>%
                                select(
                                  c("Classifier","Accuracy","95%CI",NIR = "AccuracyNull",
                                    "Kappa","Sensitivity","Specificity",
                                    "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                )
    )
    rm(tmp, cm, predictions)
}

misdiagnosed.PROVIDER.cases <- read_csv("misdiagnosed_PROVIDER_FORMATTED.csv") %>% select(-c(noteid))
#-- prepare current data to align(in terms of no.of columns) with trained data sets
diff.matrix.colnames <- setdiff(colnames(provider), colnames(misdiagnosed.PROVIDER.cases))
tmp <- data.frame(matrix( , nrow = nrow(misdiagnosed.PROVIDER.cases), 
                          ncol = length(setdiff(colnames(provider), colnames(misdiagnosed.PROVIDER.cases)))
))
colnames(tmp) <- diff.matrix.colnames
tmp[is.na(tmp)] <- 0
misdiagnosed.PROVIDER.cases <- cbind(misdiagnosed.PROVIDER.cases, tmp)
rm(tmp)
colnames(misdiagnosed.PROVIDER.cases) <- seq(1, length(misdiagnosed.PROVIDER.cases))
misdiagnosed.PROVIDER.cases$label <- 'X1'
misdiagnosed.PROVIDER.cases$label <- factor(misdiagnosed.PROVIDER.cases$label, levels = c('X0', 'X1'))
misdiagnosed.PROVIDER.cases.results <- as.data.frame(list())
for (c in 17:20) {
  predictions <- setNames(
    data.frame(
      misdiagnosed.PROVIDER.cases$label,
      predict(object = ml.train.provider[[c]], misdiagnosed.PROVIDER.cases, type = "prob"),
      predict(object = ml.train.provider[[c]], misdiagnosed.PROVIDER.cases, type = "raw")
    ),
    c("obs","X0","X1","pred")
  )
  cm <- confusionMatrix(
    reference = predictions$obs,
    data = predictions$pred,
    mode = "everything",
    positive = "X1"
  )
  tmp <- as.data.frame(t(rbind(
    fetchResults(cm$byClass, ml.train.provider[c]),                                                              # Fetch Recall,Specificity,Precision
    fetchResults(cm$overall, ml.train.provider[c]),                                                              # Fetch Accuracy,95%CI
    fetchResults(as.data.frame(cm$table)$Freq,ml.train.provider[c])                                             # Fetch TP,FP,FN,TN                                             
  )))
  misdiagnosed.PROVIDER.cases.results <- rbind(misdiagnosed.PROVIDER.cases.results,
                                             tmp%>%
                                               mutate(
                                                 "Classifier" = names(ml.train.provider[c]),
                                                 "95%CI"= paste0("(",round(AccuracyLower,3),",",round(AccuracyUpper,3),")")
                                               )%>%
                                               select(
                                                 c("Classifier","Accuracy","95%CI",NIR = "AccuracyNull",
                                                   "Kappa","Sensitivity","Specificity",
                                                   "Precision","F1","Prevalence",TP ="4",FP ="2",FN ="3",TN ="1")
                                               )
  )
  rm(tmp, cm, predictions)
}

write.xlsx(misdiagnosed.PROVIDER.cases.results,"results_misdiagnosied.xlsx", row.names = F, quote = F)



#-------------
# Extraction parameters
#-------------
#rm(list=ls(pattern="^results"))
incrementStart <- function(x){4*x-2}
incrementEnd <- function(x){4*x}
rf.parameters <- as.data.frame(list())
svm.parameters <- as.data.frame(list())
xgb.parameters <- as.data.frame(list())
models.list.20s <- list(
  'ml.train.negated.provider',
  'ml.train.negated.provider.binary',
  'ml.train.negated.triage',
  'ml.train.negated.triage.binary',
  'ml.train.provider',
  'ml.train.provider.binary',
  'ml.train.triage',
  'ml.train.triage.binary'
)
for(MLmodel in models.list.20s){
  for (i in 1:5){
    for (c in incrementStart(i):incrementEnd(i)) {
      ifelse(c %in% seq(from = 2, to = 18, by = 4), rf.parameters <- rbind(rf.parameters, get(MLmodel)[[c]]$bestTune %>%
                                                                  `row.names<-`(paste0(MLmodel, '.', names(get(MLmodel)[c])))
                                                                ),
             ifelse(c %in% seq(from = 3, to = 19, by = 4), svm.parameters <- rbind(svm.parameters, get(MLmodel)[[c]]$bestTune %>%
                                                                          `row.names<-`(paste0(MLmodel, '.', names(get(MLmodel)[c])))
                                                                        ),
                                                    xgb.parameters <- rbind(xgb.parameters, get(MLmodel)[[c]]$bestTune %>%
                                                                              `row.names<-`(paste0(MLmodel, '.', names(get(MLmodel)[c])))
                                            )
             )
      )
    }
  }
}
models.list.45s <- list(
  'ml.train.neg.triage.binary.w.controls3.simulations',
  'ml.train.neg.triage.w.controls3.simulations',
  'ml.train.triage.binary.w.controls3.simulations',
  'ml.train.triage.w.controls3.simulations'
)
incrementStart <- function(x){3*x-1}
incrementEnd <- function(x){3*x}
for(MLmodel in models.list.45s){
  for (i in 1:15){
    for (c in incrementStart(i):incrementEnd(i)) {
      ifelse(c %in% seq(from = 2, to = 44, by = 3), rf.parameters <- rbind(rf.parameters, get(MLmodel)[[c]]$bestTune %>%
                                                                  `row.names<-`(paste0(MLmodel, '.', names(get(MLmodel)[c])))
      ),
      ifelse(c %in% seq(from = 3, to = 45, by = 3), xgb.parameters <- rbind(xgb.parameters, get(MLmodel)[[c]]$bestTune %>%
                                                                              `row.names<-`(paste0(MLmodel, '.', names(get(MLmodel)[c])))
      ),
             NULL
      )
      )
    }
  }
}
