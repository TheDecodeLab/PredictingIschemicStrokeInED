#::::::::::::::
# ML Training
#:::::::::::::
library(caret)
OneHotEncode<-function(dat){
  data.encoded<-data.frame(predict(                                             # One hot encode Categorical variables
    dummyVars("~.",data=dat),
    newdata = dat)
  )%>%
    mutate(label=factor(label,                                                  # Rename the label classes to valid R variable names
                        labels = make.names(levels(as.factor(as.character(label))))))
  return(data.encoded)
}
## Split GNSIS into train and test data sets
set.seed(111)
index<-caTools::sample.split(GNSIS$label,SplitRatio=.8)
trainSet<-OneHotEncode(
  subset(GNSIS,index==TRUE)%>%                                  # Stroke:7120(17.4%) No-Stroke:33883(82.6%)
    select(-c("PT_ID")))
testSet<-OneHotEncode(
  subset(GNSIS,index==FALSE)%>%                                  # Stroke:1780(17.4%) No-Stroke:8471(82.6%)
    select(-c("PT_ID")))
paramGrid<-trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        summaryFunction = twoClassSummary,                      # Evaluate performance
                        classProbs = T,                                         # Estimate class probabilities
                        allowParallel = T,
                        search = "random")
ml.train.allVars<-list()
startTime<-list()
finishTime<-list()
classifiers<-c("glm","rf","svmRadial","xgbDART")
for (c in 1:length(classifiers)) {
  print(paste("started:",classifiers[[c]],"at",Sys.time()))
  startTime[[classifiers[[c]]]]<-Sys.time()
  ml.train.allVars[[classifiers[[c]]]]<-train(label~.,
                                              data = trainSet,
                                              method = classifiers[[c]],
                                              preProcess = c("center","scale"),
                                              metric = "ROC",
                                              trControl = paramGrid,
                                              tuneLength = 3)
  finishTime[[classifiers[[c]]]]<-Sys.time()
  print(paste("finished:",classifiers[[c]],"at",Sys.time()))
}
do.call(rbind,startTime)
#:::::::::::::::::::::::::::::
# Evaluate the trained models
#::::::::::::::::::::::::::::
library(pROC)
library(e1071)
predictions<-setNames(
  data.frame(
    testSet$label,
    predict(object = ml.train.allVars[[1]], testSet, type = "prob"),
    predict(object = ml.train.allVars[[1]], testSet, type = "raw")
  ),
  c("obs","X0","X1","pred")
)
confusionMatrix(
  reference = predictions$obs,
  data = predictions$pred,
  mode = "everything",
  positive = "X1"
)
prSummary(predictions, lev = rev(levels(predictions$obs)))
pROC::roc(
  predictor = predictions$X1,
  response = predictions$obs,
  levels = rev(levels(predictions$obs))
)
