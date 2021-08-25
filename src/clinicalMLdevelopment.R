#------------------
# Correlation
#------------------
corr<-cor(GNSIS%>%
      select("AGE_AT_INDEX","BP_SYSTOLIC","BP_DIASTOLIC","DAYS_BTW_LASTOP_INDEX","CREATININE_CLOSEST_TO_INDEX",
             "BMI_MEDIAN_BEFORE_INDEX","HB_MEDIAN_BEFORE_INDEX","HBA1C_MEDIAN_BEFORE_INDEX",
             "LDL_MEDIAN_BEFORE_INDEX","PLT_MEDIAN_BEFORE_INDEX","WBC_MEDIAN_BEFORE_INDEX"))
tiff("correlation.tiff", units = "in", width = 8, height = 8, res = 300, compression = 'jpeg')
corrplot(corr,
         type = "upper",order = "hclust",
         col=brewer.pal(n=8,name="RdYlBu"))
dev.off()
#------------------
# Feature Importance
#------------------
library(Boruta)
boruta_output <- Boruta(label~.,
                        data=GNSIS %>% 
                          select(-c("PT_ID")),
                        doTrace=1)
getSelectedAttributes(boruta_output,withTentative = TRUE)
TentativeRoughFix(boruta_output)
imps <- attStats(TentativeRoughFix(boruta_output))
imps2 <- imps[imps$decision !='Rejected',c('meanImp','decision')]
imps2[order(-imps2$meanImp),]
#------------------
# ML Training
#------------------
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
classifiers<-c("rf")
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
ml.trained<-train(label~.,
                  data=trainSet,
                  method="glm",
                  preProcess=c("center","scale"),                               # Normalize the data
                  metric="ROC",                                                 # Metric to optimize
                  trControl=paramGrid,
                  tuneLength=5
                  )
ml.train<-list()
ml.train[["logReg.withIndicators"]]<-ml.trained
#------------------
# Evaluate the trained models
#------------------
library(pROC)
library(e1071)
predictions<-setNames(
  data.frame(
    testSet$label,
    predict(object = ml.train.allVars[[4]], testSet, type = "prob"),
    predict(object = ml.train.allVars[[4]], testSet, type = "raw")
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
#------------------
# Summary Stats for paper
#------------------
library(dplyr)
library(openxlsx)
# CASES.afterImp use for cases
binaryCols<-names(Filter(function(x) x,apply(mice::complete(CONTROLS.afterImp[[3]]),2,function(x) { all(x %in% 0:1) })))
tmp<-t(mice::complete(CONTROLS.afterImp[[3]])%>%
         select(-c(binaryCols,"INDEX_INSURANCE_TYPE"))%>%
         #group_by(.)%>%
         group_by(n=n())%>% # to get length
         summarize_all(.,
                       funs(
                         missing = sum(is.na(.)),
                         mean = mean,
                         q25 = quantile(.,c(.25)),
                         q75 = quantile(.,c(.75)),
                         sd = sd
                       ),
                       na.rm=T))%>%
  round(.,3)
tmp <- data.frame(colSums(mice::complete(CONTROLS.afterImp[[3]])%>%
                 select(c(binaryCols))))

openxlsx::write.xlsx(tmp, 'tmp.xlsx', row.names = T)

tmp<-t(CONTROLS[[3]]%>%
         select(binaryCols)%>%
         group_by(n=n())%>%
         summarise_all(.,list(sum)))

# Feature Importance
feaImp_clinicalGenentech%>%
  ggplot(aes(y=Variable, x=Score, fill =Model))+
  geom_bar(stat = "identity", position = "stack")+
  facet_wrap(.~Model)+
  theme_minimal()+
  theme(legend.position = "bottom")

feaImp_clinicalGenentech%>%
  select(Variable,Score)%>%
  group_by(Variable)%>%
  summarise_all(funs(mean=mean))%>%
  arrange(desc(mean))%>%
  ggplot(aes(y=reorder(Variable, -mean), x=mean))+
  geom_bar(stat = "identity", position = "stack", fill = "steelblue")+
  theme_minimal()+
  scale_x_continuous(limits = c(0,100), expand = c(0,0))+
  labs(x = "Avg. Feature Importance Score")+
  theme(axis.title.y = element_blank())

#------------------------------
# remove PFO all time and re run the models
#----------------------------
features <- c(
  "GENDER","BP_SYSTOLIC","BP_DIASTOLIC","SMOKE_STTS","FAM_HEART_HIST","FAM_STROKE_HIST","ATRIAL_FIB_PREINDEX","ATRIAL_FLUTTER_PREINDEX","AFIB_FLUTTER_PREINDEX","HYPERTENSION_PREINDEX",
  "MI_PREINDEX","DIABETES_PREINDEX","DYSLIPIDEMIA_PREINDEX","CHF_PREINDEX","HYPERCOAG_STATES_PREINDEX","CHRONIC_LIVER_DIS_PREINDEX","CHRONIC_LUNG_DIS_PREINDEX","RHEUM_DIS_PREINDEX",
  "CHRONIC_KIDNEY_DIS_PREINDEX","NEOPLASM_PREINDEX","PERI_VASC_DIS_PREINDEX","PFO_PREINDEX","PAST_ISCHEMIC_STROKE_PREINDEX","PAST_HEMORRHAGIC_STROKE_PREINDEX",
  "ASPIRIN_PRIOR","CLOPIDOGREL_PRIOR","DIPYRIDAMOLE_PRIOR","COUMAD_WARF_PRIOR","ORAL_ANTICOAG_PRIOR","STATINS_PRIOR","ANTI_HTN_PRIOR","AGE_AT_INDEX","INDEX_INSURANCE_TYPE",
  "DAYS_BTW_LASTOP_INDEX","BMI_MEDIAN_BEFORE_INDEX","HB_MEDIAN_BEFORE_INDEX","HBA1C_MEDIAN_BEFORE_INDEX","LDL_MEDIAN_BEFORE_INDEX","PLT_MEDIAN_BEFORE_INDEX",
  "WBC_MEDIAN_BEFORE_INDEX","MIGRAINE_PREINDEX","CONVULSIONS_PREINDEX","EPILEPSY_PREINDEX","DEPRESSION_PREINDEX","MANIC_BIPOLAR_DIS_PREINDEX",
  "ANXIETY_DIS_PREINDEX","CONVERSION_DIS_PREINDEX","SYNCOPE_PREINDEX","ALCOHOL_DEP_ABUSE_PREINDEX","OPIOID_DEP_ABUSE_PREINDEX",
  "CANNABIS_DEP_ABUSE_PREINDEX","COCAINE_DEP_ABUSE_PREINDEX","OTHER_DRUG_DEP_ABUSE_PREINDEX","MULTIPLE_SCLEROSIS_PREINDEX",
  "ESRD_PREINDEX","PERIPHERAL_NEUROPATHY_PREINDEX","BRAIN_TUMOR_PREINDEX","HEPATIC_ENCEPHALOPATHY_PREINDEX",
  "CIRRHOSIS_PREINDEX","MENIERE_VERTIGO_PREINDEX","CREATININE_CLOSEST_TO_INDEX"
)
ml.train.allVars.v2 <- list()
startTime <- list()
finishTime <- list()
classifiers <- c("glm", "rf" ,"xgbDART")
for (c in 1:length(classifiers)) {
  print(paste("started:",classifiers[[c]],"at",Sys.time()))
  startTime[[classifiers[[c]]]]<-Sys.time()
  ml.train.allVars.v2[[classifiers[[c]]]]<-train(label~.,
                                              data = trainSet,
                                              method = classifiers[[c]],
                                              preProcess = c("center","scale"),
                                              metric = "ROC",
                                              trControl = paramGrid,
                                              tuneLength = 3)
  finishTime[[classifiers[[c]]]]<-Sys.time()
  print(paste("finished:",classifiers[[c]],"at",Sys.time()))
}

#--------------------
# Validation: 03/15/2021
#-------------------
validation.results.v1 <- as.data.frame(list())
for (c in 1:3) {
  predictions <- setNames(
    data.frame(
      testSet$label,
      predict(object = ml.train.allVars.v2[[c]], testSet, type = "prob"),
      predict(object = ml.train.allVars.v2[[c]], testSet, type = "raw")
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
    fetchResults(cm$byClass, ml.train.allVars.v2[c]),                                                              # Fetch Recall,Specificity,Precision
    fetchResults(cm$overall, ml.train.allVars.v2[c]),                                                              # Fetch Accuracy,95%CI
    fetchResults(as.data.frame(cm$table)$Freq,ml.train.allVars.v2[c]),                                             # Fetch TP,FP,FN,TN
    roc(predictor = predictions$X1,response = predictions$obs,levels = rev(levels(predictions$obs)))$auc,      # Calculate AUROC
    prSummary(predictions, lev = rev(levels(predictions$obs)))[1]                                              # Calculate AUPR
  )))
  validation.results.v1 <- rbind(validation.results.v1,
                          tmp%>%
                            mutate(
                              "Classifier" = names(ml.train.allVars.v2[c]),
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
#--------------------
# Validation: 07/15/2021
#-------------------
validation.cases <- read_excel('~/Genentech/v1/CASES_FOR_CLINICAL_MODEL_NO_MRN.xlsx', sheet = 'copy') %>% 
  select(-c(INDEX_DT, ID)) 
validation.cases$label <- factor(validation.cases$label, levels = c('X0', 'X1'))
validation.results <- as.data.frame(list())
for (c in 1:3) {
    predictions <- setNames(
      data.frame(
        validation.cases$label,
        predict(object = ml.train.allVars.v2[[c]], validation.cases, type = "prob"),
        predict(object = ml.train.allVars.v2[[c]], validation.cases, type = "raw")
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
      fetchResults(cm$byClass, ml.train.allVars.v2[c]),                                                              # Fetch Recall,Specificity,Precision
      fetchResults(cm$overall, ml.train.allVars.v2[c]),                                                              # Fetch Accuracy,95%CI
      fetchResults(as.data.frame(cm$table)$Freq,ml.train.allVars.v2[c])                                             # Fetch TP,FP,FN,TN
    )))
    validation.results <- rbind(validation.results,
                            tmp%>%
                              mutate(
                                "Classifier" = names(ml.train.allVars.v2[c]),
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
