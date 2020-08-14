#:::::::::::::::::::::
# Data Pre-processing
#::::::::::::::::::::
library(readxl)
library(dplyr)
library(mice)
library(caTools)
library(caret)
library(corrplot)
## Custom method to re-define levels of categorical variables
adjust_cat_levels<-function(dat){
  return(dat%>%
           mutate(
             GENDER=if_else(dat$PT_SEX=="Female",1,0),                          # Re-code Patient Sex to a binary variable
             SMOKE_STTS=if_else(dat$SMOKE_STTS=="CURRENT SMOKER",1,0),          # Dissolve Smoking Status to only 2 levels
             INDEX_INSURANCE_TYPE=if_else(dat$INDEX_INSURANCE_TYPE %in% c("Self Pay","Special Billing","VA",NA),         # Combine fewer% insurance categories into OTHERS
                                          "Others",dat$INDEX_INSURANCE_TYPE)  
           ))
}
## Selected features
features<-c(
  "GENDER","BP_SYSTOLIC","BP_DIASTOLIC","SMOKE_STTS","FAM_HEART_HIST","FAM_STROKE_HIST","ATRIAL_FIB_PREINDEX","ATRIAL_FLUTTER_PREINDEX","AFIB_FLUTTER_PREINDEX","HYPERTENSION_PREINDEX",
  "MI_PREINDEX","DIABETES_PREINDEX","DYSLIPIDEMIA_PREINDEX","CHF_PREINDEX","HYPERCOAG_STATES_PREINDEX","CHRONIC_LIVER_DIS_PREINDEX","CHRONIC_LUNG_DIS_PREINDEX","RHEUM_DIS_PREINDEX",
  "CHRONIC_KIDNEY_DIS_PREINDEX","NEOPLASM_PREINDEX","PERI_VASC_DIS_PREINDEX","PFO_PREINDEX","PAST_ISCHEMIC_STROKE_PREINDEX","PAST_HEMORRHAGIC_STROKE_PREINDEX","PFO_ALL_TIME",
  "ASPIRIN_PRIOR","CLOPIDOGREL_PRIOR","DIPYRIDAMOLE_PRIOR","COUMAD_WARF_PRIOR","ORAL_ANTICOAG_PRIOR","STATINS_PRIOR","ANTI_HTN_PRIOR","AGE_AT_INDEX","INDEX_INSURANCE_TYPE",
  "DAYS_BTW_LASTOP_INDEX","BMI_MEDIAN_BEFORE_INDEX","HB_MEDIAN_BEFORE_INDEX","HBA1C_MEDIAN_BEFORE_INDEX","LDL_MEDIAN_BEFORE_INDEX","PLT_MEDIAN_BEFORE_INDEX",
  "WBC_MEDIAN_BEFORE_INDEX","MIGRAINE_PREINDEX","CONVULSIONS_PREINDEX","EPILEPSY_PREINDEX","DEPRESSION_PREINDEX","MANIC_BIPOLAR_DIS_PREINDEX",
  "ANXIETY_DIS_PREINDEX","CONVERSION_DIS_PREINDEX","SYNCOPE_PREINDEX","ALCOHOL_DEP_ABUSE_PREINDEX","OPIOID_DEP_ABUSE_PREINDEX",
  "CANNABIS_DEP_ABUSE_PREINDEX","COCAINE_DEP_ABUSE_PREINDEX","OTHER_DRUG_DEP_ABUSE_PREINDEX","MULTIPLE_SCLEROSIS_PREINDEX",
  "ESRD_PREINDEX","PERIPHERAL_NEUROPATHY_PREINDEX","BRAIN_TUMOR_PREINDEX","HEPATIC_ENCEPHALOPATHY_PREINDEX",
  "CIRRHOSIS_PREINDEX","MENIERE_VERTIGO_PREINDEX","CREATININE_CLOSEST_TO_INDEX"
)
## Read CONTROLS 
CONTROLS<-list()
for (dat in 1:3) {
  CONTROLS[[dat]]<-adjust_cat_levels(
    readxl::read_excel(paste("CONTROLS",dat,"_DATABASE_v1.4_8.01.2020.xlsx",sep = ""),
                       na=c("NA"),col_names = T
    )%>%
      filter(ENC_TYPE %in% c("ED ONLY", "ED TO IP"))        # Consider only 'ED ONLY' and 'ED to IP' encounters
  )%>%
    select("PT_ID",features)%>%
    mutate(indicatorMissingLDL=if_else(is.na(LDL_MEDIAN_BEFORE_INDEX),1,0),     # Create a binary variable to indicate patient with missing LDL value
           indicatorMissingA1C=if_else(is.na(HBA1C_MEDIAN_BEFORE_INDEX),1,0)    # Create a binary variable to indicate patient with missing A1c value
    )
}
names(CONTROLS)<-paste0("CONTROLS",seq_along(CONTROLS))
## Impute missing values
CONTROLS.afterImp<-list()
for (dat in 1:3){                                                               # Create S3:mids MICE objects for each CONTROL set
  CONTROLS.afterImp[[dat]]<-mice(CONTROLS[[dat]][features],
                                 m=25,maxit = 25,method = "pmm",seed=99)
}
names(CONTROLS.afterImp)<-paste0("CONTROLS",seq_along(CONTROLS.afterImp))
## Read CASES
CASES<-adjust_cat_levels(
  read_excel("GNSIS_DATABASE_v7.5.1_7.15.2020.xlsx",na=c("NA"),col_names = T))%>%
  select("PT_ID",features)%>%
  mutate(indicatorMissingLDL=if_else(is.na(LDL_MEDIAN_BEFORE_INDEX),1,0),       # Create a binary variable to indicate patient with missing LDL value
         indicatorMissingA1C=if_else(is.na(HBA1C_MEDIAN_BEFORE_INDEX),1,0)      # Create a binary variable to indicate patient with missing A1c value
  )
jiang.ImpLabs<-read_excel("jiangLabs.xlsx",na=c("NA"),col_names=T)
#x<-round(colSums(is.na(CASES))/length(CASES$PT_ID)*100,0)                      # Compute missing %
CASES<-left_join(CASES,
                 jiang.ImpLabs,
                 by=c("PT_ID"))%>%
  mutate(                                                                       # Replace only the missing values in Durgesh's labs with Jiang's imputed lab values
    CREATININE_CLOSEST_TO_INDEX=coalesce(CREATININE_CLOSEST_TO_INDEX.x,CREATININE_CLOSEST_TO_INDEX.y),
    HB_MEDIAN_BEFORE_INDEX=coalesce(HB_MEDIAN_BEFORE_INDEX.x,HB_MEDIAN_BEFORE_INDEX.y),
    HBA1C_MEDIAN_BEFORE_INDEX=coalesce(HBA1C_MEDIAN_BEFORE_INDEX.x,HBA1C_MEDIAN_BEFORE_INDEX.y),
    LDL_MEDIAN_BEFORE_INDEX=coalesce(LDL_MEDIAN_BEFORE_INDEX.x,LDL_MEDIAN_BEFORE_INDEX.y),
    PLT_MEDIAN_BEFORE_INDEX=coalesce(PLT_MEDIAN_BEFORE_INDEX.x,PLT_MEDIAN_BEFORE_INDEX.y),
    WBC_MEDIAN_BEFORE_INDEX=coalesce(WBC_MEDIAN_BEFORE_INDEX.x,WBC_MEDIAN_BEFORE_INDEX.y)
  )%>%
  select("PT_ID","indicatorMissingLDL","indicatorMissingA1C",all_of(features))
#select(order(colnames(.)))                                                   # Order columns alphabetically
CASES.afterImp<-mice(CASES[features],
                     m=25,maxit = 25,method = "pmm",seed=99)
## Extract data sets from MICE objects
pull.data<-function(PT_ID,indicatorMissingLDL,indicatorMissingA1C,impDat){
  dat<-cbind(PT_ID,indicatorMissingLDL,indicatorMissingA1C,mice::complete(impDat))%>%
    filter(AGE_AT_INDEX>=18)
}
## Create the final data set to be used for ML
GNSIS<-rbind(
  # CASES imputed
  pull.data(CASES$PT_ID,CASES$indicatorMissingLDL,CASES$indicatorMissingA1C,CASES.afterImp)%>%
    mutate(label=1),
  # CONTROLS imputed
  rbind(
    pull.data(CONTROLS[[1]]$PT_ID,CONTROLS[[1]]$indicatorMissingLDL,CONTROLS[[1]]$indicatorMissingA1C,CONTROLS.afterImp[[1]]),
    pull.data(CONTROLS[[3]]$PT_ID,CONTROLS[[3]]$indicatorMissingLDL,CONTROLS[[3]]$indicatorMissingA1C,CONTROLS.afterImp[[3]])
  )%>%
    mutate(label=0)
)%>%
  select(label,everything())
colnames(GNSIS[apply(GNSIS,2,function(x){all(x %in% 0:1)})==F])
## Evaluate correlation among continuous variables 
corr<-cor(GNSIS%>%
            select("AGE_AT_INDEX","BP_SYSTOLIC","BP_DIASTOLIC","DAYS_BTW_LASTOP_INDEX","CREATININE_CLOSEST_TO_INDEX",
                   "BMI_MEDIAN_BEFORE_INDEX","HB_MEDIAN_BEFORE_INDEX","HBA1C_MEDIAN_BEFORE_INDEX",
                   "LDL_MEDIAN_BEFORE_INDEX","PLT_MEDIAN_BEFORE_INDEX","WBC_MEDIAN_BEFORE_INDEX"))
png("pearsonCorr.png",height=800, width=800)
corrplot(corr,
         type = "upper",order = "hclust",
         col=brewer.pal(n=8,name="RdYlBu"))
dev.off()