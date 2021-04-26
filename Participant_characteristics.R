suppressPackageStartupMessages(library(tidyverse))

pheno<-read.csv("KidneyPhenoFull_04252021.csv", header=TRUE, stringsAsFactors= FALSE)
pheno<-as_tibble(pheno)
pheno

#table(pheno$Race) #all British

quantpheno<-c("Age", "Creatinine", "Cystatin_C",  "BMI", 
              "SBP", "DBP", "HbA1c", "LDL", "HDL",
              "TAGs", "TC", "Waist_circumference", "eGFR1",
              "eGFR2","eGFR3","eGFR4")

catpheno<-c("Sex")

quantitative<-as.data.frame(
    pheno%>%select(quantpheno)%>%
        summarise_all( 
            funs(n = sum(!is.na(.)), 
                 min(., na.rm = TRUE),
                 max(., na.rm = TRUE),
                 mean(., na.rm = TRUE),
                 sd(., na.rm = TRUE),
                 median(., na.rm = TRUE),
                 IQR(., na.rm = TRUE),
            ))
)

categorical<-as.data.frame(
    pheno%>%select(catpheno)%>%
        summarise_all( 
            funs(n = sum(!is.na(.)), 
                 percent = (100*mean(., na.rm = TRUE))
                 
            ))
)

pheno$Sex<-pheno$Sex-1
