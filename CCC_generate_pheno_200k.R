#Mike Francis
#11-10-2020
#Prepare phenotype files for Cystatin C/ Creatinine /eGFR GWAS analysis

library(plyr)
library(dplyr)
library(tidyverse)
library(DescTools)
source('manyColsToDummy.R')

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Load data-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

source('ukb34137_loaddata.r') #generates bd (UKB dataset)

ID<-read.table("ukb23155_c1_b0_v1_s200631.fam", header=FALSE, 
            stringsAsFactors=FALSE)
ID<-ID$V1

bdkeep<-bd[bd$f.eid %in% ID,]

bd_add<-read.table("ukb42606.tab",
                   header=TRUE, sep="\t")

bd_add<-as_tibble(bd_add)

bdkeepfull<-inner_join(bdkeep, bd_add, by=c("f.eid", "f.22000.0.0"))
bdkeepfull<-as_tibble(bdkeepfull)
bdkeepfull

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Extract relevant columns-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#30700 = blood creatinine
#30720 = blood Cystatin C
#NEED TO ADD WAIST CIRCUMFERENCE COLUMN WHEN I GET IT
#f.48.0.0, "Waist_circumference",

new<-bdkeepfull%>%select(f.eid, f.31.0.0, f.21003.0.0, f.21000.0.0, 
                         f.30700.0.0, f.30720.0.0, f.21001.0.0,
                         f.54.0.0, 
                         f.4080.0.0, f.4079.0.0, 
                         f.30750.0.0, f.30780.0.0, f.30760.0.0,
                         f.30870.0.0, f.30690.0.0,
                         f.22001.0.0, f.22027.0.0, f.22019.0.0, 
                         f.22006.0.0
                         )

colnames(new)<-c("IID", "Sex", "Age",  "Race",                   
                 "Creatinine", "Cystatin_C", "BMI",
                 "Assessment_center", 
                 "SBP", "DBP",
                 "HbA1c", "LDL", "HDL",
                 "TAGs", "TC",
                 "Genetic_Sex","Outliers_for_het_or_missing", "SexchrAneuploidy",
                 "Genetic_ethnic_grouping"
                )
              

new$Age2<-new$Age^2

   
nrow(new) #Starting sample size = 200631

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Apply QC-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#Filter by white caucasian genetic race only
new<-new%>%filter(Genetic_ethnic_grouping == "Caucasian")
nrow(new) #167248. Number removed = 33383

new<-new%>%
    filter(is.na(Outliers_for_het_or_missing) | Outliers_for_het_or_missing !="Yes") 
nrow(new) #166953. Number removed = 295

new<-new%>%
    filter(is.na(SexchrAneuploidy) | SexchrAneuploidy != "Yes")
nrow(new) #166754. Number removed = 199

#If Sex does not equal genetic sex, exclude participant
new<-new[new$Sex == new$Genetic_Sex,] 
nrow(new) #166753. Number removed = 1

new$Sex<-mapvalues(as.character(new$Sex), 
                   c("Male", "Female"), c(1,2))

#Remove QC columns whose information has already been used.
new<-new%>%select(-Genetic_Sex, -Outliers_for_het_or_missing, -SexchrAneuploidy, -Genetic_ethnic_grouping)

#Remove participants who do not have measurements for Creatinine or 
#Cystatin_C (7873 missing both, 198 missing just Creatinine, 
#117 missing just Cystatin_C; 8188 total)

new2<-new[!is.na(new$Creatinine),]
new2<-new2[!is.na(new2$Cystatin_C),]
nrow(new2) #158565. Number removed = 8188

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Calculate eGFR-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#Creatinine default units umol/L 
#Formula needs mg/dL
#To convert μmol/l to mg/dL, divide by 88.4. 
new2$Scr<-new2$Creatinine/88.42

#Cystatin C default units: mg/L
#Formula needs mg/L. No conversion necessary

#eGFR1: MDRD equation
new2$eGFR1<-0 #initialize col
new2$eGFR1[new2$Sex=="1"]<-175*
    (new2$Scr[new2$Sex=="1"]^(-1.154))* 
    (new2$Age[new2$Sex=="1"]^(-0.203))

new2$eGFR1[new2$Sex=="2"]<-175*
    (new2$Scr[new2$Sex=="2"]^(-1.154))* 
    (new2$Age[new2$Sex=="2"]^(-0.203)) * 0.742 #if female

summary(new2$eGFR1)

#eGFR2: CKD-EPI Cystatin equation: 
new2$eGFR2<-0 #initialize col
new2$eGFR2[new2$Sex=="1"]<-133*
    (pmin((new2$Cystatin_C[new2$Sex=="1"])/0.8, 1))^(-0.499) * 
    (pmax(new2$Cystatin_C[new2$Sex=="1"]/0.8, 1))^(-1.328) * 
    0.996^(new2$Age[new2$Sex=="1"])

new2$eGFR2[new2$Sex=="2"]<-133*
    (pmin((new2$Cystatin_C[new2$Sex=="2"])/0.8, 1))^(-0.499) * 
    (pmax(new2$Cystatin_C[new2$Sex=="2"]/0.8, 1))^(-1.328) * 
    0.996^(new2$Age[new2$Sex=="2"]) * 0.932 #if female


#eGFR3: CKD-EPI Creatinine equation: 
new2$eGFR3<-0 #initialize col
new2$eGFR3[new2$Sex=="1"]<-141*
    (pmin((new2$Scr[new2$Sex=="1"])/0.9, 1))^(-0.411) * 
    (pmax(new2$Scr[new2$Sex=="1"]/0.9, 1))^(-1.209) * 
    0.993^(new2$Age[new2$Sex=="1"])

new2$eGFR3[new2$Sex=="2"]<-141*
    (pmin((new2$Scr[new2$Sex=="2"])/0.7, 1))^(-0.329) * 
    (pmax(new2$Scr[new2$Sex=="2"]/0.7, 1))^(-1.209) * 
    0.993^(new2$Age[new2$Sex=="2"]) *
    1.018 #if female



#eGFR4: CKD-EPI Creatinine-Cystatin equation: 
new2$eGFR4<-0 #initialize col
new2$eGFR4[new2$Sex=="1"]<-135*
    (pmin((new2$Scr[new2$Sex=="1"])/0.9, 1))^(-0.207) *
    (pmax(new2$Scr[new2$Sex=="1"]/0.9, 1))^(-0.601) * 
    (pmin((new2$Cystatin_C[new2$Sex=="1"])/0.8, 1))^(-0.375) *
    (pmax((new2$Cystatin_C[new2$Sex=="1"])/0.8, 1))^(-0.711) *
    0.995^(new2$Age[new2$Sex=="1"])

new2$eGFR4[new2$Sex=="2"]<-135*
    (pmin((new2$Scr[new2$Sex=="2"])/0.7, 1))^(-0.248) *
    (pmax(new2$Scr[new2$Sex=="2"]/0.7, 1))^(-0.601) * 
    (pmin((new2$Cystatin_C[new2$Sex=="2"])/0.8, 1))^(-0.375) *
    (pmax((new2$Cystatin_C[new2$Sex=="2"])/0.8, 1))^(-0.711) *
    0.995^(new2$Age[new2$Sex=="2"]) * 0.969 #for females


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Winsorize eGFR cols-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


new2$eGFR1<-Winsorize(new2$eGFR1, minval = 15, maxval = 200, 
                     probs = c(0,1), type=1)
new2$eGFR2<-Winsorize(new2$eGFR2, minval = 15, maxval = 200, 
                     probs = c(0,1), type=1)
new2$eGFR3<-Winsorize(new2$eGFR3, minval = 15, maxval = 200, 
                     probs = c(0,1), type=1)
new2$eGFR4<-Winsorize(new2$eGFR4, minval = 15, maxval = 200, 
                     probs = c(0,1), type=1)
new2


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Regress phenotypes, extract residuals-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#ref: https://github.com/lindgrengroup/fatdistnGWAS/commit/35ca6d65dbdb1fc3cf417432dc0575805ce80878


phenotypes<-c("Creatinine", "Cystatin_C", "eGFR1", "eGFR2", "eGFR3", 
              "eGFR4", "BMI", "SBP", "DBP", "HbA1c", "LDL", "HDL", 
              "TAGs", "TC")

resdat<-matrix(NA,nrow(new2),length(phenotypes))
colnames(resdat)<-phenotypes

for (p in 1:length(phenotypes)){
    assign(paste("lm", phenotypes[p], sep="_"), lm(new2[[phenotypes[p]]] ~ Age + Age2 + 
        Assessment_center + Sex, data=new2, na.action=na.exclude))
    lmname<-paste("lm", phenotypes[p], sep="_")
    lmobj<-get(lmname)
    resdat[,p]<-resid(lmobj)
}



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Inverse Normalization-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

resdat_inv<-matrix(NA,nrow(new2),length(phenotypes))
colnames(resdat_inv)<-phenotypes

for (i in 1:length(phenotypes)){
    resdat_inv[,i] <- qnorm((rank(resdat[,i], na.last="keep")-0.5)/sum(!is.na(resdat[,i])))
    }

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#SEX STRATIFIED --Regress phenotypes, extract residuals=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
new2M<-new2%>%filter(Sex==1)
new2F<-new2%>%filter(Sex==2)



resdatM<-matrix(NA,nrow(new2M),length(phenotypes))
resdatF<-matrix(NA,nrow(new2F),length(phenotypes))

colnames(resdatM)<-phenotypes
colnames(resdatF)<-phenotypes

for (p in 1:length(phenotypes)){
    assign(paste("lm", phenotypes[p], sep="_"), lm(new2M[[phenotypes[p]]] ~ Age + Age2 + 
                Assessment_center, data=new2M, na.action=na.exclude))
    lmname<-paste("lm", phenotypes[p], sep="_")
    lmobj<-get(lmname)
    resdatM[,p]<-resid(lmobj)
}

for (p in 1:length(phenotypes)){
    assign(paste("lm", phenotypes[p], sep="_"), lm(new2F[[phenotypes[p]]] ~ Age + Age2 + 
                Assessment_center, data=new2F, na.action=na.exclude))
    lmname<-paste("lm", phenotypes[p], sep="_")
    lmobj<-get(lmname)
    resdatF[,p]<-resid(lmobj)
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Inverse Normalization-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

resdat_invM<-matrix(NA,nrow(new2M),length(phenotypes))
resdat_invF<-matrix(NA,nrow(new2F),length(phenotypes))

colnames(resdat_invM)<-phenotypes
colnames(resdat_invF)<-phenotypes

for (i in 1:length(phenotypes)){
    resdat_invM[,i] <- qnorm((rank(resdatM[,i], na.last="keep")-0.5)/sum(!is.na(resdatM[,i])))
}

for (i in 1:length(phenotypes)){
    resdat_invF[,i] <- qnorm((rank(resdatF[,i], na.last="keep")-0.5)/sum(!is.na(resdatF[,i])))
}


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Write Tables-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

write.table(resdat_inv, "resdat_inv.txt", quote = FALSE, row.names = FALSE)
write.table(resdat_invM, "resdat_invM.txt", quote = FALSE, row.names = FALSE)
write.table(resdat_invF, "resdat_invF.txt", quote = FALSE, row.names = FALSE)

write.csv(resdat_inv, "resdat_inv.csv", quote = FALSE, row.names = FALSE)
write.csv(resdat_invM, "resdat_invM.csv", quote = FALSE, row.names = FALSE)
write.csv(resdat_invF, "resdat_invF.csv", quote = FALSE, row.names = FALSE)


