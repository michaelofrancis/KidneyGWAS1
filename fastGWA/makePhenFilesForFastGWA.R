library(tidyverse)


#Full cohort not sex-stratified

full<-read.table("/scratch/mf91122/CCC/pheno200kresid/KidneyPhenoFullpc.txt", header=TRUE, stringsAsFactors=FALSE)

full<-as_tibble(full)

phenotypes<-c("Creatinine", "Cystatin_C", "eGFR1", "eGFR2", "eGFR3", 
              "eGFR4", "BMI", "SBP", "DBP", "HbA1c", "LDL", "HDL", 
              "TAGs", "TC", "Waist_circumference")

phenotypes<-paste(phenotypes, "resinv", sep="_")

for (j in 1:length(phenotypes)){
	write.table(
		full[,c("FID", "IID", phenotypes[j])], 
		paste("/scratch/mf91122/CCC/pheno200kresid/phen/", phenotypes[j], 
			".phen", sep=""), 
		row.names=FALSE, quote=FALSE)
}



write.table(full[,c("FID", "IID", "PC1", "PC2", "PC3", "PC4",
			"PC5", "PC6", "PC7", "PC8", "PC9", "PC10")],
		"/scratch/mf91122/CCC/pheno200kresid/phen/KidneyPhenoCovarPC.txt",
		row.names=FALSE, quote=FALSE)
