library(qqman)
library(tidyverse)

pheno<-c("Creatinine", "Cystatin_C", "eGFR1", "eGFR2", "eGFR3", 
              "eGFR4", "BMI", "SBP", "DBP", "HbA1c", "LDL", "HDL", 
              "TAGs", "TC", "Waist_circumference")

pheno<-paste(pheno, "resinv", sep="_")
#pheno<-"Creatinine_resinv"

for (j in 1:length(pheno)){
	outdir=paste("/scratch/mf91122/CCC/exomeQC200k/G.fastGWA-04082021/plot-04102021/", pheno[j], "/", sep="")
	file<-read.delim(paste("/scratch/mf91122/CCC/exomeQC200k/G.fastGWA-04082021/",pheno[j],"/geno_assoc1.fastGWA", sep=""))
	file<-as_tibble(file)
	file2<-file%>%select(CHR, SNP, POS, P)
	colnames(file2)<-c("CHR", "SNP", "BP", "P")
	dir.create(outdir,recursive = TRUE)
	png(filename=paste(outdir,pheno[j],"_manhattan.png", sep=""), type="cairo")
	manhattan(file2,
        	main = paste(pheno[j], " fastGWA UK Biobank \n n= 158565 (71602 male)"),
		ylim = c(0, 40), cex = 0.6,
	        annotatePval = FALSE, annotateTop = FALSE,
	        cex.axis = 0.9, col = c('red', 'black'),
	        suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
	        chrlabs = as.character(1:22)
        )#end manhattan

	dev.off()		

	png(filename=paste(outdir,pheno[j],"_QQ.png", sep=""), type="cairo")
	qq(file2$P,
		main = paste(pheno[j], " fastGWA UK Biobank \n n= 158565 (71602 male)"))
	dev.off()

} #end pheno for loop
