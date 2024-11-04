
# saving gene expression files for Pheflux metabolic flux prediction
# NOTE: one file is generated per patient and gene names must be in the form of ensembl ID's
geneExp <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/PatientData_CPMnorm.txt', sep = '\t', header=1, row.names=1)

for (i in 1:ncol(geneExp)){
  patient_data = as.data.frame(cbind(rownames(geneExp),geneExp[,i]))
  colnames(patient_data) <- c("Gene_ID", "Expression")
  write.csv(patient_data,paste0('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/GeneExpression/',colnames(geneExp)[i],'_CPMdata.csv'),
            row.names = F)
}

