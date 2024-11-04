# Saving individual patient-specific fluxome files as a matrix of all patient fluxomes

library(limma)
setwd('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/PGK1_KD/Fluxomes')
files <- dir('.')
files <- files[grep('Homo_sapiens_',files)]

Meta <- strsplit2(read.table('../../MetaData.txt',sep = '\t',header = 1)[,1],split='_')[,1]
celltypes <- unique(strsplit2(files,split="_")[,4])


bind_files <- c()
for (i in 1:length(files)){
  open <- read.table(files[i], sep = '\t')[-1,]
  open[,1] <- gsub('^R_','',open[,1])
  if (i == 1){
    bind_files <- open
  } else {
    bind_files <- cbind(bind_files,open[,2])
  }
}
rownames(bind_files) <- bind_files[,1]
bind_files <- bind_files[,-1]
colnames(bind_files) <- rownames(annots)[1:18]

# write.table(bind_files, 'PheFlux_Results_PGK1_KD.txt', sep = '\t')
