
## Analysis and Figure Replication ##

library(dataVisEasy)
library(limma)
library(ica)
library(matrixStats)
library(umap)
library(dplyr)
library(reshape2)
library(pcaMethods)
library(DESeq2)
library(purrr)
library(rsample)
library(plyr)
library(readxl)
library(UpSetR)
library(parcutils)
library(biomaRt)
library(ggpubr)
library(stringr)
library(readxl)
library(tibble)
library(networkD3)
library(ggforce)
library(scales)


# set working directory to location of files
setwd('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper')

# read in RNAseq CPM-normalized data
geneExp <- as.matrix(read.table('PatientData_CPMnorm.txt')) 
rownames(geneExp) <- geneExp[,1]
geneExp <- geneExp[,-1]
colnames(geneExp) <- geneExp[1,]
geneExp <- geneExp[-1,]
class(geneExp) <- "numeric"

# read in MetaData
meta <- read.table('MetaData.txt',header = 1, sep = '\t')

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genenames <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                       values = rownames(geneExp), mart= mart)

# RNAseq data with gene symbols instead of ensembl ID's
geneExp_geneNames <- geneExp
rownames(geneExp_geneNames) <- genenames$hgnc_symbol[match(rownames(geneExp_geneNames),genenames$ensembl_gene_id)]

## read in Human1 model
ihuman <- as.data.frame(read_excel('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/Human-GEM.xlsx', sheet= "Sheet2"))
ihuman.subsys.table <- as.data.frame(table(ihuman$SUBSYSTEM))

metgenes = grep("ENSG",unique(as.vector(strsplit2(ihuman$`GENE ASSOCIATION`,split=" or | and |\\(|\\)| "))),value=T)

metgenes <- as.data.frame(metgenes)
metgenes$V2 <- genenames$hgnc_symbol[match(metgenes[,1],genenames$ensembl_gene_id)]
colnames(metgenes)[1] <- "V1"


## subsetting RNAseq data to contain only metabolic genes from Human1 model
metgenes_patient <- geneExp[na.omit(match(metgenes$V1,rownames(geneExp))),]
metgenes_patient_geneNames <- metgenes_patient
rownames(metgenes_patient_geneNames) <- metgenes$V2[match(rownames(metgenes_patient),metgenes$V1)]


replace_enriched <- function(match) {
  ensg_id <- match
  replacement <- genenames$hgnc_symbol[which(genenames$ensembl_gene_id == ensg_id)]
  if (length(replacement) > 0) {
    return(as.character(replacement))
  } else {
    return(ensg_id)
  }
}

## ihuman model with gene symbols instead of ensembl ID's for GPR
ihuman.geneName <- ihuman
ihuman.geneName$GeneAssociationName <- str_replace_all(ihuman.geneName$`GENE ASSOCIATION`, 'ENSG\\d+', replace_enriched)


# Setting up patient annotations for plotting
annots <- as.data.frame(cbind(strsplit2(meta$disease_state,split = '_'),meta$phenotype))
rownames(annots) <- annots$V1
annots <- annots[,-1]
colnames(annots) <- colnames(meta)
annots$explant <- "No"
annots$explant[which(annots$disease_state=="explant.AH")] <- "Yes"
annots$disease_state2 = annots$disease_state
annots$disease_state2[which(annots$disease_state=="explant.AH")] <- "severe.AH"
annots$disease_state <- factor(annots$disease_state, levels=c("NASH", "HCV", "comp.cirrhosis","healthy.control", 
                                                              "early.ASH", "nonsevere.AH", "severe.AH", "explant.AH"))

annots$disease_state2 <- factor(annots$disease_state2, levels=c("NASH", "HCV", "comp.cirrhosis","healthy.control", 
                                                              "early.ASH", "nonsevere.AH", "severe.AH"))

annots$explant <- factor(annots$explant, levels=c("Yes","No"))
annots <- as.data.frame(annots)

ann_colors <- list(disease_state = c(explant.AH = '#ea472a', severe.AH = '#ea8f59', nonsevere.AH = '#f2e32b', early.ASH = '#7aaf3e',
                                     healthy.control = '#f280ae', comp.cirrhosis = '#7580aa', HCV = '#5e8f9c', NASH = '#814284'),
                   disease_state2 = c(severe.AH = '#ea8f59', nonsevere.AH = '#f2e32b', early.ASH = '#7aaf3e',
                                     healthy.control = '#f280ae', comp.cirrhosis = '#7580aa', HCV = '#5e8f9c', NASH = '#814284'),
                   explant = c(Yes = "#ea472a",No = "white"),
                   phenotype = c('NR' = 'black', 'R' = 'red'),
                   meld = c('0'='#D3D3D3','1' = '#F8BBA5', '2'='#D58D78', '3'= '#B25E4B', '4' = '#8F301E'),
                   child = c('0'='#D3D3D3','1' = '#F8BBA5', '2'='#D58D78', '3'= '#B25E4B', '4' = '#8F301E'),
                   abic = c('0'='#D3D3D3','1' = '#F8BBA5', '2'='#D58D78', '3'= '#B25E4B', '4' = '#8F301E'),
                   clinical_score = c('NA'='#D3D3D3','Low' = '#F8BBA5', 'Moderate'='#C47662', 'High'= '#8F301E'),
                   Survival = c("Alive"="red","Deceased" = "black"))




# loading in patient clinical parameters
all_clinical_params <- as.data.frame(read_excel('../DiseaseGEMPaper-Metabolites/RNAseq-ClusterA-2019 v4_upver_anonim.xlsx',sheet = "ALL"))

clinical_params <- as.data.frame(read_excel('../DiseaseGEMPaper-Metabolites/RNAseq-ClusterA-2019 v4_upver_anonim.xlsx',sheet = "Sheet5"))
meld <- clinical_params[15,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
child <- clinical_params[14,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
abic <- clinical_params[16,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
ast <- clinical_params[7,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
alt <- clinical_params[8,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
ggt <- clinical_params[9,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
alp <- clinical_params[10,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
tot.bili <- clinical_params[11,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
inr <- clinical_params[12,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]
alb <- clinical_params[13,match(rownames(annots),strsplit2(clinical_params[2,],split = '_')[,2],)]


## setting up plotting parameters
initiate_params()
set_annotations(annots)
set_annot_cols(ann_colors)
set_annot_samps(c("disease_state"))
set_scale.range(c(-1,1))


## getting gene symbols for all metabolic genes
met.genenames <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                       values = rownames(metgenes_patient), mart= mart)

metgenes_patient_geneNames <- metgenes_patient
rownames(metgenes_patient_geneNames) <- met.genenames$hgnc_symbol[match(rownames(metgenes_patient),met.genenames$ensembl_gene_id)]


all.genenames <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                       values = rownames(geneExp), mart= mart)

geneExp_geneNames <- geneExp
rownames(geneExp_geneNames) <- all.genenames$hgnc_symbol[match(rownames(geneExp),all.genenames$ensembl_gene_id)]


# loading in the raw count RNAseq data for outlier detection 
rawct <- read.table('refseq_counts_sampleinfo.tsv',sep ='\t', header = 1)
rawct <- as.data.frame(rawct)
rownames(rawct) <- rawct[,2]
allgenes_withnames <- rawct[2:nrow(rawct),c(1:2)]
rawct <- rawct[-1,-c(1:3)]
#rawct<- rawct[match(rownames(metgenes_scale_patient),rownames(rawct)),]

# checking for outliers
rawct <- as.matrix(rawct)
class(rawct) <- "numeric"
rawct <- rawct[-which(rowSums(rawct)==0),]
results <- prcomp(t(rawct), scale = TRUE)

# Supplemental Figure 1A
pca_results <- as.data.frame(results[["x"]])
pca_results$disease_state <- c(as.character(annots$disease_state[1:32]),"explant.AH",as.character(annots$disease_state[33:89]))
pca_results$disease_state <- factor(pca_results$disease_state, levels = c("NASH","HCV","comp.cirrhosis","healthy.control","early.ASH",
                                                                         "nonsevere.AH","severe.AH","explant.AH"))
pca_results$annots <- rownames(pca_results)
ggplot(pca_results, aes(x=PC1,y=PC2,color=disease_state,label=annots)) + geom_point(size = 8)+
  scale_color_manual(values=c('#814284','#5e8f9c', '#7580aa','#f280ae','#7aaf3e','#f2e32b','#ea8f59','#ea472a'))+
  geom_text(aes(label=ifelse(PC2< -400 ,as.character(annots),''),size=20),hjust=-0.5,vjust=-0.25) + xlim(c(-200,550)) +
  theme(axis.text=element_blank(),
        axis.title=element_text(size=24,face="bold"),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=1, linetype='solid'))

apply(pca_results[,1:2], 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) ))
# D33 is an outlier

# Supplemental Figure 1B
ggplot(pca_results, aes(x=PC1,y=PC2,color=disease_state,label=annots)) + geom_point(size = 8)+
  scale_color_manual(values=c('#814284','#5e8f9c', '#7580aa','#f280ae','#7aaf3e','#f2e32b','#ea8f59','#ea472a'))+
  xlim(c(-200,550)) + 
  theme(axis.text=element_blank(),
        axis.title=element_text(size=24,face="bold"),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=1, linetype='solid')) + ylim(c(-100,100))

# removing outlier sample from rawct 
rawct <- rawct[,-33]
rawct <- as.matrix(rawct)
class(rawct) <- "numeric"



## Part of Figure 1A
p <- myHeatmapByAnnotation(metgenes_patient_geneNames[-which(rowSds(metgenes_patient_geneNames)==0),], 
                           main = "Metabolic Gene Expression - Patient", groupings = "disease_state",
                           show.rownames = T, gaps.row = T, gaps.col = T, clust.rows = T, scale = "zscore", method = "euclidean",
                           show.colnames = T, clust.cols = T)



## Using raw RNAseq data for differential gene expression
dds <- DESeqDataSetFromMatrix(countData = rawct,
                              colData = annots,
                              design = ~ disease_state)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$disease_state2 <- relevel(dds$disease_state2, ref = "healthy.control")
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficient

# DEGs for early ASH vs. healthy
res.earlyASH.vs.healthy <- as.data.frame(results(dds, contrast=c("disease_state","early.ASH","healthy.control")))
res.earlyASH.vs.healthy <- res.earlyASH.vs.healthy[which(res.earlyASH.vs.healthy$padj<0.05),]
res.earlyASH.vs.healthy <- res.earlyASH.vs.healthy[order(res.earlyASH.vs.healthy$log2FoldChange),]
earlyASHup.healthydown <- res.earlyASH.vs.healthy[which(res.earlyASH.vs.healthy$log2FoldChange>0),]
earlyASHdown.healthyup <- res.earlyASH.vs.healthy[which(res.earlyASH.vs.healthy$log2FoldChange<0),]

# DEGs for severe AH vs. healthy
res.severeAH.vs.healthy <- as.data.frame(results(dds, contrast=c("disease_state","severe.AH","healthy.control")))
res.severeAH.vs.healthy <- res.severeAH.vs.healthy[which(res.severeAH.vs.healthy$padj<0.05),]
res.severeAH.vs.healthy$Gene <- rownames(res.severeAH.vs.healthy)
res.severeAH.vs.healthy$GeneName <- genenames$hgnc_symbol[match(rownames(res.severeAH.vs.healthy),genenames$ensembl_gene_id)]

## function for UMAP plotting
plot.umap <- function(x, labels,levels, colors,
                      main="UMAP") {
  
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  }
  layout <- as.data.frame(layout)
  layout$names <- rownames(layout)
  layout$labels <- as.character(labels)
  layout$labels <- factor(layout$labels,levels = levels)
  
  ggplot(layout, aes(V1, V2, color = labels))+ 
    geom_point(size=8)+ 
    scale_color_manual(values=colors) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle(main) + 
    theme_bw() + theme(panel.grid = element_blank(), 
                       plot.title = element_text(hjust = 0.5, size = 35), 
                       axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                       legend.position = "right", axis.text.x=element_text(colour="black"),
                       axis.text.y=element_text(colour="black"))
}

## Figure 1B - UMAP of all genes
set.seed(123)
data.for.umap <- scale(t(geneExp))
data.for.umap <- data.for.umap[,-which(is.na(colSds(data.for.umap)))]
data.umap <- umap(data.for.umap) #umap(t(patientdata_metgenes_scale[,c(1:18,30:89,19:29)]))
levels <- c("explant.AH", "severe.AH", "nonsevere.AH","early.ASH","healthy.control", "comp.cirrhosis", "HCV", "NASH")
colors <- c('#ea472a', '#ea8f59', '#f2e32b', '#7aaf3e','#f280ae', '#7580aa', '#5e8f9c', '#814284')
plot.umap(data.umap,as.character(annots$disease_state),levels,colors, main = "UMAP - All genes") # code to plot umap


## Figure 1C- UMAP of metabolic genes only
set.seed(123)
data.for.umap <- scale(t(metgenes_patient))
# data.for.umap <- (data.for.umap - rowMeans(data.for.umap)) / rowSds(data.for.umap)
data.for.umap <- data.for.umap[,-which(is.na(colSds(data.for.umap)))]
data.umap <- umap(data.for.umap) #umap(t(patientdata_metgenes_scale[,c(1:18,30:89,19:29)]))
levels <- c("explant.AH", "severe.AH", "nonsevere.AH","early.ASH","healthy.control", "comp.cirrhosis", "HCV", "NASH")
colors <- c('#ea472a', '#ea8f59', '#f2e32b', '#7aaf3e','#f280ae', '#7580aa', '#5e8f9c', '#814284')
plot.umap(data.umap,as.character(annots$disease_state),levels,colors, main = "UMAP - metabolic genes") # code to plot umap


## Figure 1D - UMAP of non-metabolic genes only
set.seed(123)
patientdata_NOmetgenes <- geneExp[-na.omit(match(metgenes[,1],rownames(geneExp))),]
data.for.umap <- scale(t(patientdata_NOmetgenes))
# data.for.umap <- (data.for.umap - rowMeans(data.for.umap)) / rowSds(data.for.umap)
data.for.umap <- data.for.umap[,-which(is.na(colSds(data.for.umap)))]
data.umap <- umap(data.for.umap) #umap(t(patientdata_metgenes_scale[,c(1:18,30:89,19:29)]))
levels <- c("explant.AH", "severe.AH", "nonsevere.AH","early.ASH","healthy.control", "comp.cirrhosis", "HCV", "NASH")
colors <- c('#ea472a', '#ea8f59', '#f2e32b', '#7aaf3e','#f280ae', '#7580aa', '#5e8f9c', '#814284')
plot.umap(data.umap,as.character(annots$disease_state),levels,colors, main = "UMAP - Non-metabolic genes") # code to plot umap




## Reading in the predicted metabolic fluxes
library(readxl)
FBA_results = as.data.frame(read_excel('Pheflux/results/Pheflux_PatientSpecific_Fluxes.xlsx', sheet = 1))
FBA_mat = as.matrix(as.data.frame(FBA_results[,3:ncol(FBA_results)]))
class(FBA_mat) <- "numeric"
Human1_Rxn = FBA_results[,1]
Subsystem = FBA_results[,2]
Subsystem = as.character(Subsystem)
rownames(FBA_mat) = Human1_Rxn

## remove zero fluxes
FBA_mat_nonzero <- FBA_mat[-which(rowSums(abs(FBA_mat))==0),]


# ecdf for setting metabolic flux cutoff
flux.melt <- melt(abs(FBA_mat_nonzero))
cdf <- ecdf(log(flux.melt$value[flux.melt$value != 0],base=10))

cdf_values <- log(flux.melt$value[flux.melt$value != 0], base = 10)
ecdf_function <- ecdf(cdf_values)
x_values <- sort(cdf_values)
y_values <- ecdf_function(x_values)

cumul_data <- data.frame(x = x_values, y = y_values)

# Calculate the slopes (differences in y values)
cumul_data <- cumul_data %>%
  mutate(slope = c(0, diff(y)))

# Find the knee point (maximum difference in slopes)
cumul_data <- cumul_data %>%
  mutate(diff_slope = c(0, diff(slope)))

# Identify the index of the knee point
knee_index <- which.max(abs(cumul_data$diff_slope[which(cumul_data$x>-7.5)])) + which(cumul_data$x>-7.5)[1]

## Supplemental Figure 1C
ggplot(cumul_data[c(seq(1,nrow(cumul_data),1000)),], aes(x = x, y = y)) +
  geom_step(size=1) +
  geom_segment(aes(x = x[knee_index], y = y[knee_index],xend = x[knee_index], yend = 0),linetype = "dashed",size=0.75) +
  geom_segment(aes(x =  min(x), y = y[knee_index],xend = x[knee_index], yend = y[knee_index]),linetype = "dashed",size=0.75) +
  geom_point(data = cumul_data[knee_index, ], aes(x = x, y = y), color = "red", size = 8) +
  labs(title = "Cumulative Distribution Function",
       x = "Log(Nonzero Fluxes)",
       y = "Cumulative Probability") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5,size = 24),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 0.5),  # Axis lines
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 16)) + ylim(c(0,1))

# knee point = 0.001100236
10^cumul_data$x[knee_index]


## applying metabolic flux cutoff
incl1 = which(apply(abs(FBA_mat) > 0.001, 1, any))
FBA_mat_incl1 = FBA_mat[incl1,]
Subsystem_incl1 = Subsystem[incl1]
Human1_Rxn_incl1 = Human1_Rxn[incl1]


## Figure 2A
set.seed(123)
data.for.umap <- scale(t(abs(FBA_mat)))
data.for.umap <- data.for.umap[,-which(is.na(colSds(data.for.umap)))]
data.umap <- umap(data.for.umap)
levels <- c("explant.AH", "severe.AH", "nonsevere.AH","early.ASH","healthy.control", "comp.cirrhosis", "HCV", "NASH")
colors <- c('#ea472a', '#ea8f59', '#f2e32b', '#7aaf3e','#f280ae', '#7580aa', '#5e8f9c', '#814284')
plot(data.umap,as.character(annots$disease_state),levels,colors, main = "UMAP - Predicted Metabolic Fluxes")
data <-as.matrix(data.umap$layout)
umap_coords_res1 <- data.umap$layout

## Choosing optimal k for kmeans clustering using Silhouette Method
library(cluster)
# Compute silhouette scores for each k (from 2 to 10 clusters)
sil_width <- vector()

for (k in 2:10) {
  kmeans_result <- kmeans(data, centers = k)
  silhouette_score <- silhouette(kmeans_result$cluster, dist(data))
  sil_width[k] <- mean(silhouette_score[, 3])
}

# Plot the silhouette scores
silhouette_plot <- data.frame(
  k = 2:10,
  silhouette = sil_width[2:10]
)

## Supplemental Figure 2A
ggplot(silhouette_plot, aes(x = k, y = silhouette)) +
  geom_point(size= 4) +
  geom_line() +
  xlab("Number of Clusters (k)") +
  ylab("Average Silhouette Width") +
  ggtitle("Silhouette Method for Optimal k") + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 40), 
                     #strip.text = element_text(size = 25),
                     strip.background.x = element_blank(), 
                     #legend.position = "default", 
                     axis.title.y = element_text(size = 20), 
                     # axis.title.x = element_blank(), 
                     axis.text.x = element_text(size = 20),
                     panel.background = element_rect(fill = "white", color = NA),  # White background
                     panel.grid.major = element_blank(),  # No major grid lines
                     panel.grid.minor = element_blank(),  # No minor grid lines
                     panel.border = element_blank(),  # No panel border
                     axis.line = element_line(color = "black", size = 1),  # Axis lines
                     axis.text = element_text(color = "black", size = 16),
                     axis.title = element_text(size = 18))

## kmeans clusters for Figure 2A
kmean <- kmeans(data.umap$layout,2)
plot(data.umap,as.character(kmean$cluster),c(1,2),colors, main = "UMAP - Predicted Metabolic Fluxes") # code to plot umap

## new plotting annotations
patient_groups <- annots
patient_groups$umap_allflux <- kmean$cluster

## Figure 2B
set.seed(123)
data.for.umap <- scale(t(abs(FBA_mat_incl1)))
data.umap <- umap(data.for.umap)
levels <- c("explant.AH", "severe.AH", "nonsevere.AH","early.ASH","healthy.control", "comp.cirrhosis", "HCV", "NASH")
colors <- c('#ea472a', '#ea8f59', '#f2e32b', '#7aaf3e','#f280ae', '#7580aa', '#5e8f9c', '#814284')
plot(data.umap,as.character(annots$disease_state),levels,colors, main = "UMAP - Predicted Metabolic Fluxes") 
umap_coords_res2 <- data.umap$layout
data <-as.matrix(data.umap$layout)

## Choosing optimal k for kmeans clustering using Silhouette Method
# Compute silhouette scores for each k (from 2 to 10 clusters)
sil_width <- vector()

for (k in 2:10) {
  kmeans_result <- kmeans(data, centers = k)
  silhouette_score <- silhouette(kmeans_result$cluster, dist(data))
  sil_width[k] <- mean(silhouette_score[, 3])
}

# Plot the silhouette scores
silhouette_plot <- data.frame(
  k = 2:10,
  silhouette = sil_width[2:10]
)

## Supplemental Figure 2B
ggplot(silhouette_plot, aes(x = k, y = silhouette)) +
  geom_point(size= 4) +
  geom_line() +
  xlab("Number of Clusters (k)") +
  ylab("Average Silhouette Width") +
  ggtitle("Silhouette Method for Optimal k") + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 40), 
                     #strip.text = element_text(size = 25),
                     strip.background.x = element_blank(), 
                     #legend.position = "default", 
                     axis.title.y = element_text(size = 20), 
                     # axis.title.x = element_blank(), 
                     axis.text.x = element_text(size = 20),
                     panel.background = element_rect(fill = "white", color = NA),  # White background
                     panel.grid.major = element_blank(),  # No major grid lines
                     panel.grid.minor = element_blank(),  # No minor grid lines
                     panel.border = element_blank(),  # No panel border
                     axis.line = element_line(color = "black", size = 1),  # Axis lines
                     axis.text = element_text(color = "black", size = 16),
                     axis.title = element_text(size = 18))

## kmeans clusters for Figure 2B
set.seed(1)
kmean <- kmeans(data.umap$layout,3)
plot(data.umap,as.character(kmean$cluster),c(1,2,3),colors, main = "UMAP - Predicted Metabolic Fluxes") # code to plot umap
patient_groups$umap_gt0.1 <- kmean$cluster+2

 

## Sankey plot in Figure 2C
sankey.df <- patient_groups[,which(colnames(patient_groups)%in%c("umap_allflux","umap_gt0.1"))]
sankey.df <- as.data.frame(apply(sankey.df,2,function(x) paste0("Cluster",x)))
freq1 <- table(sankey.df) %>% melt() %>% filter(value != 0)
colnames(freq1) <- c("source","target","value")

links <- freq1
links$source <- as.character(links$source)
links$target <- as.character(links$target)

nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

my_color <- 'd3.scaleOrdinal() .domain(["Cluster1","Cluster2","Cluster3", "Cluster4", "Cluster5"]) .range(["gray","#000004FF", "#37245c","#56106EFF","#89226AFF"])'

## Figure 2C - sankey plot
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",
              fontSize = 0, nodeWidth = 300,
              colourScale = my_color)

# updating plotting annotations
patient_groups$umap_gt0.1 <- factor(patient_groups$umap_gt0.1, levels = c(3,4,5))
ann_colors[["umap_gt0.1"]] = c("3"="#37245c","4"="#56106EFF", "5"="#89226AFF")
ann_colors[["umap_allflux"]] = c("1"="gray", "2"="#000004FF")

set_annotations(patient_groups)
set_annot_samps(c("umap_gt0.1","disease_state"))
set_annot_cols(ann_colors)

## Figure 2C - heatmap
p = myHeatmapByAnnotation(abs(FBA_mat_incl1), groupings = c("disease_state","umap_gt0.1"),
              main = "Fluxes > 0.01", show.colnames = F,
              scale.rows = "zscore",show.rownames = F, method = "euclidean",
              gaps.row = FALSE, clust.rows = T, clust.cols = T, row.groups = 4,
              gap.width = 2,groupings.gaps = c(0,1))

get_row_clusters <- extractClusters(abs(FBA_mat_incl1),to.extract = "genes",nclusters = 4,heatmap = p)
get_row_clusters <- cbind(get_row_clusters, ihuman.geneName[match(rownames(get_row_clusters),ihuman.geneName$rxnRetired),])

flux.annots <- get_row_clusters

## setting up plotting annotations
initiate_params()
set_annotations(patient_groups)
set_annot_samps(c("umap_gt0.1","umap_allflux","disease_state"))
set_annot_cols(ann_colors)


## Mapping metabolic flux annotations onto UMAP of metabolic gene expression
set.seed(123)
data.for.umap <- scale(t(metgenes_patient))
data.for.umap <- data.for.umap[,-which(is.na(colSds(data.for.umap)))]
data.umap <- umap(data.for.umap)
layout <- as.data.frame(data.umap$layout)
layout$names <- rownames(layout)
layout$labels <- as.character(annots$disease_state)
layout$labels <- factor(layout$labels,levels = levels)
layout$umap_gt0.1 <- factor(patient_groups$umap_gt0.1, level = c(3,4,5))
layout$umap_allflux <- factor(patient_groups$umap_allflux, level = c(1,2))

## Supplemental Figure 2C
ggplot(layout, aes(V1, V2, color = umap_allflux,shape=umap_allflux))+ 
  geom_point(size=8) + 
  scale_color_manual(values=c("darkgray","#000004FF")) +
  scale_shape_manual(values=c(16,17)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("UMAP - Predicted Metabolic Fluxes") + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 35), 
                     axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                     legend.position = "right", axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black")) 

## Supplemental Figure 2D
ggplot(layout, aes(V1, V2, color = umap_gt0.1,shape=umap_gt0.1))+ 
  geom_point(size=8) + 
  scale_color_manual(values=c("#37245c","#56106EFF","#89226AFF")) +
  scale_shape_manual(values=c(16,17,15)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("UMAP - Predicted Metabolic Fluxes") + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 35), 
                     axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                     legend.position = "right", axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black")) 





## overrepresentation analysis for each flux group
fluxsubset<- list()
for (fluxgroup in c("A","B","C","D")){
  
  fluxsubset[[fluxgroup]] <- as.data.frame(table(flux.annots$SUBSYSTEM[which(flux.annots$Gene.Groups==fluxgroup)]))
  fluxsubset[[fluxgroup]] <- cbind(fluxsubset[[fluxgroup]],ihuman.subsys.table$Freq[match(fluxsubset[[fluxgroup]][,1],ihuman.subsys.table$Var1)])
  fluxsubset[[fluxgroup]]$Pct <- fluxsubset[[fluxgroup]][,2]/fluxsubset[[fluxgroup]][,3]
  colnames(fluxsubset[[fluxgroup]]) <- c("Subsytem",paste0("Freq.Group",fluxgroup),"Model.Freq","Subsys.Pct")
  
  fluxsubset[[fluxgroup]]$fisher.oddsratio <- NA
  fluxsubset[[fluxgroup]]$fisher.pval <- NA
  ## Overrepresentation analysis
  for (subsys in fluxsubset[[fluxgroup]][,1]){
    total = sum(fluxsubset[[fluxgroup]][,3])
    a = fluxsubset[[fluxgroup]][which(fluxsubset[[fluxgroup]][,1]==subsys),2]
    b = fluxsubset[[fluxgroup]][which(fluxsubset[[fluxgroup]][,1]==subsys),3] - a
    c = sum(fluxsubset[[fluxgroup]][,2]) - a
    d = total - a - b - c
    deTable <- matrix(c(a,b,c,d), nrow = 2, dimnames = list(DE=c("yes","no"), GeneSet = c("in","out")))
    fisher_testresult <- fisher.test(deTable, alternative = "greater")
    fluxsubset[[fluxgroup]]$fisher.oddsratio[which(fluxsubset[[fluxgroup]][,1]==subsys)] <- 
      fisher_testresult$estimate
    fluxsubset[[fluxgroup]]$fisher.pval[which(fluxsubset[[fluxgroup]][,1]==subsys)] <-
      fisher_testresult$p.value
    
  }
}


# frequency tables for metabolic subsystems
FBA_mat_incl1_subsys.table <- data.frame(Rxn=rownames(FBA_mat_incl1))
FBA_mat_incl1_subsys.table <- cbind(FBA_mat_incl1_subsys.table,ihuman.geneName[match(FBA_mat_incl1_subsys.table[,1],ihuman.geneName$rxnRetired),3:6])
FBA_mat_incl1_subsys.table <- as.data.frame(table(FBA_mat_incl1_subsys.table$SUBSYSTEM))
FBA_mat_incl1_subsys.table$total <- ihuman.subsys.table$Freq[match(FBA_mat_incl1_subsys.table$Var1,ihuman.subsys.table$Var1)]
FBA_mat_incl1_subsys.table$Pct <- (FBA_mat_incl1_subsys.table$Freq / FBA_mat_incl1_subsys.table$total) * 100


FBA_mat_nonzero_subsys.table <- data.frame(Rxn=rownames(FBA_mat_nonzero))
FBA_mat_nonzero_subsys.table <- cbind(FBA_mat_nonzero_subsys.table,ihuman.geneName[match(FBA_mat_nonzero_subsys.table[,1],ihuman.geneName$rxnRetired),3:6])
FBA_mat_nonzero_subsys.table <- as.data.frame(table(FBA_mat_nonzero_subsys.table$SUBSYSTEM))
FBA_mat_nonzero_subsys.table$total <- ihuman.subsys.table$Freq[match(FBA_mat_nonzero_subsys.table$Var1,ihuman.subsys.table$Var1)]
FBA_mat_nonzero_subsys.table$Pct <- (FBA_mat_nonzero_subsys.table$Freq / FBA_mat_nonzero_subsys.table$total) * 100



GPR <- read.table("GPR.txt", sep = '\t')
GPR_FBA_incl_1 <- GPR[match(rownames(FBA_mat_incl1),rownames(GPR)),]



## pseudotime analysis of metabolic genes for healthy and ALD patient samples
library(slingshot)
library(grDevices)
library(RColorBrewer)

## redo metabolic gene expression UMAP to get coordinates
set.seed(123)
disease_state = annots$disease_state
data.for.umap <- scale(t(metgenes_patient))
# data.for.umap <- (data.for.umap - rowMeans(data.for.umap)) / rowSds(data.for.umap)
data.for.umap <- data.for.umap[,-which(is.na(colSds(data.for.umap)))]
data.umap <- umap(data.for.umap) #umap(t(patientdata_metgenes_scale[,c(1:18,30:89,19:29)]))
levels <- c("explant.AH", "severe.AH", "nonsevere.AH","early.ASH","healthy.control", "comp.cirrhosis", "HCV", "NASH")
colors <- c('#ea472a', '#ea8f59', '#f2e32b', '#7aaf3e','#f280ae', '#7580aa', '#5e8f9c', '#814284')
plot.umap(data.umap,as.character(disease_state),levels,colors, main = "UMAP - metabolic genes") # code to plot umap


reducedDim_data <- data.umap$layout[c(1:51,80:89),]
cluster_labels <- as.character(disease_state[c(1:51,80:89)])

# Slingshot for pseudotime calculation
sds <- slingshot(reducedDim_data, cluster_labels, start.clus="healthy.control", end.clus="explant.AH")

pseudotime <- slingPseudotime(sds)
colnames(pseudotime) <- apply(pseudotime, 2,function(x) paste0("healthy.control",'to',cluster_labels[which.max(x)]))
row.names(pseudotime) <- cluster_labels

# Figure 3A
plotcol = as.vector(ann_colors$disease_state[match(annots$disease_state[c(1:51,80:89)],names(ann_colors$disease_state))])
plot(reducedDim_data, col = plotcol, pch=16, asp = 1,cex=3)
lines(SlingshotDataSet(sds), lwd=2, col='black')

# Figure 3B
colors <- colorRampPalette(brewer.pal(9, 'Greys'))(100)
plotcol <- colors[cut(pseudotime, breaks=100)]
plot(reducedDim_data, col = "black",bg = plotcol, pch=21, asp = 1,cex=3)
lines(SlingshotDataSet(sds), lwd=2, col='black')


# creating scalebar for Figure 3B
df <- data.frame(
  x = 1,
  y = seq(0, 1, length.out = 100),
  color = colors
)
  
ggplot(df, aes(x = x, y = y, fill = color)) +
  geom_tile(width = 1) +
  scale_fill_identity() +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.y = element_line(color = "black")
  ) +
  labs(y = "Values", fill = "Color")




## corelation analysis between calculated pseudotime and metabolic fluxes
all.flux.pseudo.corrs <- c()
pseudo <- as.numeric(pseudotime[,1])
for (j in rownames(FBA_mat_incl1)){
  findflux <- abs(as.numeric(FBA_mat_incl1[match(j,rownames(FBA_mat_incl1)),c(1:51,80:89)]))
  get.cor <- cor.test(pseudo,findflux, method = 'spearman')
  if (!(is.na(get.cor$estimate[['rho']]))){
    all.flux.pseudo.corrs <- rbind(all.flux.pseudo.corrs,c(j,get.cor$estimate[['rho']]))
  }
}

# subsetting for only positive correlations
all.flux.pseudo.corrs <- as.data.frame(all.flux.pseudo.corrs)
all.flux.pseudo.corrs$V2 <- as.numeric(all.flux.pseudo.corrs$V2)
all.flux.pseudo.corrs <- cbind(all.flux.pseudo.corrs,ihuman.geneName[match(all.flux.pseudo.corrs$V1,ihuman.geneName$rxnRetired),3:6])
all.flux.pseudo.corrssigcor <- all.flux.pseudo.corrs[which(all.flux.pseudo.corrs$V2< -0.5 | all.flux.pseudo.corrs$V2>0.5),]
all.flux.pseudo.corrssigcor <- as.data.frame(all.flux.pseudo.corrssigcor)
all.flux.pseudo.corrssigcor <- all.flux.pseudo.corrssigcor[order(all.flux.pseudo.corrssigcor$V2),]

# tabulating the associated subsystems
all.flux.pseudo.poscor <- all.flux.pseudo.corrs[which(all.flux.pseudo.corrs$V2>0.5),]
all.flux.pseudo.poscor <- all.flux.pseudo.poscor[!duplicated(all.flux.pseudo.poscor[,c(2,5)]),]
all.flux.pseudo.poscor.subsys.table <- as.data.frame(table(all.flux.pseudo.poscor$SUBSYSTEM))
all.flux.pseudo.poscor.subsys.table$total <- ihuman.subsys.table$Freq[match(all.flux.pseudo.poscor.subsys.table$Var1,ihuman.subsys.table$Var1)]
all.flux.pseudo.poscor.subsys.table$Pct <- (all.flux.pseudo.poscor.subsys.table$Freq / all.flux.pseudo.poscor.subsys.table$total) * 100
all.flux.pseudo.poscor.subsys.table <- all.flux.pseudo.poscor.subsys.table[which(all.flux.pseudo.poscor.subsys.table$Freq>5),]
all.flux.pseudo.poscor.subsys.table <- all.flux.pseudo.poscor.subsys.table[order(all.flux.pseudo.poscor.subsys.table$Pct, decreasing = T),]


## Figure 3C
ggplot(all.flux.pseudo.poscor.subsys.table[1:10,], aes(x = reorder(Var1, Pct), y = Pct)) +
  geom_point(aes(size = Pct, color = Freq)) +
  scale_size_continuous(name = "Percentage",range = c(4,12)) +
  scale_color_gradientn(colors = c("blue", "red"),
                        values = rescale(c(min(all.flux.pseudo.poscor.subsys.table$Freq[1:10]), max(all.flux.pseudo.poscor.subsys.table$Freq[1:10]))),
                        name = "Frequency") +
  labs(title = "Top 10 Positive Pseudotime/Flux Correlations", x = "Subsystems", y = "Percentage") + ylim(c(15,50))+
  coord_flip() +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size=24),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),  # No panel border
        axis.line = element_line(color = "black", size = 0.5),  # Axis lines
        axis.text.x = element_text(color = "black", size = 22),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title = element_text(size = 18))

## colorbar for Figure 3C
df <- data.frame(
  x = 1,
  y = seq(0, 1, length.out = 100),
  color = colorRampPalette(c("blue", "red"))(100)
)

ggplot(df, aes(x = x, y = y, fill = color)) +
  geom_tile(width = 1) +
  scale_fill_identity() +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.y = element_line(color = "black")
  ) +
  labs(y = "Values", fill = "Color")

# adding pseudotime to annotations
annots$pseudotime <- "Pseudo0"
annots$pseudotime[c(1:51,80:89)] <- paste0("Pseudo",c(1:61))
colors <- colorRampPalette(brewer.pal(9, 'Greys'))(100)
plotcol <- colors[cut(pseudotime, breaks=100)]
ann_colors[['pseudotime']] <- setNames(c("red",plotcol),c("Pseudo0",annots$pseudotime[c(1:51,80:89)]))

# subsetting for only fluxes in PPP, glycolysis, or TCA cycle metabolic subsystems
flux.groups <- data.frame(reactions = c(ihuman$rxnRetired[which(ihuman$SUBSYSTEM=="Pentose phosphate pathway")],
                                        ihuman$rxnRetired[which(ihuman$SUBSYSTEM=="Glycolysis / Gluconeogenesis")],
                                        ihuman$rxnRetired[which(ihuman$SUBSYSTEM=="Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism")]))

flux.groups$Subsystem <- c(rep("PPP",length(ihuman$rxnRetired[which(ihuman$SUBSYSTEM=="Pentose phosphate pathway")])),
                           rep("Glycolysis",length(ihuman$rxnRetired[which(ihuman$SUBSYSTEM=="Glycolysis / Gluconeogenesis")])),
                           rep("TCA",length(ihuman$rxnRetired[which(ihuman$SUBSYSTEM=="Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism")])))
flux.groups <- flux.groups[which(flux.groups$reactions%in%all.flux.pseudo.poscor$V1),]
rownames(flux.groups) <- flux.groups$reactions

initiate_params()
set_annotations(annots)
set_annot_samps(c("disease_state","pseudotime"))
set_annotations.genes(flux.groups)
set_annot_genes("Subsystem")
set_annot_cols(ann_colors)
set_scale.range(c(-1,1))

## Figure 3D
p <- myHeatmapByAnnotation(abs(FBA_mat_incl1[,c(1:51,80:89)[order(pseudotime)]]), groupings.genes = "Subsystem",
                           rownames(flux.groups),
                           scale = "zscore", method = "euclidean", show.colnames = T,clust.cols = F)



### calculating sigificant metabolic flux differences between severe AH and healthy control patients
flux_g0.1 <- as.data.frame(abs(FBA_mat_incl1))
severeAH.vs.healthy.sigflux <- c()
for(s in rownames(flux_g0.1)){
  wilcox_test = wilcox.test(as.numeric(flux_g0.1[rownames(flux_g0.1)==s,which(annots$disease_state=="severe.AH")]),
                            as.numeric(flux_g0.1[rownames(flux_g0.1)==s,which(annots$disease_state=="healthy.control")]),
                            alternative = "two.sided")
  if(!(is.na(wilcox_test$p.value)) & wilcox_test$p.value < 0.05){
    flux_test_long <- as.data.frame(cbind(as.character(annots$disease_state[which(annots$disease_state%in%c("severe.AH","healthy.control"))]),
                                          as.numeric(flux_g0.1[rownames(flux_g0.1)==s,which(annots$disease_state%in%c("severe.AH","healthy.control"))])))
    flux_test_long$V1 <- factor(flux_test_long$V1, levels = c("severe.AH","healthy.control"))
    flux_test_long$V2 <- as.numeric(flux_test_long$V2)
    mean_by_group <- aggregate(V2 ~ V1, data = flux_test_long, FUN = mean)
    if (mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]>mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]){
      FC <- -(mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]/mean_by_group$V2[which(mean_by_group$V1=="severe.AH")])
      severeAH.vs.healthy.sigflux <- rbind(severeAH.vs.healthy.sigflux,c(s,FC,wilcox_test$p.value))
    } else {
      FC <- mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]/mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]
      severeAH.vs.healthy.sigflux <- rbind(severeAH.vs.healthy.sigflux,c(s,FC,wilcox_test$p.value))
    }
  }
}

severeAH.vs.healthy.sigflux <- cbind(severeAH.vs.healthy.sigflux,ihuman$SUBSYSTEM[match(severeAH.vs.healthy.sigflux[,1],ihuman$rxnRetired)],
                                     ihuman$EQUATION[match(severeAH.vs.healthy.sigflux[,1],ihuman$rxnRetired)],
                                     ihuman$`GENE ASSOCIATION`[match(severeAH.vs.healthy.sigflux[,1],ihuman$rxnRetired)])





## update annots
update_annots <- annots
update_annots$disease_state <- as.character(update_annots$disease_state)
update_annots$disease_state[which(annots$disease_state=="explant.AH")] = "Explant AH"
update_annots$disease_state[which(annots$disease_state=="severe.AH")] = "Severe AH"
update_annots$disease_state[which(annots$disease_state=="nonsevere.AH")] = "Nonsevere AH"
update_annots$disease_state[which(annots$disease_state=="early.ASH")] = "Early ASH"
update_annots$disease_state[which(annots$disease_state=="healthy.control")] = "Healthy"
update_annots$disease_state <- factor(update_annots$disease_state, levels = c("Healthy","Early ASH","Nonsevere AH","Severe AH", "Explant AH"))

update_cols <- list(disease_state = c("Explant AH" = '#ea472a', "Severe AH" = '#ea8f59', "Nonsevere AH" = '#f2e32b', "Early ASH" = '#7aaf3e',
                                      "Healthy" = '#f280ae', comp.cirrhosis = '#7580aa', HCV = '#5e8f9c', NASH = '#814284'))

initiate_params()
set_annotations(update_annots)
set_annot_cols(update_cols)

## Figure 5B (top)
beeswarmGenes(metgenes_patient_geneNames[match(c("SLC2A2","SLC2A3","SLC2A1"),rownames(metgenes_patient_geneNames)),c(1:51,80:89)],
              rownames(metgenes_patient_geneNames),groupby.x = "disease_state",
              color.by = "disease_state",facet.wrap = T,axis.text.x.size = 18,ncols=3) +
  geom_boxplot(alpha=0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 10)) 

## Figure 5B (bottom)
beeswarmGenes(metgenes_patient_geneNames[match(c("HKDC1","HK1","HK2"),rownames(metgenes_patient_geneNames)),c(1:51,80:89)],
              list =rownames(metgenes_patient_geneNames),color.by = "disease_state",groupby.x = "disease_state",
              facet.wrap = T,axis.text.x.size = 18, ncols = 3) +
  geom_boxplot(alpha=0)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 10)) 





glucose_import <- as.data.frame(abs(FBA_mat_incl1[match(c("HMR_5029","HMR_4396","HMR_4198","HMR_4394","HMR_4373","HMR_4928","HMR_5998","HMR_4388","HMR_4281","HMR_4930"),rownames(FBA_mat_incl1)),c(1:51,80:89)]))
glucose_import <- cbind(t(glucose_import),t(metgenes_patient_geneNames[match(c("SLC2A1","SLC2A2","SLC2A3","SLC2A5","SLC2A6","SLC2A10","HKDC1","PGM1","AGXT","HK1","SLC16A7","LDHA","LDHB"),rownames(metgenes_patient_geneNames)),c(1:51,80:89)]))
glucose_import_GPR <- GPR[match(c("HMR_5029","HMR_4396","HMR_4198","HMR_4394","HMR_4373","HMR_4928","HMR_5998","HMR_4388","HMR_4281","HMR_4930"),rownames(GPR)),c(1:51,80:89)]
rownames(glucose_import_GPR) <- paste0(rownames(glucose_import_GPR),"_GPR")
glucose_import <- cbind(glucose_import,t(glucose_import_GPR))
glucose_import <- as.data.frame(glucose_import)
glucose_import$disease_state <- annots$disease_state[c(1:51,80:89)]
glucose_import$phenotype <- annots$phenotype[c(1:51,80:89)]

## Figure 5C
ggplot(glucose_import, aes(x = HMR_5029, y = HMR_5029_GPR,color = disease_state)) +
  geom_point(size = 6) + xlab("HMR_5029 Flux") + ylab("HMR_5029 NEL") + 
  scale_color_manual(values = ann_colors$disease_state) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") + theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)
  )
cor.test(as.numeric(glucose_import[,match("HMR_5029",colnames(glucose_import))]),
         as.numeric(glucose_import[,match("HMR_5029_GPR",colnames(glucose_import)),]),
         method="spearman")

## Figure 5D
ggplot(glucose_import, aes(x = HMR_4394, y = HMR_4394_GPR,color = disease_state)) +
  geom_point(size = 6) + xlab("HMR_4394 Flux") + ylab("HMR_4394 NEL") + 
  scale_color_manual(values = ann_colors$disease_state) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") + theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)
  )

cor.test(as.numeric(glucose_import[,match("HMR_4394",colnames(glucose_import))]),
         as.numeric(glucose_import[,match("HMR_4394_GPR",colnames(glucose_import)),]),
         method="spearman")

## Figure 5E
ggplot(glucose_import, aes(x = HMR_5029, y = HMR_4394,color = disease_state)) +
  geom_point(size = 6) + xlab("HMR_5029 Flux") + ylab("HMR_4394 Flux") + 
  scale_color_manual(values = ann_colors$disease_state) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") + theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)
  )
cor.test(as.numeric(glucose_import[,match("HMR_5029",colnames(glucose_import))]),
         as.numeric(glucose_import[,match("HMR_4394",colnames(glucose_import)),]),
         method="spearman")

## Figure 6B
ggplot(glucose_import, aes(x = HMR_4928, y = HMR_4928_GPR,color = disease_state)) +
  geom_point(size = 6) + xlab("HMR_4928 Flux") + ylab("HMR_4928 NEL") + 
  scale_color_manual(values = ann_colors$disease_state) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") + theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)
  )
cor.test(as.numeric(glucose_import[,match("HMR_4928",colnames(glucose_import))]),
         as.numeric(glucose_import[,match("HMR_4928_GPR",colnames(glucose_import)),]),
         method="spearman")

## Figure 6C
ggplot(glucose_import, aes(x = HMR_4930, y = HMR_4930_GPR,color = disease_state)) +
  geom_point(size = 6) + xlab("HMR_4930 Flux") + ylab("HMR_4930 NEL") + 
  scale_color_manual(values = ann_colors$disease_state) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") + theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)
  )

cor.test(as.numeric(glucose_import[,match("HMR_4930",colnames(glucose_import))]),
         as.numeric(glucose_import[,match("HMR_4930_GPR",colnames(glucose_import)),]),
         method="spearman")

## Figure 6D
ggplot(glucose_import, aes(x = HMR_4281, y = HMR_4281_GPR,color = disease_state)) +
  geom_point(size = 6) + xlab("HMR_4281 Flux") + ylab("HMR_4281 NEL") + 
  scale_color_manual(values = ann_colors$disease_state) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") + theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)
  )
cor.test(as.numeric(glucose_import[,match("HMR_4281",colnames(glucose_import))]),
         as.numeric(glucose_import[,match("HMR_4281_GPR",colnames(glucose_import)),]),
         method="spearman")


## Supplemental Figure 4

set_annotations(annots)
set_annot_samps("disease_state")
set_annot_cols(ann_colors)
beeswarmGenes(abs(FBA_mat_incl1[match(c("HMR_9048","HMR_4896"),rownames(FBA_mat_incl1)),c(1:51,80:89)]),
              rownames(FBA_mat_incl1),groupby.x = "disease_state",
              color.by = "disease_state",facet.wrap = T,axis.text.x.size = 18,ncols=3) +
  geom_boxplot(alpha=0) +
  geom_signif(comparisons = list(c("healthy.control", "severe.AH")),map_signif_level = T,
              textsize = 8, step_increase=0.1,vjust=0.5,margin_top = 0.09) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 10)) 


## TCA Cycle for Supplemental Figure 5
tca.rxns <- c("HMR_4137","HMR_4143","HMR_4145","HMR_4456","HMR_4458","HMR_4589","HMR_3957","HMR_5297",
              "HMR_4147","HMR_4152","HMR_4652","HMR_8743","HMR_4410","HMR_4141")

tca.rxn.flux <- abs(FBA_mat_incl1[match(tca.rxns,rownames(FBA_mat_incl1)),c(1:51,80:89)])
tca.rxn.flux <- t(tca.rxn.flux)
tca.rxn.flux <- as.data.frame(tca.rxn.flux)
tca.rxn.flux <- cbind(annots$disease_state[c(1:51,80:89)],tca.rxn.flux)
colnames(tca.rxn.flux)[1] <- "V1"

find_outliers <- function(df, value_col,group) {
  df %>%
    group_by(across(all_of(group))) %>%
    mutate(
      Q1 = quantile(.data[[value_col]], 0.25),
      Q3 = quantile(.data[[value_col]], 0.75),
      IQR = Q3 - Q1,
      lower_bound = Q1 - 2*IQR,
      upper_bound = Q3 + 2*IQR,
      outlier = .data[[value_col]] < lower_bound | .data[[value_col]] > upper_bound 
    ) #| .data[[value_col]] < Q1-threshold
}


for (i in 2:ncol(tca.rxn.flux)) {
  value_col <- names(tca.rxn.flux)[i]
  tca.rxn.flux[,i] <- as.numeric(tca.rxn.flux[,i])
  get.outlier <- find_outliers(tca.rxn.flux, value_col,"V1")
  tca.rxn.flux[which(get.outlier$outlier == TRUE), i] <- NA
}

set_annotations(annots)
set_annot_samps("disease_state")
set_annot_cols(ann_colors)

## Supplemental Figure 5
myHeatmapByAnnotation(abs(t(tca.rxn.flux[,2:ncol(tca.rxn.flux)])), clust.rows = F, clust.cols = T,
                      groupings = "disease_state",scale = "zscore", method = "euclidean",groupings.gaps = 2,
                      fontsize.row = 12)


tca.rxn.flux <- as.data.frame(t(tca.rxn.flux[,2:ncol(tca.rxn.flux)]))

## differential flux analysis for Supplemental Figure 5A
tca.rxn.sigflux <- c()
for(s in rownames(tca.rxn.flux)){
  wilcox_test = wilcox.test(na.omit(as.numeric(tca.rxn.flux[rownames(tca.rxn.flux)==s,1:18])),
                            na.omit(as.numeric(tca.rxn.flux[rownames(tca.rxn.flux)==s,52:61])),
                            alternative = "two.sided")
  if(!(is.na(wilcox_test$p.value)) & wilcox_test$p.value < 0.05){
    flux_test_long <- as.data.frame(cbind(as.character(annots$disease_state[which(annots$disease_state%in%c("severe.AH","healthy.control"))]),
                                          as.numeric(tca.rxn.flux[rownames(tca.rxn.flux)==s,c(1:18,52:61)])))
    if (any(is.na(flux_test_long$V2))){
      flux_test_long <- flux_test_long[-which(is.na(flux_test_long$V2)),]}
    flux_test_long$V1 <- factor(flux_test_long$V1, levels = c("severe.AH","healthy.control"))
    flux_test_long$V2 <- as.numeric(flux_test_long$V2)
    mean_by_group <- aggregate(V2 ~ V1, data = flux_test_long, FUN = mean)
    if (mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]>mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]){
      FC <- -(mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]/mean_by_group$V2[which(mean_by_group$V1=="severe.AH")])
      tca.rxn.sigflux <- rbind(tca.rxn.sigflux,c(s,FC,wilcox_test$p.value))
    } else {
      FC <- mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]/mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]
      tca.rxn.sigflux <- rbind(tca.rxn.sigflux,c(s,FC,wilcox_test$p.value))
    }
  }
}

tca.rxn.sigflux <- cbind(tca.rxn.sigflux,ihuman$SUBSYSTEM[match(tca.rxn.sigflux[,1],ihuman$rxnRetired)],
                                     ihuman$EQUATION[match(tca.rxn.sigflux[,1],ihuman$rxnRetired)],
                                     ihuman$`GENE ASSOCIATION`[match(tca.rxn.sigflux[,1],ihuman$rxnRetired)])

## OXPHOS for Supplemental Figure 6
oxphos.rxns <- c("HMR_6921","HMR_6911","HMR_6918","HMR_6914","HMR_6916")

oxphos.flux <- abs(FBA_mat_incl1[match(oxphos.rxns,rownames(FBA_mat_incl1)),c(1:51,80:89)])
oxphos.flux <- t(oxphos.flux)
oxphos.flux <- as.data.frame(oxphos.flux)
oxphos.flux <- cbind(annots$disease_state[c(1:51,80:89)],oxphos.flux)
colnames(oxphos.flux)[1] <- "V1"


for (i in 2:ncol(oxphos.flux)) {
  value_col <- names(oxphos.flux)[i]
  oxphos.flux[,i] <- as.numeric(oxphos.flux[,i])
  get.outlier <- find_outliers(oxphos.flux, value_col,"V1")
  oxphos.flux[which(get.outlier$outlier == TRUE), i] <- NA
}
## Supplemental Figure 6B
myHeatmapByAnnotation(abs(t(oxphos.flux[,2:ncol(oxphos.flux)])), clust.rows = F, clust.cols = T,
                      groupings = "disease_state",scale = "zscore", method = "euclidean",groupings.gaps = 2,
                      fontsize.row = 12)

oxphos.flux <- as.data.frame(t(oxphos.flux[,2:ncol(oxphos.flux)]))

## differential flux analysis for Supplemental Figure 6A
oxphos.sigflux <- c()
for(s in rownames(oxphos.flux)){
  wilcox_test = wilcox.test(na.omit(as.numeric(oxphos.flux[rownames(oxphos.flux)==s,1:18])),
                            na.omit(as.numeric(oxphos.flux[rownames(oxphos.flux)==s,52:61])),
                            alternative = "two.sided")
  if(!(is.na(wilcox_test$p.value)) & wilcox_test$p.value < 0.05){
    flux_test_long <- as.data.frame(cbind(as.character(annots$disease_state[which(annots$disease_state%in%c("severe.AH","healthy.control"))]),
                                          as.numeric(oxphos.flux[rownames(oxphos.flux)==s,c(1:18,52:61)])))
    if (any(is.na(flux_test_long$V2))){
      flux_test_long <- flux_test_long[-which(is.na(flux_test_long$V2)),]}
    flux_test_long$V1 <- factor(flux_test_long$V1, levels = c("severe.AH","healthy.control"))
    flux_test_long$V2 <- as.numeric(flux_test_long$V2)
    mean_by_group <- aggregate(V2 ~ V1, data = flux_test_long, FUN = mean)
    if (mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]>mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]){
      FC <- -(mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]/mean_by_group$V2[which(mean_by_group$V1=="severe.AH")])
      oxphos.sigflux <- rbind(oxphos.sigflux,c(s,FC,wilcox_test$p.value))
    } else {
      FC <- mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]/mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]
      oxphos.sigflux <- rbind(oxphos.sigflux,c(s,FC,wilcox_test$p.value))
    }
  }
}

oxphos.sigflux <- cbind(oxphos.sigflux,ihuman$SUBSYSTEM[match(oxphos.sigflux[,1],ihuman$rxnRetired)],
                         ihuman$EQUATION[match(oxphos.sigflux[,1],ihuman$rxnRetired)],
                         ihuman$`GENE ASSOCIATION`[match(oxphos.sigflux[,1],ihuman$rxnRetired)])

### PPP for Supplemental Figure 7
ppp.rxns <- c("G6PDH2c","HMR_4625","HMR_4474","HMR_4352","HMR_4477","HMR_4501","HMR_4565")

ppp.flux <- abs(FBA_mat_incl1[match(ppp.rxns,rownames(FBA_mat_incl1)),c(1:51,80:89)])
ppp.flux <- t(ppp.flux)
ppp.flux <- as.data.frame(ppp.flux)
ppp.flux <- cbind(annots$disease_state[c(1:51,80:89)],ppp.flux)
colnames(ppp.flux)[1] <- "V1"


for (i in 2:ncol(ppp.flux)) {
  value_col <- names(ppp.flux)[i]
  ppp.flux[,i] <- as.numeric(ppp.flux[,i])
  get.outlier <- find_outliers(ppp.flux, value_col,"V1")
  ppp.flux[which(get.outlier$outlier == TRUE), i] <- NA
}

annots_subset <- annots[c(1:51,80:89),]
annots_subset$disease_state <- factor(annots_subset$disease_state,levels = c("healthy.control","early.ASH","nonsevere.AH","severe.AH","explant.AH"))

ppp.order <- c()
for (disease in levels(annots_subset$disease_state)){
  if (disease == "healthy.control"){
    clust.samps <- dist(abs(ppp.flux[match(setdiff(rownames(annots_subset)[which(annots_subset$disease_state==disease)],c("I92","I95"))
                                           ,rownames(ppp.flux)),2:ncol(ppp.flux)]),method="euclidean")
    order <- hclust(clust.samps, method = "complete")$order
    clust.samp.order <- c(setdiff(rownames(annots_subset)[which(annots_subset$disease_state==disease)],c("I92","I95"))[order],
                          c("I92","I95"))
    ppp.order <- c(ppp.order,clust.samp.order)
    
  } else {
  clust.samps <- dist(abs(ppp.flux[match(rownames(annots_subset)[which(annots_subset$disease_state==disease)],rownames(ppp.flux))
                                   ,2:ncol(ppp.flux)]),method="euclidean")
  order <- hclust(clust.samps, method = "complete")$order
  clust.samp.order <- rownames(annots_subset)[which(annots_subset$disease_state==disease)][order]
  ppp.order <- c(ppp.order,clust.samp.order)
  }
}


## Supplemental Figure 7B
myHeatmap(abs(t(ppp.flux[match(ppp.order,rownames(ppp.flux)),2:ncol(ppp.flux)])), clust.rows = F, clust.cols = F,
                      gaps.col = c(10,22,33,51),scale = "zscore", method = "euclidean",
                      fontsize.row = 12)

ppp.flux <- as.data.frame(t(ppp.flux[,2:ncol(ppp.flux)]))

## differential flux analysis for Supplemental Figure 7A
ppp.sigflux <- c()
for(s in rownames(ppp.flux)){
  wilcox_test = wilcox.test(na.omit(as.numeric(ppp.flux[rownames(ppp.flux)==s,1:18])),
                            na.omit(as.numeric(ppp.flux[rownames(ppp.flux)==s,52:61])),
                            alternative = "two.sided")
  if(!(is.na(wilcox_test$p.value)) & wilcox_test$p.value < 0.05){
    flux_test_long <- as.data.frame(cbind(as.character(annots$disease_state[which(annots$disease_state%in%c("severe.AH","healthy.control"))]),
                                          as.numeric(ppp.flux[rownames(ppp.flux)==s,c(1:18,52:61)])))
    if (any(is.na(flux_test_long$V2))){
      flux_test_long <- flux_test_long[-which(is.na(flux_test_long$V2)),]}
    flux_test_long$V1 <- factor(flux_test_long$V1, levels = c("severe.AH","healthy.control"))
    flux_test_long$V2 <- as.numeric(flux_test_long$V2)
    mean_by_group <- aggregate(V2 ~ V1, data = flux_test_long, FUN = mean)
    if (mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]>mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]){
      FC <- -(mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]/mean_by_group$V2[which(mean_by_group$V1=="severe.AH")])
      ppp.sigflux <- rbind(ppp.sigflux,c(s,FC,wilcox_test$p.value))
    } else {
      FC <- mean_by_group$V2[which(mean_by_group$V1=="severe.AH")]/mean_by_group$V2[which(mean_by_group$V1=="healthy.control")]
      ppp.sigflux <- rbind(ppp.sigflux,c(s,FC,wilcox_test$p.value))
    }
  }
}

ppp.sigflux <- cbind(ppp.sigflux,ihuman$SUBSYSTEM[match(ppp.sigflux[,1],ihuman$rxnRetired)],
                        ihuman$EQUATION[match(ppp.sigflux[,1],ihuman$rxnRetired)],
                        ihuman$`GENE ASSOCIATION`[match(ppp.sigflux[,1],ihuman$rxnRetired)])

#################
# in silico knockdown experimental results

# loading in metabolic fluxes from KD experiments
HKDC1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HKDC1_KD_4June2024/Fluxomes/PheFlux_Results_HKDC1_KD.txt', sep = '\t')
FBA_mat_HKDC1_KD <- FBA_mat
FBA_mat_HKDC1_KD <- as.data.frame(cbind(HKDC1_KD,FBA_mat_HKDC1_KD[,-match(colnames(HKDC1_KD),colnames(FBA_mat))]))

HKDC1_50pctKD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HKDC1_50pctKD_try2/Fluxomes/PheFlux_Results_HKDC1_50pctKD_try2.txt', sep = '\t')
FBA_mat_HKDC1_50pctKD <- FBA_mat
FBA_mat_HKDC1_50pctKD <- as.data.frame(cbind(HKDC1_50pctKD,FBA_mat_HKDC1_50pctKD[,-match(colnames(HKDC1_50pctKD),colnames(FBA_mat))]))


PKM_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/PKM_KD/Fluxomes/PheFlux_Results_PKM_KD.txt', sep = '\t')
FBA_mat_PKM_KD <- FBA_mat
FBA_mat_PKM_KD <- as.data.frame(cbind(PKM_KD,FBA_mat_PKM_KD[,-match(colnames(PKM_KD),colnames(FBA_mat))]))


HKDC1_PKM_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HKDC1_PKM_KD/Fluxomes/PheFlux_Results_HKDC1_PKM_KD.txt', sep = '\t')
FBA_mat_HKDC1_PKM_KD <- FBA_mat
FBA_mat_HKDC1_PKM_KD <- as.data.frame(cbind(HKDC1_PKM_KD,FBA_mat_HKDC1_PKM_KD[,-match(colnames(HKDC1_PKM_KD),colnames(FBA_mat))]))


bind_matrices <- cbind(FBA_mat_HKDC1_50pctKD[,1:39],FBA_mat_HKDC1_KD[,1:39],FBA_mat_PKM_KD[,1:39],FBA_mat_HKDC1_PKM_KD[,1:39], FBA_mat[,c(1:39,80:89)])
colnames(bind_matrices)[(seq(1,39*4,39)[1]):(seq(1,39*4,39)[2]-1)] <- paste0(colnames(bind_matrices)[(seq(1,39*4,39)[1]):(seq(1,39*4,39)[2]-1)],"_HKDC1.50pctKD")
colnames(bind_matrices)[(seq(1,39*4,39)[2]):(seq(1,39*4,39)[3]-1)] <- paste0(colnames(bind_matrices)[(seq(1,39*4,39)[2]):(seq(1,39*4,39)[3]-1)],"_HKDC1.KD")
colnames(bind_matrices)[(seq(1,39*4,39)[3]):(seq(1,39*4,39)[4]-1)] <- paste0(colnames(bind_matrices)[(seq(1,39*4,39)[3]):(seq(1,39*4,39)[4]-1)],"_PKM.KD")
colnames(bind_matrices)[(seq(1,39*4,39)[4]):(39*4)] <- paste0(colnames(bind_matrices)[(seq(1,39*4,39)[4]):(39*4)],"_HKDC1.PKM.KD")

HK.KD.annots <- as.data.frame(colnames(bind_matrices))
rownames(HK.KD.annots) <- HK.KD.annots[,1]
HK.KD.annots$state <- "No.KD"
HK.KD.annots$state[grep("_HKDC1.KD",rownames(HK.KD.annots))] <- "HKDC1.KD"
HK.KD.annots$state[grep("_HKDC1.50pctKD",rownames(HK.KD.annots))] <- "HKDC1.50pctKD"
HK.KD.annots$state[grep("_PKM.KD",rownames(HK.KD.annots))] <- "PKM.KD"
HK.KD.annots$state[grep("_HKDC1.PKM.KD",rownames(HK.KD.annots))] <- "HKDC1.PKM.KD"
HK.KD.annots$state[grep("I",rownames(HK.KD.annots))] <- "Healthy"
HK.KD.annots$state <- factor(HK.KD.annots$state, levels = c("Healthy","No.KD","HKDC1.50pctKD","HKDC1.KD","PKM.KD","HKDC1.PKM.KD"))
HK.KD.annots[,1] <- strsplit2(HK.KD.annots[,1],split="_")[,1]
HK.KD.annots$disease_state <- annots$disease_state[match(HK.KD.annots[,1],rownames(annots))]
HK.KD.annots$disease_state


initiate_params()
set_annotations(HK.KD.annots)
set_annot_samps(c("disease_state","state"))
set_annot_cols(ann_colors)
set_scale.range(c(-1,1))


#glycolysis
reactions <- c("HMR_5029","HMR_4316","HMR_4394","G6PDH2c","HMR_4396","HMR_3944","HMR_4381","HMR_4379","HMR_4375","HMR_4391",
               "HMR_4373","HMR_4368","HMR_4365","HMR_4363","HMR_4358","HMR_4388","HMR_4928","HMR_5998","HMR_4198","HPYRR2x",
               "HMR_8774","HMR_8775","HMR_4930","HMR_4281","LACLt")
clust.samps <- stats::dist(t(abs(bind_matrices[match(reactions,rownames(bind_matrices)),which(HK.KD.annots$state=="No.KD")[c(1:18,30:39)]])), method = "euclidean")
a <- hclust(clust.samps, method = "complete")
AH.samp.order <-  rev(a$labels[a$order])
clust.samps <- dist(t(abs(bind_matrices[match(reactions,rownames(bind_matrices)),which(HK.KD.annots$state=="Healthy")])), method = "euclidean")
a <- hclust(clust.samps, method = "complete")
Healthy.samp.order <-  a$labels[a$order]

## Figure 7A
p = myHeatmapByAnnotation(abs(bind_matrices[match(reactions,rownames(bind_matrices)),c(match(Healthy.samp.order,colnames(bind_matrices)),match(AH.samp.order,colnames(bind_matrices)),
                                                                                       unlist(lapply(39*seq(0,1),function(x)
                                                                                         x+match(AH.samp.order,HK.KD.annots[-which(HK.KD.annots$state=="PKM.KD"),1]))))]),
                          main = "Glycolysis Reactions", gap.width = 3,fontsize.row = 14,
                          groupings = c("state"), groupings.gaps = c(1),
                          scale.rows = "zscore",show.rownames = T,
                          gaps.row = F, clust.rows = F, clust.cols = F,show.colnames = T)



# HKDC1 knockdown
scaling.factor.HKDC1 <- max(metgenes_patient[match("ENSG00000156510",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000156510",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

scaling.factor.HKDC1.explant <- max(metgenes_patient[match("ENSG00000156510",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000156510",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.HKDC1 <- metgenes_patient

metgenes_patient.HKDC1[match("ENSG00000156510",rownames(metgenes_patient.HKDC1)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.HKDC1[match("ENSG00000156510",rownames(metgenes_patient.HKDC1)),which(disease_state%in%"severe.AH")] * scaling.factor.HKDC1

metgenes_patient.HKDC1[match("ENSG00000156510",rownames(metgenes_patient.HKDC1)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.HKDC1[match("ENSG00000156510",rownames(metgenes_patient.HKDC1)),which(disease_state%in%"explant.AH")] * scaling.factor.HKDC1.explant


metgenes_patient.HKDC1.50pctKD <- metgenes_patient
metgenes_patient.HKDC1.50pctKD[match("ENSG00000156510",rownames(metgenes_patient.HKDC1.50pctKD)),which(disease_state%in%c("nonsevere.AH","severe.AH","explant.AH"))] <- 
  metgenes_patient.HKDC1.50pctKD[match("ENSG00000156510",rownames(metgenes_patient.HKDC1.50pctKD)),which(disease_state%in%c("nonsevere.AH","severe.AH","explant.AH"))]/2


combined_HKDC1 <- as.data.frame(t(c(metgenes_patient[match("ENSG00000156510",rownames(metgenes_patient)),c(1:18,30:39,80:89)],
                                    metgenes_patient.HKDC1.50pctKD[match("ENSG00000156510",rownames(metgenes_patient.HKDC1.50pctKD)),c(1:18,30:39)],
                                    metgenes_patient.HKDC1[match("ENSG00000156510",rownames(metgenes_patient.HKDC1)),c(1:18,30:39)])))
rownames(combined_HKDC1) <- "HKDC1"
colnames(combined_HKDC1) <- c(rownames(annots)[c(1:18,30:39,80:89)],paste0(rownames(annots)[c(1:18,30:39)],"_HKDC1.50pct"),
                              paste0(rownames(annots)[c(1:18,30:39)],"_HKDC1.KD"))

HKDC1.annots <- data.frame(colnames(combined_HKDC1))
rownames(HKDC1.annots) <- HKDC1.annots[,1]
HKDC1.annots$state <- "No.KD"
HKDC1.annots$state[grep("I",rownames(HKDC1.annots))] <- "Healthy"
HKDC1.annots$state[grep("_HKDC1.50pct",rownames(HKDC1.annots))] <- "HKDC1.50pct"
HKDC1.annots$state[grep("_HKDC1.KD",rownames(HKDC1.annots))] <- "HKDC1.KD"
HKDC1.annots$state <- factor(HKDC1.annots$state, levels = c("Healthy","No.KD","HKDC1.50pct","HKDC1.KD"))
HKDC1.annots$disease_state <- annots$disease_state[match(strsplit2(HKDC1.annots[,1],split = "_")[,1],rownames(annots))]

initiate_params()
set_annotations(HKDC1.annots)
set_annot_cols(ann_colors)
set_annot_samps(c("state"))
set_scale.range(c(-1,1))

## Figure 7B
beeswarmGenes(combined_HKDC1,list=rownames(combined_HKDC1),color.by = "disease_state",groupby.x = "state",facet.wrap = T) +
  geom_boxplot(alpha=0)+
  theme(axis.text.x = element_text(color = "black", size = 14,angle=45,hjust=1),
        axis.text = element_text(color = "black", size = 14),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.title = element_text(size = 16)
  ) + 
  geom_signif(comparisons = list(c("Healthy", "No.KD"),
                                 c("No.KD", "HKDC1.50pct"),
                                 c("HKDC1.50pct", "HKDC1.KD")),
              map_signif_level=T, textsize = 8, step_increase=-0.1,vjust=0.5,margin_top=0.2)


## identifying variability in patient response to Hkdc1 KD
reactions2 <- c("HMR_5029","HMR_4394","HMR_4381","HMR_4379","HMR_4375","HMR_4391",
                "HMR_4373","HMR_4368","HMR_4365","HMR_4363","HMR_4358")
to.plot <- abs(bind_matrices[match(reactions2,rownames(bind_matrices)),which(HK.KD.annots$state=="HKDC1.KD" & HK.KD.annots$disease_state%in%c("explant.AH","severe.AH"))])
to.plot.zscore <- colMeans(as.data.frame(t(scale(t(to.plot)))))
to.plot <- as.data.frame(rbind(to.plot,to.plot.zscore))
rownames(to.plot)[nrow(to.plot)] <- "score"
p = myHeatmap(to.plot[-nrow(to.plot),],
              main = "Glycolysis Reactions", method = "euclidean", col.groups = 6,
              scale.rows = "zscore",show.rownames = T,
              gaps.row = F, clust.rows = F, clust.cols = T,show.colnames = T)

HKDC1.KD.clusts <- extractClusters(to.plot[-nrow(to.plot),],
                                   heatmap = p, nclusters = 6,to.extract ="samples")


set_annotations(HK.KD.annots)
set_annot_samps(c("disease_state","state"))

clust.A <- hclust(stats::dist(t(to.plot[-nrow(to.plot),which(HKDC1.KD.clusts$Sample.Groups=="C")])), method = "complete")
clust.E <- hclust(stats::dist(t(to.plot[-nrow(to.plot),which(HKDC1.KD.clusts$Sample.Groups=="D")])), method = "complete")
clust.mid <- hclust(stats::dist(t(to.plot[-nrow(to.plot),which(!(HKDC1.KD.clusts$Sample.Groups%in%c("C","D")))])), method = "complete")

## Figure 7C
p = myHeatmap(to.plot[-nrow(to.plot),rev(match(c(rev(clust.A$labels[clust.A$order]),rev(clust.mid$labels[clust.mid$order]),clust.E$labels[clust.E$order]),
                                               colnames(to.plot)))],
              main = "Glycolysis Reactions", method = "euclidean", col.groups = 5,
              scale.rows = "zscore",show.rownames = T,gaps.col = c(7,22),#100,118,136,154),
              gaps.row = F, clust.rows = F, clust.cols = F,show.colnames = T,fontsize.row = 12,fontsize.col = 12)

# update annotations
HK.KD.annots$severeAH.annotsflux.groups <- NA
HK.KD.annots$severeAH.annotsflux.groups[grep(paste0(strsplit2(rownames(HKDC1.KD.clusts)[which(HKDC1.KD.clusts$Sample.Groups=="D")],split = "_")[,1],collapse = "|"),
                                             rownames(HK.KD.annots))] <- "Group1"
HK.KD.annots$severeAH.annotsflux.groups[grep(paste0(strsplit2(rownames(HKDC1.KD.clusts)[which(HKDC1.KD.clusts$Sample.Groups=="C")],split = "_")[,1],collapse = "|"),
                                             rownames(HK.KD.annots))] <- "Group3"
HK.KD.annots$severeAH.annotsflux.groups[which(HK.KD.annots$disease_state=="nonsevere.AH")] <- "nonsevere.AH"
HK.KD.annots$severeAH.annotsflux.groups[grep("I",rownames(HK.KD.annots))] <- "Healthy"
HK.KD.annots$severeAH.annotsflux.groups[which(is.na(HK.KD.annots$severeAH.annotsflux.groups))] <- "Group2"
HK.KD.annots$severeAH.annotsflux.groups <- factor(HK.KD.annots$severeAH.annotsflux.groups, levels = c("Healthy","nonsevere.AH","Group1","Group2","Group3","explant.AH"))

ann_colors[["severeAH.annotsflux.groups"]] <- c("Group3"="red","Group1"="blue", "Healthy"='#f280ae',"explant.AH" = "black","nonsevere.AH"="gray","Group2"="gray")

set_annotations(HK.KD.annots)
set_annot_samps("state")
set_annot_cols(ann_colors)


# PCA with pseudotime analysis - Hkdc1 KD
pca.vals <-myPCA(t(scale(t(abs(bind_matrices[match(reactions2,rownames(bind_matrices)),which(HK.KD.annots$state=="HKDC1.KD" & HK.KD.annots$disease_state%in%c("explant.AH","severe.AH"))])))),return.ggplot.input = T)
pca.dim <- pca.vals[["input_data"]][,1:2]

sds <- slingshot(pca.dim[-grep("D37|D34|B15",rownames(pca.dim)),1:2],
                      clusterLabels = rep(1,nrow(pca.dim)-3))


# colors <- colorRampPalette(brewer.pal(9, 'Greys'))(100)
# plotcol <- ann_colors$severeAH.annotsflux.groups[match(HK.KD.annots$severeAH.annotsflux.groups[match(rownames(pca.dim),rownames(HK.KD.annots))],names(ann_colors$severeAH.annotsflux.groups))]
# plot(pca.dim[,1:2], col = "black",bg = plotcol, pch=21, asp = 1,cex=3)
# lines(SlingshotDataSet(sds), lwd=2, col='black')

# Figure 7D
ggplot(pca.dim, aes(x = -PC1, y = PC2)) +
  geom_line(data = as.data.frame(SlingshotDataSet(sds)@curves[[1]]$s), aes(x = -PC1, y = PC2), color = "black", size = 1.25) +
  geom_point(color = "black", fill = ann_colors$severeAH.annotsflux.groups[match(HK.KD.annots$severeAH.annotsflux.groups[match(rownames(pca.dim),rownames(HK.KD.annots))],names(ann_colors$severeAH.annotsflux.groups))], shape = 21, size = 8) +  # Plot points
  geom_text(aes(label=ifelse(PC1< -2 & PC2 < - 0.25,as.character(rownames(annots)[c(1:18,30:39)]),'')),hjust=-0.35,vjust=-0.25,size=5) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16))


# PCA with pseudotime analysis - no KD
pca.flux.dims <- myPCA(t(scale(t(abs(FBA_mat[match(reactions2,rownames(FBA_mat)),c(1:18,30:39)])))),return.ggplot.input = T)
sds.noKD <- slingshot(pca.flux.dims$input_data[,1:2],
                      clusterLabels = rep(1,nrow(pca.flux.dims$input_data)))

# colors <- colorRampPalette(brewer.pal(9, 'Greys'))(100)
# plotcol <- ann_colors$severeAH.annotsflux.groups[match(HK.KD.annots$severeAH.annotsflux.groups[match(rownames(pca.flux.dims$input_data),rownames(HK.KD.annots))],names(ann_colors$severeAH.annotsflux.groups))]
# plot(pca.flux.dims$input_data[,1:2], col = "black",bg = plotcol, pch=21, asp = 1,cex=3)
# lines(SlingshotDataSet(sds.noKD), lwd=2, col='black')


## Figure 7E
ggplot(pca.flux.dims$input_data[,1:2], aes(x = PC1, y = PC2)) +
  geom_line(data = as.data.frame(SlingshotDataSet(sds.noKD)@curves[[1]]$s), aes(x = PC1, y = PC2), color = "black", size = 1.25) +
  geom_point(color = "black", fill = ann_colors$severeAH.annotsflux.groups[match(HK.KD.annots$severeAH.annotsflux.groups[match(rownames(pca.flux.dims$input_data),rownames(HK.KD.annots))],names(ann_colors$severeAH.annotsflux.groups))], shape = 21, size = 8) +  # Plot points
  geom_text(aes(label=ifelse(rownames(annots)[c(1:18,30:39)]%in%c("D37","D34","B15"),as.character(rownames(annots)[c(1:18,30:39)]),'')),hjust=-0.15,vjust=1.75,size=5) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)) + xlim(c(-5,7))

# legend
plot.new()
unique_groups <- c("Group1","Group2","Group3")
legend_colors <- ann_colors$severeAH.annotsflux.groups[match(unique_groups,names(ann_colors$severeAH.annotsflux.groups))]
legend("center", legend = unique_groups, pt.bg = legend_colors, pch = 21, pt.cex = 4, title = "Severe AH Groups", bty = "n", cex = 1.2,y.intersp = 1.5)


# correlation between Hkdc1 KD pseudotime vs. No KD
to.plot <- data.frame(pseudo.KD=-as.numeric(slingPseudotime(sds)), 
                      pseudo.NoKD= - (as.numeric(slingPseudotime(sds.noKD))[-grep("D37|D34|B15",rownames(pca.dim))] -
                                        max(as.numeric(slingPseudotime(sds.noKD))[-grep("D37|D34|B15",rownames(pca.dim))])),
                      group = as.character(HK.KD.annots$severeAH.annotsflux.groups[c(1:18,30:39)][-grep("D37|D34|B15",rownames(pca.dim))]))
rownames(to.plot) <- rownames(slingPseudotime(sds))

cor.test(as.numeric(to.plot$pseudo.KD),
         as.numeric(to.plot$pseudo.NoKD),
         method="spearman")

## Figure 7F
ggplot(to.plot, aes(x = pseudo.NoKD, y = pseudo.KD)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") +  # Add dotted black line , data = subset(to.plot, flux.pseudo > 130)
  geom_point(color = "black", fill =ann_colors$severeAH.annotsflux.groups[match(HK.KD.annots$severeAH.annotsflux.groups[match(rownames(to.plot),rownames(HK.KD.annots))],names(ann_colors$severeAH.annotsflux.groups))], shape = 21, size = 8) +  # Plot points
  # geom_text(aes(label=ifelse(PC1< -2 & PC2 < - 0.25,as.character(rownames(annots)[c(1:18,30:39)]),'')),hjust=-0.35,vjust=-0.25,size=5) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)) #+ ylim(c(0,9))


# update annotations
annots$severeAH.annotsflux.groups <- as.character(annots$disease_state)
annots$severeAH.annotsflux.groups[c(1:18,30:39)] <- as.character(HK.KD.annots$severeAH.annotsflux.groups[match(rownames(annots)[c(1:18,30:39)],rownames(HK.KD.annots))])
annots$severeAH.annotsflux.groups[which(is.na(annots$severeAH.annotsflux.groups))] <- "Group2"
annots$severeAH.annotsflux.groups <- factor(annots$severeAH.annotsflux.groups,levels=c("Group1","Group2","Group3"))

# DESeqs between patient groups following Hkdc1 KD
dds <- DESeqDataSetFromMatrix(countData = rawct[,c(1:18,30:39)],
                              colData = annots[c(1:18,30:39),],
                              design = ~ severeAH.annotsflux.groups)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)


severeAH.groups <- as.data.frame(results(dds, contrast=c("severeAH.annotsflux.groups","Group3","Group1")))

severeAH.groups <- severeAH.groups[which(severeAH.groups$padj<0.05 & severeAH.groups$log2FoldChange>0),]

severeAH.groups$GeneName <- genenames$hgnc_symbol[match(rownames(severeAH.groups),genenames$ensembl_gene_id)]
glycolysis.genes <- unique(as.vector(strsplit2(ihuman.geneName$GeneAssociationName[which(ihuman.geneName$SUBSYSTEM=="Glycolysis / Gluconeogenesis")],
                                               split="\\(|\\)| or | and ")))
glycolysis.genes[which(glycolysis.genes%in%severeAH.groups$GeneName)]

# Group1 vs. Group3 signif genes
#"PKM" "HK1" "HK2"

## Figure 8A
set_annotations(annots)
beeswarmGenes(metgenes_patient_geneNames[match(c("HK1","HK2","PKM"),rownames(metgenes_patient_geneNames)),c(1:18,30:39)],
              rownames(metgenes_patient_geneNames),color.by = "severeAH.annotsflux.groups",
              groupby.x = "severeAH.annotsflux.groups",facet.wrap = T,ncols = 3) + 
  geom_boxplot(alpha=0)+
  # geom_text_repel(aes(label = rownames(annots)[c(1:18,30:39)]), 
  #                 position = ggbeeswarm::position_quasirandom(),
  #                 size = 3) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 18)) +
  geom_signif(comparisons = list(c("Group1","Group3")),
              map_signif_level=T, textsize = 8, step_increase=0.1,vjust=0.5,margin_top = 0.08)



### gene knockdowns

metgenes_patient.KD <- metgenes_patient

# HKDC1 knockdown
scaling.factor.HKDC1 <- max(metgenes_patient[match("ENSG00000156510",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000156510",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000156510",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000156510",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.HKDC1

metgenes_patient.KD[match("ENSG00000156510",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000156510",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.HKDC1.explant

# HK1 knockdown 
scaling.factor.HK1 <- max(metgenes_patient[match("ENSG00000156515",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000156515",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000156515",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000156515",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.HK1

scaling.factor.HK1.explant <- max(metgenes_patient[match("ENSG00000156515",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000156515",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000156515",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000156515",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.HK1.explant

# HK2 knockdown 
scaling.factor.HK2 <- max(metgenes_patient[match("ENSG00000159399",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000159399",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000159399",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000159399",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.HK2

scaling.factor.HK2.explant <- max(metgenes_patient[match("ENSG00000159399",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000159399",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000159399",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000159399",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.HK2.explant

# GLUT1 (Slc2a1) knockdown
scaling.factor.glut1 <- max(metgenes_patient[match("ENSG00000117394",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000117394",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000117394",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000117394",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.glut1

scaling.factor.glut1.explant <- max(metgenes_patient[match("ENSG00000117394",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000117394",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000117394",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000117394",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.glut1.explant


# GLUT3 (Slc2a3) knockdown
scaling.factor.glut3 <- max(metgenes_patient[match("ENSG00000059804",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000059804",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000059804",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000059804",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.glut3

scaling.factor.glut3.explant <- max(metgenes_patient[match("ENSG00000059804",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000059804",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000059804",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000059804",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.glut3.explant


# PKM knockdown
scaling.factor.PKM <- max(metgenes_patient[match("ENSG00000067225",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000067225",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000067225",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000067225",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.PKM

scaling.factor.PKM.explant <- max(metgenes_patient[match("ENSG00000067225",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000067225",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000067225",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000067225",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.PKM.explant


# ENO1 knockdown
scaling.factor.ENO1 <- max(metgenes_patient[match("ENSG00000074800",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000074800",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000074800",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000074800",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.ENO1

scaling.factor.ENO1.explant <- max(metgenes_patient[match("ENSG00000074800",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000074800",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000074800",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000074800",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.ENO1.explant

# GAPDH knockdown
scaling.factor.GAPDH <- max(metgenes_patient[match("ENSG00000111640",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000111640",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000111640",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000111640",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.GAPDH

scaling.factor.GAPDH.explant <- max(metgenes_patient[match("ENSG00000111640",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000111640",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000111640",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000111640",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.GAPDH.explant

# LDHA knockdown
scaling.factor.LDHA <- max(metgenes_patient[match("ENSG00000134333",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000134333",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000134333",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000134333",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.LDHA

scaling.factor.LDHA.explant <- max(metgenes_patient[match("ENSG00000134333",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000134333",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000134333",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000134333",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.LDHA.explant

# LDHB knockdown
scaling.factor.LDHB <- max(metgenes_patient[match("ENSG00000111716",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000111716",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000111716",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000111716",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.LDHB

scaling.factor.LDHB.explant <- max(metgenes_patient[match("ENSG00000111716",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000111716",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000111716",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000111716",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.LDHB.explant


# PFKP knockdown
scaling.factor.PFKP <- max(metgenes_patient[match("ENSG00000067057",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000067057",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000067057",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000067057",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.PFKP

scaling.factor.PFKP.explant <- max(metgenes_patient[match("ENSG00000067057",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000067057",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000067057",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000067057",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.PFKP.explant


# PGK1 knockdown
scaling.factor.PGK1 <- max(metgenes_patient[match("ENSG00000102144",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000102144",rownames(metgenes_patient)),which(disease_state%in%"severe.AH")])

metgenes_patient.KD[match("ENSG00000102144",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] <- 
  metgenes_patient.KD[match("ENSG00000102144",rownames(metgenes_patient.KD)),which(disease_state%in%"severe.AH")] * scaling.factor.PGK1

scaling.factor.PGK1.explant <- max(metgenes_patient[match("ENSG00000102144",rownames(metgenes_patient)),which(disease_state%in%"healthy.control")])/
  max(metgenes_patient[match("ENSG00000102144",rownames(metgenes_patient)),which(disease_state%in%"explant.AH")])

metgenes_patient.KD[match("ENSG00000102144",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] <- 
  metgenes_patient.KD[match("ENSG00000102144",rownames(metgenes_patient.KD)),which(disease_state%in%"explant.AH")] * scaling.factor.PGK1.explant


metgenes_patient.KD_genenames <- metgenes_patient.KD
rownames(metgenes_patient.KD_genenames) <- make.names(genenames$hgnc_symbol[match(rownames(metgenes_patient.KD_genenames),genenames$ensembl_gene_id)],unique=T)


all.KD.geneExp <- rbind(c(metgenes_patient_geneNames[match("HKDC1",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("HKDC1",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("HK1",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("HK1",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("HK2",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("HK2",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("SLC2A1",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("SLC2A1",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("SLC2A3",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("SLC2A3",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("PKM",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("PKM",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("ENO1",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("ENO1",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("GAPDH",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("GAPDH",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("LDHA",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("LDHA",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("LDHB",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("LDHB",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("PFKP",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("PFKP",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]),
                        c(metgenes_patient_geneNames[match("PGK1",rownames(metgenes_patient_geneNames)),c(1:18,30:39,80:89)], metgenes_patient.KD_genenames[match("PGK1",rownames(metgenes_patient_geneNames)),c(1:18,30:39)]))


rownames(all.KD.geneExp) <- c("HKDC1","HK1","HK2","SLC2A1","SLC2A3","PKM","ENO1","GAPDH","LDHA","LDHB","PFKP","PGK1")
colnames(all.KD.geneExp) <- c(rownames(annots)[c(1:18,30:39,80:89)],paste0(rownames(annots)[c(1:18,30:39)],"KD"))

all.KD.annots <- data.frame(colnames(all.KD.geneExp))
rownames(all.KD.annots) <- all.KD.annots[,1]
all.KD.annots$state <- "No.KD"
all.KD.annots$state[grep("I",rownames(all.KD.annots))] <- "Healthy"
all.KD.annots$state[grep("KD",rownames(all.KD.annots))] <- "KD"
all.KD.annots$state <- factor(all.KD.annots$state, levels = c("Healthy","No.KD","KD"))
all.KD.annots$disease_state <- annots$disease_state[match(c(rownames(annots)[c(1:18,30:39,80:89)],rownames(annots)[c(1:18,30:39)]),rownames(annots))]


initiate_params()
set_annotations(all.KD.annots)
set_annot_cols(ann_colors)
set_annot_samps(c("state"))
set_scale.range(c(-1,1))


## Figure 8B
beeswarmGenes(all.KD.geneExp,list=c("HK1","HK2","PKM"),color.by = "disease_state",groupby.x = "state",facet.wrap = T, ncols=5) +
  geom_boxplot(alpha=0)+ geom_signif(comparisons = list(c("Healthy","No.KD"),
                                                        c("No.KD", "KD")), 
                                     map_signif_level=TRUE,
                                     textsize = 10, step_increase=0.1,vjust=0.5) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 18))


## Supplemental Figure 8A
beeswarmGenes(all.KD.geneExp,list=c("PFKP","PGK1","SLC2A1","SLC2A3","ENO1","GAPDH","LDHA","LDHB"),color.by = "disease_state",groupby.x = "state",facet.wrap = T, ncols=4) +
  geom_boxplot(alpha=0)+ geom_signif(comparisons = list(c("Healthy","No.KD"),
                                                        c("No.KD", "KD")), 
                                     map_signif_level=TRUE,
                                     textsize = 10, step_increase=0.1,vjust=0.5) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 18))


## Reading in the predicted metabolic fluxomes for each in silico KD experiment
HKDC1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HKDC1_KD_4June2024/Fluxomes/PheFlux_Results_HKDC1_KD.txt', sep = '\t')
FBA_mat_HKDC1_KD <- FBA_mat
FBA_mat_HKDC1_KD <- as.data.frame(cbind(HKDC1_KD,FBA_mat_HKDC1_KD[,-match(colnames(HKDC1_KD),colnames(FBA_mat))]))

HKDC1_50pctKD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HKDC1_50pctKD_try2/Fluxomes/PheFlux_Results_HKDC1_50pctKD_try2.txt', sep = '\t')
FBA_mat_HKDC1_50pctKD <- FBA_mat
FBA_mat_HKDC1_50pctKD <- as.data.frame(cbind(HKDC1_50pctKD,FBA_mat_HKDC1_50pctKD[,-match(colnames(HKDC1_50pctKD),colnames(FBA_mat))]))

HK1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HK1_KD_4June2024/Fluxomes/PheFlux_Results_HK1_KD.txt', sep = '\t')
FBA_mat_HK1_KD <- FBA_mat
FBA_mat_HK1_KD <- as.data.frame(cbind(HK1_KD,FBA_mat_HK1_KD[,-match(colnames(HK1_KD),colnames(FBA_mat))]))

HK2_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HK2_KD/Fluxomes/PheFlux_Results_HK2_KD.txt', sep = '\t')
FBA_mat_HK2_KD <- FBA_mat
FBA_mat_HK2_KD <- as.data.frame(cbind(HK2_KD,FBA_mat_HK2_KD[,-match(colnames(HK2_KD),colnames(FBA_mat))]))

HK1_HKDC1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HK1_HKDC1_KD_4June2024/Fluxomes/PheFlux_Results_HK1_HKDC1_KD.txt', sep = '\t')
FBA_mat_HK1_HKDC1_KD <- FBA_mat
FBA_mat_HK1_HKDC1_KD <- as.data.frame(cbind(HK1_HKDC1_KD,FBA_mat_HK1_HKDC1_KD[,-match(colnames(HK1_HKDC1_KD),colnames(FBA_mat_HK1_HKDC1_KD))]))

HK2_HKDC1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HK2_HKDC1/Fluxomes/PheFlux_Results_HK2_HKDC1_KD.txt', sep = '\t')
FBA_mat_HK2_HKDC1_KD <- FBA_mat
FBA_mat_HK2_HKDC1_KD <- as.data.frame(cbind(HK2_HKDC1_KD,FBA_mat_HK2_HKDC1_KD[,-match(colnames(HK2_HKDC1_KD),colnames(FBA_mat_HK2_HKDC1_KD))]))

GLUT1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/GLUT1_KD_10June2024/Fluxomes/PheFlux_Results_GLUT1_KD.txt', sep = '\t')
FBA_mat_GLUT1_KD  <- FBA_mat
FBA_mat_GLUT1_KD  <- as.data.frame(cbind(GLUT1_KD,FBA_mat_GLUT1_KD[,-match(colnames(GLUT1_KD),colnames(FBA_mat_GLUT1_KD))]))

GLUT3_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/GLUT3_KD/Fluxomes/PheFlux_Results_GLUT3_KD.txt', sep = '\t')
FBA_mat_GLUT3_KD  <- FBA_mat
FBA_mat_GLUT3_KD  <- as.data.frame(cbind(GLUT3_KD,FBA_mat_GLUT3_KD[,-match(colnames(GLUT3_KD),colnames(FBA_mat_GLUT3_KD))]))

PKM_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/PKM_KD/Fluxomes/PheFlux_Results_PKM_KD.txt', sep = '\t')
FBA_mat_PKM_KD  <- FBA_mat
FBA_mat_PKM_KD  <- as.data.frame(cbind(PKM_KD,FBA_mat_PKM_KD[,-match(colnames(PKM_KD),colnames(FBA_mat_PKM_KD))]))

ENO1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/ENO1_KD/Fluxomes/PheFlux_Results_ENO1_KD.txt', sep = '\t')
FBA_mat_ENO1_KD  <- FBA_mat
FBA_mat_ENO1_KD  <- as.data.frame(cbind(ENO1_KD,FBA_mat_ENO1_KD[,-match(colnames(ENO1_KD),colnames(FBA_mat_ENO1_KD))]))

GAPDH_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/GAPDH_KD/Fluxomes/PheFlux_Results_GAPDH_KD.txt', sep = '\t')
FBA_mat_GAPDH_KD  <- FBA_mat
FBA_mat_GAPDH_KD  <- as.data.frame(cbind(GAPDH_KD,FBA_mat_GAPDH_KD[,-match(colnames(GAPDH_KD),colnames(FBA_mat_GAPDH_KD))]))

LDHA_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/LDHA_KD/Fluxomes/PheFlux_Results_LDHA_KD.txt', sep = '\t')
FBA_mat_LDHA_KD  <- FBA_mat
FBA_mat_LDHA_KD  <- as.data.frame(cbind(LDHA_KD,FBA_mat_LDHA_KD[,-match(colnames(LDHA_KD),colnames(FBA_mat_LDHA_KD))]))

LDHB_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/LDHB_KD/Fluxomes/PheFlux_Results_LDHB_KD.txt', sep = '\t')
FBA_mat_LDHB_KD  <- FBA_mat
FBA_mat_LDHB_KD  <- as.data.frame(cbind(LDHB_KD,FBA_mat_LDHB_KD[,-match(colnames(LDHB_KD),colnames(FBA_mat_LDHB_KD))]))

HKDC1_PKM_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/HKDC1_PKM_KD/Fluxomes/PheFlux_Results_HKDC1_PKM_KD.txt', sep = '\t')
FBA_mat_HKDC1_PKM_KD  <- FBA_mat
FBA_mat_HKDC1_PKM_KD  <- as.data.frame(cbind(HKDC1_PKM_KD,FBA_mat_HKDC1_PKM_KD[,-match(colnames(HKDC1_PKM_KD),colnames(FBA_mat_HKDC1_PKM_KD))]))

# predict.new <- predict(svm.linear.cv,newdata = t(abs(FBA_mat_HK2_HKDC1_PKM_KD[match(c("HMR_5029","HMR_4394","HMR_4381"),rownames(FBA_mat_HK2_HKDC1_PKM_KD)),])))
# predict.new <- as.data.frame(t(predict.new))
# rownames(predict.new) = "predict.new"

PFKP_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/PFKP_KD/Fluxomes/PheFlux_Results_PFKP_KD.txt', sep = '\t')
FBA_mat_PFKP_KD  <- FBA_mat
FBA_mat_PFKP_KD  <- as.data.frame(cbind(PFKP_KD,FBA_mat_PFKP_KD[,-match(colnames(PFKP_KD),colnames(FBA_mat_PFKP_KD))]))
# 
PGK1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/PGK1_KD/Fluxomes/PheFlux_Results_PGK1_KD.txt', sep = '\t')
FBA_mat_PGK1_KD  <- FBA_mat
FBA_mat_PGK1_KD  <- as.data.frame(cbind(PGK1_KD,FBA_mat_PGK1_KD[,-match(colnames(PGK1_KD),colnames(FBA_mat_PGK1_KD))]))


bind_matrices <- cbind(FBA_mat_HKDC1_KD[,grep("A|B|D",colnames(FBA_mat_HKDC1_KD))],FBA_mat_HK1_KD[,grep("A|B|D",colnames(FBA_mat_HK1_KD))],FBA_mat_HK2_KD[,grep("A|B|D",colnames(FBA_mat_HK2_KD))],FBA_mat_HK1_HKDC1_KD[,grep("A|B|D",colnames(FBA_mat_HK1_HKDC1_KD))],FBA_mat_HK2_HKDC1_KD[,grep("A|B|D",colnames(FBA_mat_HK2_HKDC1_KD))],
                       FBA_mat_PKM_KD[,grep("A|B|D",colnames(FBA_mat_PKM_KD))],FBA_mat_HKDC1_PKM_KD[,grep("A|B|D",colnames(FBA_mat_HKDC1_PKM_KD))],FBA_mat_GLUT1_KD[,grep("A|B|D",colnames(FBA_mat_GLUT1_KD))],FBA_mat_GLUT3_KD[,grep("A|B|D",colnames(FBA_mat_GLUT3_KD))],
                       FBA_mat_ENO1_KD[,grep("A|B|D",colnames(FBA_mat_ENO1_KD))],FBA_mat_GAPDH_KD[,grep("A|B|D",colnames(FBA_mat_GAPDH_KD))],
                       FBA_mat_LDHA_KD[,grep("A|B|D",colnames(FBA_mat_LDHA_KD))], FBA_mat_LDHB_KD[,grep("A|B|D",colnames(FBA_mat_LDHB_KD))], FBA_mat_PFKP_KD[,grep("A|B|D",colnames(FBA_mat_PFKP_KD))],
                       FBA_mat_PGK1_KD[,grep("A|B|D",colnames(FBA_mat_PGK1_KD))],FBA_mat_HKDC1_50pctKD[,grep("A|B|D",colnames(FBA_mat_PGK1_KD))], FBA_mat[,c(1:18,30:39,80:89)])

colnames(bind_matrices)[(seq(1,476,28)[1]:(seq(1,476,28)[2]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_HKDC1.KD")
colnames(bind_matrices)[(seq(1,476,28)[2]:(seq(1,476,28)[3]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_HK1.KD")
colnames(bind_matrices)[(seq(1,476,28)[3]:(seq(1,476,28)[4]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_HK2.KD")
colnames(bind_matrices)[(seq(1,476,28)[4]:(seq(1,476,28)[5]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_HKDC1.HK1.KD")
colnames(bind_matrices)[(seq(1,476,28)[5]:(seq(1,476,28)[6]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_HKDC1.HK2.KD")
colnames(bind_matrices)[(seq(1,476,28)[6]:(seq(1,476,28)[7]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_PKM.KD")
colnames(bind_matrices)[(seq(1,476,28)[7]:(seq(1,476,28)[8]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_HKDC1.PKM.KD")
colnames(bind_matrices)[(seq(1,476,28)[8]:(seq(1,476,28)[9]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_GLUT1.KD")
colnames(bind_matrices)[(seq(1,476,28)[9]:(seq(1,476,28)[10]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_GLUT3.KD")
colnames(bind_matrices)[(seq(1,476,28)[10]:(seq(1,476,28)[11]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_ENO1.KD")
colnames(bind_matrices)[(seq(1,476,28)[11]:(seq(1,476,28)[12]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_GAPDH.KD")
colnames(bind_matrices)[(seq(1,476,28)[12]:(seq(1,476,28)[13]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_LDHA.KD")
colnames(bind_matrices)[(seq(1,476,28)[13]:(seq(1,476,28)[14]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_LDHB.KD")
colnames(bind_matrices)[(seq(1,476,28)[14]:(seq(1,476,28)[15]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_PFKP.KD")
colnames(bind_matrices)[(seq(1,476,28)[15]:(seq(1,476,28)[16]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_PGK1.KD")
colnames(bind_matrices)[(seq(1,476,28)[16]:(seq(1,476,28)[17]-1))] <- paste0(colnames(FBA_mat[,c(1:18,30:39)]),"_HKDC1.50pct.KD")


HK.KD.annots <- as.data.frame(colnames(bind_matrices))
rownames(HK.KD.annots) <- HK.KD.annots[,1]
HK.KD.annots$state <- "No.KD"
HK.KD.annots$state[grep("_HKDC1.KD",rownames(HK.KD.annots))] <- "HKDC1.KD"
HK.KD.annots$state[grep("_HK1.KD",rownames(HK.KD.annots))] <- "HK1.KD"
HK.KD.annots$state[grep("_HK2.KD",rownames(HK.KD.annots))] <- "HK2.KD"
HK.KD.annots$state[grep("_HKDC1.HK1.KD",rownames(HK.KD.annots))] <- "HKDC1.HK1.KD"
HK.KD.annots$state[grep("_HKDC1.HK2.KD",rownames(HK.KD.annots))] <- "HKDC1.HK2.KD"
HK.KD.annots$state[grep("_GLUT1.KD",rownames(HK.KD.annots))] <- "GLUT1.KD"
HK.KD.annots$state[grep("_PKM.KD",rownames(HK.KD.annots))] <- "PKM.KD"
HK.KD.annots$state[grep("_ENO1.KD",rownames(HK.KD.annots))] <- "ENO1.KD"
HK.KD.annots$state[grep("_GAPDH.KD",rownames(HK.KD.annots))] <- "GAPDH.KD"
HK.KD.annots$state[grep("_LDHA.KD",rownames(HK.KD.annots))] <- "LDHA.KD"
HK.KD.annots$state[grep("_LDHB.KD",rownames(HK.KD.annots))] <- "LDHB.KD"
HK.KD.annots$state[grep("_GLUT3.KD",rownames(HK.KD.annots))] <- "GLUT3.KD"
HK.KD.annots$state[grep("_HKDC1.PKM.KD",rownames(HK.KD.annots))] <- "HKDC1.PKM.KD"
HK.KD.annots$state[grep("_PFKP.KD",rownames(HK.KD.annots))] <- "PFKP.KD"
HK.KD.annots$state[grep("_PGK1.KD",rownames(HK.KD.annots))] <- "PGK1.KD"
HK.KD.annots$state[grep("_HKDC1.50pct.KD",rownames(HK.KD.annots))] <- "HKDC1.50pct.KD"
HK.KD.annots$state[grep("I",rownames(HK.KD.annots))] <- "Healthy"
HK.KD.annots$state <- factor(HK.KD.annots$state,
                             levels = c("Healthy","No.KD","HKDC1.50pct.KD","HKDC1.KD","HK1.KD","HK2.KD","HKDC1.HK1.KD",
                                        "HKDC1.HK2.KD","PFKP.KD","PGK1.KD",
                                        "PKM.KD","HKDC1.PKM.KD","GLUT1.KD",
                                        "GLUT3.KD","ENO1.KD","GAPDH.KD","LDHA.KD","LDHB.KD"))
HK.KD.annots[,1] <- strsplit2(HK.KD.annots[,1],split="_")[,1]


HK.KD.annots$severeAH.annotsflux.groups <- NA
HK.KD.annots$severeAH.annotsflux.groups[grep(paste0(rownames(annots)[which(annots$severeAH.annotsflux.groups=="Group1")],collapse = "|"),
                                             rownames(HK.KD.annots))] <- "Group1"
HK.KD.annots$severeAH.annotsflux.groups[grep(paste0(rownames(annots)[which(annots$severeAH.annotsflux.groups=="Group3")],collapse = "|"),
                                             rownames(HK.KD.annots))] <- "Group3"
HK.KD.annots$severeAH.annotsflux.groups[grep(paste0(rownames(annots)[which(annots$severeAH.annotsflux.groups=="Group2")],collapse = "|"),
                                             rownames(HK.KD.annots))] <- "Group2"
HK.KD.annots$severeAH.annotsflux.groups[grep("I",rownames(HK.KD.annots))] <- "Healthy"

HK.KD.annots$severeAH.annotsflux.groups <- factor(HK.KD.annots$severeAH.annotsflux.groups, levels = c("Healthy","Group1","Group2","Group3"))

ann_colors[["severeAH.annotsflux.groups"]] <- c("Healthy"='#f280ae',"Group1"="blue","Group2"="gray", "Group3"="red")

ann_colors[["state"]] <- c("Healthy"="#f280ae","No.KD"="gray","HK2.HKDC1.PKM.KD"="black","HKDC1.PKM.KD"="red")


set_annotations(HK.KD.annots)
set_annot_samps(c("severeAH.annotsflux.groups","state"))
set_annot_cols(ann_colors)

ann_colors[['state']] <- NULL
set_annot_cols(ann_colors)

reactions <- c("HMR_5029","HMR_4316","HMR_4394","G6PDH2c","HMR_4396","HMR_3944","HMR_4381","HMR_4379","HMR_4375","HMR_4391",
               "HMR_4373","HMR_4368","HMR_4365","HMR_4363","HMR_4358","HMR_4388","HMR_4928","HMR_5998","HMR_4198","HPYRR2x",
               "HMR_8774","HMR_8775","HMR_4930","HMR_4281","LACLt")
clust.samps <- dist(t(abs(bind_matrices[match(reactions,rownames(bind_matrices)),which(HK.KD.annots$state=="No.KD")])), method = "euclidean")
a <- hclust(clust.samps, method = "complete")
AH.samp.order <-  a$labels[a$order]
clust.samps <- dist(t(abs(bind_matrices[match(reactions,rownames(bind_matrices)),which(HK.KD.annots$state=="Healthy")])), method = "euclidean")
a <- hclust(clust.samps, method = "complete")
Healthy.samp.order <-  a$labels[a$order]

## Figure 8C
p = myHeatmapByAnnotation(abs(bind_matrices[match(reactions2,rownames(bind_matrices)),unlist(lapply(match(c("No.KD","HKDC1.KD","HKDC1.HK1.KD","HKDC1.HK2.KD","HKDC1.PKM.KD"),HK.KD.annots$state)-1,
                                                                                                    function(x) x+match(AH.samp.order,HK.KD.annots[which(HK.KD.annots$state%in%c("No.KD","HKDC1.KD","HKDC1.HK1.KD","HKDC1.HK2.KD","HKDC1.PKM.KD")),1])))]),
                          main = "Glycolysis Reactions", gap.width = 1,fontsize.row = 14,groupings.gaps = c(1,3),
                          scale.rows = "zscore",show.rownames = T, groupings = c("severeAH.annotsflux.groups","state"),
                          gaps.row = F, clust.rows = F, clust.cols = F,show.colnames = F)

## Supplemental Figure 8B
p = myHeatmapByAnnotation(abs(bind_matrices[match(reactions2,rownames(bind_matrices)),unlist(lapply(match(c("GLUT1.KD","GLUT3.KD","ENO1.KD","GAPDH.KD","LDHA.KD","LDHB.KD",
                                                                                                            "PFKP.KD","PGK1.KD"),HK.KD.annots$state)-1,
                                                                                                    function(x) x+match(AH.samp.order,HK.KD.annots[which(HK.KD.annots$state%in%c("HKDC1.KD","HK1.KD","HK2.KD","PKM.KD","GLUT1.KD",
                                                                                                                                                                                 "GLUT3.KD","ENO1.KD","GAPDH.KD","LDHA.KD","LDHB.KD",
                                                                                                                                                                                 "PFKP.KD","PGK1.KD")),1])))]),
                          main = "Glycolysis Reactions", gap.width = 1,fontsize.row = 18,groupings.gaps = c(1,3),
                          scale.rows = "zscore",show.rownames = T, groupings = c("severeAH.annotsflux.groups","state"),
                          gaps.row = F, clust.rows = F, clust.cols = F,show.colnames = F)

############################
## Glucose use reprogramming score

#RFE for identifying flux features to include in the model

reactions <-  c("HMR_5029","HMR_4394","HMR_4381","HMR_4379","HMR_4375","HMR_4391",
                "HMR_4373","HMR_4368","HMR_4365","HMR_4363","HMR_4358")

# training data for disease state progression
classifier.data <- as.data.frame(t(abs(FBA_mat_incl1[match(reactions,rownames(FBA_mat_incl1)),c(1:51,80:89)])))
classifier.data$disease_state <- as.character(annots$disease_state[c(1:51,80:89)])
classifier.data$disease_state[which(classifier.data$disease_state=="explant.AH")] <- "severe.AH"
classifier.data$disease_state <- factor(classifier.data$disease_state,levels=c("healthy.control","early.ASH",
                                                                               "nonsevere.AH","severe.AH"))
classifier.data$disease_state <- as.numeric(classifier.data$disease_state)



HK.KD.annots$disease_state <- unlist(lapply(HK.KD.annots[,1],function(x) as.character(annots$disease_state[match(x,rownames(annots))])))

# test data for Hkdc1 knockdown
predict.data <- as.data.frame(t(abs(bind_matrices[match(reactions,rownames(bind_matrices)),-which(HK.KD.annots$disease_state%in%c("healthy.control") |
                                                                                                    HK.KD.annots$state=="No.KD")])))

predict.data <- predict.data[grep("HKDC1.KD",rownames(predict.data)),]
groups <- as.character(HK.KD.annots$severeAH.annotsflux.groups[c(1:18,30:39)])
predict.data$groups <- groups
predict.data$groups <- factor(predict.data$groups)


# test data for Hkdc1/Pkm double knockdown
test.data2 <- as.data.frame(t(abs(bind_matrices[match(reactions,rownames(bind_matrices)),-which(HK.KD.annots$disease_state%in%c("healthy.control") |
                                                                                                  HK.KD.annots$state=="No.KD")])))

test.data2 <- test.data2[grep("HKDC1.PKM.KD",rownames(test.data2)),]
groups <- as.character(HK.KD.annots$severeAH.annotsflux.groups[c(1:18,30:39)])
test.data2$groups <- groups
test.data2$groups <- factor(test.data2$groups)


library(caret)
control <- rfeControl(functions = rfFuncs, # or other model-specific functions like lmFuncs, svmFuncs
                      method = "loocv", # or any other resampling method
                      verbose = TRUE,
                      returnResamp = "all") # to return resampling results

# function for testing performance of model using rfe
rfe_test_performance <- function(trainData, testData,testData2, target_variable) {
  # Store the results for each iteration
  results <- list()
  
  # Create a grid of features to test
  subsets <- c(1:11)  # Modify this to suit the number of features you're interested in
  
  # Perform RFE with test evaluation
  for (n_features in subsets) {
    # RFE on training data
    rfe_fit <- rfe(trainData[, -which(colnames(trainData) == target_variable)], 
                   trainData[, target_variable], 
                   sizes = n_features,
                   rfeControl = control)
    
    # Get selected features
    selected_features <- predictors(rfe_fit)
    
    if (length(selected_features) > n_features) {
      selected_features <- selected_features[1:n_features]
    }
    
    if (length(selected_features)==1){
      data = as.data.frame(trainData[, selected_features])
      rownames(data) = rownames(trainData)
      colnames(data) = selected_features
      
      
      final_model <- train(data, 
                           trainData[, target_variable],
                           method = "svmLinear")  # Change method to your chosen model
    } else {
      
      # Train the model on selected features
      final_model <- train(trainData[, selected_features], 
                           trainData[, target_variable],
                           method = "svmLinear")  # Change method to your chosen model
    }
    
    train_pred <- predict(final_model, trainData[, selected_features])
    pseudo.cor <- cor.test(as.numeric(train_pred),as.numeric(pseudotime),method = "spearman")$estimate[['rho']]
    tab <- table(round(as.numeric(train_pred)), trainData[, target_variable])
    acc <- round(sum(diag(tab)) / sum(tab),2)
    
    # Evaluate on test set
    test_pred <- predict(final_model, testData[, selected_features])
    test_performance <- wilcox.test(test_pred[which(testData$groups=="Group1")],test_pred[which(testData$groups=="Group3")],alternative = "less")$p.value
    
    
    test_pred2 <- predict(final_model, testData2[, selected_features])
    test_performance2 <- wilcox.test(train_pred[c(1:18,30:39)],test_pred2[which(testData2$groups=="Group1")],alternative = "greater")$p.value
    
    test_performance3 <- wilcox.test(test_pred,test_pred2,alternative = "greater")$p.value
    
    
    # Store performance metrics
    results[[n_features]] <- list(
      selected_features = selected_features,
      accuracy = acc,
      pseudo = pseudo.cor,
      perform1 = test_performance,
      perform2 = test_performance2,
      perform3 = test_performance3
      
    )
  }
  return(results)
}

# Call the function
results <- rfe_test_performance(classifier.data, predict.data,test.data2, "disease_state")

# View results
View(results)

## fluxes resulting in the model with the highest accuracy are "HMR_5029", "HMR_4381","HMR_4394"
reactions <- c("HMR_5029", "HMR_4381","HMR_4394")

## retrain model and use for prediction
classifier.data <- as.data.frame(t(abs(FBA_mat_incl1[match(reactions,rownames(FBA_mat_incl1)),c(1:51,80:89)])))

classifier.data$disease_state <- as.character(annots$disease_state[c(1:51,80:89)])
classifier.data$disease_state[which(classifier.data$disease_state=="explant.AH")] <- "severe.AH"
classifier.data$disease_state <- factor(classifier.data$disease_state,levels=c("healthy.control","early.ASH",
                                                                               "nonsevere.AH","severe.AH"))
classifier.data$disease_state <- as.numeric(classifier.data$disease_state)

## LOOCV
train_control_loocv <- trainControl(method = "LOOCV")
svm.linear.cv = train(disease_state ~ ., classifier.data, method = "svmLinear",trControl = train_control_loocv)


Predict <- predict(svm.linear.cv,newdata = classifier.data[,-nrow(classifier.data)] )
Predict

## plotting scores
all_clinical_params <- all_clinical_params[-33,]
plot.score <- as.data.frame(rbind(t(Predict),all_clinical_params$MELD[c(1:51,80:89)],
                                 all_clinical_params$ABIC_puntos[c(1:51,80:89)],
                                 all_clinical_params$Child_puntos[c(1:51,80:89)],
                                 as.numeric(pseudotime),
                                 all_clinical_params$`AST_U/L`[c(1:51,80:89)],
                                 all_clinical_params$`ALT_U/L`[c(1:51,80:89)]))
rownames(plot.score) <- c("score","MELD","ABIC","Child","Pseudotime","AST","ALT")
colnames(plot.score) <- rownames(annots)[c(1:51,80:89)]

annots$disease_state.group <- annots$disease_state
annots$disease_state.group[which(annots$disease_state.group=="explant.AH")]="severe.AH"

annots$severeAH.annotsflux.groups <-"Group2"
annots$severeAH.annotsflux.groups[which(rownames(annots)%in%rownames(HK.KD.annots))] <- as.character(HK.KD.annots$severeAH.annotsflux.groups[na.omit(match(rownames(annots),rownames(HK.KD.annots)))])
annots$severeAH.annotsflux.groups[-which(annots$severeAH.annotsflux.groups%in%c("Group1","Group3"))] <- "Group2"
annots$severeAH.annotsflux.groups <- factor(annots$severeAH.annotsflux.groups,levels = c("Group1","Group2","Group3"))

annots$severeAH.annotsflux.groups_disease <- as.character(annots$severeAH.annotsflux.groups)
annots$severeAH.annotsflux.groups_disease[-grep("A|B|D",rownames(annots))] <- as.character(annots$disease_state[-grep("A|B|D",rownames(annots))])
annots$severeAH.annotsflux.groups_disease <- factor(annots$severeAH.annotsflux.groups_disease,levels = c("NASH","comp.cirrhosis","HCV","healthy.control","early.ASH","nonsevere.AH",
                                                                                                         "Group1","Group2","Group3"))
## Figure 9A
set_annotations(annots)
beeswarmGenes(plot.score,c("score"),groupby.x = "severeAH.annotsflux.groups_disease",
              color.by = "disease_state",facet.wrap = T,axis.text.x.size = 12,point.size = 6) +
  geom_signif(comparisons = list(c("Group1", "Group3")),
              map_signif_level = T,
              textsize = 8, step_increase=0.1,vjust=0.5,margin_top = 0.09, test="t.test") + 
  geom_hline(yintercept = 
               mean(c(max(plot.score[1,which(annots$severeAH.annotsflux.groups=="Group1")]),
                      min(plot.score[1,which(annots$severeAH.annotsflux.groups=="Group3")]))), 
             linetype = "dashed", size = 1) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1))

mean(c(max(plot.score[1,which(annots$severeAH.annotsflux.groups=="Group1")]),
       min(plot.score[1,which(annots$severeAH.annotsflux.groups=="Group3")])))
# [1] 3.962765

plot.score.t <- as.data.frame(t(plot.score))
plot.score.t$disease_state <- as.character(annots$disease_state[c(1:51,80:89)])

## Figure 9B
ggplot(plot.score.t,
       aes(x = score, y = Pseudotime, color = disease_state)) +
  geom_smooth(method = "lm", color = "black",se=F, alpha = 0.5) +
  # stat_cor( color = "black", size = 5) +
  geom_point(size = 8, shape = 21, fill = as.character(ann_colors$disease_state[match(annots$disease_state[c(1:51,80:89)],
                                                                                      names(ann_colors$disease_state))]), color = "black") +
  # scale_fill_manual(values = ann_colors$disease_state) +
  labs(title = "Gene Score vs. HMR_5029 Flux")  + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16))

cor.test(as.numeric(plot.score.t$score),
         as.numeric(plot.score.t$Pseudotime),
         method="spearman")

# S = 5996, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.8414595 


plot.score.t$score <- as.numeric(plot.score.t$score)
plot.score.t$MELD <- as.numeric(plot.score.t$MELD)
plot.score.t$ABIC <- as.numeric(plot.score.t$ABIC)
plot.score.t$Child <- as.numeric(plot.score.t$Child)
plot.score.t$AST <- as.numeric(plot.score.t$AST)
plot.score.t$ALT <- as.numeric(plot.score.t$ALT)

## Supplemental Figure 9A
ggplot(plot.score.t[which(!(is.na(plot.score.t$MELD))),],
       aes(x = score, y = MELD, color = disease_state)) +
  geom_smooth(method = "lm", color = "black",se=F, alpha = 0.5) +
  geom_point(size = 8, shape = 21, fill = as.character(ann_colors$disease_state[match(annots$disease_state[c(1:51,80:89)[which(!(is.na(plot.score.t$MELD)))]],
                                                                                      names(ann_colors$disease_state))]), color = "black") +
  labs(title = "Gene Score vs. HMR_5029 Flux")  + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16))


cor.test(as.numeric(plot.score.t$score[which(!(is.na(plot.score.t$MELD)))]),
         as.numeric(plot.score.t$MELD[which(!(is.na(plot.score.t$MELD)))]),
         method="spearman")

# S = 5877.5, p-value = 0.06804
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.3032824 

## Supplemental Figure 9B
ggplot(plot.score.t[which(!is.na(plot.score.t$ABIC)),],
       aes(x = score, y = ABIC, color = disease_state)) +
  geom_smooth(method = "lm", color = "black",se=F, alpha = 0.5) +
  geom_point(size = 8, shape = 21, fill = as.character(ann_colors$disease_state[match(annots$disease_state[c(1:51,80:89)[which(!(is.na(plot.score.t$ABIC)))]],
                                                                                      names(ann_colors$disease_state))]), color = "black") +
  labs(title = "Gene Score vs. HMR_5029 Flux")  + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16))

cor.test(as.numeric(plot.score.t$score[which(!is.na(plot.score.t$ABIC))]),
         as.numeric(plot.score.t$ABIC[which(!is.na(plot.score.t$ABIC))]),
         method="spearman")

# S = 5930.6, p-value = 0.07425
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.2969945 
# 0.3453899 


## Supplemental Figure 9C
ggplot(plot.score.t[which(!is.na(plot.score.t$Child)),],
       aes(x = score, y = Child, color = disease_state)) +
  geom_smooth(method = "lm", color = "black",se=F, alpha = 0.5) +
  geom_point(size = 8, shape = 21, fill = as.character(ann_colors$disease_state[match(annots$disease_state[c(1:51,80:89)[which(!(is.na(plot.score.t$Child)))]],
                                                                                      names(ann_colors$disease_state))]), color = "black") +
  labs(title = "Gene Score vs. HMR_5029 Flux")  + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16))

cor.test(as.numeric(plot.score.t$score[which(!is.na(plot.score.t$Child))]),
         as.numeric(plot.score.t$Child[which(!is.na(plot.score.t$Child))]),
         method="spearman")

# S = 6597.1, p-value = 0.03876
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.3322761


## Figure 9C
plot.score.t$AST_ALT <- plot.score.t$AST/plot.score.t$ALT

ggplot(plot.score.t[which(!is.na(plot.score.t$AST_ALT)&plot.score.t$AST_ALT<30),],#&plot.score.t$disease_state%in%c("nonsevere.AH","severe.AH")),],
       aes(x = score, y = AST_ALT, color = disease_state)) +
  geom_smooth(method = "lm", color = "black",se=F, alpha = 0.5) +
  # stat_cor( color = "black", size = 5) +
  geom_point(size = 8, shape = 21, fill = as.character(ann_colors$disease_state[match(annots$disease_state[c(1:51,80:89)[which(!is.na(plot.score.t$AST_ALT)&plot.score.t$AST_ALT<30)]],
                                                                                      names(ann_colors$disease_state))]), color = "black") +
  # scale_fill_manual(values = ann_colors$disease_state) +
  labs(title = "Gene Score vs. HMR_5029 Flux")  + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16))

cor.test(as.numeric(plot.score.t$score[which(!is.na(plot.score.t$AST_ALT)&plot.score.t$AST_ALT<30)]),#&plot.score.t$disease_state%in%c("nonsevere.AH","severe.AH"))]),
         as.numeric(plot.score.t$AST_ALT[which(!is.na(plot.score.t$AST_ALT)&plot.score.t$AST_ALT<30)]),#&plot.score.t$disease_state%in%c("nonsevere.AH","severe.AH"))]),
         method="spearman")

# S = 11102, p-value = 0.0006306
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.466902


## predicting glucose use reprogramming scores following gene knockdown
predict.all.KD <- predict(svm.linear.cv,newdata = t(abs(bind_matrices[match(c("HMR_5029","HMR_4394","HMR_4381"),rownames(bind_matrices)),
                                                                      setdiff(grep("KD",HK.KD.annots$state),which(HK.KD.annots$state=="No.KD"))])) )
predict.all.KD <- as.data.frame(t(predict.all.KD))
rownames(predict.all.KD) = "predict.new"


# putting scores in a matrix
score.KD <- c()
for (i in seq(1,length(predict.all.KD),28)){
  if (i ==1){
    score.KD <- as.numeric(predict.all.KD[seq(i,i+27)])
  } else {
    score.KD <- rbind(score.KD,as.numeric(predict.all.KD[seq(i,i+27)]))
  }
}
rownames(score.KD) <- unique(strsplit2(names(predict.all.KD),split = "_")[,2])
colnames(score.KD) <- paste0(rownames(annots)[c(1:18,30:39)],"_KD")

## predicting glucose use reprogramming scores - no gene KD
predict.no.KD <- predict(svm.linear.cv,newdata = t(abs(FBA_mat_incl1[match(c("HMR_5029","HMR_4394","HMR_4381"),rownames(FBA_mat_incl1)),c(1:51,80:89)])))


score.KD <- as.data.frame(cbind(score.KD,matrix(rep(predict.no.KD, each = 1, times = nrow(score.KD)), nrow = nrow(score.KD), byrow = TRUE)))
colnames(score.KD)[-c(1:28)] <- names(Predict)

ML.annots <- data.frame(disease_state=c(rep("KD",28),as.character(annots$disease_state[c(1:51,80:89)])))
ML.annots$disease_state[which(ML.annots$disease_state=="explant.AH")] = "severe.AH"
ML.annots$disease_state <- factor(ML.annots$disease_state,levels=c("healthy.control","early.ASH","nonsevere.AH","severe.AH","KD"))
rownames(ML.annots) <- colnames(score.KD)

ML.annots$subgroups <- as.character(ML.annots$disease_state)
ML.annots$subgroups[grep(paste0(strsplit2(rownames(HKDC1.KD.clusts)[which(HKDC1.KD.clusts$Sample.Groups=="D")],split = "_")[,1],collapse = "|"),
                         rownames(ML.annots))] <- "Group1"
ML.annots$subgroups[grep(paste0(strsplit2(rownames(HKDC1.KD.clusts)[which(HKDC1.KD.clusts$Sample.Groups=="C")],split = "_")[,1],collapse = "|"),
                         rownames(ML.annots))] <- "Group3"
# ML.annots$subgroups[which(!is.na(HK.KD.annots$severeAH.annotsflux.groups[c(1:18,30:39)]))] <- na.omit(as.character(HK.KD.annots$severeAH.annotsflux.groups[c(1:18,30:39)]))
ML.annots$subgroups[which(ML.annots$subgroups%in%c("KD","severe.AH"))] <- "Group2"
ML.annots$subgroups[which(ML.annots$disease_state=="KD")] <- paste0(ML.annots$subgroups[which(ML.annots$disease_state=="KD")],"_","KD")
ML.annots$subgroups <- factor(ML.annots$subgroups,levels = c("healthy.control","early.ASH","nonsevere.AH",
                                                             "Group1","Group1_KD","Group2","Group2_KD",
                                                             "Group3","Group3_KD"))


ML.colors <- ann_colors
ML.colors[['disease_state']][["severeAH.KD"]] <- "purple"
ML.colors[['disease_state']][["explantAH.KD"]] <- "lightblue"
ML.colors[['subgroups']] <- ML.colors[["disease_state"]]
ML.colors[['subgroups']][["Group1"]] <- "blue"
ML.colors[['subgroups']][["Group3"]] <- "red"
ML.colors[['subgroups']][["Group2"]] <- "gray"
ML.colors[['subgroups']][["Group1_KD"]] <- "blue"
ML.colors[['subgroups']][["Group3_KD"]] <- "red"
ML.colors[['subgroups']][["Group2_KD"]] <- "gray"
ML.colors[['subgroups']][["explant.AH"]] <- ann_colors$disease_state[['explant.AH']]

set_annotations(ML.annots)
set_annot_cols(ML.colors)


to.plot <- as.data.frame(melt(t(score.KD[grep("HKDC1|PKM",rownames(score.KD)),])))
to.plot$subgroup <- ML.annots$subgroups[match(to.plot$Var1,rownames(ML.annots))]
to.plot$KD <- "No.KD"
to.plot$KD[grep("KD",to.plot$subgroup)] <- "KD"
to.plot$KD <- factor(to.plot$KD,levels = c("No.KD","KD"))
to.plot$Var2 <- factor(to.plot$Var2,levels = c("HKDC1.50pct.KD","HKDC1.KD","PKM.KD","HKDC1.PKM.KD"))

## Figure 9D
ggplot(to.plot[grep("HKDC1.KD|PKM.KD|HKDC1.50pct.KD",to.plot$Var2),], aes(x=subgroup,y=value,shape=KD,fill=subgroup,group=subgroup))+
  ggbeeswarm::geom_quasirandom(size=5,stroke = 0.25)+ scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = as.character(ML.colors[['subgroups']][match(levels(to.plot$subgroup),names(ML.colors[['subgroups']]))])) +
  facet_wrap(~Var2, scales="free",ncol=4)+
  geom_signif(comparisons = list(c("Group1", "Group1_KD"),
                                 c("Group2", "Group2_KD"),
                                 c("Group3", "Group3_KD"),
                                 c("early.ASH", "Group1_KD"),
                                 c("nonsevere.AH", "Group3_KD")),
              map_signif_level = function(p) {
                if(p>0.05) return("")
                if(p<0.05 & p>0.01)return("*")
                if(p<0.01 & p>0.001)return("**")
                if(p<0.001)return("***")},
              textsize = 8, step_increase=0.1,vjust=0.5,margin_top = 0.09,color = "black", test = "t.test") +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5,size=40),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16),
        panel.grid = element_blank(),
        strip.text = element_text(size = 25),
        strip.background.x = element_blank())


# significance testing
paired.testing <- c()
for (i in unique(to.plot$Var2)){
  test1 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1_KD")]),
                  alternative = "greater",paired=T)
  test2 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3_KD")]),
                  alternative = "greater",paired=T)
  test3 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group2")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group2_KD")]),
                  alternative = "greater",paired=T)
  
  paired.testing <- rbind(paired.testing,c(test1[['p.value']],test2[['p.value']],test3[['p.value']]))
}
rownames(paired.testing) <- rownames(score.KD)
colnames(paired.testing) <- c("Group1","Group3","Group2")
paired.testing



other.testing <- c()
for (i in unique(to.plot$Var2)){
  test1 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="early.ASH")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1_KD")]),
                  alternative = "less",paired=F)
  test2 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="nonsevere.AH")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3_KD")]),
                  alternative = "less",paired=F)
  other.testing <- rbind(other.testing,c(test1[['p.value']],test2[['p.value']]))
}
rownames(other.testing) <- rownames(score.KD)
colnames(other.testing) <- c("first.test","second.test")
other.testing

# plotting glucose use reprogramming scores for other KD experiments
to.plot <- as.data.frame(melt(t(score.KD[-grep("HKDC1|PKM",rownames(score.KD)),])))
to.plot$subgroup <- ML.annots$subgroups[match(to.plot$Var1,rownames(ML.annots))]
to.plot$KD <- "No.KD"
to.plot$KD[grep("KD",to.plot$subgroup)] <- "KD"
to.plot$KD <- factor(to.plot$KD,levels = c("No.KD","KD"))

## Supplemental Figure 10
ggplot(to.plot, aes(x=subgroup,y=value,shape=KD,fill=subgroup,group=subgroup))+
  # geom_point(size=5)+
  ggbeeswarm::geom_quasirandom(size=5,stroke = 0.25)+ scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = as.character(ML.colors[['subgroups']][match(levels(to.plot$subgroup),names(ML.colors[['subgroups']]))])) +
  facet_wrap(~Var2, scales="free")+
  geom_signif(comparisons = list(c("Group1", "Group1_KD"),
                                 c("Group2", "Group2_KD"),
                                 c("Group3", "Group3_KD"),
                                 c("early.ASH", "Group1_KD"),
                                 c("nonsevere.AH", "Group3_KD")),
              map_signif_level = function(p) {
                if(p>0.05) return("")
                if(p<0.05 & p>0.01)return("*")
                if(p<0.01 & p>0.001)return("**")
                if(p<0.001)return("***")},
              textsize = 8, step_increase=0.1,vjust=0.5,margin_top = 0.09,color = "black", test = "t.test") +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5,size=40),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16),
        panel.grid = element_blank(),
        strip.text = element_text(size = 25),
        strip.background.x = element_blank())


## significance testing
paired.testing <- c()
for (i in unique(to.plot$Var2)){
  test1 <- wilcox.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1")]),
                       as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1_KD")]),
                       alternative = "greater",paired=T)
  test2 <- wilcox.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3")]),
                       as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3_KD")]),
                       alternative = "greater",paired=T)
  test3 <- wilcox.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group2")]),
                       as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group2_KD")]),
                       alternative = "greater",paired=T)
  
  paired.testing <- rbind(paired.testing,c(test1[['p.value']],test2[['p.value']],test3[['p.value']]))
}
rownames(paired.testing) <- unique(to.plot$Var2)
colnames(paired.testing) <- c("Group1","Group3","Group2")
paired.testing



other.testing <- c()
for (i in unique(to.plot$Var2)){
  test1 <- wilcox.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="early.ASH")]),
                       as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1_KD")]),
                       alternative = "less",paired=F)
  test2 <- wilcox.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="nonsevere.AH")]),
                       as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3_KD")]),
                       alternative = "less",paired=F)
  other.testing <- rbind(other.testing,c(test1[['p.value']],test2[['p.value']]))
}
rownames(other.testing) <- unique(to.plot$Var2)
colnames(other.testing) <- c("first.test","second.test")
other.testing

## confirming glucose use reprogramming scoring metric in validation cohort
originalcohort_fluxes <- FBA_mat

GSE143318_fluxes <- read.table("I:/data_acquisition/public/Alex-storage/Pheflux/GSE143318/Fluxomes/PheFlux_Results_MaximumFluxes.txt",sep ='\t')
GSE143318_annots <- data.frame(Sample = colnames(GSE143318_fluxes))
rownames(GSE143318_annots) <- GSE143318_annots$Sample
GSE143318_annots$state <- "AH"
GSE143318_annots$state[grep("AC",rownames(GSE143318_annots))] <- "AC"
GSE143318_annots$state[grep("N",rownames(GSE143318_annots))] <- "Normal"
GSE143318_annots$state <- factor(GSE143318_annots$state, levels=c("Normal","AH","AC"))

# bind matrices together
bind_flux_matrices <- as.data.frame(cbind(originalcohort_fluxes[,c(1:51,80:89)],GSE143318_fluxes[,which(GSE143318_annots$state=="AH")]))

# setting up new annotations
new.annots <- data.frame(samples=c(colnames(bind_flux_matrices)),
                         disease_state = c(as.character(annots$disease_state[c(1:51,80:89)]),as.character(GSE143318_annots$state[which(GSE143318_annots$state=="AH")])))
rownames(new.annots) <- new.annots$samples

new.annots$severeAH.annotsflux.groups <- NA
new.annots$severeAH.annotsflux.groups <- c(as.character(annots$severeAH.annotsflux.groups[c(1:51,80:89)]),rep("AH",length(which(GSE143318_annots$state=="AH"))))
new.annots$severeAH.annotsflux.groups[-which(new.annots$disease_state%in%c("severe.AH","explant.AH","AH"))] <- new.annots$disease_state[-which(new.annots$disease_state%in%c("severe.AH","explant.AH","AH"))]


new.annots$severeAH.annotsflux.groups <- factor(new.annots$severeAH.annotsflux.groups, levels = c("healthy.control","early.ASH","nonsevere.AH","Group1",
                                                                                                  "Group2","Group3","AH"))

ann_colors[["severeAH.annotsflux.groups"]] <- c("Group3"="red","Group1"="blue", "healthy.control"='#f280ae',"early.ASH"="#7aaf3e","nonsevere.AH"="#f2e32b",
                                                "AH" = "black","Group2"="gray")


# batch correction between original and validation cohorts
batch_info <- rep("Batch1",nrow(new.annots))
batch_info[which(new.annots$severeAH.annotsflux.groups=="AH")] <-"Batch2"
batch_corrected_data <- removeBatchEffect(bind_flux_matrices, 
                                          batch = batch_info)


## Figure 10A
set_annotations(new.annots)
set_annot_samps("severeAH.annotsflux.groups")
set_annot_cols(ann_colors)
p = myHeatmapByAnnotation(abs(batch_corrected_data[match(reactions2,rownames(batch_corrected_data)),]),
                          main = "Glycolysis Reactions", method = "euclidean", groupings = "severeAH.annotsflux.groups",
                          scale.rows = "zscore",show.rownames = T,#gaps.col = c(10,28,46,64),#100,118,136,154),
                          gaps.row = F, clust.rows = F, clust.cols = T,show.colnames = T,fontsize.row = 12)

## Figure 10B
pca <- myPCA(t(scale(t(batch_corrected_data[match(reactions2,rownames(batch_corrected_data)),]))),
             color.by = "severeAH.annotsflux.groups",point.size = 8,return.ggplot.input = T, PCs.to.plot = c("PC1","PC2"))
pca$PCA_variability
to.plot <- pca$input_data[,1:2]
to.plot$PC1 <- -to.plot$PC1
ggplot(to.plot, aes(x = PC1, y = PC2)) +
  # geom_smooth(method = "lm", se = FALSE, color = "black", linetype ="dashed") +  # Add dotted black line , data = subset(to.plot, flux.pseudo > 130)
  geom_point(color = "black", fill =ann_colors$severeAH.annotsflux.groups[match(new.annots$severeAH.annotsflux.groups[match(rownames(to.plot),new.annots$samples)],names(ann_colors$severeAH.annotsflux.groups))], shape = 21, size = 8) +  # Plot points
  # geom_text(aes(label=ifelse(PC1< -2 & PC2 < - 0.25,as.character(rownames(annots)[c(1:18,30:39)]),'')),hjust=-0.35,vjust=-0.25,size=5) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16)) #+ ylim(c(0,9))



## combining fluxomes from original and validation cohort following gene knockdowns
bind_matrices_original <- bind_matrices

# validation cohort
HKDC1_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/GSE143318_KD/Fluxomes/PheFlux_Results_HKDC1.KDFluxes.txt', sep = '\t')
FBA_mat_HKDC1_KD <- GSE143318_fluxes
FBA_mat_HKDC1_KD <- as.data.frame(cbind(HKDC1_KD,FBA_mat_HKDC1_KD[,-match(colnames(HKDC1_KD),colnames(GSE143318_fluxes))]))

PKM_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/GSE143318_KD/Fluxomes/PheFlux_Results_PKM.KDFluxes.txt', sep = '\t')
FBA_mat_PKM_KD <- GSE143318_fluxes
FBA_mat_PKM_KD <- as.data.frame(cbind(PKM_KD,FBA_mat_PKM_KD[,-match(colnames(PKM_KD),colnames(GSE143318_fluxes))]))

HKDC1_PKM_KD <- read.table('C:/Users/arm037/OneDrive - Thomas Jefferson University/GEM_Modeling_Dec2020/PatientGEMPaper/Pheflux/GSE143318_KD/Fluxomes/PheFlux_Results_HKDC1.PKM.KDFluxes.txt', sep = '\t')
FBA_mat_HKDC1_PKM_KD <- GSE143318_fluxes
FBA_mat_HKDC1_PKM_KD <- as.data.frame(cbind(HKDC1_PKM_KD,FBA_mat_HKDC1_PKM_KD[,-match(colnames(HKDC1_PKM_KD),colnames(GSE143318_fluxes))]))


bind_matrices <- cbind(FBA_mat_HKDC1_KD[,-grep("N",colnames(GSE143318_fluxes))],FBA_mat_PKM_KD[,-grep("N",colnames(GSE143318_fluxes))],
                       FBA_mat_HKDC1_PKM_KD[,-grep("N",colnames(GSE143318_fluxes))],GSE143318_fluxes)
colnames(bind_matrices)[(seq(1,18*7,18)[1]):(seq(1,18*7,18)[2]-1)] <- paste0(colnames(bind_matrices)[(seq(1,18*7,18)[1]):(seq(1,18*7,18)[2]-1)],"_HKDC1.KD")
colnames(bind_matrices)[(seq(1,18*7,18)[2]):(seq(1,18*7,18)[3]-1)] <- paste0(colnames(bind_matrices)[(seq(1,18*7,18)[2]):(seq(1,18*7,18)[3]-1)],"_PKM.KD")
colnames(bind_matrices)[(seq(1,18*7,18)[3]):(seq(1,18*7,18)[4]-1)] <- paste0(colnames(bind_matrices)[(seq(1,18*7,18)[3]):(seq(1,18*7,18)[4]-1)],"_HKDC1.PKM.KD")


GSE142530.KD.annots <- as.data.frame(colnames(bind_matrices))
rownames(GSE142530.KD.annots) <- GSE142530.KD.annots[,1]
GSE142530.KD.annots$state <- "No.KD"
GSE142530.KD.annots$state[grep("_HKDC1.KD",rownames(GSE142530.KD.annots))] <- "HKDC1.KD"
GSE142530.KD.annots$state[grep("_PKM.KD",rownames(GSE142530.KD.annots))] <- "PKM.KD"
GSE142530.KD.annots$state[grep("_HKDC1.PKM.KD",rownames(GSE142530.KD.annots))] <- "HKDC1.PKM.KD"
GSE142530.KD.annots$state[grep("^N",GSE142530.KD.annots[,1])] <- "Normal"
GSE142530.KD.annots$state <- factor(GSE142530.KD.annots$state, levels = c("Normal","No.KD","HKDC1.KD","PKM.KD","HKDC1.PKM.KD"))
GSE142530.KD.annots[,1] <- strsplit2(GSE142530.KD.annots[,1],split="_")[,1]
GSE142530.KD.annots$disease_state <- GSE143318_annots$state[match(GSE142530.KD.annots[,1],rownames(GSE143318_annots))]


bind_matrices_valid <- bind_matrices

all_bind_matrices <- as.data.frame(cbind(bind_matrices_original[,grep("HKDC1.KD|PKM.KD",colnames(bind_matrices_original))],bind_matrices_valid))



# predicting glucose use reprogramming scores without batch correction
no.correct <- abs(all_bind_matrices[,-which(new.annots$disease_state%in%c("AC","Normal"))])

predict.all.KD <- predict(svm.linear.cv,newdata = t(abs(no.correct[match(c("HMR_5029","HMR_4394","HMR_4381"),rownames(no.correct)),])) )
predict.all.KD <- as.data.frame(t(predict.all.KD))
rownames(predict.all.KD) = "predict.new"

# put scores into matrix
predict.all.KD.subset <- predict.all.KD[-union(which(colnames(predict.all.KD)%in%rownames(annots)),which(colnames(predict.all.KD)%in%rownames(GSE143318_annots)))]
score.KD <- c()
for (i in c("_HKDC1.KD","_PKM.KD","_HKDC1.PKM.KD")){
  if (i =="_HKDC1.KD"){
    score.KD <- as.numeric(predict.all.KD.subset[grep(i,colnames(predict.all.KD.subset))])
  } else {
    score.KD <- rbind(score.KD,as.numeric(predict.all.KD.subset[grep(i,colnames(predict.all.KD.subset))]))
  }
}
rownames(score.KD) <- c("HKDC1.KD","PKM.KD","HKDC1.PKM.KD")
colnames(score.KD) <- paste0(strsplit2(colnames(predict.all.KD.subset)[grep("HKDC1.KD",colnames(predict.all.KD.subset))],"_")[,1],"_KD")

# predict scores in validation cohort prior to any gene KD
GSE143318_forML <- as.data.frame(t(abs(GSE143318_fluxes[match(reactions,rownames(GSE143318_fluxes)),])))
Predict.valid <- predict(svm.linear.cv,newdata = GSE143318_forML )

score.KD <- as.data.frame(cbind(score.KD,matrix(rep(c(predict.no.KD,Predict.valid[grep("AH",names(Predict.valid))]), each = 1, times = nrow(score.KD)), nrow = nrow(score.KD), byrow = TRUE)))
colnames(score.KD)[grep("^V",colnames(score.KD))] <- names(c(predict.no.KD,Predict.valid[grep("AH",names(Predict.valid))]))

# setting up new annotations
ML.annots <- data.frame(disease_state=c(rep("KD",length(colnames(score.KD)[grep("_KD",colnames(score.KD))])),
                                        as.character(annots$disease_state[c(1:51,80:89)]),rep("AH",length(which(GSE143318_annots$state=="AH")))))
ML.annots$disease_state[which(ML.annots$disease_state=="explant.AH")] = "severe.AH"
ML.annots$disease_state <- factor(ML.annots$disease_state,levels=c("healthy.control","early.ASH","nonsevere.AH","severe.AH","AH","KD"))
rownames(ML.annots) <- colnames(score.KD)

ML.annots$subgroups <- as.character(ML.annots$disease_state)
ML.annots$subgroups[grep(paste0(rownames(annots)[which(annots$severeAH.annotsflux.groups=="Group1")],collapse = "|"),
                         rownames(ML.annots))] <- "Group1"
ML.annots$subgroups[grep(paste0(rownames(annots)[which(annots$severeAH.annotsflux.groups=="Group3")],collapse = "|"),
                         rownames(ML.annots))] <- "Group3"
ML.annots$subgroups[which(ML.annots$subgroups%in%c("KD","severe.AH"))] <- "Group2"
ML.annots$subgroups[which(ML.annots$disease_state=="KD")] <- paste0(ML.annots$subgroups[which(ML.annots$disease_state=="KD")],"_","KD")
ML.annots$subgroups[intersect(which(ML.annots$disease_state%in%c("KD")),grep("^AH",rownames(ML.annots)))] <- "AH_KD"
ML.annots$subgroups[intersect(which(ML.annots$disease_state%in%c("AH")),grep("^AH",rownames(ML.annots)))] <- "AH"
ML.annots$subgroups <- factor(ML.annots$subgroups,levels = c("healthy.control","early.ASH","nonsevere.AH",
                                                             "Group1","Group1_KD","Group2","Group2_KD",
                                                             "Group3","Group3_KD","AH","AH_KD"))


ML.colors <- ann_colors
ML.colors[['disease_state']][["severeAH.KD"]] <- "purple"
ML.colors[['disease_state']][["explantAH.KD"]] <- "lightblue"
ML.colors[['subgroups']] <- ML.colors[["disease_state"]]
ML.colors[['subgroups']][["Group1"]] <- "blue"
ML.colors[['subgroups']][["Group3"]] <- "red"
ML.colors[['subgroups']][["Group2"]] <- "gray"
ML.colors[['subgroups']][["Group1_KD"]] <- "blue"
ML.colors[['subgroups']][["Group3_KD"]] <- "red"
ML.colors[['subgroups']][["Group2_KD"]] <- "gray"
ML.colors[['subgroups']][["AH"]] <- "black"
ML.colors[['subgroups']][["AH_KD"]] <- "black"
ML.colors[['subgroups']][["explant.AH"]] <- ann_colors$disease_state[['explant.AH']]

set_annotations(ML.annots)
set_annot_cols(ML.colors)


to.plot <- as.data.frame(melt(t(score.KD)))
to.plot$subgroup <- ML.annots$subgroups[match(to.plot$Var1,rownames(ML.annots))]
to.plot$KD <- "No.KD"
to.plot$KD[grep("KD",to.plot$subgroup)] <- "KD"
to.plot$KD <- factor(to.plot$KD,levels = c("No.KD","KD"))

## Figure 10C
ggplot(to.plot, aes(x=subgroup,y=value,shape=KD,fill=subgroup,group=subgroup))+
  geom_hline(yintercept = mean(c(max(to.plot$value[which(to.plot$subgroup=="Group1"&to.plot$Var2=="HKDC1.KD")]),
                                 min(to.plot$value[which(to.plot$subgroup=="Group3"&to.plot$Var2=="HKDC1.KD")]))),
             linetype="dashed",size=1)+
  ggbeeswarm::geom_quasirandom(size=5,stroke = 0.25)+ scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = as.character(ML.colors[['subgroups']][match(levels(to.plot$subgroup),names(ML.colors[['subgroups']]))])) +
  facet_wrap(~Var2, scales="free",ncol = 4)+
  geom_signif(comparisons = list(c("AH","AH_KD"),
                                 c("nonsevere.AH","AH"),
                                 c("early.ASH","AH_KD")),
              map_signif_level = function(p) {
                if(is.numeric(p)) return("")},
              textsize = 8, step_increase=0.1,vjust=0.5,margin_top = -0.05,color = "black") +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5,size=40),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_blank(),  # No major grid lines
        panel.grid.minor = element_blank(),  # No minor grid lines
        panel.border = element_blank(),  # No panel border
        axis.line = element_line(color = "black", size = 1),  # Axis lines
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16),
        panel.grid = element_blank(),
        strip.text = element_text(size = 25),
        strip.background.x = element_blank())


# statistical testing
paired.testing <- c()
for (i in unique(to.plot$Var2)){
  test1 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1_KD")]),
                  alternative = "greater",paired=T)
  test2 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3_KD")]),
                  alternative = "greater",paired=T)
  test3 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group2")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group2_KD")]),
                  alternative = "greater",paired=T)
  test4 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="AH")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="AH_KD")]),
                  alternative = "greater",paired=T)
  
  paired.testing <- rbind(paired.testing,c(test1[['p.value']],test2[['p.value']],test3[['p.value']],test4[['p.value']]))
}
rownames(paired.testing) <- unique(to.plot$Var2)
colnames(paired.testing) <- c("Group1","Group3","Group2","AH")
paired.testing



other.testing <- c()
for (i in unique(to.plot$Var2)){
  test1 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="early.ASH")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group1_KD")]),
                  alternative = "less",paired=F)
  test2 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="nonsevere.AH")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="Group3_KD")]),
                  alternative = "less",paired=F)
  test3 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="nonsevere.AH")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="AH")]),
                  alternative = "less",paired=F)
  test4 <- t.test(as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="early.ASH")]),
                  as.numeric(to.plot$value[which(to.plot$Var2==i & to.plot$subgroup=="AH_KD")]),
                  alternative = "less",paired=F)
  other.testing <- rbind(other.testing,c(test1[['p.value']],test2[['p.value']],test3[['p.value']],test4[['p.value']]))
}
rownames(other.testing) <- unique(to.plot$Var2)
colnames(other.testing) <- c("first.test","second.test","third.test","forth.test")
other.testing
