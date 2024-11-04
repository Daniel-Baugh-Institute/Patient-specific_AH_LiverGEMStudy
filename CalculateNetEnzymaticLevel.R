## FUNCTIONS ###

#####################################
booleanVectorRule <- function(rule, fpkmDic) {
  boolean_list <- c()
  
  vector_rule <- gsub("or |and |'|\\(|\\)", "", rule) # Removes "or", "and", "'", "(", and ")"
  vector <- strsplit2(vector_rule, " ")
  
  # g_vector <- paste0("G_", vector)
  
  for (gene in vector) {
    if (gene %in% rownames(fpkmDic)) {
      boolean_list <- c(boolean_list, 'True')
    } else {
      boolean_list <- c(boolean_list, 'False')
    }
  }
  
  return(boolean_list)
}
######################################
getG <- function(rule, fpkmDic,patientNum) {
  orList <- c()
  
  # Devide the rule by "or" and iterate over each resulting subrule
  subrules <- strsplit2(rule, " or ")
  
  for (subrule in subrules) {
    vector <- strsplit2(subrule, ' and ')
    g_vector = c()
    for (gene in vector) {
      gene <- gsub("\\(|\\)| ", "", gene) # Removes "(", ")", and spaces
      g_vector <- c(g_vector,gene)
    }
    if ((length(g_vector)==1) & (nchar(g_vector)[1]==15) & (!(g_vector[1]%in%rownames(fpkmDic)))){
      value = 0
    } else {
    value <- min(fpkmDic[na.omit(match(g_vector,rownames(fpkmDic))),patientNum]) ########### CHANGE FOR EACH PATIENT
    }
    orList <- c(orList, value) # Add the minimum to a list
  }
  
  return(sum(orList)) # Print the sum of the list
}




####################################
getEg <- function(rxns, fpkmDic,patientNum) {
  g_metab <- c()  # gene expression of reactions partaking in the metabolism
  
  for (i in 1:length(rxns)) {
    ##########################################################
    ## Gets the GPR and a boolean list with known or unknown genes
    rule <- ihuman$`GENE ASSOCIATION`[match(rxns[i],ihuman$rxnRetired)]
    rule[which(is.na(rule))] <- 0
    
    boolean_list <- booleanVectorRule(rule, fpkmDic)
    ##########################################################
    ## Gets the expression value 'g'
    if (!("False" %in% boolean_list || rule == '0')) {
      g <- getG(rule, fpkmDic,patientNum)
      g_metab <- c(g_metab, g + 1e-8)
    }
  }
  
  ##############################################################
  ## Obtains a median value
  E_g <- median(g_metab)
  return(list(E_g = E_g, g_metab = g_metab))
}

## Calculate GPR and Correlations ###

## Calculate GPR's
GPR <- matrix(NA,nrow = nrow(FBA_mat),ncol(FBA_mat))


for (j in 1:length(colnames(FBA_mat))) {
  getEg_out = getEg(rownames(FBA_mat),metgenes_patient,j) 
  for (i in 1:length(rownames(FBA_mat))){
    ##########################################################
    ## Gets the GPR and a boolean list with known or unknown genes
    rule <- ihuman$`GENE ASSOCIATION`[match(rownames(FBA_mat)[i],ihuman$rxnRetired)]
    rule[which(is.na(rule))] <- 0
    boolean_list = booleanVectorRule(rule,metgenes_patient)
    ##########################################################
    ## Gets the expression value 'g'
    # get 'g' for reaction with GPR.
    if (!('False' %in% boolean_list || rule == '0')){
      GPR[i,j] = getG(rule, metgenes_patient,j)+1e-8#*1e-8
    } else {
      GPR[i,j] = getEg_out[["E_g"]] }
  }
}

GPR <- as.data.frame(GPR)
rownames(GPR) <- rownames(FBA_mat)
colnames(GPR) <- colnames(FBA_mat)

# write.table(GPR, "GPR.txt", sep = '\t')




