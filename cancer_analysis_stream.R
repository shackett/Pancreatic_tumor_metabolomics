
setwd("/Users/Sean/Desktop/Rabinowitz/Cancer_metabolomics")
source("cancer_lib.R")
library("gplots")
library("colorRamps")
library("qvalue")
library(ggplot2)
library(data.table)

options(stringsAsFactors = F)

#types are small.molecules or lipomics
analysis.type = "small.molecules"

#either load previous bootstraps (generated locally or on cluster) or generate bootstrap null t statistics
bs.in.line = FALSE

#choose a low or high concentration based upon signal linearity (T) or just use the low concentration sample (F)
optim_DR = FALSE

#if exploratory == FALSE:
# remove 2 - 1 and 1 - 2 mappings of metabolites to peaks (for a given instrument)
# average across instruments
# average across multiple tumor samples from the same patient
exploratory = TRUE

#if thorough == TRUE, then additional plots and statistical tests will be utilized
thorough = FALSE


if(analysis.type == "small.molecules"){
  
  #small molecule metabolomics
  TNmeta <- read.delim("TNmetabolomics_6.27.13_1of2.txt", fill = TRUE, header = FALSE, colClasses = "character")
  
  ### read.delim("old
  #TNmeta_old <- read.delim("OldData/TNmetabolomics_10.6.12_1of2.txt", fill = TRUE, header = FALSE, colClasses = "character")
  #TNmeta_old[,1:3] == TNmeta[,1:3]
  
  TNmeta_sub <- read.delim("TNmetabolomics_6.27.13_2of2.txt", fill = TRUE, header = FALSE, colClasses = "character")
  TNmeta_sub <- TNmeta_sub[,apply(TNmeta_sub[1:8,] == "", 2, sum, na.rm = TRUE) != 8]
  TNmeta <- cbind(TNmeta, TNmeta_sub[,-c(1:3)])
  
  data.list = small.mol.tables(TNmeta)
  header.data <- data.list$header.data
  TNmeta <- data.list$TNmeta
  meta.data <- data.list$meta.data
  
  meta.data.raw <- meta.data
  
}



if(analysis.type == "lipomics"){
  #lipomics
  TNmeta <- read.delim("Tissue_lipid_master_list_102012.txt", fill = TRUE, header = FALSE, colClasses = "character")	
  data.list = lipomics.tables(TNmeta)
  header.data <- data.list$header.data
  TNmeta <- data.list$TNmeta
  meta.data <- data.list$meta.data
}






#read in sample covariates table

sample.covariates <- read.covariates()
clinical_data <- read.table("cancer_association/clinicalData2.tsv", header = TRUE)

#### Remove non-samples and unwanted samples ####

TNmeta <- TNmeta[,(header.data$Status %in% c("Tumor", "Normal"))]
header.data <- header.data[(header.data$Status %in% c("Tumor", "Normal")),]

invalid.sampleID <- union(c(as.character(sample.covariates$PatientID[sample.covariates$Pathology == "NPCA"][!is.na(sample.covariates$PatientID[sample.covariates$Pathology == "NPCA"])]), "P0011K", "P0013M"), clinical_data$tumorID[clinical_data$Provenance != "Adeno"][!is.na(clinical_data$tumorID[clinical_data$Provenance != "Adeno"])])

#clinical_data[clinical_data$tumorID %in% invalid.sampleID,]
#sample.covariates[sample.covariates$PatientID %in% invalid.sampleID,]

TNmeta <- TNmeta[,!(header.data$Patient %in% invalid.sampleID)]
header.data <- header.data[!(header.data$Patient %in% invalid.sampleID),]


#### Analyze the linearity of the ion counts: if 2x sample < 2*1x sample then use low concentrations

if(analysis.type == "small.molecules"){
  return_list = dilution_linearity(TNmeta, header.data)
  #either choose whether low or high conentration samples should be used on a per-peak basis
  #featureLH <- dyn_low_high(return_list)
  #or, just choose the lowest or highest concentration sample
  featureLH <- rep("L", times = nrow(TNmeta))
  
}
if(analysis.type == "lipomics"){
  featureLH <- rep("L", times = nrow(TNmeta))
}	


#cbind(meta.data, linear.comp$slope)
#par(mfrow = c(1,2))
#plot(anova.ind.tumor$patSS ~ anova.ind.tumor$disSS, main = "ANOVA: Variance explained by tumor vs ind", xlab = "tumorSS", ylab = "indSS")
#plot(anova.ind.tumor$patF ~ anova.ind.tumor$disF, main = "ANOVA: F-statistic for tumor and ind", xlab = "tumorF", ylab = "indF")
#plot(anova.ind.tumor$patMS ~ anova.ind.tumor$disMS, main = "ANOVA: Variance explained by tumor vs ind (per df)", xlab = "tumorMS", ylab = "indMS")



##### use median normalization seperately for each instrument

Instruments <- unique(meta.data$Instrument)

InormMat <- matrix(data = NA, ncol = length(TNmeta[1,]), nrow = length(TNmeta[,1]))

#TNmeta <- log(TNmeta, base = 2)

for (int in Instruments){
  
  threshold <- 100
  meta.floor.val <- threshold * 10
  meta.floor <- TNmeta[meta.data$Instrument %in% int,]
  meta.floor[meta.floor < meta.floor.val] <- NA
  
  row.median <- apply(meta.floor, 1, median, na.rm = TRUE)
  scaling.factors <- rep(NA, length(meta.floor[1,]))
  for (j in 1:length(meta.floor[1,])){
    scaling.factors[j] <- median(row.median[!is.na(meta.floor[,j])]/meta.floor[,j][!is.na(meta.floor[,j])])
  }
  
  meta.floor.val <- threshold
  meta.floor <- TNmeta[meta.data$Instrument %in% int,]
  meta.floor[meta.floor < meta.floor.val] <- NA
  
  InormMat[meta.data$Instrument %in% int,] <- t(t(meta.floor) * scaling.factors)
}

InormMat[InormMat < threshold] <- NA
InormMat <- log(InormMat, base = 2)
TNmeta.norm <- InormMat
rownames(TNmeta.norm) <- rownames(TNmeta)
colnames(TNmeta.norm) <- colnames(TNmeta)








####### renormalize by normalizing residuals, after the principal effects have been removed by regularized PCA

#not implemented





#####

if(analysis.type == "small.molecules"){
  unique.terms = c("Patient", "Status", "Dilution", "Sample Number", "exo-vivo Timecourse")
  samples <- sample.combos.fxn2(header.data, min.is.two = FALSE, no.blanks = TRUE, first.timecourse = TRUE, c(0, 10), unique.terms)
}
if(analysis.type == "lipomics"){
  unique.terms = c("Patient", "Status", "Dilution")
  samples <- sample.combos.fxn2(header.data, min.is.two = FALSE, no.blanks = TRUE, first.timecourse = FALSE, c(0, 10), unique.terms)
}


#choose features depending on featureLH and generate header and point estimate matrix.

#combine tissue samples from the same patient
if(exploratory == FALSE){unique.terms = c("Patient", "Status")}else{
  unique.terms = c("Patient", "Status", "Sample Number")
}

if(optim_DR == FALSE){
  
  TN.header = header.data[(apply(samples, 1, sum) == 1),]
  TNmeta.norm.subset = TNmeta.norm[,(apply(samples, 1, sum) == 1)]
  
  subset.names = unique(apply(TN.header[colnames(TN.header) %in% unique.terms], 1, paste, collapse = ""))
  
  TN.subset.header <- as.data.frame(matrix(data = NA, ncol = length(unique.terms), nrow = length(subset.names)))
  TN.subset.point <- matrix(data = NA, ncol = length(subset.names), nrow = length(TNmeta.norm[,1])); colnames(TN.subset.header) <- unique.terms
  
  for (i in 1:length(subset.names)){
    
    matches <- apply(TN.header[,colnames(header.data) %in% unique.terms], 1, paste, collapse = "") %in% subset.names[i]
    
    TN.subset.header[i,] <- lapply(TN.header[matches,][1,colnames(header.data) %in% unique.terms], as.character)
    
    relevant_samples <- matches & TN.header$Dilution == "Low"
    if(sum(relevant_samples) == 0){
      relevant_samples <- matches
    }
    TN.subset.point[,i] <- apply(matrix(TNmeta.norm.subset[,relevant_samples], ncol = sum(relevant_samples)), 1, mean, na.rm = TRUE)
  }
  
}else{
  
  #### double check this output if running with this option ######
  
  subset.names = unique(apply(TN.header[colnames(TN.header) %in% unique.terms], 1, paste, collapse = ""))
  
  TN.subset.header <- as.data.frame(matrix(data = NA, ncol = length(unique.terms), nrow = length(subset.names)))
  TN.subset.point <- matrix(data = NA, ncol = length(subset.names), nrow = length(TNmeta.norm[,1]))
  colnames(TN.subset.header) <- unique.terms
  
  for (i in 1:length(subset.names)){
    
    matches <- apply(header.data[,colnames(header.data) %in% unique.terms], 1, paste, collapse = "") %in% subset.names[i]
    
    TN.subset.header[i,] <- lapply(header.data[matches,][1,colnames(header.data) %in% unique.terms], as.character)
    
    TNabund.subset <- TNmeta.norm[,matches]
    
    dilution.order <- as.character(header.data[matches,]$Dilution); dilution.order[dilution.order == "High"] <- 2; dilution.order[dilution.order == "Low"] <- 1; dilution.order[dilution.order == "2x_Low"] <- 0; dilution.order <- as.numeric(dilution.order)
    
    high_samp <- dilution.order == max(dilution.order)
    low_samp <- dilution.order == min(dilution.order)
    
    TN.subset.point[,i] <- unlist(lapply(c(1:length(TNmeta.norm[,1])), function(x){
      ifelse(featureLH[x] == "H", mean(TNabund.subset[x,high_samp], na.rm = TRUE), mean(TNabund.subset[x,low_samp], na.rm = TRUE))
    }))}}

#### remove 2 - 1, 1 - 2 metabolite - peak mappings and sparse peaks, and average across instruments

if(exploratory == FALSE){
  
  #remove unknowns
  
  unks <- c("X -136.010 [90.515-91.515]", "Unk_129@10--147", "myo-inositol")
  unks <- c(unks, "D-Mannose/L-Sorbose-potential-2", "citrate-isocitrate", "2-hydroxybutyrate/1-hydroxybutanoic acid-potential-1", "2-hydroxybutyrate/1-hydroxybutanoic acid-potential-2",
            "D-sedoheptulose-7-phosphate-1", "D-sedoheptulose-7-phosphate-2")
  unk_met <- meta.data$Metabolite %in% unks
  
  degenerate_count <- data.table(meta.data)
  degenerate_count <- degenerate_count[!(Metablite_name %in% unk_met),,]
  degenerate_count <- degenerate_count[,list(npeaks = length(Peak_Number)), by = c("Metabolite", "Instrument")]
  degenerate_count <- degenerate_count[degenerate_count$npeaks != 1,]
  
  #one_two_map <- meta.data$Metabolite %in% degenerate_count$Metabolite[degenerate_count$npeaks != 1]
  
  one_two_map <- sapply(1:nrow(meta.data), function(x){
   any(meta.data$Metabolite[x] == degenerate_count[,Metabolite] & meta.data$Instrument[x] == degenerate_count[,Instrument])
  })
  
  #one_two_map <- meta.data$"Peak_Number" != ""
  
  
  #two to one mapping
  
  fraction.missing <- 0.75
  
  extra_sparse <- apply(is.nan(TN.subset.point), 1, sum) > length(TN.subset.point[1,])*fraction.missing
  
  same_peak <- sapply((1:length(TN.subset.point[,1]))[!extra_sparse], function(x){
    apply((TN.subset.point[x,] %*% t(rep(1, times = length(TN.subset.point[!extra_sparse,1]))) - t(TN.subset.point[!extra_sparse,])), 2, sum, na.rm = TRUE) == 0})
  
  two_one_map <- 1:length(TN.subset.point[,1]) %in% (1:length(TN.subset.point[,1]))[!extra_sparse][apply(same_peak, 2, sum) - (1 + sum(apply(same_peak, 2, sum) == length(TN.subset.point[!extra_sparse,1]))) != 0]
  
  #choose metabolites that weren't removed from the previously
  
  unique_matches <- c(1:length(TN.subset.point[,1]))[!extra_sparse & !one_two_map & !two_one_map & ! unk_met]
  
  #save(unique_matches, file = "unique_match.Rdata")
  
  meta.data_sub <- meta.data[unique_matches,]
  
  #table(table(meta.data_sub$Metabolite))
  
  #meta.data_sub[order(unique(meta.data_sub$Metabolite)),]
  #sort(unique(meta.data_sub$Metabolite))
  
  unique.mets <- sort(unique(meta.data_sub$Metabolite))
  
  #average rows from the same metabolite to generate TN.subset.point.simple
  
  TN.subset.point[is.nan(TN.subset.point)] <- NA
  
  TN.subset.point.simple <- matrix(NA, ncol = length(TN.subset.point[1,]), nrow = length(unique.mets))
  for(u in c(1:length(unique.mets))){
    
    matches <- unique_matches[meta.data_sub$Metabolite %in% unique.mets[u]]
    
    TN.subset.point.simple[u,] <- apply(matrix(TN.subset.point[matches,], nrow = length(matches), ncol = length(TN.subset.point[1,])) - apply(matrix(TN.subset.point[matches,], nrow = length(matches), ncol = length(TN.subset.point[1,])), 1, mean, na.rm = TRUE), 2, mean, na.rm = TRUE) + mean(apply(matrix(TN.subset.point[matches,], nrow = length(matches), ncol = length(TN.subset.point[1,])), 1, mean, na.rm = TRUE))
    
    #TN.subset.point.simple[u,] <- apply(matrix(TN.subset.point[matches,], nrow = length(matches), ncol = length(TN.subset.point[1,])), 2, mean, na.rm = TRUE)
    
  }
  TN.subset.point <- TN.subset.point.simple
  meta.data <- data.frame(Metabolite_name = unique.mets, Metabolite = unique.mets); rownames(meta.data) <- unique.mets
}

#standardize

TN.point.std <- (TN.subset.point - apply(TN.subset.point, 1, mean, na.rm = TRUE)) #/ apply(TN.subset.point, 1, sd, na.rm = TRUE)
#TN.point.std <- (TN.subset.point - apply(TN.subset.point, 1, mean, na.rm = TRUE)) / apply(TN.subset.point, 1, sd, na.rm = TRUE)
TN.point.std[is.nan(TN.point.std)] <- NA
if(exploratory == FALSE){rownames(TN.point.std) <- unique.mets}else{rownames(TN.point.std) <- rownames(TNmeta.norm)}
colnames(TN.point.std) <- subset.names




#tmp <- read.delim("~/Desktop/supplementalTable2014(1).txt")
#tmp[,1][!(tmp[,1] %in% rownames(TN.point.std))]








#write.table(TN.point.std, file = "Tables/lipidstdTNmeta.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
#write.table(TN.subset.header, file = "Tables/lipidstdheader.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

#split into tumor and normal samples

tumor.samples <- TN.subset.header[(TN.subset.header$Status == "Tumor") & (TN.subset.header$Patient %in% TN.subset.header[TN.subset.header$Status == "Normal",]$Patient),]

tumor.mat <- TN.point.std[,(TN.subset.header$Status == "Tumor") & (TN.subset.header$Patient %in% TN.subset.header[TN.subset.header$Status == "Normal",]$Patient)]

normal.mat <- matrix(data = NA, ncol = length(tumor.mat[1,]), nrow = length(tumor.mat[,1]))
for (i in 1:length(normal.mat[1,])){
  normal.mat[,i] <- TN.point.std[,(TN.subset.header[,1] ==  tumor.samples[i,1]) & (TN.subset.header$Status == "Normal")]
}




#output data, metabolite and header table
if(analysis.type == "small.molecules"){
  
  Low_dilution_patients <- unique(header.data$Patient[header.data$Dilution == "Low" & header.data$Patient %in% tumor.samples$Patient])
  High_dilution_patients <- unique(header.data$Patient[header.data$Dilution == "High" & header.data$Patient %in% tumor.samples$Patient])
  # Low dilution was used where available
  
  if(exploratory){
    
    used_raw_header <- header.data[(header.data$Patient %in% Low_dilution_patients & header.data$Dilution == "Low") | (header.data$Patient %in% High_dilution_patients[!(High_dilution_patients %in% Low_dilution_patients)] & header.data$Dilution == "High"),]
    #used_meta.data <- meta.data.raw[meta.data.raw$Metabolite %in% unique.mets & meta.data.raw$Peak_Number == "",]
    used_meta.data <- meta.data.raw
    
    write.table(TNmeta[rownames(TNmeta) %in% rownames(used_meta.data), colnames(TNmeta) %in% rownames(used_raw_header)], file = "Tables/all_totalTNmeta.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    write.table(used_raw_header[,colnames(used_raw_header) %in% c("Patient", "Sample", "Status", "Sample Number", "Replicate", "Day")], file = "Tables/all_totalheader.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    
    # add technical information for each ion - m/z, RT, ion formed
    
    ion_data <- read.delim('Tables/JKions.txt')
    
    library(data.table)
    
    new_ion_df <- ion_data[chmatch(rownames(used_meta.data), ion_data$X),]
    new_ion_df[is.na(new_ion_df$X), c(2:4)] <- used_meta.data[is.na(new_ion_df$X),c(1,2,4)]
    new_ion_df[is.na(new_ion_df$X), 1] <- rownames(used_meta.data)[is.na(new_ion_df$X)]
    rownames(new_ion_df) <- new_ion_df$X
    new_ion_df <- new_ion_df[,-1]
    
    write.table(new_ion_df, file = "Tables/all_totalmeta.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    
    write.table(TN.point.std[,(TN.subset.header$Patient %in% tumor.samples$Patient)], file = "Tables/all_stdTNmeta.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    # write.table(TN.subset.header, file = "Tables/stdheader.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    
  }else{
    
    used_raw_header <- header.data[(header.data$Patient %in% Low_dilution_patients & header.data$Dilution == "Low") | (header.data$Patient %in% High_dilution_patients[!(High_dilution_patients %in% Low_dilution_patients)] & header.data$Dilution == "High"),]
    #used_meta.data <- meta.data.raw[meta.data.raw$Metabolite %in% unique.mets & meta.data.raw$Peak_Number == "",]
    used_meta.data <- meta.data.raw[unique_matches,]
    
    
    write.table(TNmeta[rownames(TNmeta) %in% rownames(used_meta.data), colnames(TNmeta) %in% rownames(used_raw_header)], file = "Tables/totalTNmeta.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    write.table(used_raw_header[,colnames(used_raw_header) %in% c("Patient", "Sample", "Status", "Sample Number", "Replicate", "Day")], file = "Tables/totalheader.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    
    # add technical information for each ion - m/z, RT, ion formed
    
    ion_data <- read.delim('Tables/JKions.txt')
    
    library(data.table)
    
    new_ion_df <- ion_data[chmatch(rownames(used_meta.data), ion_data$X),]
    new_ion_df[is.na(new_ion_df$X), c(2:4)] <- used_meta.data[is.na(new_ion_df$X),c(1,2,4)]
    new_ion_df[is.na(new_ion_df$X), 1] <- rownames(used_meta.data)[is.na(new_ion_df$X)]
    rownames(new_ion_df) <- new_ion_df$X
    new_ion_df <- new_ion_df[,-1]
    
    write.table(new_ion_df, file = "Tables/totalmeta.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    
    write.table(TN.point.std[,(TN.subset.header$Patient %in% tumor.samples$Patient)], file = "Tables/stdTNmeta.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    # write.table(TN.subset.header, file = "Tables/stdheader.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
    
    
    
  }
  
}




#######

if(analysis.type == "small.molecules"){
  heatmap.maker.alldata(TN.point.std, TN.subset.header, "data_metabolomics/total_heatmap.pdf")
  heatmap.maker.alldata(TN.point.std, TN.subset.header, "data_metabolomics/total_heatmap_colun.pdf", reordercol = TRUE)
}
if(analysis.type == "lipomics"){
  heatmap.maker.alldata(TN.point.std, TN.subset.header, "data_lipomics/total_heatmap.pdf")
}





#######

if(bs.in.line == TRUE){
  
  measure.TNdiff <- TRUE
  measure.covariate.effect <- FALSE
  bsN <- 10^6
  
  use.covariates.forT <- FALSE
  
  tn.bs <- list()
  
  if(measure.TNdiff == TRUE){
    tn.bs[["tn"]] <- matrix(data = NA, ncol = 1, nrow = length(TN.subset.point[,1]))
    tn.bs[["tn-bs"]] <- matrix(data = NA, ncol = bsN, nrow = length(TN.subset.point[,1]))
  }
  if(measure.covariate.effect == TRUE){
    for(covar in cancer.cov){
      tn.bs[[paste(covar, "tn_p", sep = "_")]] <- matrix(data = NA, ncol = 1, nrow = length(TN.subset.point[,1]))
      tn.bs[[paste(covar, "tn_p-bs", sep = "_")]] <- matrix(data = NA, ncol = bsN, nrow = length(TN.subset.point[,1]))
    }}	
  
  
  if(measure.TNdiff == TRUE & use.covariates.forT == FALSE){
    for (i in 1:length(TN.subset.point[,1])){
      n.val <- normal.mat[i,!is.na(normal.mat[i,]) & !is.na(tumor.mat[i,])]
      t.val <- tumor.mat[i,!is.na(normal.mat[i,]) & !is.na(tumor.mat[i,])]
      
      if(length(n.val) > 5){
        
        sample.sqrt <- sqrt(length(n.val))
        test.stat <- mean(n.val - t.val)/(sd(n.val - t.val) / sample.sqrt)
        tn.bs[["tn"]][i,1] <- test.stat
        
        n.val.adj <- n.val - mean(n.val-t.val)
        cent.val <- mean(c(n.val.adj, t.val))
        n.val.resid <- n.val.adj - cent.val
        t.val.resid <- t.val - cent.val
        n.val.resid <- n.val.resid * sqrt(length(n.val) / (length(n.val) - 1))
        t.val.resid <- t.val.resid * sqrt(length(n.val) / (length(n.val) - 1))
        
        for (j in 1:bsN){
          
          resid.samples <- sample(c(1:length(n.val)), length(n.val), replace = TRUE)
          tn.bs[["tn-bs"]][i,j] <- mean(n.val.resid[resid.samples] - t.val.resid[resid.samples])/(sd(n.val.resid[resid.samples] - t.val.resid[resid.samples]) / sample.sqrt)
          
        }}
      print(paste("metabolite", i, "done", collapse = " "))
    }
    
    if(analysis.type == "small.molecules"){
      if(exploratory == FALSE){
        save(tn.bs, TN.point.std, meta.data, header.data, file = "data_metabolomics/tndiff-tn-bs.R")
      }else{
        save(tn.bs, TN.point.std, meta.data, header.data, file = "data_metabolomics/exp_tndiff-tn-bs.R")
      }
    }
    if(analysis.type == "lipomics"){
      if(exploratory == FALSE){
        save(tn.bs, file = "data_lipomics/tn-bs-lipids.R")
      }else{
        save(tn.bs, file = "data_lipomics/exp_tn-bs-lipids.R")
      }
    }	
    
  }
  
}else{
  
  if(analysis.type == "small.molecules"){
    if(exploratory == FALSE){
      load("data_metabolomics/tndiff-tn-bs.R")}else{
        load("data_metabolomics/exp_tndiff-tn-bs.R")
      }
  }
  if(analysis.type == "lipomics"){
    if(exploratory == FALSE){
      load("data_lipomics/tn-bs-lipids.R")}else{
        load("data_lipomics/exp_tn-bs-lipids.R")
      }
  }	
}


if(thorough == TRUE){
  #test for normality of the T-N differences
  KS_test_TN_diff(normal.mat, tumor.mat)
}

#FDR correction

TN_FDR_stats <- data.frame(fractionNull = rep(NA, times = 1), testStatCO = rep(NA, times = 1), FDRattained = rep(NA, times = 1), discoveries = rep(NA, times = 1))
rownames(TN_FDR_stats) = "TNcomp"

test.pvalue <- tn.bs[["tn"]]
bs.pvalue <- tn.bs[["tn-bs"]]
test.bs.p <- rep(NA, times = length(test.pvalue))

for (i in 1:length(test.pvalue)){
  valid.bs <- bs.pvalue[i,(!is.na(bs.pvalue[i,]) & !is.nan(bs.pvalue[i,]))]
  
  test.bs.p[i] <- 1 - (sum(abs(test.pvalue[i]) > abs(valid.bs))-1)/length(valid.bs)
}
test.bs.p <- sapply(test.bs.p, function(x){min(x, 1)}) #1 pseuo-count added to limit resolution of empirical null


nbs <- length(bs.pvalue[1,])

rm(tn.bs)

##### Determine the proportion of true positives (pi)

tn_fdr <- qvalue(test.bs.p)
plot(tn_fdr)

TN_FDR_stats$fractionNull <- tn_fdr$pi0


###### Optimize a value of the p-value cutoff that yields a given FDR
defined.mets <- !is.na(test.bs.p)
num.p <- sum(defined.mets)

TN_FDR_stats$testStatCO <- seq(0.0001, 1, by = 0.0001)[which.min(unlist(lapply(seq(0.0001, 1, by = 0.0001), FDR.eval, pi = TN_FDR_stats$fractionNull, num.p = num.p, test.stats = test.bs.p, FDRlev = 0.05, use.unif = TRUE)))]
TN_FDR_stats$FDRattained <- FDR.calc(pi = TN_FDR_stats$fractionNull, lambda = TN_FDR_stats$testStatCO, num.p = num.p, test.stats = test.bs.p[defined.mets], use.unif = TRUE)


Discovery_pvalues <- test.bs.p[defined.mets][test.bs.p[defined.mets] < TN_FDR_stats$testStatCO]
Discovery_qvalues <- rep(NA, times = length(Discovery_pvalues))
for (i in 1:length(Discovery_pvalues)){
  Discovery_qvalues[i] <- (TN_FDR_stats$fractionNull*Discovery_pvalues[i]*num.p) / ((sum(test.bs.p[defined.mets] <= Discovery_pvalues[i])))
}
Discovery_qvalues <- signif(Discovery_qvalues, digits = 5)
all_qvalues <- unlist(lapply(c(1:length(test.bs.p)), function(x){
  (TN_FDR_stats$fractionNull*test.bs.p[x]*num.p) / (sum(test.bs.p[defined.mets] <= test.bs.p[x]))
}))

Discoveries <- data.frame(Discoveries = rownames(TN.point.std)[defined.mets][test.bs.p[defined.mets] < TN_FDR_stats$testStatCO], p_values = Discovery_pvalues, q_values = Discovery_qvalues)

TN_FDR_stats$discoveries <- length(Discoveries[,1])

########## FINAL OUTPUT PLOTS ################

if(analysis.type == "lipomics"){
  write.table(Discoveries, file = "data_lipomics/TNdiscoveries.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  barplot.maker(Discoveries, TN.subset.header, meta.data, TN.point.std, "data_lipomics/lipid_significance_barplot.pdf", "data_lipomics/TN_lipid_fold_change.csv")
  heatmap.makerTN(Discoveries, TN.subset.header, TN.point.std, "data_lipomics/TN_lipid_sig_diff.pdf")
  discovery.diff.plot(meta.data, Discoveries, tumor.mat, normal.mat, "data_lipomics/diffplot.pdf")
  boxplot.prep(Discoveries, tumor.mat, normal.mat, "data_lipomics/boxplot_prep.R")
  load("data_lipomics/boxplot_prep.R")
  
  pathway_col <- data.frame(Pathway = c("non-Lysophospholipid", "Lysophospholipid"), Colorz = c("blue", "chartreuse"), stringsAsFactors = FALSE)
  pathway.sort.order <- c("non-Lysophospholipid", "Lysophospholipid")
  pathway_col <- pathway_col[sapply(c(1:length(pathway.sort.order)), function(x){c(1:length(pathway.sort.order))[pathway_col$Pathway %in% pathway.sort.order[x]]}),]
  
  all_disc_annotate <- data.frame(Discoveries = Discoveries$Discoveries, Palatable.Name = sapply(names(save.sig.diff), function(x){unlist(strsplit(x, split = "-[[:digit:]]-"))[1]}), Pathway = ifelse(c(1:length(names(save.sig.diff))) %in% grep("^L", names(save.sig.diff)), "Lysophospholipid", "non-Lysophospholipid"), stringsAsFactors = FALSE)
  
  subset.data <- data.frame(internal_name = names(save.sig.diff), palatable_name = all_disc_annotate$Palatable.Name, color = unlist(lapply(all_disc_annotate$Pathway, function(x){pathway_col$Colorz[pathway_col$Pathway == x]})), stringsAsFactors = FALSE)
  
  Discoveries[,1] <- as.character(Discoveries[,1])
  
  median_diffs <- unlist(lapply(c(1:length(Discoveries$Discoveries)), function(x){median(save.sig.diff[[Discoveries$Discoveries[x]]], na.rm = TRUE)}))
  
  reordering_path <- NULL
  for(k in 1:length(pathway.sort.order)){
    reordering_path <- c(reordering_path, c(1:length(all_disc_annotate$Pathway))[all_disc_annotate$Pathway %in% pathway.sort.order[k]][order(median_diffs[c(1:length(all_disc_annotate$Pathway))[all_disc_annotate$Pathway %in% pathway.sort.order[k]]])])
  }
  
  
  TN.boxplot(save.sig.diff, subset.data, "data_lipomics/discBP1.pdf", usepval = TRUE, pvals = Discoveries$q_values, legend_spec = list(text = pathway_col$Pathway, color = pathway_col$Colorz), subset.index = order(median_diffs))
  
  
  TN.boxplot(save.sig.diff, subset.data, "data_lipomics/discBP2.pdf", usepval = TRUE, pvals = Discoveries$q_values, legend_spec = list(text = pathway_col$Pathway, color = pathway_col$Colorz), subset.index = reordering_path)
  
  
  
  
  
}





if(analysis.type == "small.molecules"){
  
  #### make lots of summary plots comparing matched tumor and benign samples ####  
  
  write.table(Discoveries, file = "data_metabolomics/TNdiscoveries.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  barplot.maker(Discoveries, TN.subset.header, meta.data, TN.point.std, "data_metabolomics/met_significance_barplot.pdf", "data_metabolomics/TN_fold_change.csv")
  heatmap.makerTN(Discoveries, TN.subset.header, TN.point.std, "data_metabolomics/TN_sig_diff.pdf")
  discovery.diff.plot(meta.data, Discoveries, tumor.mat, normal.mat, "data_metabolomics/diffplot.pdf")
  boxplot.prep(Discoveries, tumor.mat, normal.mat, "data_metabolomics/boxplot_prep.R")
  
  load("data_metabolomics/boxplot_prep.R")
  
  ##### Boxplot all discoveries #######
  
  #boxplot.prep(Discoveries, tumor.mat, normal.mat, file.name = "data_metabolomics/discbp")
  #load("data_metabolomics/discbp")
  
  all_disc_annotate <- read.delim("data_metabolomics/disc_colors.txt", stringsAsFactors = FALSE)
  
  Discoveries$Discoveries[!(Discoveries$Discoveries %in% all_disc_annotate$Discoveries)]
  
  Discoveries_renamed <- Discoveries
  Discoveries_renamed$Discoveries <- sapply(Discoveries_renamed$Discoveries, function(x){all_disc_annotate$Palatable.Name[all_disc_annotate$Discoveries == x]})
  write.table(Discoveries_renamed, file = "data_metabolomics/TNdiscoveriesRenamed.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  all_disc_annotate <- all_disc_annotate[all_disc_annotate$Discoveries %in% Discoveries$Discoveries,]
  
  #pathway_col <- data.frame(Pathway = as.character(unique(all_disc_annotate$Pathway)), Colorz = c("cornflowerblue", "green3", "khaki3", "cadetblue2", "plum1", "firebrick1"), stringsAsFactors = FALSE)
  
  pathway.sort.order <- c("Non-essential AA", "Essential AA", "AA related", "Nucleotide related", "Central carbon metabolite", "Misc") 
  pathway.colors <- c("green3", "firebrick1", "khaki3", "cadetblue2", "plum1", "cornflowerblue")
  pathway_col <- data.frame(Pathway = pathway.sort.order, Colorz = pathway.colors, stringsAsFactors = FALSE)
  
  ###### rewrite boxplot code #####
  
  pathway_col$nmembers <- sapply(pathway_col$Pathway, function(x){sum(all_disc_annotate$Pathway == x)})
  pathway_col$xoffset <- (c(1:nrow(pathway_col)) - 1)*1
  pathway_col$xcoord <- cumsum(pathway_col$nmembers) - pathway_col$nmembers/2 + pathway_col$xoffset + 1
  
  # split pathway labeling position so that long names are in two pieces
  
  pathway_colPlot <- pathway_col
  pathway_colPlot$Display <- pathway_colPlot$Pathway
  pathway_colPlot$y <- 4
  
  for(pw in c("Non-essential AA", "Nucleotide related", "Central carbon metabolite")){
    words <- strsplit(pathway_colPlot$Pathway[pathway_colPlot$Pathway == pw], split = ' ')[[1]]
    l1 <- pathway_colPlot[pathway_colPlot$Pathway == pw,]
    l2 <- pathway_colPlot[pathway_colPlot$Pathway == pw,]
    l1$Display <-  paste(words[-length(words)], collapse = " "); l1$y = 4
    l2$Display <-  words[length(words)]; l2$y = 3.8
    
    pathway_colPlot[pathway_colPlot$Pathway == pw,] <- l1
    pathway_colPlot <- rbind(pathway_colPlot, l2)
  }
  
  
  
  
  boxplotStack <- NULL
  
  for(i in 1:nrow(Discoveries)){
    
    boxplotStack <- rbind(boxplotStack, data.frame(metabolite = Discoveries$Discoveries[i], Change = save.sig.diff[[i]]))
    
  }
  
  boxplotStack <- as.data.table(boxplotStack[!is.na(boxplotStack$Change),])
  boxplotStack[,medianChange := median(Change), by = metabolite]
  boxplotStack <- boxplotStack[order(medianChange),]
  
  boxplotStack$renamedMet <- sapply(boxplotStack$metabolite, function(x){all_disc_annotate$Palatable.Name[all_disc_annotate$Discoveries == x]})
  boxplotStack$Pathway_string <- sapply(boxplotStack$metabolite, function(x){all_disc_annotate$Pathway[all_disc_annotate$Discoveries == x]})
  boxplotStack$Pathway <- factor(boxplotStack$Pathway_string, levels = pathway_col$Pathway)
  
  boxplotStack$FILL <- sapply(boxplotStack$Pathway, function(x){pathway_col$Colorz[pathway_col$Pathway == x]})
  boxplotStack <- boxplotStack[order(boxplotStack$Pathway, boxplotStack$medianChange),]
  
  boxplotStack$renamedMet_string <- boxplotStack$renamedMet
  boxplotStack$renamedMet <- factor(boxplotStack$renamedMet, levels = unique(boxplotStack$renamedMet))
  
  boxplotStack[,xoffset := pathway_col$xoffset[pathway_col$Pathway == as.character(Pathway_string)], by = Pathway_string]
  boxplotStack[,xmet := c(1:nrow(Discoveries))[renamedMet_string == levels(renamedMet)], by = renamedMet_string]
  boxplotStack[,xpos := xoffset + xmet,]
  
  
 
  
  
  boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), legend.position = "top",
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 14, angle = 90, hjust = 1, face = "plain"), axis.line = element_blank(),
                         axis.text = element_text(color = "black"))
  
  
  TN_boxplot <- ggplot()
  TN_boxplot + geom_boxplot(data = boxplotStack, aes(x = xpos, y = Change, group = renamedMet, fill = FILL), outlier.size = 0) + geom_hline(yintercept = 0, size = 1, col = "darkgray") +
    scale_y_continuous(expression(Tumor / "Benign Adjacent"), breaks = seq(-4,4), labels = 2^seq(-4,4), lim = c(-4,4)) + scale_x_discrete("Metabolites", breaks = unique(boxplotStack$xpos), labels = unique(boxplotStack$renamedMet), expand = c(0.02,0.02)) +
    boxplot_theme + geom_text(data = pathway_colPlot, aes(label = Display, color = Colorz, x = xcoord, y = y), fontface = "bold", size = 6.5) + scale_color_identity() + scale_fill_identity()
  
  ggsave("data_metabolomics/TN_boxplot.pdf", width = 14, height = 12)
  
  boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), legend.position = "top",
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 13, angle = 70, hjust = 1, face = "plain"), axis.line = element_blank(),
                         axis.text = element_text(color = "black"))
  
  TN_boxplot <- ggplot()
  TN_boxplot + geom_boxplot(data = boxplotStack, aes(x = xpos, y = Change, group = renamedMet, fill = FILL), outlier.size = 0) + geom_hline(yintercept = 0, size = 1, col = "darkgray") +
    scale_y_continuous(expression(Tumor / "Benign Adjacent"), breaks = seq(-4,4), labels = 2^seq(-4,4), lim = c(-4,4)) + scale_x_discrete("Metabolites", breaks = unique(boxplotStack$xpos), labels = unique(boxplotStack$renamedMet), expand = c(0.02,0.02)) +
    boxplot_theme + geom_text(data = pathway_colPlot, aes(label = Display, color = Colorz, x = xcoord, y = y), fontface = "bold", size = 6.5) + scale_color_identity() + scale_fill_identity()
  ggsave("data_metabolomics/TN_boxplot_angle.pdf", width = 14, height = 12)
  
  
  
  ##########
  
  #pathway_col <- pathway_col[sapply(c(1:length(pathway.sort.order)), function(x){c(1:length(pathway.sort.order))[pathway_col$Pathway %in% pathway.sort.order[x]]}),]
  
  
  subset.data <- data.frame(internal_name = names(save.sig.diff), palatable_name = all_disc_annotate$Palatable.Name, color = unlist(lapply(all_disc_annotate$Pathway, function(x){pathway_col$Colorz[pathway_col$Pathway == x]})), stringsAsFactors = FALSE)
  
  #subset.data <- data.frame(internal_name = names(save.sig.diff), palatable_name = c("Glutamine", "Glycine", "Kynurenine"), color = "orange", stringsAsFactors = FALSE)
  #TN.boxplot(save.sig.diff, subset.data, "data_metabolomics/discBP1.pdf", usepval = TRUE, pvals = Discoveries$q_values)
  
  #all discoveries, sorted alphabetically
  
  TN.boxplot(save.sig.diff, subset.data, "data_metabolomics/discBP1.pdf", usepval = TRUE, pvals = Discoveries$q_values, legend_spec = list(text = pathway_col$Pathway, color = pathway_col$Colorz))
  
  #all discoveries, sorted overall by median
  
  median_diffs <- unlist(lapply(c(1:length(Discoveries$Discoveries)), function(x){median(save.sig.diff[[Discoveries$Discoveries[x]]], na.rm = TRUE)}))
  
  TN.boxplot(save.sig.diff, subset.data, "data_metabolomics/discBP2.pdf", usepval = TRUE, pvals = Discoveries$q_values, legend_spec = list(text = pathway_col$Pathway, color = pathway_col$Colorz), subset.index = order(median_diffs))
  
  #all discoveries, sorted first by category, then by median
  
  reordering_path <- NULL
  for(k in 1:length(pathway.sort.order)){
    reordering_path <- c(reordering_path, c(1:length(all_disc_annotate$Pathway))[all_disc_annotate$Pathway %in% pathway.sort.order[k]][order(median_diffs[c(1:length(all_disc_annotate$Pathway))[all_disc_annotate$Pathway %in% pathway.sort.order[k]]])])
  }
  
  TN.boxplot(save.sig.diff, subset.data, "data_metabolomics/discBP3.pdf", usepval = TRUE, pvals = Discoveries$q_values, legend_spec = list(text = pathway_col$Pathway, color = pathway_col$Colorz), subset.index = reordering_path)
  
  TN.boxplot(save.sig.diff, subset.data, "data_metabolomics/discBP3_minimal.pdf", usepval = TRUE, subset.index = reordering_path)
  
  Discoveries_out <- Discoveries; Discoveries_out$Discoveries <- all_disc_annotate$Palatable.Name; colnames(Discoveries_out) <- c("Discoveries", "p-value", "q-value")
  library(xtable)
  print(xtable(Discoveries_out[reordering_path,], digits = 5), include.rownames= FALSE)
  
  ##### determine 
  
  fuconate = tumor.mat[meta.data[,1] == "D-Mannose/L-Sorbose-potential-2",] - normal.mat[meta.data[,1] == "D-Mannose/L-Sorbose-potential-2",]
  fuconate <- fuconate[!is.na(fuconate)]
  gln = tumor.mat[meta.data[,1] == "glutamine",] - normal.mat[meta.data[,1] == "glutamine",]
  gln <- gln[!is.na(gln)]
  
  library(ggplot2)
  library(scales) 
  
  scatter_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), legend.position = "right", 
                         panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan"),
                         axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  
  plot_breaks <- c(2^seq(-5, 5))
  plot_breaks_labels <- sapply(plot_breaks, function(x){as.character(x)})
  fuconate_plot = ggplot(data = data.frame(relative_abund = 2^fuconate), aes(y = relative_abund, x = 1))
  fuconate_plot + geom_boxplot(fill = "coral1") + scale_y_continuous("Tumor:Normal", trans = log2_trans(), limits = c(1/20, 20), breaks = plot_breaks, labels = plot_breaks_labels) + scatter_theme
  ggsave("data_metabolomics/fuconateBP.pdf", width = 6, height = 10)
  fuconate_plot + geom_violin(fill = "coral1") + scale_y_continuous("Tumor:Normal", trans = log2_trans(), limits = c(1/20, 20), breaks = plot_breaks, labels = plot_breaks_labels) + scatter_theme
  ggsave("data_metabolomics/fuconateVP.pdf", width = 6, height = 10)
  
  
  names(sort(gln))[round(length(gln)*seq(0.1,0.9, by = 0.1))]
  names(sort(fuconate))[round(length(fuconate)*seq(0.1,0.9, by = 0.1))]
  
}


