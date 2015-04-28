


bs.with.missing.data <- function(fitted.resid, sample.thresh = 7){

	# sample residuals from a matrix with NaNs.  Only sample rows with > sample.thresh non-NaN members and do so with a probability proportional to the number of samples in a row.
	residuals <- fitted.resid
	non.nan.points <- length(fitted.resid[!is.nan(fitted.resid)])
	
	num.valid.samples <- ifelse(num.samples >= 7, num.samples, 0)
	prob.row <- num.valid.samples/sum(num.valid.samples)

	row.samp <- sample(c(1:length(fitted.resid[,1])), non.nan.points, replace = TRUE, prob = prob.row)
	
	residuals[!is.nan(residuals)] <- apply(fitted.resid[row.samp,], 1, sample.col)
	residuals
	}


sample.col <- function(a.row){
		sample(a.row[!is.nan(a.row)], 1) * length(a.row[!is.nan(a.row)])/(length(a.row[!is.nan(a.row)]) - 1)
		
		}

half.life.optim.nomiss <- function(tau, vec.tI, vec.tF, time.decay){

#tau is actually lambda (1/tau) for conventional half-life exponential decay.

#optimiztion function for half-life (tau).  minimizing the difference between the final abundance and the initial abundance after decaying for time.decay with half-life tau.
		
#taking the sum of square diff rather than the sum of absolute values resulted in overvaluing of single measurements		
		
	sum((vec.tF - vec.tI*exp(-time.decay*tau))^2)
	
	}







half.life.optim <- function(tau, vec.tI, vec.tF, missing, time.decay){

#tau is actually lambda (1/tau) for conventional half-life exponential decay.

#optimization function for half-life (tau).  minimizing the difference between the final abundance and the initial abundance after decaying for time.decay with half-life tau.
		
	sum((vec.tF[!missing] - vec.tI[!missing]*exp(-time.decay[!missing]*tau))^2)
	
	}



all.sample.corr <- function(tall.data){
	
#outputs spearman and pearson correlations for all pairwise comparisons of columns for a dataset as well as the number of comparisons.  Extract with $Spearman, $Pearson and $SampleSize
	
	mat.dim <- length(tall.data[1,])
	header.names <- colnames(tall.data)
	
	pearson.corr.mat <- matrix(data = NA, ncol = mat.dim, nrow = mat.dim)
	colnames(pearson.corr.mat) <- header.names
	rownames(pearson.corr.mat) <- header.names
	spearman.corr.mat <- matrix(data = NA, ncol = mat.dim, nrow = mat.dim)
	colnames(spearman.corr.mat) <- header.names
	rownames(spearman.corr.mat) <- header.names
	sampleSize.mat <- matrix(data = NA, ncol = mat.dim, nrow = mat.dim)
	colnames(sampleSize.mat) <- header.names
	rownames(sampleSize.mat) <- header.names
	
	for (i in 1:mat.dim){
	for (j in 1:mat.dim){	
		
		comp.mat <- tall.data[,c(i,j)]
		comp.mat <- comp.mat[apply(is.na(comp.mat), 1, sum) == 0,]
		pearson.corr.mat[i,j] <- cor(comp.mat[,1], comp.mat[,2], method = "pearson")
		spearman.corr.mat[i,j] <- cor(comp.mat[,1], comp.mat[,2], method = "spearman")
		sampleSize.mat[i,j] <- length(comp.mat[,1])
		}}
		
		corr.out <- list()
		corr.out$Pearson <- pearson.corr.mat
		corr.out$Spearman <- spearman.corr.mat
		corr.out$SampleSize <- sampleSize.mat
		
		corr.out
		
		}


take.the.first <- function(vec){
	order(vec == 1, decreasing = TRUE)[1]
	}





calculate.CV.curve <- function(TNmeta, header.data, meta.data, threshold){
	
	#create a list relating the signal/noise of each MS at different peak intensities
	
	sample.combos <- sample.combos.fxn(header.data, min.is.two = TRUE, no.blanks = TRUE)
	
	CV.list <- list(Exactive = NULL, Ultra = NULL, Max = NULL)
	
	for (j in 1:length(sample.combos[1,])){
		rep.pairs <- TNmeta[,sample.combos[,j] == 1]
		nreps <- length(rep.pairs[1,])
		#take sets of replicates with no NAs
		m.data.pairs <- meta.data[apply(is.na(rep.pairs), 1, sum) == 0,]
		rep.pairs <- rep.pairs[apply(is.na(rep.pairs), 1, sum) == 0,]
		
		row.var <- apply(rep.pairs, 1, var)*(nreps/(nreps-1))
		CV <- apply(rep.pairs, 1, mean)/sqrt(row.var)
		mean.cv <- cbind(apply(rep.pairs, 1, mean), CV)
		colnames(mean.cv) <- c("IC", "CV")
		
		CV.list$Exactive <- rbind(CV.list$Exactive, mean.cv[m.data.pairs$Instrument == "Exactive",])
		CV.list$Ultra <- rbind(CV.list$Ultra, mean.cv[m.data.pairs$Instrument == "Ultra",])
		CV.list$Max <- rbind(CV.list$Max, mean.cv[m.data.pairs$Instrument == "Max",])
		
		}
	
	CV.list}
	
	
if.match.write <- function(string, rem.terms){
	for (te in rem.terms){
	if(length(grep(te, string, perl = TRUE)) != 0){
		string <- unlist(strsplit(string, te))[1]
		}}
		string}
		
		
cv.calc <- function(meta, j, poly.fit.e, fit.degree){
#predict the coefficient of variation (CV) as a function of peak abundance based upon the variance of replicates of different sizes of peaks 

apply((matrix(log(meta[,j][!is.na(meta[,j])], base = 2), ncol = fit.degree +1, nrow = length(meta[,j][!is.na(meta[,j])]))^matrix(0:fit.degree,  ncol = fit.degree +1, nrow = length(meta[,j][!is.na(meta[,j])]), byrow = TRUE))*matrix(poly.fit.e,  ncol = fit.degree +1, nrow = length(meta[,j][!is.na(meta[,j])]), byrow = TRUE), 1, sum)

}




fit.to.mean.abund <- function(meta.floor, poly.fit.e, fit.degree, just.median = FALSE){
row.means <- apply(meta.floor, 1, mean, na.rm = TRUE)
row.median <- apply(meta.floor, 1, median, na.rm = TRUE)
scaling.factors <- rep(NA, length(meta.floor[1,]))
for (j in 1:length(meta.floor[1,])){
	if(just.median == TRUE){
#the median of the fractional difference between a samples abundance and the mean abundance for metabolites over the threshold
	scaling.factors[j] <- median(row.median[!is.na(meta.floor[,j])]/meta.floor[,j][!is.na(meta.floor[,j])])
	} else { 
#fractional difference between a samples abundance and the mean abundance for metabolites over the threshold weighted by the fitted signal:noise (the CV)
	cv.j <- 2^cv.calc(meta.floor[!is.na(meta.floor[,j]),], j, poly.fit.e, fit.degree)
		
	scaling.factors[j] <- sum((row.means[!is.na(meta.floor[,j])]/meta.floor[,j][!is.na(meta.floor[,j])])*cv.j)/sum(cv.j)
	}}
	scaling.factors
	}		




sample.combos.fxn <- function(header.data, min.is.two = FALSE, no.blanks = TRUE){
	if(no.blanks == FALSE){
		header.data.red <- header.data
		} else {
		header.data.red <- header.data[header.data$Status %in% c("Normal", "Tumor"),]		}
	
	unique.samples <- unique(paste(header.data.red$Patient, header.data.red$Status, header.data.red$Dilution, sep = ""))

	sample.combos <- matrix(data = NA, ncol = length(unique.samples), nrow = length(header.data[,1]))
	for (i in 1:length(unique.samples)){
		sample.combos[,i] <- ifelse(paste(header.data$Patient, header.data$Status, header.data$Dilution, sep = "") %in% unique.samples[i], 1, 0)
		}
		if(min.is.two == TRUE){
		sample.combos <- sample.combos[,apply(sample.combos, 2, sum) > 1]
		}
		sample.combos
	}

sample.combos.fxn2 <- function(header.data, min.is.two = FALSE, no.blanks = TRUE, first.timecourse = FALSE, timecourse.valid = NA, unique.terms){
	if(no.blanks == FALSE){
		header.data.red <- header.data
		} else {
		header.data.red <- header.data[header.data$Status %in% c("Normal", "Tumor"),]		}
			if("exo-vivo Timecourse" %in% unique.terms){
			header.data.red$"exo-vivo Timecourse" <- as.character(header.data.red$"exo-vivo Timecourse")}
		
	if(first.timecourse == TRUE){
		header.data.red <- header.data.red[is.na(header.data.red$"exo-vivo Timecourse") | (header.data.red$"exo-vivo Timecourse" %in% timecourse.valid),]
		}
	
	unique.samples <- unique(apply(header.data.red[colnames(header.data) %in% unique.terms], 1, paste, collapse = ""))

	sample.combos <- matrix(data = NA, ncol = length(unique.samples), nrow = length(header.data[,1]))
	for (i in 1:length(unique.samples)){
		sample.combos[,i] <- ifelse(apply(header.data[,colnames(header.data) %in% unique.terms], 1, paste, collapse = "") %in% unique.samples[i], 1, 0)
		}
		if(min.is.two == TRUE){
		sample.combos <- sample.combos[,apply(sample.combos, 2, sum) > 1]
		}
		sample.combos
	}


FDR.eval <- function(pi, lambda, test.stats, num.p, FDRlev = 0.1, use.unif = FALSE, use.perms = NULL){
	
	if(use.unif == TRUE){
		n.fp <- pi*lambda*num.p
		}
	if(use.unif == FALSE){
		n.fp <- pi*sum(apply(use.perms <= lambda, 1, sum)) / (length(use.perms[1,]))
		}	
		abs((n.fp / length(test.stats[test.stats <= lambda])) - FDRlev)}




FDR.calc <- function(pi, lambda, test.stats, num.p, use.unif = FALSE, use.perms = NULL){
	
	if(use.unif == TRUE){
		n.fp <- pi*lambda*num.p
		}
	if(use.unif == FALSE){
		n.fp <- pi*sum(apply(use.perms <= lambda, 1, sum)) / (length(use.perms[1,]))
		}	
		n.fp / length(test.stats[test.stats <= lambda])}
	
brok.stick.var <- function(cumvar){
### Determine the expected variance explained by principle components based upon the length of a PCA eigenvalue vector 
	
	p = length(cumvar)
	brokstick <- data.frame(prcomp = c(1:p), lk = rep(NA, times = p))
	for(i in 1:p){
	jsum <- 0
	for(j in i:p){
	jsum <- jsum + (1/j)
	}
	brokstick[i,2] <- jsum*(1/p)
	}
	brokstick
	}
		
		
small.mol.tables = function(TNmeta){
####### Small molecule table input #########

data.list = list()	
########## collapsing the sample information
header.data <- as.data.frame(t(TNmeta[1:9,]))
colnames(header.data) <- unlist(lapply(header.data[1,], as.character))
header.data <- header.data[-c(1:3),]
header.data[header.data == ""] <- NA
header.data$Patient <- sub(22210, 2210, header.data$Patient)
header.data$Patient <- sub(2210, 22210, header.data$Patient)
header.data$Patient <- sub("22210_", "22210-", header.data$Patient)
rownames(header.data) <- paste(header.data$Patient, header.data$Sample, sep = "-")


########### collapsing instrument/metabolite information

meta.data <- TNmeta[-c(2:9), -c(4:length(TNmeta[1,]))]
meta.data <- cbind(meta.data, NA)
colnames(meta.data) <- c("Metablite_name", "Instrument", "Peak_Number", "Metabolite")
meta.data <- meta.data[-1,]

rem.terms <- c("_Xact", "_QQQ+", "_QQQ-", "-nega", "_nega", "_posi", "-posi")

for (i in 1:length(meta.data[,1])){
	meta.data$Metabolite[i] <- if.match.write(meta.data[i,1], rem.terms)}

meta.name <- rep(NA, times = length(meta.data$Metabolite))
for (i in 1:length(meta.data$Metabolite)){
	if(meta.data$Peak_Number[i] == ""){meta.name[i] <- paste(meta.data$Metabolite[i], meta.data$Instrument[i], sep = "_")}else{
		meta.name[i] <- paste(meta.data$Metabolite[i], meta.data$Instrument[i], meta.data$Peak_Number[i], sep = "_")
		}}
rownames(meta.data) <- meta.name
	
######### Retag data table

colnames(TNmeta) <- c(TNmeta[1,1:3], rownames(header.data))
rownames(TNmeta) <- c(TNmeta[1:9,1], rownames(meta.data))
TNmeta <- TNmeta[-c(1:9),]
TNmeta <- TNmeta[,-c(1:3)]
TNmeta[TNmeta == ""] <- NA
TNmeta <- apply(TNmeta, c(1,2), as.numeric)

data.list$header.data <- header.data
data.list$meta.data <- meta.data
data.list$TNmeta <- TNmeta

data.list

	}
	
lipomics.tables <- function(TNmeta){

data.list <- list()
	
header.data <- as.data.frame(t(TNmeta[1:5,]), stringtofactors = FALSE)
colnames(header.data) <- c("Patient", "Dilution", "exo-vivo Timecourse", "Status", "Sample")	
header.data <- header.data[-1,]	
header.data$Patient <- sub(22210, 2210, header.data$Patient)
header.data$Patient <- sub(2210, 22210, header.data$Patient)
header.data$Dilution <- sub("L", "Low", header.data$Dilution); header.data$Dilution <- sub("H", "High", header.data$Dilution); header.data$Dilution[header.data$Dilution == ""] <- "Low"
header.data$"exo-vivo Timecourse"[header.data$"exo-vivo Timecourse"==""] <- NA 	
header.data$Status <- sub("N", "Normal", header.data$Status); header.data$Status <- sub("T", "Tumor", header.data$Status)
header.data <- apply(header.data, 2, as.character)

rownames(header.data) <- apply(cbind(header.data[,colnames(header.data) == "Patient"], header.data[,colnames(header.data) == "Status"], header.data[,colnames(header.data) == "Sample"]), 1, function(x){paste(x, collapse = "")})
	
	
	
metabolites <- TNmeta[-c(1:5), -c(2:length(TNmeta[1,]))]

meta.data <- data.frame("Metablite_name" = metabolites, "Instrument" = rep("E", times = length(metabolites)), "Peak_Number" = rep(NA, times = length(metabolites)), "Metabolite" = unlist(lapply(metabolites, function(x){paste(c(unlist(strsplit(x, ")"))[1], ")"), collapse ="")})), stringsAsFactors = FALSE)

unique.met <- unique(meta.data$"Metabolite")
for(i in 1:length(unique.met)){
	
	sub.met <- meta.data[meta.data$Metabolite == unique.met[i],]
	multi.peak <- grep('-[0-9]+', sub.met$"Metablite_name")
	peak.assignment <- rep(NA, times = length(sub.met[,1]))
	
	peak.assignment[multi.peak] <- unlist(lapply(sub.met[multi.peak]$"Metablite_name", function(x){unlist(strsplit(x, '-'))[length(unlist(strsplit(x, '-')))]}))
	peak.assignment <- as.numeric(peak.assignment)
	if(length(unique(peak.assignment)) != length(peak.assignment)){peak.assignment[!is.na(peak.assignment)] <- NA}
	
	peak.assignment[is.na(peak.assignment)] <- c(1:length(sub.met[,1]))[!(c(1:length(sub.met[,1])) %in% peak.assignment)]
	
	meta.data$"Peak_Number"[meta.data$Metabolite == unique.met[i]] <- peak.assignment
	}
	
rownames(meta.data) <- apply(meta.data[,c(4,3,2)], 1, function(x){paste(x, collapse = "-")}) 


colnames(TNmeta) <- c(TNmeta[1,1], rownames(header.data))
rownames(TNmeta) <- c(TNmeta[1:5,1], rownames(meta.data))
TNmeta <- TNmeta[-c(1:5),]
TNmeta <- TNmeta[,-1]
TNmeta[TNmeta == ""] <- NA
TNmeta <- apply(TNmeta, c(1,2), as.numeric)
header.data <- data.frame(header.data, stringsAsFactors = FALSE)

data.list$header.data <- header.data
data.list$meta.data <- meta.data
data.list$TNmeta <- TNmeta

data.list
	
	}
	
read.covariates <- function(){
	
######### Read in epidemiological data

epi.data <- read.delim("patientInfo/epi_risk.txt")
epi2.data <- read.delim("patientInfo/diseaseinfo.txt")
surgery.data <- read.delim("patientInfo/sampleAq.txt")
ca19levels <- read.delim("patientInfo/tumor_markers_10.5.11.txt")

weight <- epi.data$WT
names(weight) <- epi.data$SPONSORSUBJECTID

age <- epi2.data$age
names(age) <- epi2.data$SPONSORSUBJECTID
#get an approximate age from birthdate, where surgery date is not provided
ageTOdate <- rep(NA, times = length(age))
default.date <- as.Date("11/1/2011", "%m/%d/%Y")
for(i in 1:length(age)){
	d.m.y <- as.numeric(unlist(strsplit(as.character(epi2.data$BIRTHDT[i]), split = "/")))
	d.m.y[3] <- paste("19", d.m.y[3], sep ="")
	b.day <- as.Date(paste(d.m.y, collapse = "/"), "%m/%d/%Y")
	ageTOdate[i] <- round(difftime(default.date, b.day, units = "days")/365.25, digits = 2)
	}
age[is.na(age)] <- ageTOdate[is.na(age)]

sex <- as.character(epi2.data$SEX)
names(sex) <- epi2.data$SPONSORSUBJECTID

path_class <- as.character(epi2.data$PATHRESULTS)
names(path_class) <- epi2.data$SPONSORSUBJECTID


sample.covariates <- cbind(data.frame(unique(header.data[,colnames(header.data) %in% "Patient"]), stringsAsFactors = FALSE), NA, NA, NA, NA)
colnames(sample.covariates) <- c("PatientID", "Weight", "Age", "Sex", "Pathology")

for (i in 1:length(sample.covariates[,1])){
	sample.weight <- weight[names(weight) %in% sample.covariates$"PatientID"[i]]
	sample.age <- age[names(age) %in% sample.covariates$PatientID[i]]
	sample.sex <- sex[names(sex) %in% sample.covariates$PatientID[i]]
	sample.path <- path_class[names(path_class) %in% sample.covariates$PatientID[i]]
	
	sample.covariates$Weight[i]  <-ifelse(length(sample.weight) == 0, NA, sample.weight)
	sample.covariates$Age[i] <- ifelse(length(sample.age) == 0, NA, sample.age)
	sample.covariates$Sex[i] <- ifelse(length(sample.sex) == 0, NA, sample.sex)
	sample.covariates$Pathology[i] <- ifelse(length(sample.path) == 0, NA, sample.path)
	}
sample.covariates[sample.covariates == ""] <- NA

sample.covariates

}
		
		
TN.subsetfxn <- function(TN.header, unique.terms){

valid.s <- rep(NA, times = length(TN.header[,1]))
#only take non-Low sample concentrations when other samples aren't available" 
reduced.header <- apply(TN.header[,colnames(TN.header) %in% unique.terms], 1, function(x){paste(x, collapse="")})
for(i in 1:length(unique(reduced.header))){
	sampleZ <- reduced.header %in% unique(reduced.header)[i]
	
	if("Low" %in% TN.header[sampleZ,]$Dilution){
		valid.s[sampleZ] <- ifelse(TN.header[sampleZ,]$Dilution == "Low", TRUE, FALSE)
		}else{
			valid.s[sampleZ] <- TRUE
			}}
		valid.s
	}
	
	

heatmap.maker.alldata <- function(TN.point.std, TN.subset.header, filename, reordercol = FALSE){

TN.subset <- TN.point.std[apply(!is.na(TN.point.std), 1, sum) > length(TN.point.std[1,])*0.6,]
colnames(TN.subset) <- TN.subset.header$Patient

pdf(file = filename)
heatmap.2(TN.subset[,order(TN.subset.header$Status)], Colv = reordercol, dendrogram = ifelse(reordercol, "both", "row"), cexRow = 0.3, col = blue2yellow(50), trace = "n", ColSideColors = ifelse(TN.subset.header[order(TN.subset.header$Status),]$Status == "Normal", "BLUE", "RED"), na.color = "darkgrey", symkey = TRUE)
legend("topright", c("Normal", "Tumor"), text.col = c("BLUE", "RED"))
dev.off()
	}




heatmap.makerTN <- function(Discoveries, TN.subset.header, TN.point.std, filename){

reordered.header <- TN.subset.header[order(TN.subset.header$Status),]
reordered.data <- TN.point.std[,order(TN.subset.header$Status)]
reordered.data <- reordered.data[rownames(reordered.data) %in% Discoveries$Discoveries,]
colnames(reordered.data) <- paste(reordered.header$Patient)
colcolor <- ifelse(reordered.header$Status == "Normal", "BLUE", "RED")

pdf(file = filename, width = 13, height = 10)
heatmap.2(reordered.data, Colv = "NULL", trace = "none", colsep = c(1:length(reordered.header$Status)-1)[reordered.header$Status == "Tumor"][1], col = blue2yellow(50), ColSideColors = colcolor, cexRow = 1, na.color = "darkgrey", symkey = TRUE)
legend("topright", c("Normal", "Tumor"), text.col = c("BLUE", "RED"))
dev.off()
	}


barplot.maker <- function(Discoveries, TN.subset.header, meta.data, TN.point.std, filename1, filename2){
#make a barplot of each discovery where the T-N values are averaged across tumors and values are ordered by difference.  Also output a .csv file of the TNdifferences.

pdf(file = filename1)

fold.change <- as.data.frame(matrix(NA, nrow = length(Discoveries[,1])+2, ncol = length(unique(TN.subset.header$Patient))*2))
rownames(fold.change) <- c("Tumor/Normal", "Sample", as.character(Discoveries[,1]))
fold.change[1,] <- rep(c("Normal", "Tumor"), each = length(unique(TN.subset.header$Patient)))
fold.change[2,] <- rep(unique(TN.subset.header$Patient), times = 2)



for (element in c(1:length(Discoveries[,1]))){

met.matches <- rownames(meta.data) %in% Discoveries[element,1]
data.subset <- TN.point.std[met.matches,]
patients <- unique(TN.subset.header$Patient)

cancer.abund <- rep(NA, times = length(patients))
normal.abund <- rep(NA, times = length(patients))

for (i in 1:length(patients)){
	cancer.abund[i] <- mean(data.subset[(TN.subset.header$Patient %in% patients[i]) & TN.subset.header$Status == "Tumor"])
	normal.abund[i] <- mean(data.subset[(TN.subset.header$Patient %in% patients[i]) & TN.subset.header$Status == "Normal"])
	}
	barplot(sort(cancer.abund - normal.abund), main = Discoveries[element,1], col = "RED")
	fold.change[element+2,] <- c(normal.abund, cancer.abund)

}

dev.off()

write.table(fold.change, file = filename2, col.names = FALSE, row.names = TRUE, sep = ",")
}




discovery.diff.plot <- function(meta.data, Discoveries, tumor.mat, normal.mat, file.name){

#make a difference plot for each discovery - looking at all peaks for a given metabolite that was significant in at least instrument/peak

pdf(file = file.name)

unique.disc <- unique(meta.data$Metabolite[rownames(meta.data) %in% Discoveries$Discoveries])

for (el in unique.disc){

met.matches <- meta.data$Metabolite %in% el
n.match <- sum(met.matches) 
met.color <- rainbow(n.match, start = 0.6, 0.95)

sig.matches <- rownames(tumor.mat)[met.matches] %in% Discoveries$Discoveries

q.val.match <- rep(NA, times = n.match)
for(i in 1:n.match){
	if(rownames(tumor.mat)[met.matches][i] %in% Discoveries$Discoveries){
		q.val.match[i] <- Discoveries$q_values[Discoveries$Discoveries %in% rownames(tumor.mat)[met.matches][i]]}}


cancer.match = matrix(tumor.mat[met.matches,], nrow = n.match)
normal.match = matrix(normal.mat[met.matches,], nrow = n.match)

test.stat.dif = apply(matrix(cancer.match - normal.match, nrow = n.match), 1, mean, na.rm = TRUE)
test.stat.se = apply(matrix(cancer.match - normal.match, nrow = n.match), 1, sd, na.rm = TRUE)/sqrt(apply(matrix((!is.na(cancer.match) & !is.na(normal.match)), nrow = n.match), 1, sum, na.rm = TRUE))


data.range <- range(c(cancer.match, normal.match), na.rm = TRUE)

npatients <- length(unique(tumor.samples$Patient))
patient_num <- rep(NA, times = length(tumor.samples[,1]))
for(pat in 1:length(patient_num)){
	patient_num[pat] <- c(1:npatients)[unique(tumor.samples$Patient) %in% tumor.samples$Patient[pat]]
	}

for(comp in 1:n.match){

non.na.subset <- !(is.na(cancer.match[comp,]) & is.na(normal.match[comp,]))

xpos <- patient_num - 0.45 + 0.9*(comp/(n.match + 1))

if(comp == 1){
plot(cancer.match[comp,non.na.subset] ~ xpos[non.na.subset], col = "RED", ylim = data.range, xlim = c(0, npatients + 2), main = paste("Metabolite abundances for ", el, sep = ""), xlab = "Patient", ylab = "Standardized log normalized ion-count")
points(normal.match[comp,non.na.subset] ~ xpos[non.na.subset], col = "BLUE")

}else{

points(cancer.match[comp,non.na.subset] ~ xpos[non.na.subset], col = "RED", ylim = data.range, xlim = c(0, npatients + 2), main = paste("Metabolite abundances for ", el, sep = ""))
points(normal.match[comp,non.na.subset] ~ xpos[non.na.subset], col = "BLUE")

}
	
for (i in 1:npatients){
	cancer.val <- cancer.match[comp,][patient_num == i][!is.na(cancer.match[comp,][patient_num == i])]
	normal.val <- normal.match[comp,][patient_num == i][!is.na(normal.match[comp,][patient_num == i])]

	

	if((length(cancer.val) >= 1) & length(normal.val) >= 1){
	if(sum(patient_num[non.na.subset] == i) > 1){
		nc.combos <- expand.grid(c(1:length(cancer.val)), c(1:length(normal.val)))
	
	diff.fact <- cancer.val[nc.combos[,1]] - normal.val[nc.combos[,2]]
	
	for(seg in 1:length(nc.combos[,1])){ 
	segments(xpos[patient_num == i][1], cancer.val[nc.combos[seg,1]], xpos[patient_num == i][1], normal.val[nc.combos[seg,2]], col = met.color[comp])
	}}
	segments(xpos[patient_num == i][1], cancer.val, xpos[patient_num == i][1], normal.val, col = met.color[comp])
	}}
	
	segments(npatients + 1 - 0.45 + 0.9*(comp/(n.match + 1)), 0, npatients + 1 - 0.45 + 0.9*(comp/(n.match + 1)), test.stat.dif[comp], col = met.color[comp])
	segments(rep(npatients + 1 - 0.45 + 0.9*(comp/(n.match + 1)), times = 2), test.stat.dif[comp], c(rep(npatients + 1 - 0.45 + 0.9*(comp/(n.match + 1)), times = 2) + c(-0.2,0.2)), test.stat.dif[comp] + ifelse(test.stat.dif[comp] < 0, 1, -1)*0.2, col = met.color[comp])
	segments(npatients + 1 - 0.45 + 0.9*(comp/(n.match + 1)) - 0.2, ifelse(test.stat.dif[comp] < 0, -1, 1)*test.stat.se[comp], npatients + 1 - 0.45 + 0.9*(comp/(n.match + 1)) + 0.2, ifelse(test.stat.dif[comp] < 0, -1, 1)*test.stat.se[comp])
		
	}
	
	legend("topright", legend = ifelse(sig.matches, paste(meta.data$Instrument[met.matches], ", q-value: ", q.val.match, sep = ""), meta.data$Instrument[met.matches]), text.col = c(met.color), fill = c(sig.matches), cex = 1)
	
	legend("topleft", legend = c("tumor", "normal"), text.col = c("RED", "BLUE"))
	

}

dev.off()
}
	
	
TN.boxplot <- function(save.sig.diff, subset.data, file.name, usepval = FALSE, pvals = NULL, legend_spec = NA, subset.index = c(1:length(names(save.sig.diff))), suppress.whiskers = FALSE){

#make a boxplot showing the fold difference between tumor and normal samples for each feature provided as elements of the list generated by boxplot.prep.  No outliers are allowed for the boxplot and the axis is converted from log differences to the postive reals by exponentiation.



subset.diffs <- save.sig.diff[subset.index]

logbox <- boxplot(subset.diffs, range = 0, plot = FALSE)

logbox$stats <- 2^logbox$stats
logbox$conf <- 2^logbox$conf

if(suppress.whiskers == TRUE){
	logbox$stats[1,] <- logbox$stats[2,]
	logbox$stats[5,] <- logbox$stats[4,]}

pdf(file = file.name, width = 10.6, height = 8)

bxp(logbox, log = "y", axes = FALSE, frame.plot = TRUE, boxfill = subset.data$color[subset.index], ylab = "Tumor:Normal", varwidth = TRUE)
abline(h = 1, lwd = 3, col = "gray48")
axis(1, at=seq(1,length(subset.diffs),by=1), labels= subset.data$palatable_name[subset.index], las = 2)
axis(2, at=c(1/8, 1/4, 1/2, 1, 2, 4, 8), labels= c(1/8, 1/4, 1/2, 1, 2, 4, 8), las = 2)
if(usepval == TRUE){
	if(!is.na(legend_spec)){
	legend("topleft", legend_spec$text , fill = legend_spec$color)}
	text(c(1:length(subset.index)), 0.11 + 0.125*odd(c(1:length(subset.index)))*0.12, labels = pvals[subset.index], cex = 12/length(subset.data[,1])*0.7 + 0.25)
	}


dev.off()
	
}



boxplot.prep <- function(Discoveries, tumor.mat, normal.mat, file.name){

#take a set of discoveries which you want to view as a boxplot, measure the T - N difference for each non NA feature and put this vector into a list which is then saved to the root directory.

save.sig.diff <- list()
save.sig.diff.q <- NULL
save.sig.diff.name <- NULL

for (i in c(1:length(Discoveries[,1]))){
cancer.match <- tumor.mat[rownames(tumor.mat) %in% Discoveries$Discoveries[i],]
normal.match <- normal.mat[rownames(tumor.mat) %in% Discoveries$Discoveries[i],]
cancer.diff <- cancer.match - normal.match

save.sig.diff$"x" <- cancer.diff
names(save.sig.diff)[names(save.sig.diff) == "x"] <- as.character(Discoveries$Discoveries[i])

save.sig.diff.q <- c(save.sig.diff.q, Discoveries$q_values[i])
save.sig.diff.name <- c(save.sig.diff.name, meta.data$Metabolite[rownames(meta.data) %in% Discoveries$Discoveries[i]])

}

save(save.sig.diff, save.sig.diff.q, save.sig.diff.name, file = file.name)

}



dilution_linearity <- function(TNmeta, header.data){

#analyze the linearity of ion counts versus sample dilutions using regression and anova

linear.comp <- list(nested.comp = rep(NA, times = length(TNmeta[,1])), slope.comp = rep(NA, times = length(TNmeta[,1])), slope = rep(NA, times = length(TNmeta[,1])))
anova.ind.tumor <- data.frame(patSS = rep(NA, times = length(TNmeta[,1])), patF = rep(NA, times = length(TNmeta[,1])), disSS = rep(NA, times = length(TNmeta[,1])), disF = rep(NA, times = length(TNmeta[,1])), patMS = rep(NA, times = length(TNmeta[,1])), disMS = rep(NA, times = length(TNmeta[,1])))

for(i in 1:length(TNmeta[,1])){
non.missing.val <- !(is.na(TNmeta[i,]) | TNmeta[i,] == 0 | (header.data$"exo-vivo Timecourse" %in% "30"))

if(!(sum(!non.missing.val)/length(non.missing.val) > 0.6)){

red.head <- header.data[non.missing.val,]
red.head$Patient <- as.factor(red.head$Patient)
red.head$"Sample Number" <- as.factor(red.head$"Sample Number") 
red.head$Dilution <- as.factor(red.head$Dilution)
red.head$Status <- as.factor(red.head$Status)

numeric.dil <- red.head$"Dilution"; numeric.dil <- as.character(numeric.dil)
numeric.dil[numeric.dil == "Low"] <- 1
numeric.dil[numeric.dil == "2x_Low"] <- 0
numeric.dil[numeric.dil == "High"] <- 2
numeric.dil  <- as.numeric(numeric.dil)

dil.fact.model <- lm(log(TNmeta[i, non.missing.val], base = 2) ~ red.head$Patient + red.head$"Patient"*red.head$"Status" + 
red.head$"Patient" + red.head$"Dilution")

dil.numeric.model <- lm(log(TNmeta[i, non.missing.val], base = 2) ~ red.head$"Patient" + red.head$"Patient"*red.head$"Status" + 
red.head$"Patient" + numeric.dil)

anova.ind.tumor$patSS[i] <- anova(dil.fact.model)$"Sum Sq"[rownames(anova(dil.fact.model)) == "red.head$Patient"]
anova.ind.tumor$patF[i] <- anova(dil.fact.model)$"F value"[rownames(anova(dil.fact.model)) == "red.head$Patient"]
anova.ind.tumor$disSS[i] <- anova(dil.fact.model)$"Sum Sq"[rownames(anova(dil.fact.model)) == "red.head$Status"]
anova.ind.tumor$disF[i] <- anova(dil.fact.model)$"F value"[rownames(anova(dil.fact.model)) == "red.head$Status"]
anova.ind.tumor$patMS[i] <- anova(dil.fact.model)$"Mean Sq"[rownames(anova(dil.fact.model)) == "red.head$Patient"]
anova.ind.tumor$disMS[i] <- anova(dil.fact.model)$"Mean Sq"[rownames(anova(dil.fact.model)) == "red.head$Status"]

linear.comp$nested.comp[i] <- anova(dil.fact.model, dil.numeric.model, test = "Chisq")[2,5]

linear.comp$slope.comp[i] <- (1 - pt(abs((summary(dil.numeric.model)$coef[rownames(summary(dil.numeric.model)$coef) == "numeric.dil",][1] - 1)/summary(dil.numeric.model)$coef[rownames(summary(dil.numeric.model)$coef) == "numeric.dil",][2]), dil.numeric.model$df))*2

linear.comp$slope[i] <- summary(dil.numeric.model)$coef[rownames(summary(dil.numeric.model)$coef) == "numeric.dil",][1]
	}}

return_list = list("linear.comp" = linear.comp, "anova.ind.tumor" = anova.ind.tumor)
return_list

	}



dyn_low_high <- function(return_list){

### determine whether the low or high concentration samples are more reliable for a feature
#determine whether the low or high concentration sample should be used based on the regression from dilution_linearity
	
	featureLH <- rep(NA, times = length(TNmeta[,1]))

for(i in 1:length(TNmeta[,1])){
	if(is.na(return_list[["linear.comp"]]$slope.comp[i]) & is.na(return_list[["linear.comp"]]$nested.comp[i])){
		featureLH[i] <- "H"
		}else{
	if(is.na(return_list[["linear.comp"]]$slope.comp[i])){
		if(return_list[["linear.comp"]]$nested.comp[i] > 0.05){
			featureLH[i] <- "H"
			}else{
				featureLH[i] <- "L"
				}
		}else{
			if(return_list[["linear.comp"]]$slope.comp[i] > 0.05){featureLH[i] <- "H"}else{
			featureLH[i] <- "L"
			}
		}
	}
	}
	
	featureLH
	}



KS_test_TN_diff <- function(normal.mat, tumor.mat){

#use kolmogrov-smirnoff tests to investigate the normality of the residuals and T-N differences: an assumption used by a parameteric t-test's.

KS_test <- rep(NA, times = length(normal.mat[,1]))
KS_test2 <- rep(NA, times = length(normal.mat[,1]))

for (i in 1:length(TN.subset.point[,1])){
  	n.val <- normal.mat[i,!is.na(normal.mat[i,]) & !is.na(tumor.mat[i,])]
 	t.val <- tumor.mat[i,!is.na(normal.mat[i,]) & !is.na(tumor.mat[i,])]
	KS_test[i] <- ks.test(n.val - t.val, pt, df = length(n.val)-1)$p.value
 	KS_test2[i] <- ks.test(c(n.val - mean(n.val - t.val), t.val), pt, df = length(n.val)-1)$p.value
	#KS_test[i] <- ks.test(n.val - t.val, pnorm)$p.value
 	#KS_test2[i] <- ks.test(c(n.val - mean(n.val - t.val), t.val), pnorm)$p.value

 } 
	hist(KS_test, main = "KS-test of T-N normality", xlab = "pvalue", ylab = "frequency", breaks = 20)
	hist(KS_test2, main = "KS-test of normality of residuals", xlab = "pvalue", ylab = "frequency", breaks = 20)
	ks.test(KS_test2, punif)
	}

	
					