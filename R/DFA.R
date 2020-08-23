
DFA <- function(data, groups, variables, normtests='yes', priorprob='SIZES', predictive='yes') {

cat('\n\n\n\nLinear Discriminant Function Analysis:\n')


donnes <- as.data.frame(data[,c(groups,variables)])

if (anyNA(donnes) == TRUE) {
	donnes <- na.omit(donnes)
	cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
}


grpnames <- as.vector(as.matrix(donnes[,1])) # group names, in the same order as in the data matrix
grpnames <- unique(grpnames)
grpnums  <- seq(1:length(grpnames))
ngroups  <- length(grpnames)
nDVs <- ncol(donnes) - 1
N <- nrow(donnes)

Ncases <- nrow(donnes)


if (is.factor(donnes[,1]) == F)  donnes[,1] <- factor(donnes[,1], ordered = FALSE, labels=grpnames)

donnes <- as.data.frame(donnes)

grpFreqs <- as.matrix(table(donnes[,1]))


# descriptive statistics & tests of univarite & multivariate normality -- from the MVN package
if (is.null(normtests) | normtests == 'YES' | normtests == 'yes') {

	# whole sample / all groups combined
	cat('\n\nWhole-sample statistics\n\n')
	umvn(donnes[,2:ncol(donnes)])

	# separate stats for each group
	cat('\n\nGroup statistics\n\n')
	for (lupeg in 1:ngroups) { 
		dum <- subset(donnes, donnes[,1] == grpnames[lupeg])
		cat('\n\nGroup', paste(grpnames[lupeg]),'\n')
		umvn(dum[,2:ncol(dum)]) 
	}
}
cat('\n\n')


# run the tests for Homogeneity of Variances & Covariances function, print the results, & get the sscp matrices
sscps <- homovarcovar(donnes)

# LDA from MASS package


#  the lda function from MASS produces different raw, lda coefficients when different priors are used
#  but SPSS produces the same coefficients regardless of the priors that are used
#  to produce the SPSS results, use priors based on the group sizes, as in prior=(grpFreqs/sum(grpFreqs))

#  from the Details for the lda function in MASS:
#  Specifying the prior will affect the classification unless over-ridden in predict.lda. Unlike in 
#  most statistical packages, it will also affect the rotation of the linear discriminants within their space,
#  as a weighted between-groups covariance matrix is used. Thus the first few linear discriminants emphasize 
#  the differences between groups with the weights given by the prior, which may differ from their prevalence in the dataset. 


# SPSS options for priors are "All groups equal" or "Compute from group sizes"
if (is.null(priorprob)) priorprob = 'SIZES'
if (priorprob == 'EQUAL') priors = matrix((1/ngroups), 1, ngroups) 
if (priorprob == 'SIZES') priors = grpFreqs/sum(grpFreqs) 

ldaoutput <- MASS::lda(x = as.matrix(donnes[,2:ncol(donnes)]), grouping=donnes[,1], prior=priors)

lda.values <- stats::predict(ldaoutput, donnes[,2:ncol(donnes)]) # obtain scores on the DFs



# eigenvalues, canonical correlations, & one-way anovas on the DFs
dfc <- data.frame(donnes[,1], lda.values$x)
evals <- matrix(-9999, ncol(ldaoutput$scaling), 4)
evals[,1] <- 1:ncol(ldaoutput$scaling)
anovaDFoutput <- matrix(-9999, ncol(ldaoutput$scaling), 5)
ttestDFoutput <- lapply(1:ncol(ldaoutput$scaling), function(x) matrix(-9999, nrow=choose(ngroups,2), ncol=13))
names(ttestDFoutput) = c(paste("Discriminant Function ", 1:ncol(ldaoutput$scaling), sep=""))
for (luper in 1:nrow(evals)) {
	dd <- data.frame(dfc[,1], dfc[luper+1])
	colnames(dd) <- c('grp','dv')
	fit <- stats::lm(dv ~ as.factor(grp), data = dd)
	betwss <- stats::anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- stats::anova(fit)["Residuals", "Sum Sq"]
 	anovaDFoutput[luper,1] <- as.numeric(summary(fit)["r.squared"])
	anovaDFoutput[luper,2] <- stats::anova(fit)["as.factor(grp)","F value"]
	anovaDFoutput[luper,3] <- stats::anova(fit)["as.factor(grp)","Df"]
	anovaDFoutput[luper,4] <- fit$df.residual
	anovaDFoutput[luper,5] <- stats::anova(fit)["as.factor(grp)","Pr(>F)"]

	ttestDFoutput[[luper]] <- ttestboc(dd, varest=FALSE)			

	evals[luper,2] <- betwss / withss # eigenvalue
	evals[luper,4] <- sqrt(betwss / (betwss+withss)) # canonical correlation
}
sv <- ldaoutput$svd;  svproprotions <- sv^2/sum(sv^2) # % of variance
evals[,3] <- svproprotions
cat('\n\n\nEigenvalues & canonical correlations:\n\n')
dimnames(evals) <- list(rep("", dim(evals)[1]))
colnames(evals) <- c('Function','  eigenvalue','     proportion of variance','     canonical correlation')
print(round(evals,3))



cat('\n\n\nMultivariate peel-down significance tests:\n\n')  # using p.asym from the CCP package

# July, 2019: for CCA, the CCP::p.asym function uses the #s of vars in set1 and set2
# as the values for p & q; but to get correct results for DFA using this function,
# p = the # of DVs and q = the # of DFs -- at least that what seems to yield correct values

rho <- evals[,4]
NCVs <- nrow(evals)

mvFmat <- Wilks(rho=rho, Ncases=Ncases, p = NCVs, q = nDVs)
colnames(mvFmat) <- c('            Wilks Lambda', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')

mvFmat <- Pillai(rho=rho, Ncases=Ncases, p = NCVs, q = nDVs)
colnames(mvFmat) <- c('  Pillai-Bartlett Trace', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')

mvFmat <- Hotelling(rho=rho, Ncases=Ncases, p = NCVs, q = nDVs)
colnames(mvFmat) <- c('   Hotelling-Lawley Trace', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')

mvFmat <- RoyRoot(rho=rho, Ncases=Ncases, p = NCVs, q = nDVs)
colnames(mvFmat) <- c('      Roy\'s Largest Root', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')



cat('\n\n\n\nCanonical Discriminant Function (raw) Coefficients:\n')
colnames(ldaoutput$scaling) <-  c(paste("Function ", 1:ncol(ldaoutput$scaling), sep=""))
print(round(ldaoutput$scaling,3))

# centering each variable within groups
group.center <- function(var,grp) { return(var-tapply(var,grp,mean,na.rm=T)[grp]) }
cdonnes <- matrix(-9999,nrow(donnes),(ncol(donnes)-1))
dfc <- cbind(donnes[,1], lda.values$x)
cdfc <- matrix(-9999,nrow(dfc),(ncol(dfc)-1))
for (lupec in 1:(ncol(donnes)-1)) { cdonnes[,lupec] <- group.center(donnes[,(lupec+1)], donnes[,1]) }
for (lupec in 1:(ncol(dfc)-1)) { cdfc[,lupec] <- group.center(dfc[,(lupec+1)], dfc[,1]) }
cdonnes <- cbind(donnes[,1], cdonnes) # placing the grouping variable back in the centered matrix
cdonnesdf <- data.frame(cdonnes)
structCoef <- stats::cor(x = cdonnes[,2:ncol(cdonnes)], y = cdfc) ; # round(structCoef,2)
rownames(structCoef) <- rownames(ldaoutput$scaling)
colnames(structCoef) <- colnames(ldaoutput$scaling)
cat('\n\nStructure Coefficients:\n')
colnames(structCoef) <-  c(paste("Function ", 1:ncol(structCoef), sep=""))
print(round(structCoef,3))

# standardized coefficients
pooledSDs <- as.matrix(apply(cdonnes, 2, FUN = sd)) # cdonnes contains the mean-centered data
standCoef <- (pooledSDs[2:nrow(pooledSDs),]) * ldaoutput$scaling
cat('\n\nStandardized Coefficients:\n')
colnames(standCoef) <-  c(paste("Function ", 1:ncol(standCoef), sep=""))
print(round(standCoef,3))



sscpwith <- sscps$sscpwith
sscpbetw <- sscps$sscpbetw

# provides the standardized coefficients from SPSS:
poolwith <- sscpwith * (1/(nrow(donnes)-ngroups))
pooledSDs <- sqrt(diag(poolwith)) # pooled SDs for SPSS results
standCoefSPSS <- pooledSDs * ldaoutput$scaling
cat('\n\n\nStandardized Coefficients from SPSS:\n')
colnames(standCoefSPSS) <- c(paste("Function ", 1:ncol(standCoefSPSS), sep=""))
print(round(standCoefSPSS,3))

	
# group means & SDs on the raw ldfs
ldfscores  <- data.frame(donnes[,1],lda.values$x)
ldfscoresZ <- data.frame(donnes[,1],scale(lda.values$x))
centroids <- centroidSDs <- centroidsZ <- centroidSDsZ <- matrix(-9999,ngroups,(ncol(ldfscores)-1))
for (lupec in 2:ncol(ldfscores)) {
	aggM  <- stats::aggregate(x = ldfscores[,lupec], by= list(ldfscores[,1]), FUN = mean)
	aggSD <- stats::aggregate(x = ldfscores[,lupec], by= list(ldfscores[,1]), FUN = sd)
	centroids[,(lupec-1)] <- aggM[,2]
	centroidSDs[,(lupec-1)] <- aggSD[,2]

	aggMZ  <- stats::aggregate(x = ldfscoresZ[,lupec], by= list(ldfscoresZ[,1]), FUN = mean)
	aggSDZ <- stats::aggregate(x = ldfscoresZ[,lupec], by= list(ldfscoresZ[,1]), FUN = sd)
	centroidsZ[,(lupec-1)] <- aggMZ[,2]
	centroidSDsZ[,(lupec-1)] <- aggSDZ[,2]
}
cat('\n\n\nFunctions at Group Centroids\n')
cat('\nUnstandardized canonical discriminant functions evaluated at group means:\n\n')
rownames(centroids) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroids) <- c(paste("Function ", 1:(ncol(centroids)), sep=""))
print(round(centroids,3))

cat('\n\nGroup Standard Deviations on the unstandardized functions:\n\n')
rownames(centroidSDs) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroidSDs) <- c(paste("Function ", 1:(ncol(centroidSDs)), sep=""))
print(round(centroidSDs,3))

cat('\n\nStandardized canonical discriminant functions evaluated at group means:\n\n')
rownames(centroidsZ) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroidsZ) <- c(paste("Function ", 1:(ncol(centroidsZ)), sep=""))
print(round(centroidsZ,3))

cat('\n\nGroup Standard Deviations on the standardized functions:\n\n')
rownames(centroidSDsZ) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroidSDsZ) <- c(paste("Function ", 1:(ncol(centroidSDsZ)), sep=""))
print(round(centroidSDsZ,3))


cat('\n\n\nOne-way ANOVAs using the scores on a discriminant function as the DV:\n\n')
dimnames(anovaDFoutput) <-list(rep("", dim(anovaDFoutput)[1]))
colnames(anovaDFoutput) <- c('Eta-squared','          F','    df','    df','        p')
rownames(anovaDFoutput) <- colnames(ldaoutput$scaling)
anovaDFoutput[,1:4] <- round(anovaDFoutput[,1:4],2)
anovaDFoutput[,5] <- round(anovaDFoutput[,5],4)
print(anovaDFoutput)


# one-way anovas & t-tests on the DVs
anovaDVoutput <- matrix(-9999, length(variables), 5)
ttestDVoutput <- lapply(1:length(variables), function(x) matrix(-9999, nrow=choose(ngroups,2), ncol=13))
names(ttestDVoutput)=variables 
for (lupec in 1:length(variables)) {
	ddd <- data.frame(donnes[,1], donnes[,(lupec+1)])
	colnames(ddd) <- c('grp','dv')
	fit <- stats::lm(dv ~ as.factor(grp), data = ddd)
	betwss <- stats::anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- stats::anova(fit)["Residuals", "Sum Sq"]
 	anovaDVoutput[lupec,1] <- as.numeric(summary(fit)["r.squared"])
	anovaDVoutput[lupec,2] <- stats::anova(fit)["as.factor(grp)","F value"]
	anovaDVoutput[lupec,3] <- stats::anova(fit)["as.factor(grp)","Df"]
	anovaDVoutput[lupec,4] <- fit$df.residual
	anovaDVoutput[lupec,5] <- stats::anova(fit)["as.factor(grp)","Pr(>F)"]
	
	ttestDVoutput[[lupec]] <- ttestboc(ddd, varest=FALSE)			
}
cat('\n\n\n\nOne-way ANOVAs on the original DVs\n') 
cat('\n(provided for comparisons with the ANOVAs on the discriminant functions):\n\n')
dimnames(anovaDVoutput) <-list(rep("", dim(anovaDVoutput)[1]))
colnames(anovaDVoutput) <- c('Eta-squared','          F','    df','    df','        p')
rownames(anovaDVoutput) <- variables
anovaDVoutput[,1:4] <- round(anovaDVoutput[,1:4],2)
anovaDVoutput[,5] <- round(anovaDVoutput[,5],4)
print(anovaDVoutput)


cat('\n\n\n\nt-tests and effect sizes for group differences on the discriminant functions:\n\n\n')
#print(ttestDFoutput)
lsnoms <- names(ttestDFoutput)
for (lupe in 1:length(ttestDFoutput)) {
	cat(lsnoms[lupe],'\n\n')
	print(ttestDFoutput[[lupe]], row.names = F); cat('\n\n')
}


cat('\n\n\n\nt-tests and effect sizes for group differences on the original DVs\n')
cat('\n(provided for comparisons with the t-tests on the discriminant functions):\n\n\n')
#print(ttestDVoutput)
lsnoms <- names(ttestDVoutput)
for (lupe in 1:length(ttestDVoutput)) {
	cat(lsnoms[lupe],'\n\n')
	print(ttestDVoutput[[lupe]], row.names = F); cat('\n\n')
}





# 2014 Bird - Controlling the Maximum Familywise Type I Error Rate in Analyses of Multivariate Experiments

cat('\n\nAppropriate univariate significance tests (ignoring the discriminant functions)')
cat('\nfor controlling the maximum familywise Type I error rate (MFWER), based on:')
cat('\n\nBird, K. D., & Hadzi-Pavlovic, D. (2013). Controlling the maximum familywise Type I ')
cat('\nerror rate in analyses of multivariate experiments. Psychological Methods, 19(2), 265-280.\n')

lsnoms <- names(ttestDVoutput)


if (ngroups == 2) {

	# Two Groups 
	# *** Critical values for analyses controlling MFWER when k = 2 ***. 
	N <- 40
	P <- 4
	ALPHA <- 0.05 
	# *** Critical t and p values and confidence level for t* analysis ***. 
	TCRIT <- qt(1-(1-(1-ALPHA)^(1/P))/2,N-2)
	PCRIT <- 1-(1-ALPHA)^(1/P)
	CLEVEL <- 1-PCRIT 
	# *** Critical t and p values for Mt** analysis ***. 
	TCRIT2 <- qt(1-(1-(1-ALPHA)**(1/(P-1)))/2,N-2) 
	PCRIT2 <- 1-(1-ALPHA)^(1/(P-1))
	# # output 
	# RESULT <- c("N", "P", "ALPHA", "TCRIT", "PCRIT", "CLEVEL", "TCRIT2", "PCRIT2")
	# RESULT <- sapply(RESULT,get)
	# print(RESULT,digits=5)

	cat('\nWhen there are only two groups and the researcher decides to control the MFWER') 
	cat('\nat .05 by (only) carrying out multiple t tests, then the MANOVA test is')
	cat('\nirrelevant and the critical t value is the t1-critical value below.')
	cat('\nThe t value for the data must be greater than t1-critical value for a pairwise')
	cat('\ncomparison to be significant.\n\n\n')
		
	for (lupe in 1:length(ttestDVoutput)) {
		dddd <- ttestDVoutput[[lupe]]
		BirdOut <- dddd[,c(1,5,9)]
		BirdOut[,4] <- round(TCRIT,4)
		BirdOut[,5] <- ifelse( abs(BirdOut[,3]) > TCRIT, 'significant', 'not significant')
		colnames(BirdOut)[4:5] <- c('    t1-critical','         t1 decision')
#		rownames(BirdOut) <- ""
		cat(lsnoms[lupe],'\n\n')
		print(BirdOut, row.names = F); cat('\n\n')
	}

	cat('\n\nWhen there are only two groups and the researcher prefers the two-stage')
	cat('\napproach to controling the MFWER at .05, then the MANOVA F test must be')
	cat('\nsignificant and the t value for each comparison must be greater than') 
	cat('\nt2-critical value below for a pairwise comparison to be significant.\n\n\n')

	for (lupe in 1:length(ttestDVoutput)) {
		dddd <- anovaDVoutput[[lupe]]
		BirdOut <- BirdOut[,1:3]
		BirdOut[,4] <- round(TCRIT2,4)
		BirdOut[,5] <- ifelse( abs(BirdOut[,3]) > TCRIT2, 'significant', 'not significant')
		colnames(BirdOut)[4:5] <- c('    t2-critical','         t2 decision')
#		rownames(BirdOut) <- ""
		cat(lsnoms[lupe],'\n\n')
		print(BirdOut, row.names = F); cat('\n\n')
	}
}


if (ngroups > 2) {

	# More Than Two Groups  
	# *** Critical values for analyses controlling MFWER when k >= 3 ***. 
	# N <- N
	K <- ngroups
	P <- nDVs 
	ALPHA <- .05
	# *** Critical t and p values and confidence level for t' analysis ***.
	TCRIT <- qtukey((1-ALPHA)**(1/P),K,N-K)/sqrt(2)
	PCRIT <- 2*(1-pt(TCRIT,N-K))
	CLEVEL <- 1-PCRIT
	# *** Critical F and critical t and p values for F't'' analysis ***.
	FCRIT <- qf((1-ALPHA)^(1/P),K-1,N-K)
	TCRIT2 <- qtukey((1-ALPHA)^(1/P),K-1,N-K)/sqrt(2)
	PCRIT2 <- 2*(1-pt(TCRIT2,N-K))
	# #output
	# RESULT <- c("N", "K", "P", "ALPHA", "TCRIT", "PCRIT", "CLEVEL", "FCRIT", "TCRIT2", "PCRIT2") 
	# RESULT <- sapply(RESULT,get)
	# print(RESULT,digits=5) 
	
	cat('\n\nFor designs with more than two groups, Bird et al. (2014) could find no')
	cat('\njustification for the inclusion of an initial MANOVA test in MCPs designed to') 
	cat('\ncontrol the MFWER with protected t tests.\n')
	
	cat('\nWhen a researcher decides to control the MFWER at .05 by soley carrying out') 
	cat('\nmultiple t tests, then the MANOVA and ANOVA tests are irrelevant and the')
	cat('\ncritical t value is the t1-critical value below. The t value for the data must')
	cat('\nbe greater than t1-critical value for a pairwise comparison to be significant.\n\n\n')
		
	for (lupe in 1:length(ttestDVoutput)) {
		dddd <- ttestDVoutput[[lupe]]
		BirdOut <- dddd[,c(1,5,9)]
		BirdOut[,4] <- round(TCRIT,4)
		BirdOut[,5] <- ifelse( abs(BirdOut[,3]) > TCRIT, 'significant', 'not significant')
		colnames(BirdOut)[4:5] <- c('    t1-critical','         t1 decision')
		cat(lsnoms[lupe],'\n\n')
		print(BirdOut, row.names = F); cat('\n\n')
	}

	cat('\n\nWhen a researcher prefers the two-stage approach to controling the MFWER at .05,')
	cat('\nthen the ANOVA F value for the data must be greater than the F-critical value below,')
	cat('\nand the t value for each comparison must be greater than t2-critical value') 
	cat('\nbelow for a pairwise comparison to be significant.\n\n\n')

	for (lupe in 1:length(ttestDVoutput)) {
		cat(lsnoms[lupe],'\n')

		decisionF <- ifelse( anovaDVoutput[lupe,2] > FCRIT, 'significant', 'not significant')
		cat('\nF =', anovaDVoutput[lupe,2],'   F-critical =', round(FCRIT,4),'    decision =', decisionF,'\n\n')

		dddd <- ttestDVoutput[[lupe]]
		BirdOut <- dddd[,c(1,5,9)]
		BirdOut[,4] <- round(TCRIT2,4)
		BirdOut[,5] <- ifelse( abs(BirdOut[,3]) > TCRIT2, 'significant', 'not significant')
		colnames(BirdOut)[4:5] <- c('    t2-critical','         t2 decision')
		print(BirdOut, row.names = F); cat('\n\n')
	}
}



# plot
colnames(centroidsZ) <-  c(paste("DF ", 1:ncol(centroidSDsZ), sep=""))

graphics::matplot(1:length(grpnames), centroidsZ, type = "l", lty=1, lwd=3, 
        xaxt='n', xlab=groups, cex.axis=1.2, cex.lab = 1.3,
        ylab='Discriminant Function z Scores', ylim = c(-3, 3), cex.axis=1.2        )
graphics::axis(side=1, at=grpnums, labels=rownames(centroidsZ), xlab="groups")
graphics::title(main='Mean Standardized Discriminant Function Scores for the Groups')
graphics::legend("topright", legend = colnames(centroidsZ), bty="n", lwd=2, col=1:ncol(centroidsZ))



if (predictive == 'YES' | predictive == 'yes' | is.null(predictive)) {

cat('\n\n\n\nPREDICTIVE DISCRIMINANT ANALYSIS:\n')

cat('\n\nPrior Probabilities for Groups:\n')
print(round(ldaoutput$prior,3))

# Frequencies: Original vs Predicted

# freqs_op <- data.frame(cbind(donnes[,1], lda.values$class))
# colnames(freqs_op) <-  c("Original", "Predicted") 
#freqs <- table(freqs_op) # doesn't produce a square table if the are no values for a factor level

squareTable <- function(x,y) {
    Original  <- factor(x, levels = grpnames)
    Predicted <- factor(y, levels = grpnames)
    table(Original, Predicted)
}

freqs <- squareTable(donnes[,1], lda.values$class)

cat('\n\nCross-Tabulation of the Original and Predicted Group Memberships:\n\n')
print(freqs)	

cat('\n\nProportion of original grouped cases correctly classified:  ',
    round((sum(diag(freqs)) / sum(freqs)),3))

cat('\n\n\nChi-square test of independence:\n')
print(summary(freqs)) # chi-square test of indepedence

cat('\n\nRow Frequencies:\n\n')
print(margin.table(freqs, 1)) # A frequencies (summed over B) 
cat('\n\nColumn Frequencies:\n\n')
print(margin.table(freqs, 2)) # B frequencies (summed over A)

cat('\n\nCell Proportions:\n\n')
print(round(prop.table(freqs),2)) 
cat('\n\nRow-Based Proportions:\n\n')
print(round(prop.table(freqs, 1),2)) 
cat('\n\nColumn-Based Proportions:\n\n')
print(round(prop.table(freqs, 2),2)) 

# 
cat('\n\nAgreement (kappas) between the Predicted and Original Group Memberships:\n\n')
kappas_cvo <- kappas(stats::na.omit(cbind(lda.values$class, donnes[,1])))
print(kappas_cvo,3)


# Frequencies: Original vs Cross-Validated (leave-one-out cross-validation)

# classifications from leave-one-out cross-validation
ldaoutputCV <- MASS::lda(x = as.matrix(donnes[,c(2:ncol(donnes))]),
   grouping=donnes[,1], prior=priors, CV = TRUE)

freqs_cvp <- data.frame(cbind(ldaoutputCV$class, lda.values$class))
colnames(freqs_cvp) <-  c("Cross-Validated", "Predicted") 
# freqsCVP <- table(freqs_cvp) # doesn't produce a square table if the are no values for a factor level
freqsCVP <-  squareTable(ldaoutputCV$class, lda.values$class)
# colnames(freqsCVP) <-  paste("Group ", 1:ngroups, sep="") 
# rownames(freqsCVP) <-  paste("Group ", 1:ngroups, sep="") 
colnames(freqsCVP) <- paste(grpnames)
rownames(freqsCVP) <- paste(grpnames)
cat('\n\n\n\nCross-Tabulation of the Cross-Validated and Predicted Group Memberships:\n\n')
print(freqsCVP)	

cat('\n\nProportion of cross-validated grouped cases correctly classified:  ',
    round((sum(diag(freqsCVP)) / sum(freqsCVP)),3))

cat('\n\nChi-square test of indepedence:\n')
print(summary(freqsCVP)) # chi-square test of indepedence

cat('\n\nRow Frequencies:\n\n')
print(margin.table(freqsCVP, 1)) # A frequencies (summed over B) 
cat('\n\nColumn Frequencies:\n\n')
print(margin.table(freqsCVP, 2)) # B frequencies (summed over A)

cat('\n\nCell Proportions:\n\n')
print(round(prop.table(freqsCVP),2)) 
cat('\n\nRow-Based Proportions:\n\n')
print(round(prop.table(freqsCVP, 1),2))  
cat('\n\nColumn-Based Proportions:\n\n')
print(round(prop.table(freqsCVP, 2),2)) 

# kappas
cat('\n\nAgreement (kappas) between the Cross-Validated and Original Group Memberships:\n\n')
kappas_cvo <- kappas(stats::na.omit(cbind(ldaoutputCV$class, donnes[,1])))
print(kappas_cvo)

# kappas
cat('\n\nAgreement (kappas) between the Cross-Validated and Predicted Group Memberships:\n\n')
kappas_cvp <- kappas(stats::na.omit(cbind(ldaoutputCV$class, lda.values$class)))
print(kappas_cvp,3)
cat('\n\n\n')

}

DFAoutput <- list(  
   rawCoef=ldaoutput$scaling,
   structCoef=structCoef,
   standCoef=standCoef,
   standCoefSPSS=standCoefSPSS,
   centroids=centroids,
   centroidSDs=centroidSDs,
   centroidsZ=centroidsZ,
   centroidSDsZ=centroidSDsZ,
   DFAscores=lda.values$x )

return(invisible(DFAoutput))

}



