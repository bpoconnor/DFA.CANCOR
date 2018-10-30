

DFA <- function( data, groups, variables, normtests='yes', priorprob='SIZES', predictive='yes') {

cat('\n\n\n\nLinear Discriminant Function Analysis:\n')

donnes <- cbind(data[,groups],data[,variables])

#donnes[,1] <- as.numeric((donnes[,1]))
grpnames <- unique(donnes[,1])
ngroups <- length(grpnames)

if (is.factor(donnes[,1]) == F)  donnes[,1] <- factor( donnes[,1], ordered = FALSE, labels=grpnames)

donnes <- as.data.frame(donnes)

grpFreqs <- as.matrix(table(donnes[,1]))


# descriptive statistics & tests of univarite & multivariate normality -- from the MVN package
if (is.null(normtests) | normtests == 'YES' | normtests == 'yes' ) {

	# whole sample / all groups combined
	cat('\n\nWhole-sample statistics\n\n')
	umvn(donnes[,2:ncol(donnes)])

	# separate stats for each group
	cat('\n\nGroup statistics\n\n')
	for (lupeg in 1:ngroups) { 
#		dum <- subset(donnes, donnes[,1] == lupeg)
		dum <- subset(donnes, donnes[,1] == grpnames[lupeg] )
		cat('\n\nGroup ', paste(grpnames[lupeg]),':\n')
		umvn(dum[,2:ncol(dum)]) 
	}
}
cat('\n\n')


# run the tests for Homogeneity of Variances & Covariances function, print the results, & get the sscp matrices
sscps <- homovarcovar( donnes )

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

ldaoutput <- lda(x = as.matrix(donnes[,2:ncol(donnes)]), grouping=donnes[,1], prior=priors)
lda.values <- predict(ldaoutput, donnes[,2:ncol(donnes)]) # obtain scores on the DFs


# eigenvalues, canonical correlations, & one-way anovas on the DFs
dfc <- data.frame(donnes[,1], lda.values$x)
evals <- matrix(-9999, ncol(ldaoutput$scaling), 4)
evals[,1] <- 1:ncol(ldaoutput$scaling)
anovaDFoutput <- matrix(-9999, ncol(ldaoutput$scaling), 5)
ttestDFoutput <- lapply(1:ncol(ldaoutput$scaling), function(x) matrix(-9999, nrow=choose(ngroups,2), ncol=13))
names(ttestDFoutput) = c(paste("Function ", 1:ncol(ldaoutput$scaling), sep="") )
for (luper in 1:nrow(evals)) {
	dd <- data.frame( dfc[,1], dfc[luper+1] )
	colnames(dd) <- c('grp','dv')
	fit <- lm(dv ~ as.factor(grp), data = dd)
	betwss <- anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- anova(fit)["Residuals", "Sum Sq"]
 	anovaDFoutput[luper,1] <- as.numeric(summary(fit)["r.squared"])
	anovaDFoutput[luper,2] <- anova(fit)["as.factor(grp)","F value"]
	anovaDFoutput[luper,3] <- anova(fit)["as.factor(grp)","Df"]
	anovaDFoutput[luper,4] <- fit$df.residual
	anovaDFoutput[luper,5] <- anova(fit)["as.factor(grp)","Pr(>F)"]

	ttestDFoutput[[luper]] <- ttestboc(dd, varest=FALSE)			

	evals[luper,2] <- betwss / withss # eigenvalue
	evals[luper,4] <- sqrt( betwss / (betwss+withss)) # canonical correlation
}
sv <- ldaoutput$svd;  svproprotions <- sv^2/sum(sv^2) # % of variance
evals[,3] <- svproprotions
cat('\n\n\nEigenvalues & canonical correlations:\n\n')
dimnames(evals) <- list(rep("", dim(evals)[1]))
colnames(evals) <- c('Function','  eigenvalue','     proportion of variance','     canonical correlation')
print(round(evals,3))


cat('\n\n\nMultivariate peel-down significance tests:\n\n')  # using p.asym from the CCP package

rho <- evals[,4]
N <- nrow(donnes)
p <- length(variables)
q <- length(evals[,4])


invisible(capture.output( pW <-  ((p.asym(rho, N, p, q, tstat = "Wilks"))) ))
dmat <- cbind(pW$stat, pW$approx, pW$df1, pW$df2, pW$p.value )
colnames(dmat) <- c('            Wilks Lambda', '     F-approx. ', '     df1', '          df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')

invisible(capture.output( pH <- p.asym(rho, N, p, q, tstat = "Hotelling") ))
dmat <- cbind(pH$stat, pH$approx, pH$df1, pH$df2, pH$p.value )
colnames(dmat) <- c('  Hotelling-Lawley Trace', '     F-approx. ', '     df1', '          df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')

invisible(capture.output( pP <- p.asym(rho, N, p, q, tstat = "Pillai") ))
dmat <- cbind(pP$stat, pP$approx, pP$df1, pP$df2, pP$p.value )
colnames(dmat) <- c('   Pillai-Bartlett Trace', '     F-approx. ', '     df1', '          df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')

invisible(capture.output( pR <- p.asym(rho, N, p, q, tstat = "Roy") ))
dmat <- cbind(pR$stat, pR$approx, pR$df1, pR$df2, pR$p.value )
colnames(dmat) <- c('      Roy\'s Largest Root', '     F-approx. ', '     df1', '          df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')

cat('\n\n\n\nCanonical Discriminant Function (raw) Coefficients:\n')colnames(ldaoutput$scaling) <-  c(paste("Function ", 1:ncol(ldaoutput$scaling), sep="") )
print(round(ldaoutput$scaling,3))# centering each variable within groupsgroup.center <- function(var,grp) { return(var-tapply(var,grp,mean,na.rm=T)[grp]) }cdonnes <- matrix(-9999,nrow(donnes),(ncol(donnes)-1))dfc <- cbind(  donnes[,1], lda.values$x   )cdfc <- matrix(-9999,nrow(dfc),(ncol(dfc)-1))for (lupec in 1:(ncol(donnes)-1)) { cdonnes[,lupec] <- group.center( donnes[,(lupec+1)], donnes[,1] ) }for (lupec in 1:(ncol(dfc)-1)) { cdfc[,lupec] <- group.center( dfc[,(lupec+1)], dfc[,1] ) }cdonnes <- cbind( donnes[,1], cdonnes) # placing the grouping variable back in the centered matrixcdonnesdf <- data.frame(cdonnes)structCoef <- cor( x = cdonnes[,2:ncol(cdonnes)], y = cdfc) ; # round(structCoef,2)rownames(structCoef) <- rownames(ldaoutput$scaling)colnames(structCoef) <- colnames(ldaoutput$scaling)cat('\n\nStructure Coefficients:\n')colnames(structCoef) <-  c(paste("Function ", 1:ncol(structCoef), sep="") )
print(round(structCoef,3))# standardized coefficientspooledSDs <- as.matrix(apply(cdonnes, 2, FUN = sd) ) # cdonnes contains the mean-centered datastandCoef <- (pooledSDs[2:nrow(pooledSDs),]) * ldaoutput$scalingcat('\n\nStandardized Coefficients:\n')colnames(standCoef) <-  c(paste("Function ", 1:ncol(standCoef), sep="") )
print(round(standCoef,3))


sscpwith <- sscps$sscpwith
sscpbetw <- sscps$sscpbetw

# provides the standardized coefficients from SPSS:
poolwith <- sscpwith * (1/(nrow(donnes)-ngroups))
pooledSDs <- sqrt(diag(poolwith)) # pooled SDs for SPSS results
standCoefSPSS <- pooledSDs * ldaoutput$scaling
cat('\n\n\nStandardized Coefficients from SPSS:\n')
colnames(standCoefSPSS) <-  c(paste("Function ", 1:ncol(standCoefSPSS), sep="") )
print(round(standCoefSPSS,3))


# group means & SDs on the raw ldfs
ldfscores <- as.matrix(cbind(donnes[,1],lda.values$x))
centroids <- centroidSDs <- matrix(-9999,ngroups,(ncol(ldfscores)-1))

ldfscoresZ <- as.matrix(cbind(donnes[,1],scale(lda.values$x)))
centroidsZ <- centroidSDsZ <- matrix(-9999,ngroups,(ncol(ldfscores)-1))

for (lupec in 2:ncol(ldfscores)) {
	aggM  <- aggregate( formula = ldfscores[,lupec] ~ ldfscores[,1], data = ldfscores, FUN = mean )
	aggSD <- aggregate( formula = ldfscores[,lupec] ~ ldfscores[,1], data = ldfscores, FUN = sd )
	centroids[,(lupec-1)] <- aggM[,2]
	centroidSDs[,(lupec-1)] <- aggSD[,2]

	aggMZ  <- aggregate( formula = ldfscoresZ[,lupec] ~ ldfscoresZ[,1], data = ldfscoresZ, FUN = mean )
	aggSDZ <- aggregate( formula = ldfscoresZ[,lupec] ~ ldfscoresZ[,1], data = ldfscoresZ, FUN = sd )
	centroidsZ[,(lupec-1)] <- aggMZ[,2]
	centroidSDsZ[,(lupec-1)] <- aggSDZ[,2]
}
cat('\n\n\nFunctions at Group Centroids:\n')
cat('\nUnstandardized canonical discriminant functions evaluated at group means:\n\n')
rownames(centroids) <-  c(paste("Group ", grpnames, sep="") ) 
colnames(centroids) <-  c(paste("Function ", 1:(ncol(centroids)), sep="") )
print(round(centroids,3))

cat('\n\nGroup Standard Deviations on the unstandardized functions:\n\n')
rownames(centroidSDs) <-  c(paste("Group ", grpnames, sep="") ) 
colnames(centroidSDs) <-  c(paste("Function ", 1:(ncol(centroidSDs)), sep="") )
print(round(centroidSDs,3))

cat('\n\nStandardized canonical discriminant functions evaluated at group means:\n\n')
rownames(centroidsZ) <-  c(paste("Group ", grpnames, sep="") ) 
colnames(centroidsZ) <-  c(paste("Function ", 1:(ncol(centroidsZ)), sep="") )
print(round(centroidsZ,3))

cat('\n\nGroup Standard Deviations on the standardized functions:\n\n')
rownames(centroidSDsZ) <-  c(paste("Group ", grpnames, sep="") ) 
colnames(centroidSDsZ) <-  c(paste("Function ", 1:(ncol(centroidSDsZ)), sep="") )
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
	dd <- data.frame(  donnes[,1], donnes[,(lupec+1)]   )
	colnames(dd) <- c('grp','dv')
	fit <- lm(dv ~ as.factor(grp), data = dd)
	betwss <- anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- anova(fit)["Residuals", "Sum Sq"]
 	anovaDVoutput[lupec,1] <- as.numeric(summary(fit)["r.squared"])
	anovaDVoutput[lupec,2] <- anova(fit)["as.factor(grp)","F value"]
	anovaDVoutput[lupec,3] <- anova(fit)["as.factor(grp)","Df"]
	anovaDVoutput[lupec,4] <- fit$df.residual
	anovaDVoutput[lupec,5] <- anova(fit)["as.factor(grp)","Pr(>F)"]
	
	ttestDVoutput[[lupec]] <- ttestboc(dd, varest=FALSE)			
}
cat('\n\n\n\nOne-way ANOVAs on the DVs (not part of the DFA; provided for comparisons with the ANOVAs on the DFs):\n\n')
dimnames(anovaDVoutput) <-list(rep("", dim(anovaDVoutput)[1]))
colnames(anovaDVoutput) <- c('Eta-squared','          F','    df','    df','        p')
rownames(anovaDVoutput) <- variables
anovaDVoutput[,1:4] <- round(anovaDVoutput[,1:4],2)
anovaDVoutput[,5] <- round(anovaDVoutput[,5],4)
print(anovaDVoutput)


cat('\n\n\n\nT-tests and effect sizes for group differences on the DFs:\n\n')
print(ttestDFoutput)


cat('\n\n\n\nT-tests and effect sizes for group differences on the DVs:\n')
cat('\n(not part of the DFA; provided for comparisons with the T-tests on the DFs)\n\n')
print(ttestDVoutput)


# plot
#colnames(centroidsZ) <-  c( "", paste("DF ", 1:(ncol(centroidSDsZ)-1), sep="") )
colnames(centroidsZ) <-  c(paste("DF ", 1:ncol(centroidSDsZ), sep="") )

matplot( 1:length(grpnames), centroidsZ[,1:2], type = "l", lty=1, lwd=3, 
        xaxt='n', xlab=groups, cex.axis=1.2, cex.lab = 1.3,
        ylab='Discriminant Function z Scores', ylim = c( -3, 3 ), cex.axis=1.2         )
axis(side=1, at=grpnames, labels=grpnames, xlab="groups")


# matplot(centroidsZ[,1], centroidsZ[,2:3], type = "l", lty=1, lwd=3, 
        # xaxt='n', xlab=groups, cex.axis=1.2, cex.lab = 1.3,
        # ylab='Discriminant Function z Scores', ylim = c( -3, 3 ), cex.axis=1.2         )
#axis(side = 1, at = x,labels = T, xlab='Group')
#axis(side=1, at=centroidsZ[,1], labels=centroidsZ[,1], xlab="groups")

title(main='Mean Standardized Discriminant Function Scores for the Groups')
#legend("topright", legend = colnames(centroidsZ[,2:3]), bty="n", lwd=2, col=1:ncol(centroidsZ[,2:3]) )
legend("topright", legend = colnames(centroidsZ), bty="n", lwd=2, col=1:ncol(centroidsZ) )



if (predictive == 'YES' | predictive == 'yes' | is.null(predictive) ) {

cat('\n\n\n\nPREDICTIVE DISCRIMINANT ANALYSIS:\n')

cat('\n\nPrior Probabilities for Groups:\n')
print(round(ldaoutput$prior,3))

# Frequencies: Original vs Predicted

# freqs_op <- data.frame( cbind( donnes[,1], lda.values$class ) )
# colnames(freqs_op) <-  c("Original", "Predicted") 
#freqs <- table(freqs_op) # doesn't produce a square table if the are no values for a factor level

squareTable <- function(x,y) {
    # x <- factor(x)
    # y <- factor(y)
#    commonLevels <- 1:max(nlevels(x), nlevels(y))
    # commonLevels <- 1:max(nlevels(x), nlevels(y))
    # x <- factor(x, levels = commonLevels)
    # y <- factor(y, levels = commonLevels)

    Original <- factor(x, levels = grpnames)
    Predicted <- factor(y, levels = grpnames)

    table(Original, Predicted)
}

freqs <- squareTable( donnes[,1], lda.values$class  )

cat('\n\nCross-Tabulation of the Original and Predicted Group Memberships:\n\n')
print(freqs)	

cat('\n\nProportion of original grouped cases correctly classified:  ',  round((sum(diag(freqs)) / sum(freqs)),3) )

cat('\n\n\nChi-square test of indepedence:\n')
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

# kappas
cat('\n\nAgreement (kappas) between the Predicted and Original Group Memberships:\n\n')
kappas_cvo <- kappas(na.omit(cbind( lda.values$class, donnes[,1] )) )
print(round(kappas_cvo,3))


# Frequencies: Original vs Cross-Validated (leave-one-out cross-validation)

# classifications from leave-one-out cross-validation
ldaoutputCV <- lda(x = as.matrix(donnes[,c(2:ncol(donnes))]), grouping=donnes[,1], prior=priors, CV = TRUE)

freqs_cvp <- data.frame( cbind( ldaoutputCV$class, lda.values$class  ) )
colnames(freqs_cvp) <-  c("Cross-Validated", "Predicted") 
# freqsCVP <- table(freqs_cvp) # doesn't produce a square table if the are no values for a factor level
freqsCVP <-  squareTable( ldaoutputCV$class, lda.values$class  )
# colnames(freqsCVP) <-  paste("Group ", 1:ngroups, sep="") 
# rownames(freqsCVP) <-  paste("Group ", 1:ngroups, sep="") 
colnames(freqsCVP) <- paste(grpnames)
rownames(freqsCVP) <- paste(grpnames)
cat('\n\n\n\nCross-Tabulation of the Cross-Validated and Predicted Group Memberships:\n\n')
print(freqsCVP)	

cat('\n\nProportion of cross-validated grouped cases correctly classified:  ',  round((sum(diag(freqsCVP)) / sum(freqsCVP)),3) )

cat('\n\nChi-square test of indepedence:\n')
print(summary(freqsCVP)) # chi-square test of indepedence

#mytable <- table(freqs_cvp) # A will be rows, B will be columns 
# mytable <- squareTable(freqs_cvp[,1], freqs_cvp[,2]) # A will be rows, B will be columns 
# colnames(mytable) <- paste(grpnames)
# rownames(mytable) <- paste(grpnames)
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
kappas_cvo <- kappas(na.omit(cbind( ldaoutputCV$class, donnes[,1] )) )
print(round(kappas_cvo,3))

# kappas
cat('\n\nAgreement (kappas) between the Cross-Validated and Predicted Group Memberships:\n\n')
kappas_cvp <- kappas(na.omit(cbind( ldaoutputCV$class, lda.values$class )) )
print(round(kappas_cvp,3))
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
   centroidSDsZ=centroidSDsZ  )

return(invisible(DFAoutput))

}
