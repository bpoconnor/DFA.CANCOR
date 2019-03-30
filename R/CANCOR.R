
CANCOR <- function(data, set1, set2, plot='yes', plotCV=1 ) {

cat('\n\nCanonical Correlation Analysis:\n\n')

data <- as.data.frame(data[,c(set1,set2)])
 
if (anyNA(data) == TRUE) {
	data <- na.omit(data)
	cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
}

set1data <- as.data.frame(data[,set1])
set2data <- as.data.frame(data[,set2])


# descriptive statistics & tests of univarite & multivariate normality -- from the MVN package
umvn(cbind(set1data,set2data))


# Pearson correlations
cat('\n\nPearson correlations for Set 1:\n\n');print(round(stats::cor(set1data),2))
cat('\n\nPearson correlations for Set 2:\n\n');print(round(stats::cor(set2data),2))
cat('\n\nPearson correlations between Set 1 & Set 2:\n\n');print(round(stats::cor(set2data,set1data),2))


CCoutput <- CCA::cc(set1data, set2data) # cc function from the CCA package


cat('\n\n\nMultivariate peel-down significance tests:\n\n')  # using p.asym from the CCP package

rho <- CCoutput$cor
N <- nrow(set1data)
p <- ncol(set1data)
q <- ncol(set2data)


invisible(utils::capture.output( pW <-  ((CCP::p.asym(rho, N, p, q, tstat = "Wilks"))) ))
dmat <- cbind(pW$stat, pW$approx, pW$df1, pW$df2, pW$p.value )
colnames(dmat) <- c('            Wilks Lambda', '   F-approx. ', '  df1', '    df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')

invisible(utils::capture.output( pH <- CCP::p.asym(rho, N, p, q, tstat = "Hotelling") ))
dmat <- cbind(pH$stat, pH$approx, pH$df1, pH$df2, pH$p.value )
colnames(dmat) <- c('  Hotelling-Lawley Trace', '   F-approx. ', '  df1', '    df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')

invisible(utils::capture.output( pP <- CCP::p.asym(rho, N, p, q, tstat = "Pillai") ))
dmat <- cbind(pP$stat, pP$approx, pP$df1, pP$df2, pP$p.value )
colnames(dmat) <- c('   Pillai-Bartlett Trace', '   F-approx. ', '  df1', '    df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')

invisible(utils::capture.output( pR <- CCP::p.asym(rho, N, p, q, tstat = "Roy") ))
dmat <- cbind(pR$stat, pR$approx, pR$df1, pR$df2, pR$p.value )
colnames(dmat) <- c('      Roy\'s Largest Root', '   F-approx. ', '  df1', '    df2', '         p')
rownames(dmat) <- paste(1:nrow(dmat), paste("through ", nrow(dmat), sep = ""))
print(round(dmat,4)); cat('\n\n')


# bivariate correlations for the canonical variates
bicortests <- matrix(-9999,ncol(CCoutput$scores$xscores),6)
for (nfunction in 1:ncol(CCoutput$scores$xscores)) {
	correl <- stats::cor(CCoutput$scores$xscores[,nfunction],CCoutput$scores$yscores[,nfunction] )		
	eigval <- correl**2 / (1 - correl**2)	
	dd <- stats::cor.test( CCoutput$scores$xscores[,nfunction],CCoutput$scores$yscores[,nfunction]) 	
	bicortests[nfunction,] <- cbind(eigval, correl, correl**2, dd$statistic, dd$parameter, dd$p.value)
}
cat('\n\n\nCanonical correlations:\n\n')
colnames(bicortests) <- c('  eigenvalue', '  canonical r', '   canonical r sq.','         t','     df','   p value')
rownames(bicortests) <- paste(" Canonical function ", 1:nrow(bicortests),sep = "")
print(round(bicortests,3))
cat('\nThe above t tests are for the Pearson correlations between the canonical variate scores')
cat('\nfor each function, i.e., they are not the multivariate significance tests.\n\n')


# raw canonical coefficients
cat('\n\n\n\nRaw canonical coefficients for Set 1:\n\n')
raw1 <- CCoutput$xcoef
colnames(raw1) <- paste("     CV", 1:ncol(raw1),sep = "")
print(round(raw1,2))
cat('\n\nRaw canonical coefficients for Set 2:\n\n')
raw2 <- CCoutput$ycoef
colnames(raw2) <- paste("     CV", 1:ncol(raw2),sep = "")
print(round(raw2,2))

# structure coefficients (canonical loadings)
cat('\n\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 1 variates:\n\n')
struct11 <- CCoutput$scores$corr.X.xscores
colnames(struct11) <- paste("     CV", 1:ncol(struct11),sep = "")
print(round(struct11,2))
cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 1 variates:\n\n')
struct21 <- CCoutput$scores$corr.Y.xscores
colnames(struct21) <- paste("     CV", 1:ncol(struct21),sep = "")
print(round(struct21,2))
cat('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 2 variates:\n\n')
struct12 <- CCoutput$scores$corr.X.yscores
colnames(struct12) <- paste("     CV", 1:ncol(struct12),sep = "")
print(round(struct12,2))
cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 2 variates:\n\n')
struct22 <- CCoutput$scores$corr.Y.yscores
colnames(struct22) <- paste("     CV", 1:ncol(struct22),sep = "")
print(round(struct22,2))

# standardized canonical coefficients for Set 1
s1 <- diag(sqrt(diag(stats::cov(set1data)))) #  diagonal matrix of sd's
rownames(s1) <- colnames(set1data)
cat('\n\n\nStandardized coefficients for Set 1 variables:\n\n')
stand1 <- s1 %*% CCoutput$xcoef
colnames(stand1) <- paste("     CV", 1:ncol(stand1),sep = "")
print(round(stand1,2))

# standardized canonical coefficients for Set 2
s2 <- diag(sqrt(diag(stats::cov(set2data)))) #  diagonal matrix of sd's
rownames(s2) <- colnames(set2data)
cat('\n\nStandardized coefficients for Set 2 variables:\n\n')
stand2 <- s2 %*% CCoutput$ycoef
colnames(stand2) <- paste("     CV", 1:ncol(stand2),sep = "")
print(round(stand2,2))


# helio plot -- from the yacca package
if (plot == 'yes' | plot == 'YES' | is.null(plot)) {
	if (is.null(plot)) plotCV = 1
#	cca.fit <- yacca::cca(set1data, set2data)
#	yacca::helio.plot(cca.fit, x.name="Set 1", y.name="Set 2", cv=plotCV)
	boc.fit <- list( xstructcorr=struct11, ystructcorr=struct22, xstructcorrsq=struct11**2,
                     ystructcorrsq=struct22**2, xlab=rownames(struct11), ylab=rownames(struct22) )
	yacca::helio.plot(boc.fit, x.name="Set 1", y.name="Set 2", cv=plotCV)
	cat('\n\n\nThe plot is provided by the yacca package. Helio plots display data in')
	cat('\nradial bars, with larger values pointing outward from a base reference circle')
	cat('\nand smaller (more negative) values pointing inward.')
}


CANCORoutput <- list(  
   cancors = bicortests,
   rawCoefSet1 = raw1,
   rawCoefSet2 = raw2,
   structCoef11 = struct11,
   structCoef21 = struct21,
   structCoef12 = struct12,
   structCoef22 = struct22,
   standCoefSet1 = stand1,
   standCoefSet2 = stand2  )

return(invisible(CANCORoutput))

cat('\n\n\n\n')
}

