

umvn <- function(ContinuousVariables) {

# descriptive statistics & tests of univarite & multivariate normality -- from the MVN package


# # uvn1 <- uniNorm(ContinuousVariables, type = "SW", desc = TRUE)  # Shapiro-Wilk test of univariate normality
# # uvn2 <- uniNorm(ContinuousVariables, type = "CVM", desc = FALSE)  # Cramer- von Mises test of univariate normality
# # uvn3 <- uniNorm(ContinuousVariables, type = "Lillie", desc = FALSE)  # Lilliefors (Kolmogorov-Smirnov) test of univariate normality
# # uvn4 <- uniNorm(ContinuousVariables, type = "SF", desc = FALSE)  # Shapiro-Francia test of univariate normality
# # uvn5 <- uniNorm(ContinuousVariables, type = "AD", desc = FALSE)  # Anderson-Darling test of univariate normality

# cat('\n\nDescriptive Statistics:\n\n')
# print(uvn1$"Descriptive Statistics")

# cat('\n\n\nShapiro-Wilk tests of univariate normality:\n\n')
# print(uvn1$"Shapiro-Wilk's Normality Test")

# mvn1 <- hzTest(ContinuousVariables,      qqplot = F)  # Henze-Zirkler's multivariate normality test 
# mvn2 <- mardiaTest(ContinuousVariables,  qqplot = F)  # Mardia's multivariate normality test 
# mvn3 <- roystonTest(ContinuousVariables, qqplot = F)  # Royston's multivariate normality test 

# cat('\n\n\n'); print(mvn1)
# cat('\n\n');   print(mvn2)
# cat('\n\n');   print(mvn3); cat('\n\n')


res1 = MVN::mvn(data = ContinuousVariables,  mvnTest = "mardia", univariateTest = "SW") 
res2 = MVN::mvn(data = ContinuousVariables,  mvnTest = "hz")
res3 = MVN::mvn(data = ContinuousVariables,  mvnTest = "royston")
res4 = MVN::mvn(data = ContinuousVariables,  mvnTest = "dh",  desc = FALSE)    


cat('\n\nDescriptive Statistics:\n\n')
descstats <- res1$"Descriptives"
descstats <- descstats[,-c(7,8)]
print(round(descstats,3))

cat('\n\n\nShapiro-Wilk tests of univariate normality:\n\n')
print(res1$"univariateNormality")

cat('\n\n\nTests of multivariate normality:\n\n')
print(res1$"multivariateNormality"[1:2,], row.names = F); cat('\n\n')
print(res2$"multivariateNormality", row.names = F); cat('\n\n')
print(res3$"multivariateNormality", row.names = F); cat('\n\n')
print(res4$"multivariateNormality", row.names = F); cat('\n\n')

}






kappa.cohen <- function (kapdon) {

	# the data for this function (kapdon) are the category values for each of 2 columns,
	# but the analyses are performed on the contingency table	

	kapdonCT <- table(kapdon[,1],kapdon[,2])  

	# based on Valiquette 1994 BRM, Table 1
	n <- sum(kapdonCT)  # Sum of Matrix elements
	kapdonP <- kapdonCT / n  # table of proportions
	po <- sum(diag(kapdonP))
	c <- rowSums(kapdonP)
	r <- colSums(kapdonP)
	pe <- r %*% c
	num <- po - pe
	den <- 1 - pe
	kappa <- num / den

	# SE and variance of kappa from Cohen 1968
	sek <- sqrt((po*(1-po))/(n*(1-pe)^2))
		
	if (n < 100) var <- pe/(n*(1-pe))  # kappa variance as reported by Cohen in his original work
	
	if (n >= 100) {
		s <- t(matrix(colSums(kapdonP))) # columns sum
		var <- (pe+pe^2-sum(diag(r%*%s) * (r+t(s)))) / (n*(1-pe)^2)
		# asymptotic kappa variance as reported by 
		# Fleiss, J. L., Lee, J. C. M., & Landis, J. R. (1979).
	    # The large sample variance of kappa in the case of different sets of raters. 
	    # Psychological Bulletin, 86, 974-977
	}

	zkappa <- kappa / sqrt(var)
	
	sig <- round(pnorm(abs(zkappa),lower.tail = F),5) * 2 # 2-tailed test

#	print( c(kappa, sek, var, zkappa, sig))

	return(invisible(c(kappa, zkappa, sig)))

	# # based on kappa.m	
	# n <- sum(kapdon) # Sum of Matrix elements	
	# kapdonP <- kapdon/n # proportion		
	# r <- matrix(rowSums(kapdonP)) # rows sum
	# s <- t(matrix(colSums(kapdonP))) # columns sum	
	# Ex <- (r) %*% s # expected proportion for random agree
	# f <- diag(1,3)	
	# # pom <- apply(rbind(t(r),s),2,min)	
	# po <- sum(sum(kapdonP * f))  # sum(sum(x.*f))
	# pe <- sum(sum(Ex * f))
	# k <- (po-pe)/(1-pe)
	# # km <- (pom-pe)/(1-pe) # maximum possible kappa, given the observed marginal frequencies
	# # ratio <- k/km # observed as proportion of maximum possible	
	# # kappa standard error for confidence interval as reported by Cohen in his original work
	# sek <- sqrt((po*(1-po))/(n*(1-pe)^2))	
	# var <- pe/(n*(1-pe))  # kappa variance as reported by Cohen in his original work
	# # var <- (pe+pe^2-sum(diag(r*s).*(r+s')))/(n*(1-pe)^2)  # for N > 100	
	# zk  <-  k/sqrt(var)

}



# Fleiss's kappa

# source: https://en.wikipedia.org/wiki/Fleiss%27_kappa

# Fleiss's kappa is a generalisation of Scott's pi statistic, a
# statistical measure of inter-rater reliability. It is also related to
# Cohen's kappa statistic. Whereas Scott's pi and Cohen's kappa work for
# only two raters, Fleiss's kappa works for any number of raters giving
# categorical ratings (see nominal data), to a fixed number of items. It
# can be interpreted as expressing the extent to which the observed amount
# of agreement among raters exceeds what would be expected if all raters
# made their ratings completely randomly. Agreement can be thought of as
# follows, if a fixed number of people assign numerical ratings to a number
# of items then the kappa will give a measure for how consistent the
# ratings are. The scoring range is between 0 and 1. 

# Conger, A.J. (1980). Integration and generalisation of Kappas for multiple raters. Psychological Bulletin, 88, 322-328. 
# Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bul- letin, 76, 378-382. 
# Fleiss, J.L., Levin, B., & Paik, M.C. (2003). Statistical Methods for Rates and Proportions, 3rd Edition. New York: John Wiley & Sons. 

kappa.fleiss <- function(kapdon) {
	
	# the data for this function (kapdon) are the category values for each column,
	# but the analyses are performed on a count matrix (not a contin table) = the fleissmat below	
	fleissmat <- matrix(0,nrow(kapdon),max(kapdon))
	for (luper in 1:nrow(kapdon)) {
		for (lupec in 1:ncol(kapdon)) {
			fleissmat[luper,kapdon[luper,lupec]] <- fleissmat[luper,kapdon[luper,lupec]] + 1				
		}
	}	
	n <- nrow(fleissmat) 
	m <- sum(fleissmat[1,]) 	
	a <- n * m		
	pj <- colSums(fleissmat) / a 
	b <- pj * (1-pj)
	c <- a*(m-1)
	d <- sum(b)
	kj <- 1-(colSums((fleissmat * (m-fleissmat))) / (c %*% b)) # the value of kappa for the j-th category
	# sekj <- sqrt(2/c) 
	# zkj <- kj / sekj
	# pkj <- round(pnorm(abs(zkj),lower.tail = F),5) * 2  # 2-tailed test
	k <- sum(b*kj) / d  # Fleiss's (overall) kappa
	sek <- sqrt(2*(d^2-sum(b *(1-2 *pj))))/sum(b *sqrt(c)) 	
	#ci <- k+(c(-1,1) * (pnorm(zkj)) * sek) 	
	zk <- k / sek  # normalized kappa
	sig <- round(pnorm(abs(zk),lower.tail = F),5) * 2  # 2-tailed test
	
#	print( c(k, sek, zk, sig))
	
	return(invisible(c(k, zk, sig)))
}



kappas <- function(grpdat) {
	
	kc <- kappa.cohen(grpdat)  
	kf <- kappa.fleiss(grpdat)
	kappasOUT <- rbind(kc,kf)
	kappasOUT[,1:2] <- round(kappasOUT[,1:2],3)
	kappasOUT[,3] <- round(kappasOUT[,3],5)	
	rownames(kappasOUT) <- c( "Cohen's kappa", "Fleiss's kappa")
	colnames(kappasOUT) <- c( "    kappa", "        z", "         p" )
		
	# k2 <- irr::kappa2( grpdat )  # Unweighted Kappa for categorical data without a logical order
	# kf <- irr::kappam.fleiss( grpdat, exact = FALSE, detail = TRUE)
	# kl <- irr::kappam.light( grpdat )
	# kappasOUT <- matrix( -9999, 3, 4)
	# kappasOUT[1,] <- cbind( k2$subjects, k2$value, k2$statistic, k2$p.value)
	# kappasOUT[2,] <- cbind( kl$subjects, kl$value, kl$statistic, kl$p.value)
	# kappasOUT[3,] <- cbind( kf$subjects, kf$value, kf$statistic, kf$p.value)
	# rownames(kappasOUT) <- c( "Cohen's kappa", "Light's kappa", "Fleiss's kappa")
	# colnames(kappasOUT) <- c( "   N", "    kappa", "         z", "      p" )

	return (kappasOUT)
}









############################# T tests  ############################################################


"ttestboc" <-  function (donnesT, varest=FALSE) {
	
# reads raw data; the groups variable is in 1st column; the DV(s) are in the subsequent columns
# the groups variable can be categorical
# the function compares the lowest & highest values of the group variable

# var.equal -- a logical variable indicating whether to treat the two variances as being equal. 
# If TRUE then the pooled variance is used to estimate the variance otherwise the  
# Welch (or Satterthwaite) approximation to the degrees of freedom is used.

# p.adjust.methods options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

#cat("\n\nGroups t-test:\n")

grpnames <- as.vector(as.matrix(donnesT[,1])) # group names, in the same order as in the data matrix
grpnames <- unique(grpnames)
grpnums  <- seq(1:length(grpnames))

donnesT[,1] <- as.numeric(donnesT[,1])

resultsM <- matrix(-9999,1,13)
ngroups <- max(donnesT[,1])

for (lupe1 in 1:(ngroups-1)) {
	for (lupe2 in (lupe1+1):ngroups) {

		dum <- subset(donnesT, donnesT[,1] == lupe1 | donnesT[,1] == lupe2  )
		
		newdon = data.frame(dum[,1], dum[,2])


#		newdon <- stats::na.omit(cbind( as.numeric(donnesT[,1]), donnesT[,2] ))

		groupmin <- min(newdon[,1])
		groupmax <- max(newdon[,1])

		mgrp1 <- mean(subset(newdon[,2],newdon[,1]==groupmin))
		mgrp2 <- mean(subset(newdon[,2],newdon[,1]==groupmax))

		sdgrp1 <- stats::sd(subset(newdon[,2],newdon[,1]==groupmin))
		sdgrp2 <- stats::sd(subset(newdon[,2],newdon[,1]==groupmax))

		N1 <- nrow(subset(newdon,newdon[,1]==groupmin))
		N2 <- nrow(subset(newdon,newdon[,1]==groupmax))

		SE1 <- sdgrp1 / sqrt(N1)
		SE2 <- sdgrp2 / sqrt(N2)

		tresults <- stats::t.test(newdon[,2]~newdon[,1],data=newdon, var.equal=varest) 
		tgroups  <- tresults$statistic
		dfgroups <- tresults$parameter
		plevel   <- tresults$p.value

		# r effect size
		reffsiz =  sqrt( tgroups**2 / (tgroups**2 + dfgroups)  )

		# d effect size -- from R&R p 303 General Formula  best, because covers = & not = Ns
		deffsiz = (tgroups * (N1+N2)) / ( sqrt(dfgroups) * sqrt(N1*N2) )

		results <- cbind( groupmin, N1, round(mgrp1,2), round(sdgrp1,2), groupmax, N2, 
						  round(mgrp2,2), round(sdgrp2,2), round(tgroups,2), round(dfgroups,2), 
						  round(plevel,5), round(reffsiz,2), round(abs(deffsiz),2) )

		results <- as.matrix(cbind(results))
		resultsM <- rbind( resultsM, results)
	}  	
}

resultsM2 <- data.frame(resultsM)

for (lupe in 2:nrow(resultsM)) {
	resultsM2[lupe,1] <- grpnames[resultsM[lupe,1]]
	resultsM2[lupe,5] <- grpnames[resultsM[lupe,5]]
}
rownames(resultsM2) <- c()
colnames(resultsM2) <- c("Group"," N1"," Mean1","  SD1"," Group"," N2"," Mean2","  SD2",
                         "     t","    df","        p","  r effsize","  d effsize")

resultsM2 <- resultsM2[-c(1),]

return(invisible(as.data.frame(resultsM2)))
}








homovarcovar <- function( donnes ) {

cat('\n\nTests for Homogeneity of Variances & Covariances:\n')

grpnames <- as.vector(as.matrix(donnes[,1])) # group names, in the same order as in the data matrix
grpnames <- unique(grpnames)
ngroups  <- length(grpnames)
nDVs <- ncol(donnes) - 1
N <- nrow(donnes)

if (is.factor(donnes[,1]) == F)  donnes[,1] <- factor( donnes[,1], ordered = FALSE, labels=grpnames)

grpFreqs <- as.matrix(table(donnes[,c(1)]))

logdetgrps <- 0 # Box's M test
logdets <- matrix(-9999,ngroups,1) # for Box's M test
for (lupeg in 1:ngroups) {
	dum <- subset(donnes, donnes[,1] == grpnames[lupeg] )
	cat('\nCovariance matrix for Group', paste(grpnames[lupeg]),'\n\n')
	print(round(stats::cov(dum[,2:ncol(dum)]),2))
	logdetgrps <- logdetgrps + (nrow(dum) - 1) * log(det(stats::cov(dum[,2:ncol(dum)]))) # for Box's M test
	logdets[lupeg,1] <- log(det(stats::cov(dum[,2:ncol(dum)]))) # for Box's M test
}

# Homogeneity of Variances
# Bartlett Test of Homogeneity of variance-covariance matrices (parametric, for K samples)
bb <- stats::bartlett.test(x=(donnes[,c(2:ncol(donnes))]), g=donnes[,1], data=donnes)
cat('\n\nBartlett Test of Homogeneity of Variances (parametric):\n')
cat('\nBartlett,s K-squared =', round(bb$statistic,3), '  df =', bb$parameter, '  p value =', round(bb$p.value,5) )

# Figner-Killeen Test of Homogeneity of Variances (non parametric, for K samples)
ff <- stats::fligner.test(donnes[,c(2:ncol(donnes))], donnes[,1], data=donnes)
cat('\n\n\nFigner-Killeen Test of Homogeneity of Variances (non parametric):\n')
cat('\nchi-squared =', round(ff$statistic,3), '  df =', ff$parameter, '  p value =', round(ff$p.value,5) )

# for the var-covar matrices from SPSS - requires Type III sums of squares
# using MANOVA to obtain the sums of squares and cross-products matrix for error, &
# the sums of squares and cross-products matrix for group/IV
# www.webpages.uidaho.edu/~kirk/calendar/R/MANOVA.doc
MV2 <-manova( as.matrix(donnes[,2:ncol(donnes)]) ~ donnes[,1], data=donnes)
sscpwith <- (N-1)*cov(MV2$residuals) # E
sscpbetw <- (N-1)*cov(MV2$fitted.values) # H

# # using Anova from the car package to obtain the sums of squares and cross-products matrices
# outcome <- as.matrix(donnes[,c(2:ncol(donnes))])
# fit <- stats::lm(outcome ~ as.factor(donnes[,c(1)]), data = donnes)
# sscp <- car::Anova(fit, type="III") #summary(sscp)
# sscpwith <- sscp$SSPE
# sscpbetw <- sscp$SSP$'as.factor(donnes[, c(1)])'

# # another approach
# options(contrasts = c('contr.sum','contr.poly'))
# Next, store the model:
# model <- lm(time ~ topic * sys, data=search)
# Finally, call the drop1 function on each model component:
# drop1(model, .~., test='F')
# The results give the type III SS, including the p-values from an F-test.



# the pooled within groups covariance matrix from SPSS:
poolwith <- sscpwith * (1/(N-ngroups))
cat('\n\n\nPooled Within Groups Covariance Matrix from SPSS:\n')
print(round(poolwith,3))

# the pooled within groups correlation matrix from SPSS:
bigDv <- diag(suppressWarnings(sqrt(1 / poolwith)))
scalemat <- diag(bigDv)
cormat <- t(scalemat) %*% poolwith %*% scalemat
rownames(cormat) <- rownames(poolwith)
colnames(cormat) <- colnames(poolwith)
cat('\n\nPooled Within Groups Correlation Matrix from SPSS:\n')
print(round(cormat,3))

# Box's M test
# formulas from http://www.real-statistics.com/multivariate-statistics/boxs-test-equality-covariance-matrices/boxs-test-basic-concepts/
# IBM SPSS Statistics Algorithms 20
logdetpw <- log(det(poolwith)) 
BoxM <- ( (N - ngroups) * logdetpw ) - logdetgrps
df1 <- ( (ngroups-1) * nDVs * (nDVs+1) ) / 2
c <- (((2*(nDVs**2))+(3*nDVs)-1) / (6*(nDVs+1)*(ngroups-1))) *
     (sum( (1/(grpFreqs-1))) - (1 / (N-ngroups)))   
c2 <- (((nDVs-1)*(nDVs+2)) / (6*(ngroups-1)))  * 
      (sum( (1/((grpFreqs-1)^2)) )  - (1 / ((N-ngroups)^2))) 
df2 <- (df1 + 2) / abs(c2-c**2)
aplus <- df1 / ( 1 - c - (df1/df2))
Fplus <- BoxM / aplus
aminus <- df2 / (1 - c + (2/df2))
Fminus <- (df2 * BoxM) / (df1 * (aminus - BoxM))
if (c2 > c**2)  bigF <- Fplus
if (c2 < c**2)  bigF <- Fminus
pbigF <- stats::pf(bigF, df1, df2, lower.tail=FALSE) # p level
cat('\n\nBox Test of Equality of Covariance Matrices:\n')
cat('\nLog Determinants:\n')
#rownames(logdets) <- paste("Group ", 1:ngroups, sep="") 
rownames(logdets) <- paste(grpnames) 
colnames(logdets) <- 'Log Determinant'
Pooled <- logdetpw
logdets <- rbind( logdets, Pooled)
print(round(logdets,3))
cat('\n\nM =', round(BoxM,3), '  F =', round(bigF,3), '  df1 =', df1, '  df2 =', df2, '  p = ', pbigF, '\n\n\n')

sscpOutput <- list(sscpwith=sscpwith, sscpbetw=sscpbetw)

return(invisible(sscpOutput))

}







# # sources for the following multivariate significance tests:

# Manly, B. F. J., & Alberto, J. A. (2017). Multivariate Statistical 
# Methods: A Primer (4th Edition). Chapman & Hall/CRC, Boca Raton, FL.

# https://rpubs.com/aaronsc32/manova-test-statistics

# http://www.statpower.net/Software.html

# Gatignon, H. (2014). Statistical Analysis of Management Data. New York, NY: Springer.




Wilks <- function(rho, Ncases, p, q) {
	
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			
	
	m <- Ncases - 3/2 - (p + q)/2
	
	wilks <- Fwilks <- df1 <- df2 <- vector("numeric", NCVs) 

	for (lupe in 1:NCVs) {
		
		wilks[lupe] <- prod(1 / (1 + eigval[lupe:NCVs]))	
		
		t <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
		
		df1[lupe] <- p * q
		df2[lupe] <- m * t - p * q/2 + 1
		
		Fwilks[lupe] <- (1 - wilks[lupe]^(1/t)) / wilks[lupe]^(1/t) * df2[lupe] / df1[lupe]
		
		p <- p - 1
		q <- q - 1
	}
	
	pvalue <- pf(Fwilks, df1, df2, lower.tail = FALSE)
	
	WilksOutput <- cbind(wilks,Fwilks,df1,df2,pvalue)
	return(invisible(WilksOutput))
}






Pillai <- function(rho, Ncases, p, q) {
		
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	pManly  <- max(p,q) # the # of variables in the larger set
	pManly2 <- p + q  # the total # of variables

	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			

	# 2017 Manly p 69 -- works 
	m = 1 # number of samples
	m2 = 1 # number of samples
	d = max(pManly,(m-1)) # the greater of pManly and m - 1
	
	Pillai <- Fpillai <- df1 <- df2 <- vector("numeric", NCVs)
	s2 <- NCVs
	
	for (lupe in 1:NCVs) {
		
		Pillai[lupe] <- sum(eigval[lupe:NCVs] / (1 + eigval[lupe:NCVs]))  
		# https://rpubs.com/aaronsc32/manova-test-statistics	
		
		s <- NCVs + 1 - lupe  # number of positive eigenvalues
		
		d = max(pManly,(m-1)) # update d
		df1[lupe] <- s * d
		pManly = pManly - 1
		
		df2[lupe] <- s2 * (Ncases - m2 - pManly2 + s2)
	
		Fpillai[lupe] <- ((Ncases - m2 - pManly2 + s2) * Pillai[lupe]) / ( d * (s - Pillai[lupe]))
		
		m2 = m2 - 2
	}
	
	pvalue <- pf(Fpillai, df1, df2, lower.tail = FALSE)
		
	PillaiOutput <- cbind(Pillai,Fpillai,df1,df2,pvalue)
	return(invisible(PillaiOutput))
}





Hotelling <- function(rho, Ncases, p, q) {
			
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	pManly  <- abs(p - q) # p  = 2 # number of variables  8 - 6?
	pManly2 <- p + q  # the total # of variables

	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			

	# 2017 Manly p 69 
	
	Hotelling <- Fhotelling <- df1 <- df2 <- vector("numeric", NCVs)
	s2 <- NCVs
	m  <- 1 # number of samples
	m2 <- 1 # number of samples

	for (lupe in 1:NCVs) {
		
		Hotelling[lupe] <- sum(eigval[lupe:NCVs])  # https://rpubs.com/aaronsc32/manova-test-statistics
			
		s <- NCVs + 1 - lupe  # number of positive eigenvalues
				
		A <- ( abs(m - pManly - 1) - 1) / 2
		df1[lupe] <- s * (2 * A + s + 1)
				
		B <- (Ncases - m2 - pManly2 - 1) / 2
		df2[lupe] <- 2 * (s2 * B + 1)
		m2 = m2 - 2
		
		Fhotelling[lupe] <- df2[lupe] * (Hotelling[lupe] / (s2 * df1[lupe]))
		}

	pvalue <- pf(Fhotelling, df1, df2, lower.tail = FALSE)
	HotellingOutput <- cbind(Hotelling,Fhotelling,df1,df2,pvalue)	
	return(invisible(HotellingOutput))
}




RoyRoot <- function(rho, Ncases, p, q) {

	# rho = the first canonical correlation (or all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	NDVs <- p + q

	eval <- rho[1]**2 / (1 - rho[1]**2)  # convert rho (canonical r) to eigenvalue
	# based on 2017 Manly p 69 
	m = 1
	d = min(p,q) # should be d <- max(p,(m-1)), but m (# IVs) is 1 for 1-way MANOVA
	df1 = d
	df2 = Ncases - m - d   # - 1
	Froy = df2 / df1 * eval
	pvalue <- pf(Froy, df1, df2, lower.tail = FALSE)
	RoyRootOutput <- cbind(rho[1]**2,Froy,df1,df2,pvalue)
	return(invisible(RoyRootOutput))
}







# Bartlett's V test (peel-down) of the significance of canonical correlations (for CCA, not DFA)

BartlettV <- function(rho, Ncases, p, q) {
	
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
		
	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			
		
	wilks <- X2 <- df <- pvalue <- vector("numeric", NCVs) 

	for (lupe in 1:NCVs) {
		
		wilks[lupe] <- prod(1 / (1 + eigval[lupe:NCVs]))	
		
		X2[lupe] <- -( (Ncases - 1) - (p + q + 1)/2) * log(wilks[lupe])
		          
		df[lupe]  <- (p - lupe +1) * (q - lupe + 1) 
		          
		pvalue[lupe] <- pchisq(X2[lupe], df[lupe], ncp=0, lower.tail = F) 
	}
		
#	print(cbind(wilks,X2,df,pvalue))
	BartlettVOutput <- cbind(round(wilks,2),round(X2,2),df,pvalue)
	return(invisible(BartlettVOutput))
}





# Rao's test (peel-down) of the significance of canonical correlations (for CCA, not DFA)

Rao <- function(rho, Ncases, p, q) {
	
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
		
	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			
		
	wilks <- Frao <- df1 <- df2 <- pvalue <- vector("numeric", NCVs) 

	for (lupe in 1:NCVs) {		
		wilks[lupe] <- prod(1 / (1 + eigval[lupe:NCVs]))			
		pnew <- p - lupe + 1
		qnew <- q - lupe + 1
		t <- (Ncases - 1) - (pnew + qnew + 1) / 2       
		s <- ifelse((pnew^2 + qnew^2) <= 5, 1, sqrt((pnew^2 * qnew^2 -4) / (pnew^2 + qnew^2 -5)))              
		df1[lupe] <- pnew * qnew 		
		df2[lupe] <- (1 + t*s - pnew*qnew/2)		
		Frao[lupe] <- ((1 - wilks[lupe]^(1/s)) / wilks[lupe]^(1/s)) * df2[lupe] / df1[lupe] 
		pvalue[lupe] <- suppressWarnings(pf(Frao[lupe], df1[lupe], df2[lupe], ncp=0, lower.tail = FALSE))
    }
#	print(cbind(wilks,Frao,df1,df2,pvalue))
	RaoOutput <- cbind(round(wilks,2),round(Frao,2),df1,round(df2,2),round(pvalue,5))
	return(invisible(RaoOutput))
}
  
   
   
  
  