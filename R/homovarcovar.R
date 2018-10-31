


homovarcovar <- function( donnes ) {

cat('\n\nTests for Homogeneity of Variances & Covariances:\n')

grpnames <- unique(donnes[,1])

donnes[,1] <- as.numeric((donnes[,1]))
ngroups <- max(donnes[,1])
donnes[,1] <- factor( donnes[,1], ordered = FALSE)
grpFreqs <- as.matrix(table(donnes[,c(1)]))

logdetgrps <- 0 # Box's M test
logdets <- matrix(-9999,ngroups,1) # for Box's M test
for (lupeg in 1:ngroups) {
	dum <- subset(donnes, donnes[,1] == lupeg)	
	cat('\nCovariance matrix for Group', paste(grpnames[lupeg]),'\n\n')
	print(round(stats::cov(dum[,2:ncol(dum)]),2))
	logdetgrps <- logdetgrps + (nrow(dum) - 1) * log(det(stats::cov(dum[,2:ncol(dum)]))) # for Box's M test
	logdets[lupeg,1] <- log(det(stats::cov(dum[,2:ncol(dum)]))) # for Box's M test
}

# Homogeneity of Variances
# Bartlett Test of Homogeneity of Variances (parametric, for K samples)
bb <- stats::bartlett.test(x=(donnes[,c(2:ncol(donnes))]), g=donnes[,1], data=donnes)
cat('\n\nBartlett Test of Homogeneity of Variances (parametric):\n')
cat('\nBartlett,s K-squared =', round(bb$statistic,3), '  df =', bb$parameter, '  p value =', round(bb$p.value,5) )

# Figner-Killeen Test of Homogeneity of Variances (non parametric, for K samples)
ff <- stats::fligner.test(donnes[,c(2:ncol(donnes))], donnes[,1], data=donnes)
cat('\n\n\nFigner-Killeen Test of Homogeneity of Variances (non parametric):\n')
cat('\nchi-squared =', round(ff$statistic,3), '  df =', ff$parameter, '  p value =', round(ff$p.value,5) )

# for the var-covar matrices from SPSS - requires Type III sums of squares
outcome <- as.matrix(donnes[,c(2:ncol(donnes))])
fit <- stats::lm(outcome ~ as.factor(donnes[,c(1)]), data = donnes)
sscp <- car::Anova(fit, type="III") #summary(sscp)
sscpwith <- sscp$SSPE
sscpbetw <- sscp$SSP$'as.factor(donnes[, c(1)])'

# the pooled within groups covariance matrix from SPSS:
poolwith <- sscpwith * (1/(nrow(donnes)-ngroups))
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
bigN <- nrow(donnes)
logdetpw <- log(det(poolwith)) 
BoxM <- ( (bigN - ngroups) * logdetpw ) - logdetgrps
ndvs <- ncol(outcome)
df1 <- ( (ngroups-1) * ndvs * (ndvs+1) ) / 2
c <- (((2*(ndvs**2))+(3*ndvs)-1) / (6*(ndvs+1)*(ngroups-1))) *
     (sum( (1/(grpFreqs-1))) - (1 / (bigN-ngroups)))   
c2 <- (((ndvs-1)*(ndvs+2)) / (6*(ngroups-1)))  * 
      (sum( (1/((grpFreqs-1)^2)) )  - (1 / ((bigN-ngroups)^2))) 
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


