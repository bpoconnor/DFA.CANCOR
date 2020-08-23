
# # CANCOR <- function(obj) {
# UseMethod("CANCOR")
# }



CANCOR <- function(data, set1, set2, plot='yes', plotCV=1, plotcoefs='structure' ) {

#UseMethod("CANCOR")


cat('\n\nCanonical Correlation Analysis:\n\n')

data <- as.data.frame(data[,c(set1,set2)])
 
if (anyNA(data) == TRUE) {
	data <- na.omit(data)
	cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
}

set1data <- as.data.frame(data[,set1])
set2data <- as.data.frame(data[,set2])

Ncases <- nrow(set1data)
NVset1 <- ncol(set1data)
NVset2 <- ncol(set2data)

# descriptive statistics & tests of univariate & multivariate normality
umvn(cbind(set1data,set2data))


# Pearson correlations
cat('\n\nPearson correlations for Set 1:\n\n');print(round(stats::cor(set1data),2))
cat('\n\nPearson correlations for Set 2:\n\n');print(round(stats::cor(set2data),2))
cat('\n\nPearson correlations between Set 1 & Set 2:\n\n');print(round(stats::cor(set2data,set1data),2))

# the CCA
output <- canonical.cor(set1data, set2data) # source("http://www.statpower.net/R312/CanCorr.r")


cat('\n\n\nMultivariate peel-down significance tests:\n\n')

#print(round(output$'Canonical Correlations',4))

cancorrels <- output$cancorrels

mvFmat <- Wilks(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mvFmat) <- c('            Wilks Lambda', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')

mvFmat <- Pillai(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mvFmat) <- c('  Pillai-Bartlett Trace', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')

mvFmat <- Hotelling(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mvFmat) <- c('   Hotelling-Lawley Trace', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')

mvFmat <- RoyRoot(rho = cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mvFmat) <- c('      Roy\'s Largest Root', '   F-approx. ', '  df1', '    df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')

mvX2mat <- BartlettV(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
cat('\n\n\nBartlett\'s V test:\n')
colnames(mvX2mat) <- c('            Wilks Lambda', '   F-approx. ', '  df', '         p')
rownames(mvX2mat) <- paste(1:nrow(mvX2mat), paste("through ", nrow(mvX2mat), sep = ""))
print(round(mvX2mat,4)); cat('\n\n')

mvFmat <- Rao(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
cat('\n\n\nRao\'s V test:\n')
colnames(mvFmat) <- c('   Wilks Lambda', '   F-approx. ', '    df1', '       df2', '         p')
rownames(mvFmat) <- paste(1:nrow(mvFmat), paste("through ", nrow(mvFmat), sep = ""))
print(round(mvFmat,4)); cat('\n\n')



# bivariate correlations for the canonical variates
cat('\n\n\nCanonical correlations:\n\n')
colnames(output$cancorrels) <- c('  Eigenvalue', '  Canonical r', '   Canonical r sq.','         t','     df','   p value')
rownames(output$cancorrels) <- paste(" Canonical function ", 1:nrow(output$cancorrels),sep = "")
print(round(output$cancorrels,3))
cat('\nThe above t tests are for the Pearson correlations between the canonical variate scores')
cat('\nfor each function, i.e., they are not the multivariate significance tests.\n\n')

  
# raw canonical coefficients
cat('\n\n\n\nRaw canonical coefficients for Set 1:\n\n')
colnames(output$raw1) <- paste("     CV", 1:ncol(output$raw1),sep = "")
print(round(output$raw1,2))

cat('\n\nRaw canonical coefficients for Set 2:\n\n')
colnames(output$raw2) <- paste("     CV", 1:ncol(output$raw2),sep = "")
print(round(output$raw2,2))


# structure coefficients (canonical loadings)
cat('\n\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 1 variates:\n\n')
colnames(output$struct11) <- paste("     CV", 1:ncol(output$struct11),sep = "")
print(round(output$struct11,2))

cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 1 variates:\n\n')
colnames(output$struct21) <- paste("     CV", 1:ncol(output$struct21),sep = "")
print(round(output$struct21,2))

cat('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 2 variates:\n\n')
colnames(output$struct12) <- paste("     CV", 1:ncol(output$struct12),sep = "")
print(round(output$struct12,2))

cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 2 variates:\n\n')
colnames(output$struct22) <- paste("     CV", 1:ncol(output$struct22),sep = "")
print(round(output$struct22,2))


# standardized canonical coefficients 
cat('\n\n\nStandardized coefficients for Set 1 variables:\n\n')
colnames(output$stand1) <- paste("     CV", 1:ncol(output$stand1),sep = "")
print(round(output$stand1,2))

cat('\n\nStandardized coefficients for Set 2 variables:\n\n')
colnames(output$stand2) <- paste("     CV", 1:ncol(output$stand2),sep = "")
print(round(output$stand2,2))



if (plot == 'yes' | plot == 'YES' | is.null(plot)) {
	if (is.null(plot)) plotCV = 1

	# bar plots - structure coefficients
	if (plotcoefs == 'structure' | is.null(plotcoefs)) {
		par(mfrow=c(1,2), mar = c(9,5,3,3)) # set the margin on all sides
		barplot(output$struct11[,plotCV], ylim=c(-1,1), ylab='CV1', 
		        main='Set 1 Structure Coefficients', col='blue', las=2)
		box()	
		barplot(output$struct22[,plotCV], ylim=c(-1,1), ylab='CV1', 
		        main='Set 2 Structure Coefficients', col='blue', las=2)
		box()
	}
		
	# bar plots - standardized coefficients
	if (plotcoefs == 'standardized') {
		par(mfrow=c(1,2), mar = c(9,5,3,3)) # set the margin on all sides
		barplot(output$stand1[,plotCV], ylim=c(-1.3,1.3), ylab='CV1', 
		        main='Set 1 Standardized Coefficients', col='blue', las=2)
		box()	
		barplot(output$stand2[,plotCV], ylim=c(-1.3,1.3), ylab='CV1', 
		        main='Set 2 Standardized Coefficients', col='blue', las=2)
		box()
	}
	#dev.off()  

	# # helio plot -- from the yacca package
	# cca.fit <- yacca::cca(set1data, set2data)
	# yacca::helio.plot(cca.fit, x.name="Set 1", y.name="Set 2", cv=plotCV)
	# boc.fit <- list( xstructcorr=struct11, ystructcorr=struct22, xstructcorrsq=struct11**2,
                     # ystructcorrsq=struct22**2, xlab=rownames(struct11), ylab=rownames(struct22) )
	# yacca::helio.plot(boc.fit, x.name="Set 1", y.name="Set 2", cv=plotCV)
	# cat('\n\n\nThe plot is provided by the yacca package. Helio plots display data in')
	# cat('\nradial bars, with larger values pointing outward from a base reference circle')
	# cat('\nand smaller (more negative) values pointing inward.')
}


CANCORoutput <- list(  
   cancorrels = cancorrels,
   rawCoefSet1 = output$raw1,
   rawCoefSet2 = output$raw2,
   structCoef11 = output$struct11,
   structCoef21 = output$struct21,
   structCoef12 = output$struct12,
   structCoef22 = output$struct22,
   standCoefSet1 = output$stand1,
   standCoefSet2 = output$stand2  )

#return(invisible(CANCORoutput))

class(CANCORoutput) <- "CANCORout"


CANCORoutput


cat('\n\n\n\n')
}





# # CANCOR.CANCORout <-
# function (x,  ...) {
	
    # if (!inherits(x, "CANCORout"))
        # stop("Use only with 'CANCORout' objects.\n")

	# cat('\n\nStandardized coefficients for Set 2 variables:\n\n')
# #	colnames(output$stand2) <- paste("     CV", 1:ncol(output$stand2),sep = "")
	# print(round(x$standCoefSet2,6))
	# print(x)

	
# # # 	
    # # if (!inherits(x, "cronbachAlpha"))
        # # stop("Use only with 'cronbachAlpha' objects.\n")
    # # if (x$standardized)
        # # cat("\nStandardized Cronbach's alpha for the", paste("'", x$name, "'", sep = ""), "data-set\n")
    # # else
        # # cat("\nCronbach's alpha for the", paste("'", x$name, "'", sep = ""), "data-set\n")
    # # cat("\nItems:", x$p)
    # # cat("\nSample units:", x$n)
    # # cat("\nalpha:", round(x$alpha, digits))
    # # if (!is.null(x$ci)) {
        # # cat("\n\nBootstrap ", 100 * diff(x$probs), "% CI based on ", x$B, " samples\n", sep = "")
        # # print(round(x$ci, digits))
    # # }
    # # cat("\n\n")
    # invisible(x)
# }




