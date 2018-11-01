


############################# T tests  ############################################################


"ttestboc" <-  function (donnesT, varest=FALSE) {
	
# reads raw data; the groups variable is in 1st column; the DV(s) are in the subsequent columns
# the groups variable can be categorical
# the function compares the lowest & highest values of the group variable

# var.equal -- a logical variable indicating whether to treat the two variances as being equal. 
# If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) 
# approximation to the degrees of freedom is used.

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

		tresults <- stats::t.test(newdon[,2]~newdon[,1],data=newdon, var.equal=varest); # t-test
		tgroups  <- tresults$statistic
		dfgroups <- tresults$parameter
		plevel   <- tresults$p.value

		# r effect size
		reffsiz =  sqrt( tgroups**2 / (tgroups**2 + dfgroups)  )

		# d effect size -- from R&R p 303 General Formula  best, because covers = & not = Ns
		deffsiz = (tgroups * (N1+N2)) / ( sqrt(dfgroups) * sqrt(N1*N2) )

		results <- cbind( groupmin, N1, round(mgrp1,2), round(sdgrp1,2), groupmax, N2, round(mgrp2,2), round(sdgrp2,2),
		           round(tgroups,2), round(dfgroups,2), round(plevel,5), round(reffsiz,2), round(deffsiz,2) )

		results <- as.matrix(cbind(results))
		resultsM <- rbind( resultsM, results)
	}  	
}
resultsM <- resultsM[-1,]

resultsM2 <- data.frame(resultsM)
for (lupe in 1:nrow(resultsM)) {
	resultsM2[lupe,1] <- grpnames[resultsM[lupe,1]]
	resultsM2[lupe,5] <- grpnames[resultsM[lupe,5]]
}
rownames(resultsM2) <- c()
colnames(resultsM2) <- c("Group","  N1"," Mean1","  SD1","   Group","  N2"," Mean2","  SD2","         t","      df","        p","   r effsize","  d effsize")

return(invisible(resultsM2))
}