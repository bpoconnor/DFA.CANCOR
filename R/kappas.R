


kappas <- function(grpdat) {
k2 <- kappa2( grpdat )  # Unweighted Kappa for categorical data without a logical order
kf <- kappam.fleiss( grpdat, exact = FALSE, detail = TRUE)
kl <- kappam.light( grpdat )
kappasOUT <- matrix( -9999, 3, 4)
kappasOUT[1,] <- cbind( k2$subjects, k2$value, k2$statistic, k2$p.value)
kappasOUT[2,] <- cbind( kl$subjects, kl$value, kl$statistic, kl$p.value)
kappasOUT[3,] <- cbind( kf$subjects, kf$value, kf$statistic, kf$p.value)
rownames(kappasOUT) <- c( "Cohen's kappa", "Light's kappa", "Fleiss's kappa")
colnames(kappasOUT) <- c( "N", "kappa", "z", "p" )
return (kappasOUT)
}


