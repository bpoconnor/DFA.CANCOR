

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