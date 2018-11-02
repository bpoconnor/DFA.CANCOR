\name{DFA}
\alias{DFA}
\title{Discriminant function analysis}
\description{This function provides SPSS- and SAS-like output for linear discriminant 
	function analysis. It uses functions from the car, CPP, MASS, MVN, and irr packages.}
\usage{DFA(data, groups, variables, normtests, priorprob, predictive)}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables.}

  \item{groups}{
  \code{}The name of the groups variable in the dataframe, \cr e.g., groups = 'Group'.}
  	    
  \item{variables}{ 
  \code{}The names of the continuous variables in the dataframe that will be used in the DFA, 
         e.g., variables = c('varA', 'varB', 'varC').}
  	    
  \item{normtests}{ 
  \code{}Should tests of univariate and multivariate normality be conducted? \cr The options are: 'yes' or 'no'.}
  	    
  \item{priorprob}{
  \code{}How should the prior probabilities of the group sizes be computed? \cr The options are:\cr 
  	    'EQUAL' for equal group sizes; or\cr 
  	    'SIZES' for the group sizes to be based on the sizes of the groups in the dataframe.}
                  
  \item{predictive}{
  \code{}Should a predictive DFA be conducted? \cr The options are: 'yes' or 'no'.}
  
%  \item{CV}{
%  \code{}Should cross-validation (leave-one-out cross-validation) be conducted for the predictive DFA? The options are: 'YES' or 'NO'.}

}
\value{The displayed output includes descriptive statistics for the groups, 
	   tests of univariate and multivariate normality, 
	   the results of tests of the homogeneity of the group variance-covariance matrices, 
	   eigenvalues & canonical correlations, Wilks lambda & peel-down statistics,
	   one-way ANOVA tests of group differences in scores on each discriminant function,
	   raw and standardized discriminant function coefficients, 
	   structure coefficients, functions at group centroids, 
	   a plot of the group means on the standardized discriminant functions,
	   and extensive output from predictive discriminant function analyses (if requested).
	   		   	   
	   The returned output is a list with elements
  \item{rawCoef}{canonical discriminant function coefficients}
  \item{structCoef}{structure coefficients}
  \item{standCoef}{standardized coefficients}
  \item{standCoefSPSS}{standardized coefficients from SPSS}
  \item{centroids}{unstandardized canonical discriminant functions evaluated at the group means}
  \item{centroidSDs}{group standard deviations on the unstandardized functions}
  \item{centroidsZ}{standardized canonical discriminant functions evaluated at the group means}
  \item{centroidSDsZ}{group standard deviations on the standardized functions}
}
\references{
       {Sherry, A. (2006). Discriminant analysis in counseling research. Counseling Psychologist, 34, 661-683.}
\cr\cr {Tabachnik, B. G., & Fidell, L. S. (2013). Using multivariate statistics (6th ed.). New York, NY: Pearson.}
}
\author{Brian P. O'Connor }
\examples{
DFA(data = data_DFA_Field, 
    groups = 'Group', 
    variables = c('Actions','Thoughts'),
    normtests = 'yes', priorprob = 'SIZES', predictive = 'yes')

\dontrun{
DFA(data = data_DFA_Sherry, 
    groups = 'Group',
    variables = c('Neuroticism','Extroversion','Openness',
                  'Agreeableness','Conscientiousness'),
    normtests = 'yes', priorprob = 'SIZES', predictive = 'yes')
}
}
