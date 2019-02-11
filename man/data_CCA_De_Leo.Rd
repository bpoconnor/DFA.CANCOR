\name{data_CCA_De_Leo}
\alias{data_CCA_De_Leo}
\docType{data}
\title{data_CCA_De_Leo}
\description{A data frame with scores on 14 variables that have the same correlational 
	structure, and which produce the same canonical correlation analysis results, 
	as those reported in De Leo and Wulfert (2013).}
\usage{data(data_CCA_De_Leo)}

\source{De Leo, J. A., & Wulfert, E. (2013). Problematic internet use and other risky 
	behaviors in college students: An application of problem-behavior theory. Psychology 
	of Addictive Behaviors, 27(1), 133-141. }

\examples{ 
\dontrun{
head(data_CCA_De_Leo)
     
# data that simulate those from De Leo & Wulfert (2013)
CANCOR(data = data_CCA_De_Leo,
       set1 = c("yrtob","yralc","yrdrugTOT","yrgambTOT","yrunsex","CiasTOT"),
       set2 = c("iveimpul","siasTOT","poms_dep","mdspss_TOT","aidTOT","fesmorals",
                "fesconflict","gpaspecified"),
       plot = 'yes', plotCV = 1)
}
}