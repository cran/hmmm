\name{getnames}
\alias{getnames}
\title{
vector of frequencies from a data frame
}
\description{
Function to extract the vector of frequencies associated with the category names from a data frame.
}
\usage{
getnames(dat, st = 3, sep = " ")
}
\arguments{
  \item{dat}{A contincency table in data frame format}
  \item{st}{Length of the string for every category name}
  \item{sep}{Separator of category names}
}
\value{The function returns a column vector of the frequencies of every combination of categories
of the involved variables in the data frame \code{dat} where each column corresponds to a variable,
each row to a combination of categories and the last column reports the frequencies. The variables
are arranged so that the farther to the left the 
column is the faster the category changes. Every frequency of each combination of categories is associated with a string of short
category names. The length of the names is determined by setting \code{st} and
consecutive names are separated by the symbol declared by \code{sep}.
}
\examples{
data(madsen)
y<-getnames(madsen,st=3,sep=";")}
\keyword{package}

