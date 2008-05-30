\name{USecon}
\docType{data}
\alias{USecon}
\title{US macroeconomic time series}
\description{
  US macroeconomic data.
}
\usage{data(USecon)}
\format{
  The format is:
 mts [1:40, 1:2]  0.1364  0.0778 -0.3117 -0.5478 -1.2636 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "M1" "GNP"
 - attr(*, "tsp")= num [1:3] 1978 1988    4
 - attr(*, "class")= chr [1:2] "mts" "ts"
}
\details{
  The series are \code{100*diff(log())} of seasonally adjusted real 
  U.S. money 'M1' and GNP from 1978 to 1987.
}
\source{
  The complete data set is available in package \code{tseries}.
}

\examples{
data(USecon)
plot(USecon)
}
\keyword{datasets}

