\name{NelPlo}
\docType{data}
\alias{NelPlo}
\title{Nelson-Plosser macroeconomic time series}
\description{
  A subset of Nelson-Plosser data.
}
\usage{data(NelPlo)}
\format{
  The format is:
 mts [1:43, 1:2] -4.39  3.12  1.08 -1.50  3.91 ...
 - attr(*, "tsp")= num [1:3] 1946 1988    1
 - attr(*, "class")= chr [1:2] "mts" "ts"
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "ip" "stock.prices"
}
\details{
  The series are \code{100*diff(log())} of 
  industrial production and stock prices (S&P500) from 1946 to 1988.
}
\source{
  The complete data set is available in package \code{tseries}.
}
\references{}

\examples{
data(NelPlo)
plot(NelPlo)
}
\keyword{datasets}

