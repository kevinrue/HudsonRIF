\name{calculateEiAB}
\alias{calculateEiAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the average expression of gene features in each of two conditions.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateEiAB(eSet, classCol, A, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eSet}{
%%     ~~Describe \code{eSet} here~~
}
  \item{classCol}{
%%     ~~Describe \code{classCol} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kevin Rue-Albrecht
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (eSet, classCol, A, B) 
{
    EiAB = calculateEiA(eSet, classCol, A)
    EiAB = cbind(EiAB, calculateEiA(eSet, classCol, B))
    colnames(EiAB) = c(A, B)
    return(EiAB)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
