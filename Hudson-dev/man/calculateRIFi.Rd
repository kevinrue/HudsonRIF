\name{calculateRIFi}
\alias{calculateRIFi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the Regulatory Impact Factor of gene features from PIF and dCij.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateRIFi(PIFi, dCij, DElist, abs.PIF = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PIFi}{
%%     ~~Describe \code{PIFi} here~~
}
  \item{dCij}{
%%     ~~Describe \code{dCij} here~~
}
  \item{DElist}{
%%     ~~Describe \code{DElist} here~~
}
  \item{abs.PIF}{
%%     ~~Describe \code{abs.PIF} here~~
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
function (PIFi, dCij, DElist, abs.PIF = T) 
{
    DE.PIFi = PIFi[DElist]
    dCij.coefs = dCij^2
    if (abs.PIF == T) {
        return(apply(X = ((dCij.coefs \%*\% abs(DE.PIFi))/length(DElist)), 
            MARGIN = 1, FUN = sum))
    }
    else {
        return(apply(X = ((dCij.coefs \%*\% DE.PIFi)/length(DElist)), 
            MARGIN = 1, FUN = sum))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
