\name{Hudson}
\alias{Hudson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main wrapper to compute Regulatory Impact Factor from gene expression set.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Hudson(eSet, contrast, classCol, DElist, abs.PIF = TRUE, regulator.list = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eSet}{
%%     ~~Describe \code{eSet} here~~
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
}
  \item{classCol}{
%%     ~~Describe \code{classCol} here~~
}
  \item{DElist}{
%%     ~~Describe \code{DElist} here~~
}
  \item{abs.PIF}{
%%     ~~Describe \code{abs.PIF} here~~
}
  \item{regulator.list}{
%%     ~~Describe \code{regulator.list} here~~
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


## The function is currently defined as
function (eSet, contrast, classCol, DElist, abs.PIF = TRUE, regulator.list = "") 
{
    parameterCheck(eSet, contrast, classCol, DElist, abs.PIF, 
        regulator.list)
    conds = list(A = as.character(contrast[[2]]), B = as.character(contrast[[3]]))
    if (regulator.list == "") {
        regulator.list = rownames(exprs(farms_informative.65))[!rownames(exprs(farms_informative.65)) \%in\% 
            probes]
    }
    EiAB = calculateEiAB(eSet = eSet, classCol = classCol, conds$A, 
        conds$B)
    Ai = (EiAB[, conds$A] + EiAB[, conds$B])/2
    dEi = EiAB[, conds$A] - EiAB[, conds$B]
    PIFi = Ai * dEi
    rAij = cor(x = t(exprs(eSet[, which(pData(eSet)[, classCol] == 
        conds$A)])), method = "spearman")
    rAij = rAij[!rownames(rAij) \%in\% DElist, DElist]
    rBij = cor(x = t(exprs(eSet[, which(pData(eSet)[, classCol] == 
        conds$B)])), method = "spearman")
    rBij = rBij[!rownames(rBij) \%in\% DElist, DElist]
    dCij = rAij - rBij
    RIFi = calculateRIFi(PIFi = PIFi, dCij = dCij, DElist = DElist, 
        abs.PIF = abs.PIF)
    output = list(eSet = eSet, DElist = DElist, classCol = classCol, 
        conds = conds, EiAB = EiAB, Ai = Ai, dEi = dEi, PIFi = PIFi, 
        rAij = rAij, rBij = rBij, dCij = dCij, RIFi = RIFi)
    class(output) = "Hudson"
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
