\name{calculateEiA}
\alias{calculateEiA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the average expression of gene features in a condition.
}
\description{
\code{calculateEiA} is called twice by \code{calculateEiAB} to compute the 
average expression of each gene feature in each of the two conditions
compared.
}
\usage{
calculateEiA(eSet, classCol, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eSet}{
An object of class \code{ExpressionSet} containing the gene expression data.
}
  \item{classCol}{
The column name in \code{pData(eSet)} which contains the label 
of the condition to subset)
}
  \item{A}{
The label of the condition to subset samples from).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Array of average expression value for each gene feature in
\code{rownames(exprs(eSet))}
}
\references{
{
Hudson, N. J., B. P. Dalrymple, et al. (2012). "Beyond differential expression:
the quest for causal mutations and effector molecules." BMC Genomics 13: 356.
\href{http://www.ncbi.nlm.nih.gov/pubmed/22849396}{Pubmed}\cr
Hudson, N. J., A. Reverter, et al. (2009). "A differential wiring analysis of
expression data correctly identifies the gene containing the causal mutation." PLoS Comput Biol 5(5): e1000382.
\href{http://www.ncbi.nlm.nih.gov/pubmed/19412532}{Pubmed}
}
}
\author{
Kevin Rue-Albrecht
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Related methods \code{\link{calculateEiAB}()}, \code{\link{exprs}()}, \code{\link{pData}()} and class \code{\link{eSet}}. 
}
\examples{
calculateEiA(eSet, "Treatment", "control")

## The function is currently defined as
function (eSet, classCol, A) 
{
    return(apply(X = exprs(eSet[, which(pData(eSet)[, classCol] == 
        A)]), MARGIN = 1, FUN = "mean"))
  }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Hudson }
\keyword{ RIF }
