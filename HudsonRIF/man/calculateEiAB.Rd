\name{calculateEiAB}

\alias{calculateEiAB}

\title{
Computes the average expression of gene features in each of two conditions.
}

\description{
\code{calculateEiAB} calls twice \code{calculateEiA} to compute the 
average expression of each gene feature in each of the two conditions
compared.
}

\usage{
calculateEiAB(eSet, classCol, A, B)
}

\arguments{
  \item{eSet}{
An object of class \code{ExpressionSet} containing the gene expression data.
}
  \item{classCol}{
The column name in \code{pData(eSet)} which contains the label 
of the condition to subset.
}
  \item{A}{
The label of the first condition to calculate the average expression for.
}
  \item{B}{
The label of the second condition to calculate the average expression for.
}
}

\value{
Matrix of average expression value for each gene feature (rows) present in
\code{rownames(exprs(eSet))} for each of the two conditions considered (columns).
}

\references{
\itemize{
\item{Hudson, N. J., B. P. Dalrymple, et al. (2012). "Beyond differential expression:
the quest for causal mutations and effector molecules." BMC Genomics 13: 356.
{\href{http://www.ncbi.nlm.nih.gov/pubmed/22849396}{Pubmed}}}
\item{Hudson, N. J., A. Reverter, et al. (2009). "A differential wiring analysis of
expression data correctly identifies the gene containing the causal mutation."
PLoS Comput Biol 5(5): e1000382.
\href{http://www.ncbi.nlm.nih.gov/pubmed/19412532}{Pubmed}}
}
}

\author{
Kevin Rue-Albrecht <\href{mailto:kevin.rue@ucdconnect.ie}{kevin DOT rue AT ucdconnect DOT ie}>
and Paul McGettigan.
}

\seealso{
Related methods \code{\link{calculateEiA}}, \code{\link{exprs}}, \code{\link{pData}} and class \code{\link{eSet}}. 
}

\examples{

## The function is currently defined as
function (eSet, classCol, A, B) 
{
    EiAB = calculateEiA(eSet, classCol, A)
    EiAB = cbind(EiAB, calculateEiA(eSet, classCol, B))
    colnames(EiAB) = c(A, B)
    return(EiAB)
}
}

\keyword{ Hudson }
\keyword{ RIF }