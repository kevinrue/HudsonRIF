\name{Hudson}
\alias{Hudson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main wrapper to compute Regulatory Impact Factor from gene expression set.
}
\description{
\code{Hudson} is the main wrapper to use for the Regulatory Impact Factor 
analysis on a gene expression dataset. It expects an ExpressionSet object
from the Biobase package, with phenotypic data including a column describing
the conditions to compare with at least two factor levels.\cr\cr
A list of suspected regulator genes can be specified, to restrict the RIF 
analysis to those.\cr\cr
We sugest the use of absolute(PIF) values for the calculation of RIF values
(see "Details" section of \link{calculateRIFi})
}
\usage{
Hudson(eSet, contrast, classCol, DElist, abs.PIF = TRUE, regulator.list = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eSet}{
An object of class \code{ExpressionSet} containing the gene expression data.
}
  \item{contrast}{
A \code{formula}} of the form \code{A~B} describing the two conditions to
compare. Positive differential expression means higher expression in 
condition \code{A} compared to condition \code{B}.
  \item{classCol}{
The column name in \code{pData(eSet)} which contains the label 
of the condition to subset)}
  \item{DElist}{
List of genes differentially expressed between the two conditions considered.
}
  \item{abs.PIF}{
\code{logical}. If \code{FALSE}, uses PIF values as published by Hudson et al.
See "Details" section of \code{\link{calculateRIFi}}.
}
  \item{regulator.list}{
List of suspected regulator genes to compute RIF values for.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{list} of variables computed during the RIF analysis. Class of the list
is set to "Hudson", associated with a \code{print.Hudson()} function.
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
Kevin Rue-Albrecht <\href{mailto:kevin.rue@ucdconnect.ie}{email}>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Related methods \code{\link{parameterCheck}}, \code{\link{calculateEiAB}},
 \code{\link{exprs}}, \code{\link{pData}},
\code{\link{calculateRIFi}}, \code{\link{cor}} and class \code{\link{ExpressionSet}},



}
\examples{
# Simplified example using Affymetrix Bovine Genome Array probe sets.
Hudson(eSet, contrast = infected~control, classCol = "Treatment",
       DElist = c("Bt.12756.1.S1_at", "Bt.188.1.S1_at", "Bt.523.1.S1_at"), abs.PIF = TRUE)

## The function is currently defined as
function (eSet, contrast, classCol, DElist, abs.PIF = TRUE, regulator.list = "") 
{
    parameterCheck(eSet, contrast, classCol, DElist, abs.PIF, 
        regulator.list)
    conds = list(A = as.character(contrast[[2]]), B = as.character(contrast[[3]]))
    if (regulator.list == "") {
        regulator.list = rownames(exprs(eSet))[!rownames(exprs(eSet)) \%in\% 
            DElist]
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
\keyword{ Hudson }
\keyword{ RIF }