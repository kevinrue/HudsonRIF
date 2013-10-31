\name{Hudson}

\alias{Hudson}

\title{
Main wrapper to compute Regulatory Impact Factor from gene expression set.
}

\description{
'\code{Hudson}' is the main wrapper to use for the Regulatory Impact Factor 
analysis for a gene expression dataset. It takes as input an ExpressionSet object
from the Biobase package, with phenotypic data including a column describing
the conditions to compare (two factor levels).\cr\cr
A list of suspected regulator genes can be specified (such as a list of transcription
factors), otherwise all genes in the expressionset that are not listed in the list of
differentially expressed genes are used.\cr\cr
We sugest the use of absolute PIF values for the calculation of RIF values
(see "Details" section of \link{calculateRIFi})
}

\usage{
Hudson(eSet, contrast, classCol, DElist, abs.PIF = TRUE, regulator.list = "")
}

\arguments{
  \item{eSet}{
An object of class \code{ExpressionSet} containing the gene expression data.
}
  \item{contrast}{
A \code{formula} of the form \code{A~B} describing the two conditions to
compare. Positive differential expression means higher expression in 
condition \code{A} compared to condition \code{B}.}
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

\value{
A \code{list} of variables computed during the RIF analysis. Class of the list
is set to "Hudson", associated with a \code{print.Hudson()} function.
}

\references{
\itemize{
\item{Hudson, N. J., B. P. Dalrymple, et al. (2012). "Beyond differential expression:
the quest for causal mutations and effector molecules." BMC Genomics 13: 356.
PMID:\href{http://www.ncbi.nlm.nih.gov/pubmed/22849396}{22849396}}
\item{Reverter, A., N. J. Hudson, et al. (2010). "Regulatory impact factors: unraveling the transcriptional regulation of complex traits from expression data." Bioinformatics 26(7): 896-904.
PMID:\href{http://www.ncbi.nlm.nih.gov/pubmed/20144946}{20144946}}
\item{Hudson, N. J., A. Reverter, et al. (2009). "A differential wiring analysis of
expression data correctly identifies the gene containing the causal mutation."
PLoS Comput Biol 5(5): e1000382.
PMID:\href{http://www.ncbi.nlm.nih.gov/pubmed/19412532}{19412532}}
}
}

\author{
Kevin Rue-Albrecht <\href{mailto:kevin.rue@ucdconnect.ie}{kevin DOT rue AT ucdconnect DOT ie}>
and Paul McGettigan.
}

\seealso{
Related methods \code{\link{paramCheck.Hudson}}, \code{\link{calculateEiAB}},
 \code{\link{exprs}}, \code{\link{pData}},
\code{\link{calculateRIFi}}, \code{\link{cor}}, class \code{\link{ExpressionSet}},
and dedicated print statement \code{\link{print.Hudson}}.
}

\examples{
# Simulate a result of the Hudson function
hud = list(
  eSet = new("ExpressionSet"),
  DElist = paste("gene", 6:9),
  classCol = "Treatment",
  conds = list(A="infected", B="control"),
  EiAB = matrix(data=rnorm(n=10, mean=10, sd=4), ncol=2, dimnames=list(paste("gene", 1:5), c("S1","S2"))), 
  Ai = rnorm(n=5, mean=10, sd=4),
  dEi = rnorm(n=5, mean=0, sd=1),
  PIFi = rnorm(n=4, mean=15, sd=10),
  rAij = matrix(data = runif(n = 20, min = -1, max = 1), ncol = 4, dimnames = list(paste("gene", 1:5), paste("gene", 6:9))),
  rBij = matrix(data = runif(n = 20, min = -1, max = 1), ncol = 4, dimnames = list(paste("gene", 1:5), paste("gene", 6:9))),
  dCij = matrix(data = runif(n = 20, min = -2, max = 2), ncol = 4, dimnames = list(paste("gene", 1:5), paste("gene", 6:9))),
  RIFi = rnorm(n=5, mean=2, sd=1))
# Named vectors are required for most plots
names(hud$PIFi) = paste("gene", 6:9)
# The Hudson function changes the class of its output to "Hudson" for 
class(hud) = "Hudson"

## The function is currently defined as
function (eSet, contrast, classCol, DElist, abs.PIF = TRUE, regulator.list = "")
{
    paramCheck.Hudson(eSet, contrast, classCol, DElist, abs.PIF, 
        regulator.list)
    conds = list(A = as.character(contrast[[2]]), B = as.character(contrast[[3]]))
    if (any(regulator.list == "")) {
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
    rAij = rAij[regulator.list, DElist]
    rBij = cor(x = t(exprs(eSet[, which(pData(eSet)[, classCol] == 
        conds$B)])), method = "spearman")
    rBij = rBij[regulator.list, DElist]
    dCij = rAij - rBij
    RIFi = calculateRIFi(PIFi = PIFi, dCij = dCij, DElist = DElist, 
        abs.PIF = abs.PIF)
    output = list(eSet = eSet, DElist = DElist, classCol = classCol, 
        conds = conds, regulator.list = regulator.list, EiAB = EiAB, 
        Ai = Ai, dEi = dEi, PIFi = PIFi, rAij = rAij, rBij = rBij, 
        dCij = dCij, RIFi = RIFi)
    class(output) = "Hudson"
    return(output)
}
}

\keyword{ Hudson }
\keyword{ RIF }