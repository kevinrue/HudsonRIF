\name{diffExpr.plot}

\alias{diffExpr.plot}

\title{
Comparative background plot of coexpression in two conditions.
}

\description{
\code{diffExpr.plot} plots the expression of two given regulatory genes
in each sample in the two conditions. Each condition is represented in
a different color, with its corresponding linear regression line.
}

\usage{
diffExpr.plot(hudson, formula, xlim = 0, ylim = 0) 
}

\arguments{
  \item{hudson}{
The output of the main wrapper \code{Hudson} function.
}
  \item{formula}{
A \code{formula} of the form \code{A~B} describing the two genes to
compare. Gene \code{A} will be plotted in Y and gene \code{B} in X.
}
  \item{xlim}{
User-defined X-axis range. If not specified, the default is the range of values
+/- 0.5.
}
  \item{ylim}{
User-defined X-axis range. If not specified, the default is the range of values
+/- 0.5.
}
}

\details{
Treatment \code{A} will be colored in red and treatment \code{A} in black.
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
Related methods \code{\link{Hudson}}, \code{\link{plot}},  \code{\link{abline}}.
}

\examples{
# This function requires a valid ExpressionSet with attached phenodata.
 
## The function is currently defined as
function (hudson, formula, xlim = 0, ylim = 0) 
{
    probeX = as.character(formula[[3]])
    probeY = as.character(formula[[2]])
    if (xlim == 0) {
        xlim = range(exprs(hudson$eSet)[probeX, which(pData(hudson$eSet)[, 
            hudson$classCol] \%in\% unlist(hudson$conds))]) + c(-0.5, 
            0.5)
    }
    if (ylim == 0) {
        ylim = range(exprs(hudson$eSet)[probeY, which(pData(hudson$eSet)[, 
            hudson$classCol] \%in\% unlist(hudson$conds))]) + c(-0.5, 
            0.5)
    }
    plot(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$B)]))), 
        main = paste(hudson$conds$B, ": Black\n", hudson$conds$A, 
            ": Red"), xlim = xlim, ylim = ylim, pch = 16)
    abline(lm(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$B)])))))
    points(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$A)]))), 
        col = "red", pch = 17)
    abline(lm(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$A)])))), 
        col = "red")
}
}

\keyword{ Hudson }
\keyword{ RIF }