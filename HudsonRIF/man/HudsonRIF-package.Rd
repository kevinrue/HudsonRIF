\name{HudsonRIF-package}
\alias{HudsonRIF-package}
\alias{HudsonRIF}
\docType{package}
\title{
Calculation of Regulatory Impact Factor from gene expression data.
}
\description{
Based on the publications by Hudson et al. (see \code{References}), this
package provides tools to compare two phenotypes, identify genes features
with high Phenotype Impact Factor (PIF) between them, and gene features with
high Regulatory Impact Factor (RIF) between them.

On the one hand, high absolute PIF indicates high level of expression and
differential expression between the two conditions.
On the other hand, high RIF indicates high differential co-expression
between a gene feature and a list of user-defined genes (typically
differentially expressed genes). 
}
\details{
\tabular{ll}{
Package: \tab HudsonRIF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-24\cr
License: \tab GPL-2\cr
}
Important function: \link{Hudson}()\cr

The main output is given a list of all intermediate objects useful
for subsequent plots.
}
\author{
Kevin Rue-Albrecht

Maintainer: Kevin Rue-Albrecht <\href{mailto:kevin.rue@ucdconnect.ie}{email}>
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
\keyword{ package }
\seealso{
\code{\link[Biobase:Biobase-package]{Biobase-package}}, \code{\link[BiocGenerics:BiocGenerics-package]{BiocGenerics-package}}
}
\examples{
# Main wrapper
\code{output = Hudson(eSet, contrast=treated~control, classCol="Treatment",
                DElist=DE.probesets, abs.PIF = TRUE)}
# Summary of output list object
\code{output}
# Access to output variables
\code{output$RIFi}
}
