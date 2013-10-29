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
Version: \tab 0.2\cr
Date: \tab 2013-10-24\cr
License: \tab GPL-2\cr
}
Wrapper function: \link{Hudson}\cr

Plotting functions: \link{diffCoexPIF.plot}

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
\code{\link[Biobase:Biobase-package]{Biobase-package}},
\code{\link[BiocGenerics:BiocGenerics-package]{BiocGenerics-package}}
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

# Use of the dedicated print function for class "Hudson"
hud

# plots the differential coexpression between all pairs of regulatory and DE genes
# WARNING: Struggles with massive data sets
diffCoex.bgplot(hudson=hud)

# plot the differential coexpression of a regulatory gene against all DE genes weighted by PIF value
diffCoexPIF.plot(hudson=hud, probe="gene 1", pch.max=15)

# plot the differential expression against the mean abundance of all genes
MA.plot(hudson=hud, symmetric=TRUE)
}
