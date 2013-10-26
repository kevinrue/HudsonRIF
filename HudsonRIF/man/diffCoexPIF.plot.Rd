\name{diffCoexPIF.plot}

\alias{diffCoexPIF.plot}

\title{
Comparative plot of coexpression in two conditions.
}

\description{
\code{diffCoexPIF.plot} plots the coexpression of the specified gene 
with each of the DE.

}
\usage{
diffCoexPIF.plot(hudson, probe, pch.max = 0)
}

\arguments{
  \item{hudson}{
The output of the main wrapper \code{Hudson} function.
}
  \item{probe}{
The regulatory probe set to plot coexpression values of.
}
  \item{pch.max}{
A scaling factor defining the size of the data point associated with the
maximal PIF value.
}
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

\seealso{
Related methods \code{\link{Hudson}}, \code{\link{plot}},  \code{\link{abline}}.
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

# plot the differential coexpression of a regulatory gene against all DE genes weighted by PIF value
diffCoexPIF.plot(hudson=hud, probe="gene 1", pch.max=15)
 
## The function is currently defined as
function (hudson, probe, pch.max = 0) 
{
    if (pch.max == 0) {
        pch.max = max(hudson$PIFi)
    }
    PIF.DE = hudson$PIFi[hudson$DElist]/max(hudson$PIFi) * pch.max
    plot(x = hudson$rBij[probe, ], y = hudson$rAij[probe, ], 
        xlim = c(-1, 1), ylim = c(-1, 1), cex = PIF.DE, xlab = hudson$conds$B, 
        ylab = hudson$conds$A, pch = 20, main = probe)
    abline(coef = c(0, 1), col = "green", lwd = 2)
}
}

\keyword{ Hudson }
\keyword{ RIF }