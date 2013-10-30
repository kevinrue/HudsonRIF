\name{MA.plot}

\alias{MA.plot}

\title{
Differential expression of all genes plotted against the average abundance.
}

\description{
\code{MA.plot} plots the differential expression of each gene in the 
expression set between the two conditions against its average abundance
across those two conditions.
}

\usage{
MA.plot(hudson, symmetric = TRUE) 
}

\arguments{
  \item{hudson}{
The output of the main wrapper \code{Hudson} function.
}
  \item{symmetric}{
\code{logical}. Whether the Y-axis (differential expression) should have
a symmetrical range.
}
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
Related methods \code{\link{Hudson}}, \code{\link{plot}}.
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

# plot the differential expression against the mean abundance of all genes
MA.plot(hudson=hud, symmetric=TRUE)
 
## The function is currently defined as
function (hudson, symmetric = TRUE) 
{
    if (symmetric) {
        plot(x = hudson$Ai, y = hudson$dEi, ylim = rep(max(abs(hudson$dEi)), 
            2) * c(-1, 1), xlab = "A", ylab = "M")
    }
    else {
        plot(x = hudson$Ai, y = hudson$dEi)
    }
    abline(h = 0, col = "red")
}
}

\keyword{ Hudson }
\keyword{ RIF }