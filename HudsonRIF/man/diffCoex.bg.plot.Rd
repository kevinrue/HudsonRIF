\name{diffCoex.bg.plot}
\alias{diffCoex.bg.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparative background plot of coexpression in two conditions.
}
\description{
\code{diffCoexPIF.plot} plots the coexpression of all regulatory genes
with each of the DE.
}
\usage{
diffCoexPIF.bg.plot(hudson)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hudson}{
The output of the main wrapper \code{Hudson} function.
}
}
\details{
This function struggles with very large expression datasets.
}
\value{
% Whatever abline() returns
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
Related methods \code{\link{Hudson}}, \code{\link{plot}},  \code{\link{abline}}.
}
\examples{
# Simplified example using Affymetrix Bovine Genome Array probe sets.
output = Hudson(eSet, contrast = infected~control, classCol = "Treatment",
       DElist = c("Bt.12756.1.S1_at", "Bt.188.1.S1_at", "Bt.523.1.S1_at"), abs.PIF = TRUE)

# Plots the coexpression of probe set "Bt.6947.1.S1_at" against all DE genes in the two
# treatments.
diffCoex.bg.plot(hudson=output) 
 
## The function is currently defined as
function (hudson) 
{
    plot(x = hudson$rBij, y = hudson$rAij, xlim = c(-1, 1), ylim = c(-1, 
        1), xlab = hudson$conds$B, ylab = hudson$conds$A, pch = 20, 
        col = "grey")
    abline(coef = c(0, 1), col = "green", lwd = 2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Hudson }
\keyword{ RIF }