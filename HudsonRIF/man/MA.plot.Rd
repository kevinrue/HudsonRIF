\name{MA.plot}
\alias{MA.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hudson}{
The output of the main wrapper \code{Hudson} function.
}
  \item{symmetric}{
\code{logical}. Whether the Y-axis (differential expression) should have
a symmetrical range.
}
}
\details{

}
\value{
% Whatever plot() returns
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
Related methods \code{\link{Hudson}}, \code{\link{plot}}.
}
\examples{
# Simplified example using Affymetrix Bovine Genome Array probe sets.
output = Hudson(eSet, contrast = infected~control, classCol = "Treatment",
       DElist = c("Bt.12756.1.S1_at", "Bt.188.1.S1_at", "Bt.523.1.S1_at"), abs.PIF = TRUE)

# Plots the coexpression of probe set "Bt.6947.1.S1_at" against all DE genes in the two
# treatments.
MA.plot(hudson=output) 
 
## The function is currently defined as
function (hudson, formula, xlim = 0, ylim = 0) 
{
    if (symmetric) {
        plot(x = hudson$Ai, y = hudson$dEi, ylim = rep(max(abs(hudson$dEi)), 
            2) * c(-1, 1))
    }
    else {
        plot(x = hudson$Ai, y = hudson$dEi)
    }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Hudson }
\keyword{ RIF }