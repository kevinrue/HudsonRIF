\name{paramCheck.diffCoexPIF.plot}

\alias{paramCheck.diffCoexPIF.plot}

\title{
Checks the user input for invalid values.
}

\description{
Returns an error if one on the user-defined inputs does not match the 
requirements of the package.
}

\usage{
paramCheck.diffCoexPIF.plot(hudson, probe, pch.max)
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

\examples{

## The function is currently defined as
function (hudson, probe, pch.max) 
{
    if (class(hudson) != "Hudson") {
        warning("\"hudson\" (", class(hudson), ") is not an Hudson object.", 
            call. = F)
    }
    if (any(!c("conds", "DElist", "PIFi", "rAij", "rBij", "RIFi") \%in\% 
        names(hudson))) {
        stop("\"hudson\" (", names(hudson), ") is not a valid Hudson object.\n         One item of $conds, $DElist, $PIFi, $rAij, $rBij or $RIFi is missing.", 
            call. = F)
    }
    if (!is.character(probe)) {
        stop("\"probe\" (", probe, ") is not a character object.", 
            call. = F)
    }
    if (!probe \%in\% hudson$rAij[probe, ]) {
        stop("\"probe\" (", probe, ") is not a valid probe ID in rAij rows.", 
            call. = F)
    }
    if (!probe \%in\% hudson$rBij[probe, ]) {
        stop("\"probe\" (", probe, ") is not a valid probe ID in rBij rows.", 
            call. = F)
    }
    if (!is.numeric(pch.max)) {
        stop("\"pch.max\" (", pch.max, ") is not a numeric object.", 
            call. = F)
    }
    if (pch.max < 0) {
        stop("\"pch.max\" (", pch.max, ") should be a positive integer.", 
            call. = F)
    }
}
}


\keyword{ Hudson }
\keyword{ RIF }