\name{print.Hudson}

\alias{print.Hudson}

\title{
Print function for the Hudson class.
}

\description{
Print a summary of each component of the Hudson object.\cr\cr
}

\usage{
\method{print}{Hudson}(x, \dots)
}

\arguments{
  \item{x}{
The output of the main wrapper \code{Hudson} function.
}
  \item{\dots}{
Additional arguments to be passed on to the function.
}
}

\value{
A \code{list} of variables computed during the Regulatory Impact Factor analysis.\cr

  \item{eSet}{
  The input \code{ExpressionSet}.}
  \item{DElist}{The input list of differentially expressed genes.}
  \item{classCol}{The column name in \code{pData(eSet)} which contains the label 
of the conditions to compare)}
  \item{conds}{A \code{list} with elements \code{A} and \code{B} corresponding
to the labels of the conditions to compare.}
  \item{EiAB}{Computed matrix of average expression value for each gene feature (rows)
present in \code{rownames(exprs(eSet))} for each of the two conditions considered
(columns).}
  \item{Ai}{Array of average abundance of genes across the conditions to compare}
  \item{dEi}{Array of differential expression of genes between the conditions to
compare for all genes in \code{eSet}.}
  \item{PIFi}{Array of Phenotype Impact Factor for all genes in \code{eSet}.}
  \item{rAij}{Pairwise coexpression value between regulatory and phenotypic genes
in condition \code{A}.}
  \item{rBij}{Pairwise coexpression value between regulatory and phenotypic genes
in condition \code{B}.}
  \item{dCij}{Pairwise differential coexpression value between regulatory and
phenotypic genes between the two conditions.}
\item{RIFi}{Array of Regulatory Impact Factor values for genes in \code{DElist}.}
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
Generic method \link{print}.
}

\examples{
# Simulate a result of the Hudson function
hud = list(eSet = new("ExpressionSet"), DElist = c("DEgene1"), classCol = "Treatment",
  conds = list(A="infected", B="control"),
  EiAB = matrix(data=c(1,2,3,4), nrow=2, dimnames=list(c("gene1","gene2"), c("S1","S2"))), 
  Ai = c(1,2,3,4), dEi = c(-2,-1,1,2), PIFi = c(-2,-1,1,2),
  rAij = matrix(data=c(1,2,3,4), nrow=2, dimnames=list(c("gene1","gene2"), c("S1","S2"))),
  rBij = matrix(data=c(1,2,3,4), nrow=2, dimnames=list(c("gene1","gene2"), c("S1","S2"))),
  dCij = matrix(data=c(1,2,3,4), nrow=2, dimnames=list(c("gene1","gene2"), c("S1","S2"))),
  RIFi = c(1,2,3,4))
  
# The Hudson function changes the class of its output to "Hudson" for 
class(hud) = "Hudson"

# Use of the dedicated print function for class "Hudson"
hud

## The function is currently defined as
function (x, ...) 
{
    cat("$eSet\n")
    print(x$eSet)
    cat("\n")
    cat("$DElist\n")
    print(head(x$DElist))
    if (length(x$DElist) > 6) {
        cat("...\n")
    }
    cat(" (Total:", length(x$DElist), "values)\n\n")
    cat("$conds\n")
    print(unlist(x$conds))
    cat(" (", length(x$conds), " conditions)\n\n", sep = "")
    cat("$EiAB\n")
    print(head(x$EiAB))
    if (nrow(x$EiAB) > 6) {
        cat("...\n")
    }
    cat(" (", paste(dim(x$EiAB), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$Ai\n")
    print(head(x$Ai))
    if (length(x$Ai) > 6) {
        cat("...\n")
    }
    cat(" (", paste(length(x$Ai)), " numeric)\n\n", sep = "")
    cat("$dEi\n")
    print(head(x$dEi))
    if (length(x$dEi) > 6) {
        cat("...\n")
    }
    cat(" (", paste(length(x$dEi)), " numeric)\n\n", sep = "")
    cat("$PIFi\n")
    print(head(x$PIFi))
    if (length(x$PIFi) > 6) {
        cat("...\n")
    }
    cat(" (", paste(length(x$PIFi)), " numeric)\n\n", sep = "")
    cat("$rAij\n")
    print(x$rAij[1:min(5, nrow(x$rAij)), 1:min(5, ncol(x$rAij))])
    if (nrow(x$rAij) > 6 | ncol(x$rAij) > 6) {
        cat("...\n")
    }
    cat(" (", paste(dim(x$rAij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$rBij\n")
    print(x$rBij[1:min(5, nrow(x$rBij)), 1:min(5, ncol(x$rBij))])
    if (nrow(x$rBij) > 6 | ncol(x$rBij) > 6) {
        cat("...\n")
    }
    cat(" (", paste(dim(x$rBij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$dCij\n")
    if (nrow(x$dCij) > 6 | ncol(x$dCij) > 6) {
        cat("...\n")
    }
    print(x$dCij[1:min(5, nrow(x$dCij)), 1:min(5, ncol(x$dCij))])
    cat(" (", paste(dim(x$dCij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$RIFi\n")
    print(head(x$RIFi))
    if (length(x$Ai) > 6) {
        cat("...\n")
    }
    cat(" (", paste(length(x$RIFi)), " numeric)\n\n", sep = "")
}
}

\keyword{ Hudson }
\keyword{ RIF }