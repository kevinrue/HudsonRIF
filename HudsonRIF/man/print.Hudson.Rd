\name{print.Hudson}
\alias{print.Hudson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print function for the Hudson class.
}
\description{
Print a summary of each component of the Hudson object.\cr\cr
}
\usage{
print.Hudson(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The output of the main wrapper \code{Hudson} function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Generic method \link{print}.
}
\examples{
# Simplified example using Affymetrix Bovine Genome Array probe sets.
output = Hudson(eSet, contrast = infected~control, classCol = "Treatment",
       DElist = c("Bt.12756.1.S1_at", "Bt.188.1.S1_at", "Bt.523.1.S1_at"), abs.PIF = TRUE)

# Overview of output
output

## The function is currently defined as
function (x, ...) 
{
    cat("$eSet\n")
    print(x$eSet)
    cat("\n")
    cat("$DElist\n")
    cat(paste(cat(x$DElist[1:5]), "\n...\n (Total:", length(x$DElist), 
        "values)\n\n"))
    cat("$conds\n")
    print(unlist(x$conds))
    cat(" (", length(x$conds), " conditions)\n\n", sep = "")
    cat("$EiAB\n")
    print(head(x$EiAB))
    cat("...\n (", paste(dim(x$EiAB), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$Ai\n")
    print(x$Ai[1:5])
    cat("...\n (", paste(length(x$Ai)), " numeric)\n\n", sep = "")
    cat("$dEi\n")
    print(x$dEi[1:5])
    cat("...\n (", paste(length(x$dEi)), " numeric)\n\n", sep = "")
    cat("$PIFi\n")
    print(x$PIFi[1:5])
    cat("...\n (", paste(length(x$PIFi)), " numeric)\n\n", sep = "")
    cat("$rAij\n")
    print(x$rAij[1:5, 1:5])
    cat("...\n (", paste(dim(x$rAij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$rBij\n")
    print(x$rBij[1:5, 1:5])
    cat("...\n (", paste(dim(x$rBij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$dCij\n")
    print(x$dCij[1:5, 1:5])
    cat("...\n (", paste(dim(x$dCij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$RIFi\n")
    print(x$RIFi[1:5])
    cat("...\n (", paste(length(x$RIFi)), " numeric)\n\n", sep = "")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Hudson }
\keyword{ RIF }