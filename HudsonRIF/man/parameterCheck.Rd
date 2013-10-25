\name{parameterCheck}
\alias{parameterCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks the user input for invalid values.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parameterCheck(eSet, contrast, classCol, DElist, abs.PIF, regulator.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eSet}{
An object of class \code{ExpressionSet} containing the gene expression data.
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
}
  \item{classCol}{
The column name in \code{pData(eSet)} which contains the label 
of the condition to subset)}
  \item{DElist}{
%%     ~~Describe \code{DElist} here~~
}
  \item{abs.PIF}{
%%     ~~Describe \code{abs.PIF} here~~
}
  \item{regulator.list}{
%%     ~~Describe \code{regulator.list} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (eSet, contrast, classCol, DElist, abs.PIF, regulator.list) 
{
    if (class(eSet) != "ExpressionSet") {
        stop("\"eSet\" (", class(eSet), ") is not an ExpressionSet", 
            call. = F)
    }
    if (class(contrast) != "formula") {
        stop("\"contrast\" (", class(contrast), ") is not a formula", 
            call. = F)
    }
    if (length(contrast[[2]]) > 1 | length(contrast[[3]]) > 1) {
        stop("\"contrast\" (", paste(as.character(contrast[c(2, 
            3)]), collapse = " ~ "), ") contains more than one term on each side", 
            call. = F)
    }
    if (!classCol \%in\% names(pData(eSet))) {
        stop("\"classCol\" (", classCol, ") is not a valid column name in \"pData(eSet)\"", 
            call. = F)
    }
    conds = list(A = as.character(contrast[[2]]), B = as.character(contrast[[3]]))
    if (!conds$A \%in\% levels(as.factor(pData(eSet)[, classCol]))) {
        stop("(", conds$A, ") is not a valid class level in \"pData(eSet)[,", 
            classCol, "]\"", call. = F)
    }
    if (!conds$B \%in\% levels(as.factor(pData(eSet)[, classCol]))) {
        stop("(", conds$B, ") is not a valid class level in \"pData(eSet)[,", 
            classCol, "]\"", call. = F)
    }
    if (sum(DElist \%in\% rownames(exprs(eSet))) != length(DElist)) {
        stop("\"DElist contains ", sum(!DElist \%in\% rownames(exprs(eSet))), 
            " feature(s) absent from \"rownames(exprs(eSet))\"", 
            call. = F)
    }
    if (!is.logical(abs.PIF)) {
        stop("\"abs.PIF\" (", class(abs.PIF), ") is not a logical", 
            call. = F)
    }
    if (any(regulator.list != "")) {
        if (sum(!regulator.list \%in\% rownames(exprs(eSet))) > 
            0) {
            stop("\"regulator.list contains ", sum(!regulator.list \%in\% 
                rownames(exprs(eSet))), " feature(s) absent from \"rownames(exprs(eSet))\"", 
                call. = F)
        }
        if (sum(DElist \%in\% regulator.list) > 0) {
            stop("\"DElist and regulator.list overlap by (", 
                sum(DElist \%in\% regulator.list), " feature(s) ", 
                call. = F)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Hudson }
\keyword{ RIF }