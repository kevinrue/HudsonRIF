\name{diffExpr.plot}
\alias{diffExpr.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparative background plot of coexpression in two conditions.
}
\description{
\code{diffExpr.plot} plots the expression of two given regulatory genes
in each sample in the two conditions. Each condition is represented in
a different color, with its corresponding linear regression line.
}
\usage{
diffExpr.plot(hudson, formula, xlim = 0, ylim = 0) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hudson}{
The output of the main wrapper \code{Hudson} function.
}
  \item{formula}{
A \code{formula} of the form \code{A~B} describing the two genes to
compare. Gene \code{A} will be plotted in Y and gene \code{B} in X.
}
  \item{xlim}{
User-defined X-axis range. If not specified, the default is the range of values
+/- 0.5.
}
  \item{ylim}{
User-defined X-axis range. If not specified, the default is the range of values
+/- 0.5.
}
}
\details{
Treatment \code{A} will be colored in red and treatment \code{A} in black.
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
diffExpr.plot(hudson=output, formula=Bt.12756.1.S1_at~Bt.278.1.S1_at) 
 
## The function is currently defined as
function (hudson, formula, xlim = 0, ylim = 0) 
{
    probeX = as.character(formula[[3]])
    probeY = as.character(formula[[2]])
    if (xlim == 0) {
        xlim = range(exprs(hudson$eSet)[probeX, which(pData(hudson$eSet)[, 
            hudson$classCol] \%in\% unlist(output$conds))]) + c(-0.5, 
            0.5)
    }
    if (ylim == 0) {
        ylim = range(exprs(hudson$eSet)[probeY, which(pData(hudson$eSet)[, 
            hudson$classCol] \%in\% unlist(output$conds))]) + c(-0.5, 
            0.5)
    }
    plot(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$B)]))), 
        main = paste(hudson$conds$B, ": Black\n", hudson$conds$A, 
            ": Red"), xlim = xlim, ylim = ylim, pch = 16)
    abline(lm(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$B)])))))
    points(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$A)]))), 
        col = "red", pch = 17)
    abline(lm(formula = formula, data = as.data.frame(t(exprs(hudson$eSet[, 
        which(pData(hudson$eSet)[, hudson$classCol] == hudson$conds$A)])))), 
        col = "red")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Hudson }
\keyword{ RIF }