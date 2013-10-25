\name{calculateRIFi}
\alias{calculateRIFi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the Regulatory Impact Factor of gene features from PIFi and dCij.
}
\description{
\code{calculateRIFi} calls uses the Phenotypic Impact Factor (PIF)
differential coexpression (dC) to compute the Regulatory Impact
Factor (RIF).\cr
A variation from the original Hudson published formula allows to use the 
absolute PIF value to compute the RIF.
}
\usage{
calculateRIFi(PIFi, dCij, DElist, abs.PIF = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PIFi}{
The Phenotype Impact Factor of all gene features in \code{DElist}.
}
  \item{dCij}{
The differential coexpression value for each pair of gene features in
\code{DElist}
and \code{regulator.list}.
}
  \item{DElist}{
List of genes differentially expressed between the two conditions considered.
}
  \item{abs.PIF}{
\code{logical}. If \code{FALSE}, uses PIF values as published by Hudson et al.
See section "Details".
}
}
\details{
\code{abs.PIF} allows to use absolute PIF values for genes in \code{DElist}.\cr
RIF values are the  sum of squared differential coexpression
(positive) weighted by PIF values (signed).
The use of signed PIF values may lead large positive and negative PIF genes
to cancel each other out, even in the situation of a regulator gene differentially
coexpressed with a large positive PIF gene and a large negative PIF gene.\cr
PIF values describe the abundance of genes across two conditions and their 
differential expression between those two conditions. Therefore, large
positive values indicate abundant genes with higher expression in condition A
compared to B, and vice versa.



}
\value{
Array of RIF values for each gene feature present in
\code{regulator.list} if specified; otherwise all gene features in
\code{rownames(exprs(eSet))} except \code{DElist}
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
Related methods \code{\link{Hudson}}, \code{\link{calculateEiAB}}. 
}
\examples{
# Adapt to your variables
classCol ="Treatment"
conds = list(A = "infected" , B = "control")

# Calculate the average expression for each condition
EiAB = calculateEiAB(eSet = eSet, classCol = classCol, conds$A, 
        conds$B)
# Calculate the average abundance across conditions
Ai = (EiAB[, conds$A] + EiAB[, conds$B])/2
# Calculate the average differential expression between conditions
dEi = EiAB[, conds$A] - EiAB[, conds$B]
# Calculate the PIF values
PIFi = Ai * dEi
# Calculate pairwise coexpression values in condition A
rAij = cor(x = t(exprs(eSet[, which(pData(eSet)[, classCol] == 
           conds$A)])), method = "spearman")
# Subset to the coexpression values to pairs between regulators and DE genes
rAij = rAij[!rownames(rAij) %in% DElist, DElist]
# Calculate pairwise coexpression values in condition B
rBij = cor(x = t(exprs(eSet[, which(pData(eSet)[, classCol] == 
           conds$B)])), method = "spearman")
# Subset to the coexpression values to pairs between regulators and DE genes           
rBij = rBij[!rownames(rBij) %in% DElist, DElist]
# Calculate the differential coexpression
dCij = rAij - rBij
# Finally: Calculate the RIF values
RIFi = calculateRIFi(PIFi = PIFi, dCij = dCij, DElist = DElist, 
                     abs.PIF = abs.PIF)
                     
## The function is currently defined as
function (PIFi, dCij, DElist, abs.PIF = T) 
{
    DE.PIFi = PIFi[DElist]
    dCij.coefs = dCij^2
    if (abs.PIF == T) {
        return(apply(X = ((dCij.coefs \%*\% abs(DE.PIFi))/length(DElist)), 
            MARGIN = 1, FUN = sum))
    }
    else {
        return(apply(X = ((dCij.coefs \%*\% DE.PIFi)/length(DElist)), 
            MARGIN = 1, FUN = sum))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Hudson }
\keyword{ RIF }