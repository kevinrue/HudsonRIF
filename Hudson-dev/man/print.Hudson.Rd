\name{print.Hudson}
\alias{print.Hudson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print function for the Hudson class.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.Hudson(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
Kevin Rue-Albrecht
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
function (x, ...) 
{
    cat("$eSet\n")
    print(x$eSet)
    cat("\n")
    cat("$DElist\n")
    cat(paste(cat(x$DElist[1:5]), "\n...\n (Total:", length(x$DElist), 
        "values)\n\n"))
    cat("$conds\n")
    print(unlist(t$conds))
    cat(" (", length(x$conds), " conditions)\n\n", sep = "")
    cat("$EiAB\n")
    print(head(x$EiAB))
    cat("...\n (", paste(dim(t$EiAB), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$Ai\n")
    print(x$Ai[1:5])
    cat("...\n (", paste(length(t$Ai)), " numeric)\n\n", sep = "")
    cat("$dEi\n")
    print(x$dEi[1:5])
    cat("...\n (", paste(length(t$dEi)), " numeric)\n\n", sep = "")
    cat("$PIFi\n")
    print(x$PIFi[1:5])
    cat("...\n (", paste(length(t$PIFi)), " numeric)\n\n", sep = "")
    cat("$rAij\n")
    print(t$rAij[1:5, 1:5])
    cat("...\n (", paste(dim(t$rAij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$rBij\n")
    print(t$rBij[1:5, 1:5])
    cat("...\n (", paste(dim(t$rBij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$dCij\n")
    print(t$dCij[1:5, 1:5])
    cat("...\n (", paste(dim(t$dCij), collapse = "x"), " matrix)\n\n", 
        sep = "")
    cat("$RIFi\n")
    print(x$RIFi[1:5])
    cat("...\n (", paste(length(t$RIFi)), " numeric)\n\n", sep = "")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
