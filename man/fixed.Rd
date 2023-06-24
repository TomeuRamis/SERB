\name{fixed}
\alias{fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Identify the Fixed of a Rule
}
\description{
Identifies the Fixed part of a rule. The Fixed part explains the parameters that are most relevant to the decision chose. So this function returns the Most Relevant attributes of a set of evidence.
}
\usage{
fixed(kbm2l, evidence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kbm2l}{KBM2L structure to query.}
  \item{evidence}{Must be a vector containing evidence for all the atributes of a KB, including the Decision and Prob columns.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the list of the s that comprise the Fixed Part of a rule.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Create the associated KB
KB <- createKB(bn, "tub")
#Optimize KB
optimalKB <- genetic(KB, 20, 0.05, 0.9, 100)
#Create the KBM2L
KBM2L <- createKBM2L(optimalKB)

#Identify the fixed part of the first item of a KB
fix1 <- fixed(KBM2L, optimalKB[1,])

#Identify the fixed part of the second rule (this time using the KBM2L structure)
fix2 <- fixed(KBM2L, KBM2L[2,])
}
