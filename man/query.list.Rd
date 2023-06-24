\name{query.list}
\alias{query.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Query the KBM2L with partial evidence
}
\description{
Make a query to the KBM2L structure in order to obtain all possible instances of a partial (or complete) evidence. Given some partial evidence, this function return all possible alternatives with their associated decision and probability.
Using complete evidence is a particular case where the there is only one possible decision.
}
\usage{
query.list(KBM2L, evidence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KBM2L}{KBM2L structure to query.}
  \item{evidence}{ Vector of partial evidence to search in the KB structure.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of the possible instances of the evidence given with their associated Decision and Probability.
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

#Create some partial evidence
evidence <- setNames(c("yes", "yes", "no"), c("Tub", "Asia", "Smoke"))

query.list(KBM2L, evidence)
}
