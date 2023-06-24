\name{genetic}
\alias{genetic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Genetic Algorithm
}
\description{
Genetic algorithm that searches for the Base that optimizes the number of rules of a KB.
}
\usage{
genetic(KBo, np, m, q, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KBo}{Table structure of the KB to optimize.}
  \item{np}{(Integer) Number of individuals in the population.}
  \item{m}{(Double) Mutation chance for each new individual. Should be between 0 and 1.}
  \item{q}{(Double) Selection chance for a individual. Should be between 0 and 1.}
  \item{t}{(Integer) Number of generations to iterate through.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the KB that has been found that has the less number of rules, and therefore is the more space optimal.
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

optimalKB <- genetic(KB, 20, 0.05, 0.9, 100)
}
