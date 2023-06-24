\name{sima}
\alias{sima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulated Annealing algorithm
}
\description{
Simulated Annealing search algorithm that searches for the Base that optimizes the number of rules of a KB.
}
\usage{
sima(KB, T0, k, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KB}{Table structure of the KB to optimize.}
  \item{T0}{(Integer) Initial temperature of the system.}
  \item{k}{(double) Random chance of choosing a non-optimal partial solution. Should be between 0 and 1.}
  \item{c}{(double) Rate of cooling of the system. Should be between 0 and 1.}
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

optimalKB <- sima(KB, 1000, 0.1, 0.9)
}
