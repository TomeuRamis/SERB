\name{tabu}
\alias{tabu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Taboo Algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Taboo search algorithm that searches for the Base that optimizes the number of rules of a KB.
}
\usage{
tabu(KB, m, maxtabu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KB}{ Table structure of the KB to optimize.}
  \item{m}{ (Integer) Number of iterations to search for.}
  \item{maxtabu}{ (Integer) Max length of the Taboo list. The longer the list is made, the more memory it has to avoid previous solutions.
  This parameter has a maximum length defined by the number of attributes of a KB. If the taboo list is too long, at some point all neighbors of a base will be considered taboo and the program will throw out an error. }
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

optimalKB <- tabu(KB, 100, 10)
}
