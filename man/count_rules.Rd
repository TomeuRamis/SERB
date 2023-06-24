\name{count_rules}
\alias{count_rules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Count Rules of Knowledge Base
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Count the number of independent "Rules" that define the Knowledge Base.
}
\usage{
count_rules(KB, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KB}{Table structure of a Knowledge Base. The two last columns must be Decision and Prob.}
  \item{verbose}{Boolean value. If TRUE prints out the ruleset in indices notation, if FALSE does not print any information.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A Rule of a KB is defined as the last item of a adjacent set of evidences with the same Decision. The fewer rules a KB has, the more optimal its space optimization is.
}
\value{
Returns an the number of rules that define the KB.
}
\references{
Fernández del Pozo, J.A., Bielza, C. and Gómez, M. (2003) \emph{A list-based compact representation for large decision tables management}, European Journal of Operational Research.
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
count_rules(KB, FALSE)
}
