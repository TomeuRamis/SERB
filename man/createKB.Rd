\name{createKB}
\alias{createKB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Create Knowledge Base
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Create a Knowledge Base from a given Bayesian Network and a Target node from that same network. The table structure created contains two new columns: Decision and Prob.
}
\usage{
createKB(BN, node)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BN}{Must be a Discrete Bayesian Network created with the \dQuote{bnlearn} package.}
  \item{node}{String name of the Target Node we want to study.}
}
\details{
This function creates a table with all possible combinations of input values for the target node. From this table we query the model repeatedly, saving the Decision given and the Probability of said Decision.
Each Decision and Probability is saved alongside the evidence that generated it, and it is return in the structure of a Table, called Knowledge Base.
}
\value{
A Table with all possible evidence combinations and the value and probability of the Target node for each instance.
}
\references{
%% ~put references to the literature/web site here ~
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
#Create the KB from a preloaded Bayesian Network (bn)
KB <- createKB(bn, "tub")
}
