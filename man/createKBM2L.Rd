\name{createKBM2L}
\alias{createKBM2L}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Creation of Knowledge Base Matrix To List structure
}
\description{
Creates a KBM2L list from a KB structure. A KBM2L is a list that reduces the space needed to store a KB structure with, almost, no information loss.}
\usage{
createKBM2L(KB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KB}{Knowledge Base table from which to create the KBM2L.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list that contains the rules that form the KB entered as input.
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

KBM2L <- createKBM2L(optimalKB)
}
