\name{readKB}
\alias{readKB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Read Knowledge Base
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Load the Knowledge Base structure from a text file.
}
\usage{
readKB(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ String name of the file that contains the KB.
%%     ~~Describe \code{x} here~~
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
Returns the Table structure saved in the file.
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
  #File should be a KB previously saved
  #KB <- readKB("output.txt")
}
