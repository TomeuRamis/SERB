\name{writeKB}
\alias{writeKB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Write Knowledge Base
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Write a Knowledge Base structure to an external file.
}
\usage{
writeKB(file, KB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KB}{Knowledge Base to write to file.}
  \item{file}{String name of the file name where you want to store the KB. If the file does   not exist it creates it.}
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
Does not return anything.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Bartomeu Ramis Tarragó
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Create a example KB
KB <- data.frame(var1=c(1, 3, 3, 4, 5),
               var2=c(7, 7, 8, 3, 2),
               var3=c(3, 3, 6, 6, 8),
               var4=c(1, 1, 2, 8, 9))
#Save it to file_output.txt
writeKB("file_output.txt", KB)
}
