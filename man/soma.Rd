\name{soma}
\alias{soma}
\title{The Self-Organising Migrating Algorithm}
\description{
The Self-Organising Migrating Algorithm is a general-purpose, stochastic optimisation algorithm. The approach is similar to that of genetic algorithms, although it is based on the idea of a series of ``migrations'' by a fixed set of individuals, rather than the development of successive generations. It can be applied to any cost-minimisation problem with a finite parameter space, and is robust to local minima.
}
\usage{
soma(costFunction, bounds, options = list(), strategy = "all2one", ...)
}
\arguments{
  \item{costFunction}{A cost function which takes a numeric vector of parameters as its first argument, and returns a numeric scalar representing the associated cost value.}
  \item{bounds}{A list with elements \code{"min"} and \code{"max"}, each a numeric vector giving the upper and lower bounds for each parameter, respectively.}
  \item{options}{A list of options for the SOMA algorithm itself. See Details.}
  \item{strategy}{The strategy type to use. Currently \code{"all2one"} is the only supported value.}
  \item{\dots}{Additional parameters to \code{costFunction}.}
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
}
\references{
I. Zelinka (2004). SOMA - self-organizing migrating algorithm. In G.C. Onwubolu & B.V. Babu, eds, New optimization techniques in engineering. Volume 141 of ``Studies in Fuzziness and Soft Computing'', pp. 167-217. Springer.
}
\author{R port by Jon Clayden <jon.clayden@gmail.com>}
\note{
%%  ~~further notes~~
}
\seealso{\code{\link{optim}}}
\examples{
}
