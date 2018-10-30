%% File Name: pbvnorm.Rd
%% File Version: 0.09

\name{pbvnorm}
\alias{pbvnorm}
\alias{dbvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Probabilities for Bivariate Normal Distribution
}

\description{
The function \code{pbvnorm} computes probabilities \eqn{\Phi_2(x,y,\rho)} for 
the standardized bivariate normal distribution (Drezner & Wesolowsky, 1990;
West, 2004).

The function \code{dbvnorm} computes the corresponding density \eqn{\phi_2(x,y,\rho)}.
}

\usage{
pbvnorm(x, y, rho)

dbvnorm(x, y, rho, log=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of first ordinate
}
  \item{y}{
Vector of second ordinate
}
  \item{rho}{
Vector of correlations
}
\item{log}{Logical indicating whether logarithm of the density
should be calculated}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Vector
}

\references{
Drezner, Z., & Wesolowsky, G. O. (1990). On the computation of the bivariate normal 
integral. \emph{Journal of Statistical Computation and Simulation, 35}(1-2), 101-107.

Genz, A. (1992). Numerical computation of multivariate normal probabilities. 
\emph{Journal of Computational and Graphical Statistics, 1}(2), 141-149. 

West, G. (2005). Better approximations to cumulative normal functions. 
\emph{Wilmott Magazine, 9}, 70-76.
}

%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link[pbivnorm:pbivnorm]{pbivnorm::pbivnorm}} and
\code{mnormt::biv.nt.prob} for alternative implementations (Genz, 1992)
which are not (yet) completely vectorized.
}



\examples{
#############################################################################
# EXAMPLE 1: Comparison with alternative implementations
#############################################################################

#*** simulate different values of ordinates and correlations
set.seed(9898)
N <- 3000
x <- stats::runif(N,-3,3)
y <- stats::runif(N,-3,3)
rho <- stats::runif(N,-.95,.95)

#*** compute probabilities
res1 <- pbv::pbvnorm(x=x,y=y,rho=rho)

#-- compare results with pbivnorm package
library(pbivnorm)
res2 <- rep(NA, N)
for (ii in 1:N){
    res2[ii] <- pbivnorm::pbivnorm(x=x[ii], y=y[ii], rho=rho[ii])
}

summary(abs(res1-res2))

#*** compute density values
log <- TRUE    # logical indicating whether log density should be evaluated
res1 <- pbv::dbvnorm(x=x, y=y, rho=rho, log=log )

\dontrun{
#-- compare results with mvtnorm package
library(mvtnorm)
res2 <- rep(NA, N)
sigma <- diag(2)
for (ii in 1:N){
    sigma[1,2] <- sigma[2,1] <- rho[ii]
    res2[ii] <- mvtnorm::dmvnorm(x=c(x[ii],y[ii]), sigma=sigma, log=log)
}
summary(abs(res1-res2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
