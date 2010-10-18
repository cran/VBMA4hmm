\name{rKnorm}
\title{Log-Density for Multinomial Distribution}
\alias{rKnorm}
\usage{
rKnorm(n,Mat.transition,init,mean,sd,eta)
}

\description{
Simulate a 2-states HMM with a gaussian distribution(S=0) and a mixture of gaussian(S=1) for the emission probability.
}
\arguments{
\item{n}{an integer corresponding to the number of simulated data.}
\item{Mat.transition}{a 2x2 matrix containing the transition probabilities between the two states.}
\item{init}{a numeric vector specifying the initial distribution of the 2-states HMM. The sum must be equal to 1 and the default is c(0.5,0.5). }
\item{mean}{vector of means. The first corresponds to S=0.}
\item{sd}{vector of standard deviations. The first corresponds to S=0.}
\item{eta}{vector of mixture proportion for S=1. The sum must be equal to 1.}
}

\value{
\code{rKnorm} returns a list with items:

  \item{data}{a numeric vector of simulated data.}
   \item{class}{a vector containing the theoretical class of each observation.}
}


\examples{
n = 1000
Mat.transition = matrix(c(0.7,0.3,0.2,0.8),2,2,byrow=TRUE)
init = rep(1/2,2)
mean = c(0,10,20)
sd   = c(1,1,1)
eta  = c(0.7,0.3)

## Simulate the data
simul = rKnorm(n,Mat.transition,init,mean,sd,eta)
data = simul[[1]]

}


\keyword{internal}
