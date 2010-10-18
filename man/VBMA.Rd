\name{VBMA}
\alias{VBMA}
\title{Variational Bayesian Models Averaging}

\description{ Based on optimal variational weights, it computes the averaged estimator of the posterior probabilities.}

\usage{
VBMA(data,kmin,kmax,dir.trans,gamma,normal,FixedDist,eqvar,
IterationMax,Threshold)
}


\arguments{
  \item{data}{a numeric vector of observations}
  \item{kmin}{an integer corresponding to the minimal number of mixture
          components. The default is 2.}  
  \item{kmax}{an integer corresponding to the maximal number of mixture
          components. The default is 5.}  
  \item{dir.trans}{a numeric vector specifying the prior on the transition probabilities.
          The default is c(1,1).}
\item{gamma}{a numeric vector specifying the prior on the precision. It corresponds to the hyperparameters of the gamma distribution. The default is c(0.01,0.01).}
\item{normal}{a numeric vector specifying the prior on the mean. It corresponds to the hyperparameters of the normal distribution. The default is c(-4,0.01).}
\item{FixedDist}{a numeric vector containing the parameters of the fixed gaussian distribution. The default corresponds to the standard gaussian distribution, i.e c(0,1).}
\item{eqvar}{if equals TRUE, the Gaussian distributions for the alternative are supposed to have the same variance. The default is TRUE.}
\item{IterationMax}{an integer that corresponds to the maximum number of iterations allowed for the EM algorithm. The default is 1000.}
\item{Threshold}{a positive real number corresponding to the threshold of convergence for the VBEM algorithm. The default is 1e-3.}
}

\value{

\code{VBMA} returns a list with items:

  \item{weight}{vector of the optimal variational weights.}
  \item{VBEM.estimation}{a list of objects dealing with the inference parameters for each model m.}
  \item{Converge}{a vector containing boolean values which indicates whether the algorithm converged or not.}
  \item{TauAverage}{aggregated estimator of the posterior probabilities for the cluster of interest.}
}



\examples{
## Define the parameters 
n = 1000
Mat.transition = matrix(c(0.7,0.3,0.2,0.8),2,2,byrow=TRUE)
init = rep(1/2,2)
mean = c(0,10,20)
sd   = c(1,1,1)
eta  = c(0.7,0.3)

## Simulate the data
simul = rKnorm(n,Mat.transition,init,mean,sd,eta)
data = simul[[1]]

result = VBMA(data,kmax=4)

}

\keyword{internal}


