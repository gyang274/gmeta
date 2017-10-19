\name{np.gmeta}
\alias{np.gmeta}

\title{
	Make inference for the study-level parameters with no distributional assumption
}
\description{
This method synthesizes information from multiple studies and make inference that is not dependent on any distributional assumption for the study-level parameters. Specifically, the study-level parameters are assumed to be unknown, fixed parameters, it draws inferences about the quantiles of this set of parameters using study-specific summary statistics
}
\usage{
np.gmeta(Thetahat, se, alpha = c(0.025, 0.975), n, m, 
band_pwr = 0.5, resample = 200, B = 40, len = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Thetahat}{
input, point estimate of the true parameter for all K studies. For the kth study, the kth \code{Thetahat} is consistent estimator for the true parameter, the kth \code{Theta}.
}


\item{se}{
input,  standard error estimate of th e true parameter estimator for all K studies. For the kth study, there is a standard error estimate, denoted as the kth \code{se}, which is associated with the kth \code{Thetahat}..
}

\item{alpha}{
quantile vector by default.
}
  
\item{n}{
total sample size of all K studies.
}
  
\item{m}{
mth ordered parameter.
}
  
\item{band_pwr}{
a constant in (0,1).
}
  
\item{resample}{
tuning parameter: R realization of CD-random variables for a confidence distribution for the mth ordered parameter.

}
  
\item{B}{
tuning parameter: process with B new "observed" data.
}
  
\item{len}{
tuning parameter: grid search for possible pairs.
}
}

\details{
	The function produces point or quantile estimation for the parameter whether tie or near tie condition exists or not.}


\value{
  An object of class "gmeta.interval", which is a list of following elements:
  
\code{percentils}       
A three dimensional array containing the min.unif, max.smooth, mean.smooth.
\code{min.unif} contains quantiles using minimum pair.
\code{max.smooth} contains quantiles using maximum pair.
\code{mean.smooth} contains quantiles using mean of all satisfied pairs. 

\code{shrink}  The shrinkage used to obtain presumed "true values".

\code{smoohlist} A sequence.

\code{distance} The ten by ten dimensional array containing results of loss function. 

\code{elig.ind} A fourteen by two dimensional array containing counts that how many pairs are smaller than threshold. 

}
\author{
Wei Qian <wq24@rutgers.edu>, Jerry Q. Cheng <jcheng1@rwjms.rutgers.edu>
}
\note{
Revised on 2017/10/01.
}

\references{
  Brian Claggett, Minge Xie & Lu Tian(2014)
Meta-Analysis With Fixed, Unknown, Study-Specific Parameters.
\emph{Journal of the American Statistical Association}, \bold{109} 1660-1671.
}

\seealso{
	\code{\link{gmeta}}
}

\examples{
 Thetahat = c(-0.102158346,  0.020746333,  0.009118397,  0.163044549, -0.098892204,
              0.161200470,  0.237464018,  0.162000380, -0.088128443,  0.337433537,
              0.032277816,  0.142151631, -0.225430197, -0.057114409,  0.100874063,
              0.177078003, -0.752891602,  0.069141934,  0.114787217, 0.239802656, -0.029858223)
 se = c(0.1742206, 0.1342359, 0.1660455, 0.1372563, 0.1756401, 
 0.1288152, 0.2001390, 0.1603925, 0.1276070, 0.1470815, 
 0.1666182, 0.1694948, 0.1602534, 0.1688223, 0.1512014, 
 0.1857019, 0.1236781, 0.1510481, 0.1733524, 0.1627342, 0.1905231) 
 np.gmeta(Thetahat = Thetahat, se = se, m=10, 
 n =rep(40,21), band_pwr = 0.5, resample=200, B=40, len=10)
}
\keyword{meta analysis}
