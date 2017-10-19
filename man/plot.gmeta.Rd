\name{plot.gmeta}
\alias{plot.gmeta}
\alias{plot.gmeta.m}
\alias{plot.gmeta.e}
\title{
	Plot For Model-Based and 2x2 Table-Based Meta-Analysis Result
}
\description{
	Produce a forest-like plot for object of \code{gmeta} class, the result of model-based or 2x2-table-based meta-analysis using \code{gmeta} function.
}
\usage{
\method{plot}{gmeta}(gmo, studies=NULL, 
	plot.option=c('confidence-density',
		'confidence-curve', 'cv',
		'confidence-distribution', 'cdf'), 
	type='l', xlab='x', ylab='density', xlim=NULL, ylim=NULL, ...)
}
\arguments{
	\item{gmo}{an object of \code{gmeta} holding results of model-based or 2x2-table-based meta-analysis from \code{gmeta} function.}
	\item{studies}{a vector specifying the index of individual studies shown on the plot beside the combined one.}
	\item{plot.option}{an option for plot.  The default is \code{'confidence-density'} for CD density. Other choices are \code{'confidence-curve'} or \code{'cv'} for confidence curve, \cr and \code{'confidence-distribution'} or \code{'cdf'} for CD.}
	\item{type}{as the \code{type} argument in \code{plot} with the default of \code{'l'}.}
	\item{xlab}{as the \code{xlab} argument in \code{plot} with the default of \code{'x'}.} 
	\item{ylab}{as the \code{ylab} argument in \code{plot} with the default of \code{'confidence density'} when \code{plot.option='confidence-density'}. Other choices are \code{'confidence curve'} when \code{plot.option} is \code{'confidence-curve'} or \code{'cv'}, and \code{'confidence distribution'} when \code{plot.option} is \code{'confidence-distribution'} or \code{'cdf'}.}
	\item{xlim}{as the \code{xlim} argument in \code{plot}.}
	\item{ylim}{as the \code{ylim} argument in \code{plot}, scaled if \code{studies} is specified.}
	\item{...}{other arguments that can be specified in \code{plot}.}
}
\value{
  A figure of combined and inidividual confidence densities, curves, or distributions.
}
\details{
	This function produces a forest-like plot for model-based or 2x2-table-based meta-analysis from the result of \code{gmeta} function call. The argument \code{plot.option} provides a choice of displaying confidence densities, curves, or distributions.
}
\author{
	Guang Yang <gyang.rutgers@gmail.com>, Jerry Q. Cheng <jcheng1@rwjms.rutgers.edu> and Minge Xie <mxie@stat.rutgers.edu>
}
\note{
	Revised on 2014/12/10.
}

\references{
	Xie, M. and Singh, K. (2013)
	Confidence distribution, the frequentist distribution estimator of a parameter (with discussions). 
	\emph{International Statistical Review}, \bold{81} 3-39.

	Xie, M., Singh, K., and Strawderman, W. E. (2011).
	Confidencedence distributions and a unifying framework for meta-analysis.
	\emph{Journal of the American Statistical Association}, \bold{106} 320-333.
	
	Yang, G., Liu, D., Wang, J. and Xie, M. (2016). 
	Meta-analysis framework for exact inferences with application to the analysis of rare events.
	\emph{Biometrics}, \bold{72} 1378-1386.
}
\seealso{
	\code{\link{gmeta}}
}
\examples{
#### gmeta: generalized meta-analysis approach ####

data(ulcer)
ulcer.o <- as.matrix(ulcer)

# impute 0.5
ulcer <- ifelse(ulcer.o == 0, 0.5, ulcer.o)
# summary statistics 
ulcer.theta <- log( (ulcer[,1]*ulcer[,4]) / (ulcer[,2]*ulcer[,3]) ) 
ulcer.sigma <- sqrt(1/ulcer[,1] + 1/ulcer[,2] + 1/ulcer[,3] + 1/ulcer[,4]) 
ulcer.pivots = data.frame(mns=ulcer.theta, sds=ulcer.sigma)

# fixed-effect model
gmo.mdlfx <- gmeta(ulcer.pivots, method='fixed-mle', gmo.xgrid=seq(from=-10,to=10,by=0.01))
summary(gmo.mdlfx)

# random-effects model, method of moments
gmo.mdlrm <- gmeta(ulcer.pivots, method='random-tau2', weight=rep(1,41), tau2=2, 
	gmo.xgrid=seq(from=-10,to=10,by=0.01)) 
summary(gmo.mdlrm)

# Forest plot of CDs
plot(gmo.mdlrm, studies=c(4,8,15,16,23,41)) # default: confidence-distribution-density
plot(gmo.mdlrm, studies=c(4,8,15,16,23,41), plot.option='cv') # using confidence-curve

# 2x2 table-based meta-analysis 
ulcer.2x2 <- cbind(ulcer[,1], ulcer[,1]+ulcer[,2], ulcer[,3], ulcer[,3]+ulcer[,4])
# Mantel-Haenszel odd-ratio
gmo.2x2MH <- gmeta(ulcer.2x2, gmi.type='2x2', method='MH', gmo.xgrid=seq(-5,5,by=0.001))
summary(gmo.2x2MH)
plot(gmo.2x2MH, studies=c(4,8,15,16,23,41))

# Peto's log-odds-ratio
gmo.2x2Pt <- gmeta(ulcer.2x2, gmi.type='2x2', method='Peto', gmo.xgrid=seq(-5,5,by=0.001))
summary(gmo.2x2Pt)
plot(gmo.2x2Pt, studies=c(4,8,15,16,23,41))

# Exact meta-analysis on log-odds-ratio based on Liu et al (2012) 
ulcer.exact <- cbind(ulcer.o[,1], ulcer.o[,1]+ulcer.o[,2], ulcer.o[,3], ulcer.o[,3]+ulcer.o[,4])
#gmo.exact1 <- gmeta(ulcer.exact, gmi.type='2x2', method='exact1',
#	gmo.xgrid=seq(-5,5,by=0.001), report.error=TRUE) 
#summary(gmo.exact1)
#plot(gmo.exact1, studies=c(4,8,15,16,23,41))

# Exact meta-analysis on risk difference based on Tian et al (2009)
#gmo.exact2 <- gmeta(ulcer.exact, gmi.type='2x2', method='exact2', 
#	gmo.xgrid=seq(-1,1,by=0.001), report.error=TRUE) 
#summary(gmo.exact2)
#plot(gmo.exact2, studies=c(4,8,15,16,23,41), plot.option='cv')
}
\keyword{meta-analysis}
