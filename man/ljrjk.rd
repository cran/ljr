\name{ljrjk}
\alias{ljrjk}
\title{Perform test of j vs k joinpoints.}
\description{
 This function tests the null hypothesis of j joinpoint(s) versus the alternative of k joinpoint(s) based on the likelihood ratio test statistic.  The p-value is determined by a Monte Carlo method.
}
\usage{
ljrjk(j,k,y,n,tm,X,ofst,R=1000,alpha=.05)
}
\arguments{
   \item{j,k}{pre-specified number of joinpoints in the null and alternative hpyotheses (the smaller is used for the null).}
   \item{y}{the vector of Binomial responses.}
   \item{n}{the vector of sizes for the Binomial random variables.}
   \item{tm}{the vector of ordered observation times.}
   \item{X}{a design matrix containing other covariates.}
   \item{ofst}{a vector of known offsets for the logit of the response.}
   \item{R}{number of Monte Carlo simulations.}
   \item{alpha}{significance level of the test.}
}
\value{
   \item{pval}{The estimate of the p-value via simulation.}
   \item{Coef}{A table of coefficient estimates.}
   \item{Joinpoint}{The estimates of the joinpoint, if it is significant.}
   \item{wlik}{The maximum value of the re-weighted log-likelihood.}
}
\details{
 The re-weighted log-likelihood is the log-likelihood divided by the largest component of n. 
}
\author{
The authors are Michal Czajkowski, Ryan Gill, and Greg Rempala.
The software is maintained by Ryan Gill \email{rsgill01@louisville.edu}.
}
\references{ 
 Czajkowski, M., Gill, R. and Rempala, G. (2008). Model selection in logistic joinpoint regression with applications to analyzing cohort mortality patterns. {\emph Statistics in Medicine} 27, 1508-1526.
}
\seealso{
 \code{\link{ljrk}}
}
\examples{
 data(kcm)
 attach(kcm)
 set.seed(12345)
 ljrjk(0,1,Count,Population,Year+.5,R=20)
}
\keyword{nonlinear}
