\name{ljrkk}
\alias{ljrkk}
\title{Test coefficients conditioned on K=k joinpoint.}
\description{
 This function performs the likelihood ratio tests to find p-values in testing the significance of each of the coefficients as well as the intercept and ordered observation times.  The p-values are determined by a Monte Carlo method.
}
\usage{
ljrkk(k,y,n,tm,X,ofst,R=1000)
}
\arguments{
   \item{k}{the pre-specified number of joinpoints (with unknown locations).}
   \item{y}{the vector of Binomial responses.}
   \item{n}{the vector of sizes for the Binomial random variables.}
   \item{tm}{the vector of ordered observation times.}
   \item{X}{a design matrix containing other covariates.}
   \item{ofst}{a vector of known offsets for the logit of the response.}
   \item{R}{number of Monte Carlo simulations.}
}
\value{
   \item{pvals}{The estimates of the p-values via simulation.}
}
\details{
 The re-weighted log-likelihood is the log-likelihood divided by the largest component of n. 
}
\author{
The authors are Michal Czajkowski, Ryan Gill, and Greg Rempala.
The software is maintained by Ryan Gill \email{rsgill01@louisville.edu}.
}
\references{ 
 Czajkowski, M., Gill, R. and Rempala, G. (2008). Model selection in logistic joinpoint regression with applications to analyzing cohort mortality patterns. \emph{Statistics in Medicine} 27, 1508-1526.
}
\seealso{
 \code{\link{ljrk}}
}
\examples{
 data(kcm)
 attach(kcm) 
 set.seed(12345)
## Not run: ljrkk(1,Count,Population,Year+.5,R=20) 
}
\keyword{nonlinear}
