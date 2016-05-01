\name{ljr01}
\alias{ljr01}
\title{Perform test of 0 vs 1 joinpoints.}
\description{
 This function tests the null hypothesis of 0 joinpoints versus the alternative of one joinpoint based on the likelihood ratio test statistic.  The p-value is determined by a Monte Carlo method.
}
\usage{
ljr01(y,n,tm,X,ofst,R=1000,alpha=.05)
}
\arguments{
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
 Czajkowski, M., Gill, R. and Rempala, G. (2008). Model selection in logistic joinpoint regression with applications to analyzing cohort mortality patterns. \emph{Statistics in Medicine} 27, 1508-1526.
}
\seealso{
 \code{\link{ljr0},\link{ljr1}}
}
\examples{
 data(kcm)
 attach(kcm)
 set.seed(12345)
## Not run: ljr01(Count,Population,Year+.5,R=20)
}
\keyword{nonlinear}
