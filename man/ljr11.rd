\name{ljr11}
\alias{ljr11}
\title{Test coefficients conditioned on K=1 joinpoint.}
\description{
 This function performs the likelihood ratio tests to find p-values in testing the significance of each of the coefficients as well as the intercept and ordered observation times.  The p-values are determined by a Monte Carlo method.
}
\usage{
ljr11(y,n,tm,X,ofst,R=1000)
}
\arguments{
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
\references{ Czajkowski, M., Gill, R. and Rempala, G. (2007). Model selection in logistic joinpoint regression with applications to analyzing cohort mortality patterns. To appear.
}
\seealso{
 \code{\link{ljr1},\link{ljr12}}
}
\examples{
 N=20
 m=2
 k=1
 beta=c(0.1,0.1,-0.05)
 gamma=c(0.1,-0.05)
 tau=c(5)
 ofst=runif(N,-2.5,-1.5)
 x1=round(runif(N,-0.5,9.5))
 x2=round(runif(N,-0.5,9.5))
 X=cbind(x1,x2)
 n=rep(1e9,N)
 tm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
 eta=ofst+beta[1]+gamma[1]*tm
 if (m>0)
 for (i in 1:m)
  eta=eta+beta[i+1]*X[,i]
 if (k>0)
  for (i in 1:k) 
   eta=eta+gamma[i+1]*pmax(tm-tau[i],0) 
 y=rbinom(N,size=n,prob=exp(eta)/(1+exp(eta)))
 temp.ljr=ljr11(y,n,tm,X,ofst,R=1000)
}
\keyword{nonlinear}
