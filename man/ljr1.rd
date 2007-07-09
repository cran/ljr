\name{ljr1}
\alias{ljr1}
\title{MLE with 1 joinpoint}
\description{
   Determines the maximum likelihood estimates of model coefficients
in the logistic joinpoint regression model with one joinpoint.
}
\usage{
ljr1(y,n,tm,X,ofst,summ=TRUE)
}
\arguments{
   \item{y}{the vector of Binomial responses.}
   \item{n}{the vector of sizes for the Binomial random variables.}
   \item{tm}{the vector of ordered observation times.}
   \item{X}{a design matrix containing other covariates.}
   \item{ofst}{a vector of known offsets for the logit of the response.}
   \item{summ}{a boolean indicator of whether summary tables should be returned.}
}
\value{
   \item{Coef}{A table of coefficient estimates.}   
   \item{Joinpoint}{The estimate of the joinpoint.}
   \item{wlik}{The maximum value of the re-weighted log-likelihood.}
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
 \code{\link{ljrb2},\link{ljrf2},\link{ljr01},link{ljr12},\link{ljrb},\link{ljrf}}
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
 temp.ljr=ljr1(y,n,tm,X,ofst)
}
\keyword{nonlinear}
