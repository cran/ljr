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
\references{ 
 Czajkowski, M., Gill, R. and Rempala, G. (2008). Model selection in logistic joinpoint regression with applications to analyzing cohort mortality patterns. {\emph Statistics in Medicine} 27, 1508-1526.
}
\seealso{
 \code{\link{ljr01},\link{ljrb},\link{ljrf}}
}
\examples{
 data(kcm)
 attach(kcm)
 ljr1(Count,Population,Year+.5)
}
\keyword{nonlinear}
