ljr22 <- function(y,n,tm,X=NULL,ofst=0,R=1000)
{
 N=length(y)
 m=ncol(X)
 if (is.null(m)) 
  m=0
 if (length(ofst)==1)
  ofst=as.double(rep(ofst,N))
 else
  ofst=as.double(ofst)
 out=.C("ljr22",as.double(y),as.double(n),as.double(tm),as.double(X),ofst,N,m,as.integer(R),p=double(m+2),PACKAGE="ljr")
 cat('Model:\n')
 cat('y~Binom(n,p) where p=invlogit(eta)\n')
 if (m>0){
  if (is.null(dimnames(X)[[2]])==FALSE)
   m.variables=dimnames(X)[[2]]
  else
   m.variables=paste("X",1:ncol(X),sep="")
 }
 else
  m.variables=NULL
 m.variables=c('Intercept',m.variables)
 t.variables=c('t','max(t-tau1,0)','max(t-tau2,0)')
 if ((ofst[1]==0)&(length(ofst)==1))
  cat('eta=b0')
 else
  cat('eta=ofst+b0')
 if (length(m.variables)>1)
  for (i in 2:length(m.variables))
   cat(paste(paste('+b',i-1,sep=""),'*',m.variables[i],sep=""))
 for (i in 0:(length(t.variables)-1))
  cat(paste(paste('+g',i,sep=""),'*',t.variables[i+1],sep=""))
 cat("\n\n")
 if (length(m.variables)>1)
  f1=data.frame(Variables=c(m.variables[-1],'Intercept','tm'),p.values=out$p)
 else
  f1=data.frame(Variables=c('Intercept','tm'),p.values=out$p)
 print(f1)
 return(list(pvals=out$p))
}
