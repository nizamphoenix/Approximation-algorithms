normalmm.Gibbs<-function(iter,X,y,burnin,beta_0.vec,K,tau_0,a,b){
  n   <-length(y) #no. observations
  p   <-dim(X)[2] #no of fixed effect predictors.
  tau<-tau_0
  library(mvtnorm)
  
  #storing results.
  par <-matrix(0,iter,p+1)  #p beta coefficients, and 1 precision coefficient.
  Tyrep <<- Ty <<- 0

  
  for(i in 1:iter){
    #Conditional posteriors.
   
    
    XTX <- t(X)%*%X
    K.inv <- solve(K)
    mean.vector <- solve(XTX + K.inv) %*% (t(X) %*% y   +  K.inv %*% beta_0.vec)
    #dimensions:         pxp    pxp        pxn     nx1       pxp       px1
    #dimensions:            pxp                px1                px1
    #dimensions:            pxp                        px1       
    #dim(mean.vector):                      px1 
    
    cov.matrix <-  solve(XTX + K.inv)/tau
    #dim(cov.matrix):    pxp
   
    beta.vec  <-rmvnorm(1,mean=mean.vector,sigma=cov.matrix) #sample beta vector
    beta.vec <-as.numeric(beta.vec)
    
    err      <- y-X%*%beta.vec
    
    alpha <- a+(n+p)/2
    gamma <- b+(sum(err^2)+t(beta.vec - beta_0.vec)%*%K.inv%*%(beta.vec - beta_0.vec))/2
    tau   <-rgamma(1,alpha,gamma) #sample tau
    par[i,]<-c(beta.vec,1/tau)
    
    #Posterior predective checking.
    X <- #THE PREDICTORS
    y <- #TARGET
    n <-length(y)

    err.full<-y-X%*%beta.vec#residual for the entire data set
    Ty[i] <<- sum(err.full^2)#The original test quantity
    yrep <- rnorm(n, mean = X%*%beta.vec, sd =1/sqrt(tau))#generating the replicated data
    err.rep <- yrep - X%*%beta.vec
    Tyrep[i] <<- sum(err.rep^2)#The replicated(synthetic) test quantity
    
 }
par <-par[-c(1:burnin),] #removing 1000 initial iterations
colnames(par)<-c('beta.intercept','beta.time','sigma.square')
thin.seq <- seq(2, dim(par)[1], 2)#Performing thinning by discarding every 2nd sample
par <-par[-thin.seq,]     
return(par) 
}


X0<-cbind(1,hiron.data[hiron.data$homc282y==0,]$time)
y0<-hiron.data[hiron.data$homc282y==0,]$logsf
n0<-length(y0)
p   <-dim(X0)[2]

sample.var <- (summary(lm.model)$sigma)^2#as in 2a
err1     <- y1-X1%*%betahat1
#s^2 = sum(err1^2)/(n1-p) = sample.var
round(sum(err1^2)/(n1-p),5)==round(sample.var,5)


chain1<-normalmm.Gibbs(iter=10000,X=X0,y=y0,burnin=1000,beta_0.vec=betahat1,K= XTXinv1,tau_0=0.1,a=(n1-p)/2,b=(n1-p)*sample.var/2)

chain2<-normalmm.Gibbs(iter=10000,X=X0,y=y0,burnin=1000,beta_0.vec=betahat1,K= XTXinv1,tau_0=1.0,a=(n1-p)/2,b=(n1-p)*sample.var/2)

plot(chain1[,1],type='l',col='green',ylab=expression(beta[0]),main='Trace plot of the intercept coefficient')
lines(chain2[,1],type='l',col='yellow',ylab=expression(beta[0]))
legend('topright',legend=c('chain1','chain2'),col=c('green','yellow'),lty=1,bty='n')

plot(chain1[,2],type='l',col='red',ylab=expression(beta[1]),main='Trace plot of the time coefficient')
lines(chain2[,2],type='l',col='black',ylab=expression(beta[1]))
legend('topright',legend=c('chain1','chain2'),col=c('red','black'),lty=1,bty='n')

plot(chain1[,3],type='l',col='darkblue',ylab=expression(sigma^2),main='Trace plot of the variance')
lines(chain2[,3],type='l',col='gray',ylab=expression(sigma))
legend('topright',legend=c('chain1','chain2'),col=c('darkblue','gray'),lty=1,bty='n')
