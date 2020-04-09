Metropolis.fn<-function(y,n,X,c,Sigma,iter,burnin,chain){ 
#Inputs:
#y: vector of responses
#n: vector (or scalar) of trial sizes. 
#X: predictor matrix including intercept.
#c: rescaling for variance-covariance matrix, scalar J(theta*|theta(t-1)) = N(theta(t-1), c^2*Sigma)
#iter: number of iterations
#burnin: number of initial iterations to throw out.

ar=0#acceptance rate
p <-dim(X)[2]   #number of parameters
library(mvtnorm)#library for multivariate gaussian, used for sampling
theta0<-rnorm(p) #initial values.
theta.sim<-matrix(0,iter,p) #matrix to store iterations
theta.sim[1,]<-theta0   
for(i in 1:(iter-1)){
theta.cand <-rmvnorm(1,mean=theta.sim[i,],sigma=c^2*Sigma)
theta.cand <-as.numeric(theta.cand)  
xbc        <-X%*%theta.cand      
lambda.c   <-exp(xbc)
xb         <-X%*%theta.sim[i,]
lambda.b   <-exp(xb)
#difference of log joint distributions.
r <- exp(sum(dpois(y,lambda = lambda.c,log=TRUE)) - sum(dpois(y,lambda=lambda.b,log=TRUE)))
#Draw an indicator whether to accept/reject candidate
ind <- min(r,1)
if(ind==1){
  ar=ar+1
  }
# jumping with probability ind
theta.sim[i+1,]<- ind*theta.cand + (1-ind)*theta.sim[i,]
}
ar<-ar/10000
cat("Acceptance rate of",chain,":", ar)
results<-theta.sim[-c(1:burnin),]#Removing the iterations in burnin phase
names(results)<-c('beta_intercept','beta_woolB','beta_tensionL','beta_tensionM') #column names
return(results)
}

poisson.model<-glm(breaks~wool+tension, data=warpbreaks, family = poisson(link = "log"))
poisson.model$coefficients

X<-model.matrix(poisson.model)
y<-warpbreaks$breaks
n<-length(y)
sigma<-vcov(poisson.model)
p <-dim(X)[2]

#vary burnin and iterations, also change c=2.4/sqrt(p) to c=1.6/sqrt(p),c=3.2/sqrt(p)
chain<-Metropolis.fn(y=y,n=n,X=X,c = /sqrt(p),Sigma=sigma,iter=10000,burnin=5000,'chain')
#Posterior means
colMeans(chain)
#Posterior standard deviations
apply(chain,2,FUN=sd)
#95 % credible intervals  
apply(chain,2,FUN=function(x) quantile(x,c(0.025,0.975)))

library(coda)
#Estimating Gelman -Rubin diagnostics.
ml<-as.mcmc.list(as.mcmc((chain)))

estml<-c(ml)

#Gelman-Rubin diagnostic.
gelman.diag(estml)[[1]]#values should be close to 1.01

effectiveSize(estml)#number of independent samples, the more the better
plot(ml)#plotting diagnostic plots
