#Step1)Initialization: Determining the Initial GMM Parameters
set.seed(1)
#Function to generate Gaussian mixture for given parameters
get.mixture.sample<-function(samples.mu,samples.sd,samples.size){
k <- length(samples.size)
mixture <- c()
for(i in seq(1:k)){
#generating normal random samples for mixture
mixture<-append(mixture,rnorm(n=samples.size[i],mean=samples.mu[i],sd=samples.sd[i]))
}
mixture <- sample(mixture)#shuffling
return(mixture)
}

#Function to obtain information from k-means
get.kmeans.results<-function(x,k){
identity <- kmeans(x,k)$cluster
mu.vec<-c()
sd.vec<-c()
lambda.vec<-c()
for(i in seq(1:k)){
mu.vec[i] <- mean(x[identity==i])
sd.vec[i] <- sd(x[identity==i])
lambda.vec[i] <- sum(identity==i)/length(identity)
}
res<-data.frame(mean=mu.vec,sd=sd.vec,lambda=lambda.vec)
return(res)
}


#Step2): Expect-Maximize
expect.maximize.GMM<-function(x,lambda.vec.init,mu.vec,sd.vec,epsilon){
k <- length(mu.vec)
lambda.matrix<-matrix(NA,nrow=k,ncol=length(x))#to store distribtuions of each lambda
Q     <- c()
Q[1]  <- 0
Q[2]  <- 0
for(i in seq(1:k)){
Q[2]   <- Q[2] + sum(log(lambda.vec.init[i] * dnorm(x, mu.vec[i], sd.vec[i])))
}
j <- 2

comp <- matrix(NA,nrow=k,ncol=length(x))

while (abs(Q[j]-Q[j-1])>=epsilon) { #checking convergence
  # E step
  for(i in seq(1:k)){
    comp[i,] <- lambda.vec.init[i] * dnorm(x, mu.vec[i], sd.vec[i])
  }
  #cat('\ndim(comp):',dim(comp))
  comp.sum  <- colSums(comp)
  #cat('\nlength(comp.sum):',length(comp.sum))
  
  for(i in seq(1:k)){
    lambda.matrix[i,] <- comp[i,]/comp.sum
  }
  #cat('\ndim(lambda.matrix):',dim(lambda.matrix))

  # M step
  for(i in seq(1:k)){
    lambda.vec.init[i] <- sum(lambda.matrix[i,])/length(x)
    #cat('\n',i,':',lambda.vec.init)
    mu.vec[i]          <- sum(lambda.matrix[i,]*x)/sum(lambda.matrix[i,])
    sd.vec[i]          <- sqrt(sum(lambda.matrix[i,]*(x-mu.vec[i])^2)/sum(lambda.matrix[i,]))
    lambda.matrix[i,]  <- lambda.vec.init[i]
    }
  #cat('\nlength(comp.sum),length(x)',length(comp.sum),length(x))
  
  j <- j + 1
  Q[j] <- sum(log(comp.sum))
  #cat('\nQ[j]=',Q[j],' j=',j)
}
result<-data.frame(mean=mu.vec,sd=sd.vec,lambda=lambda.matrix[,1],likelihood=Q[j])
return(result)
}

##Step3)Running the EM and stopping at convergence
#generate a gaussian mixture with k=3
mus<-c(10,15,20)
sds<-c(1,1,1)
sizes<-c(100,200,300)
x<-get.mixture.sample(mus,sds,sizes)
plot(density(x),main = 'Simulated GMM')
########################################
k <- length(mus)
#random initializations
init.values<-data.frame(rand.lambda=c(0.1,0.8,0.1),rand.means=c(0,0,0),rand.sd=c(10,15,20))
final.result<-expect.maximize.GMM(x,lambda.vec.init=init.values$rand.lambda,
                                  init.values$rand.means,init.values$rand.sd,epsilon=1e-6)
init.values#initializations of parameters with random values 
final.result[,1:3]
#These are the estimated parameters by EM with random initializations of lambdas as below

kmeans.result<-get.kmeans.results(x,k)
final.result2<-expect.maximize.GMM(x,lambda.vec.init=kmeans.result$lambda,mu.vec=kmeans.result$mean,
                                   sd.vec=kmeans.result$sd,epsilon=1e-6)
#initializations of parameters with values obtained via k=3-means clustering 
kmeans.result
final.result2[,1:3]
