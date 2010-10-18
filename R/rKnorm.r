  rKnorm <- function(n,Mat.transition,init,mean,sd,eta)
  {

  ## Check
  stopifnot(setequal(rowSums(Mat.transition),c(1,1)))
  stopifnot(min(sd)>0)
  stopifnot(min(Mat.transition)>0) 
  stopifnot(min(eta)>0)

  Markov.chain	= rmultinom(1,1,init)
  class		= which(Markov.chain[,1]==1)
  
  if(class[1]==1)
  {
  data  = rnorm(1,mean=mean[1],sd=sd[1]) 
  }
  if(class[1]==2)
  {
  group.mix 	= rmultinom(1,1,eta)
  data  	= rnorm(1,mean=mean[which(group.mix==1)+1],sd=sd[which(group.mix==1)+1])
  }

  for(i in 2:n)
  {
    Group.neigh = class[i-1]
    Z 		= rmultinom(1,1,Mat.transition[Group.neigh,])
    if(Z[1,1]==1)
    {
    data  = c(data,rnorm(1,mean=mean[1],sd=sd[1]))
    }
    if(Z[2,1]==1)
    {
    group.mix	= rmultinom(1,1,eta)
    data  	= c(data,rnorm(1,mean=mean[which(group.mix==1)+1],sd=sd[which(group.mix==1)+1]))
    }
  
    class = c(class,which(Z[,1]==1))
  }

  return(list(data=data,class=class))
  }



