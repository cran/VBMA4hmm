  VBMA<-function(data,kmin=2,kmax=5,dir.trans=c(1,1),gamma=c(0.01,0.01),normal=c(-4,0.01),FixedDist=c(0,1),eqvar=TRUE,IterationMax=1000,Threshold=1e-3)
  {
  
  ## Check
  stopifnot(kmin<=kmax)
  stopifnot(IterationMax>=1)
  stopifnot(is.vector(data))
  stopifnot(min(gamma)>0)
  stopifnot(normal[2]>0)
  stopifnot(FixedDist[2]>0)
  stopifnot(length(FixedDist)==2)

  ## Declaration and initialisation of local variables
  result      = numeric()
  resultEM    = vector(mode="list",length=abs((kmax-kmin)+1))
  nbInd       = length(data)
  Converge    = TRUE
  
  for(k in kmin : kmax){

  if(!eqvar)
  {
  CalculEst	= .VBHMM_hetero(data,dir=rep(1,k-1),dir.trans,gamma=gamma,normal=normal,FixedDist=FixedDist,IterationMax=IterationMax,Threshold=Threshold)
  a		= CalculEst$hyperparameters$a
  b		= CalculEst$hyperparameters$b
  t		= CalculEst$hyperparameters$t 
  m		= CalculEst$hyperparameters$m
  p		= CalculEst$hyperparameters$p
  p.eta		= CalculEst$hyperparameters$p.eta
  p.trans	= CalculEst$hyperparameters$p.trans
  PhiBarre	= CalculEst$vectors$PhiBarre
  loglik	= CalculEst$loglik
  Test		= CalculEst$Test
  Converge	= c(Converge,Test)

  UBarre  	     	= CalculEst$vectors$UBarre$UBarre
  PhiBarre		= PhiBarre[-c((3),(k+3),(2*k+3),(3*k+3))]
  result[(k-kmin)+1] 	= sum(PhiBarre*UBarre) + .Calcul.constant(a,b,m,t,p,p.trans,p.eta,FixedDist=FixedDist)- .Calcul.constant(a=rep(gamma[1],k-1),b=rep(gamma[2],k-1),m=rep(normal[1],k-1),t=rep(normal[2],k-1),p=rep(1,2),p.trans=rep(1,4),p.eta=rep(1,k-1),FixedDist=FixedDist) -loglik 
  }

  if(eqvar){
  CalculEst	= .VBHMM_homo(data,dir=rep(1,k-1),dir.trans,gamma=gamma,normal=normal,FixedDist=FixedDist,IterationMax=IterationMax,Threshold=Threshold)
  a		= CalculEst$hyperparameters$a
  b		= CalculEst$hyperparameters$b
  t		= CalculEst$hyperparameters$t 
  m		= CalculEst$hyperparameters$m
  p		= CalculEst$hyperparameters$p
  p.eta		= CalculEst$hyperparameters$p.eta
  p.trans	= CalculEst$hyperparameters$p.trans
  PhiBarre	= CalculEst$vectors$PhiBarre
  loglik	= CalculEst$loglik
  Test		= CalculEst$Test
  Converge	= c(Converge,Test)

  UBarre  		= CalculEst$vectors$UBarre$UBarre
  PhiBarre		= PhiBarre[-c((3),(5),(k+5),(2*k+5))]
  result[(k-kmin)+1]	= sum(PhiBarre*UBarre)+.Calcul.constant(a,b,m,t,p,p.trans,p.eta,FixedDist=FixedDist)-.Calcul.constant(a=gamma[1],b=gamma[2],m=rep(normal[1],k-1),t=rep(normal[2],k-1),p=rep(1,2),p.trans=rep(1,4),p.eta=rep(1,k-1),FixedDist=FixedDist)-loglik
  }

  resultEM[[k]] = CalculEst
  }

  Converge	= Converge[-1]
  weight	= rep(0,length(Converge))
  result	= result[Converge]

  if(length(result)==0)
  {
  warnings("The VBEM algorithm has not converged, weights can not been estimated")
  cat("\n")
  stop()
  }                  

  res			= exp(-(result-min(result)))
  weight[Converge] 	= res/sum(res)
  TauAverage 		= .AverageEstimator(resultEM,weight,kmin,kmax)

  return(list(weight=weight,VBEM.estimation=resultEM,Converge=Converge,TauAverage=TauAverage))
  }
  



  
  
