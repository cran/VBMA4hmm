  ## __________________________________________________________
  ##
  ## .ForwardR
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Call to the C function which computes the forward loop	
  ## __________________________________________________________
  ##

  .ForwardR = function(Phi.vect,init,Mat.trans.norm.vect,ParamB.vect)
  {

  K		= length(init)
  nbParamB 	= K+5
  nbInd		= length(Phi.vect)/K

  .C("Forward",as.double(Phi.vect),as.double(init),as.double(Mat.trans.norm.vect),as.double(ParamB.vect),as.integer(nbInd),as.integer(K),as.integer(nbParamB),F=double(nbInd*K),Lambda=double(nbInd),LambdaParamB=double(nbInd), PACKAGE="VBMA4hmm")
  }


  ## __________________________________________________________
  ##
  ## .BackwardR
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Call to the C function which computes the backward loop	
  ## __________________________________________________________
  ##

  .BackwardR = function(F,Mat.transition,K)
  {
  nbInd = length(F)/K
  .C("Backward", as.double(F),as.double(Mat.transition),as.integer(nbInd),as.integer(K),tau=double(nbInd*K),G=double(nbInd*K), PACKAGE="VBMA4hmm")
  }



  ## __________________________________________________________
  ##
  ## .InitialisationHMM 
  ##
  ## Author : S. Volant
  ##
  ##		
  ## DESCRIPTION : Initialise the HMM model parameters
  ## __________________________________________________________
  ##

  .InitialisationHMM <- function(Tau)
  {

  nbInd	= dim(Tau)[1]
  K 	= dim(Tau)[2]

  ## Trasition Matrix
  Tau.tmp1 		= Tau[-nbInd,]
  Tau.tmp2 		= Tau[-1,]
  Mat.transition.tmp	= t(Tau.tmp1)%*%Tau.tmp2

  ## Normalise the transition matrix, it must be a stochastic matrix
  Mat.transition = Mat.transition.tmp/rowSums(Mat.transition.tmp)

  ## Initial distribution
  val.propre	= round(eigen(t(Mat.transition))$values,3)
  pos		= which(val.propre == 1.000)
  init		= eigen(t(Mat.transition))$vectors[,pos]
  init		= init / sum(init)
  init		= as.numeric(init)

  return(list(init=init,Mat.transition=Mat.transition))
  }




  ## __________________________________________________________
  ##
  ## .MatToVect 
  ##
  ## Author : S. Volant
  ##
  ##		
  ## DESCRIPTION : Transform a matrix into a vector
  ## __________________________________________________________
  ##

  .MatToVect <- function(Mat)
  {

  Mat	= as.matrix(Mat)
  nrow 	= dim(Mat)[1]
  ncol 	= dim(Mat)[2]
  vect	= Mat[,1]

  if(ncol>=2)
  {
  for( i in 2:ncol){vect = c(vect,Mat[,i])}
  }
  
  return(vect)
  }



  ## __________________________________________________________
  ##
  ## Name : .NuVectorHMM_hetero
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ##  Description : Function which allows one to calculate the vector Nu
  ##		    in the heterogeneous case
  ##
  ## __________________________________________________________


  .NuVectorHMM_hetero<-function(dir,dir.trans,gamma,normal,K)
  {

  t 		= normal[2]
  m 		= normal[1]
  a 		= gamma[1]
  b 		= gamma[2]
  p.eta 	= dir
  p.trans 	= dir.trans
  Nu1 		= dir.trans-1
  Nu1.trans 	= dir.trans-1
  Nu1.eta 	= dir-1
  Nu2 		= rep(-(t*m*m*0.5 +b),K-1) 
  Nu3 		= rep(t*m,K-1)
  Nu4 		= rep(-t*0.5,K-1)
  Nu5 		= rep(a-0.5,K-1)
  Nu		=c(Nu1,Nu2,Nu3,Nu4,Nu5,rep(Nu1.trans,2),Nu1.eta)

  return(list(Nu1=Nu1,Nu2=Nu2,Nu3=Nu3,Nu4=Nu4,Nu5=Nu5,Nu=Nu,Nu1.trans=rep(Nu1.trans,2),Nu1.eta=Nu1.eta))
  }


  ## __________________________________________________________
  ##
  ## Name : .NuVectorHMM_homo
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ##  Description : Function which allows one to calculate the vector Nu
  ##		    in the homogeneous case
  ##
  ## __________________________________________________________

  .NuVectorHMM_homo<-function(dir,dir.trans,gamma,normal,K)
  {

  ## Declaration and initialisation of local variables
  t		= normal[2]
  m 		= normal[1]
  a 		= gamma[1]
  b 		= gamma[2]
  p.eta 	= dir
  p.trans 	= dir.trans
  Nu1 		= dir.trans-1
  Nu1.trans 	= dir.trans-1
  Nu1.eta 	= dir-1
  Nu2		= -t*m*m*0.5 -b
  Nu3		= rep(t*m,K-1)
  Nu4		= rep(-t*0.5,K-1)
  Nu5		= a-1+0.5*(K-1)
  Nu		= c(Nu1,Nu2,Nu3,Nu4,Nu5,rep(Nu1.trans,2),Nu1.eta)

  return(list(Nu1=Nu1,Nu2=Nu2,Nu3=Nu3,Nu4=Nu4,Nu5=Nu5,Nu=Nu,Nu1.trans=rep(Nu1.trans,2),Nu1.eta=Nu1.eta))
  }



  ## __________________________________________________________
  ##
  ## Name : PhiBarreVector_hetero
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Calculate the expected value of vector Phi
  ##               in the heterogeneous case
  ##
  ## __________________________________________________________

  .PhiBarreVectorHMM_hetero<-function(Nu,K,FixedDist=NULL)
  {

  #################### Hyper Parameters #######################	
  Nu1		= Nu[[1]]
  Nu2		= Nu[[2]]
  Nu3		= Nu[[3]]
  Nu4		= Nu[[4]]
  Nu5		= Nu[[5]]
  Nu1.trans	= Nu[[7]]
  Nu1.eta	= Nu[[8]]
  p		= Nu1 +1
  p.trans	= matrix(Nu1.trans,ncol=2,nrow=2,byrow=TRUE) +1
  p.eta		= Nu1.eta+1
  a 		= Nu5 + 0.5
  b 		= -Nu2  + 0.25*Nu3*Nu3/Nu4 
  t 		= -2 * Nu4
  m 		= -0.5*Nu3/Nu4

  ## Expectations
  Row.sum.p.trans  	= rowSums(p.trans)
  pi.tmp1.trans    	= apply(as.matrix(Row.sum.p.trans),1,digamma)
  pi.tmp2.trans    	= apply(as.matrix(p.trans),2,digamma)
  tmp3.trans       	= pi.tmp2.trans - pi.tmp1.trans
  p.trans 		= .MatToVect.row(p.trans)
  LogPikBarre		= digamma(p)-digamma(sum(p))
  LogPikBarre.eta	= digamma(p.eta)-digamma(sum(p.eta))
  LogPikBarre.trans	= .MatToVect.row(tmp3.trans)
  LambdaBarre    	= (a/b)
  LambdaMuBarre  	= (m*a/b)
  LambdaMu2Barre 	= (m*m*(a/b))+1/t
  LogLambdaBarre 	= digamma(a)-log(b) 
  PhiBarre		= c(LogPikBarre,LambdaBarre,LambdaMuBarre,LambdaMu2Barre,LogLambdaBarre,LogPikBarre.trans,LogPikBarre.eta)
  
  return(list(LogPikBarre=LogPikBarre,LambdaBarre=LambdaBarre,LambdaMuBarre=LambdaMuBarre,LambdaMu2Barre=LambdaMu2Barre,LogLambdaBarre=LogLambdaBarre,PhiBarre=PhiBarre,
  a=a,b=b,t=t,m=m,p=p,p.trans=p.trans,p.eta=p.eta,LogPikBarre.mat=tmp3.trans,LogPikBarre.eta=LogPikBarre.eta,LogPikBarre.trans=LogPikBarre.trans))
  }



  ## __________________________________________________________
  ##
  ## Name : PhiBarreVector_homo
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Calculate the expected value of vector Phi
  ##               in the homogeneous case
  ##
  ## __________________________________________________________

  .PhiBarreVectorHMM_homo<-function(Nu,K,FixedDist=NULL)
  {

  #################### Hyper Parameters #######################
  Nu1		= Nu[[1]]
  Nu2		= Nu[[2]]
  Nu3		= Nu[[3]]
  Nu4		= Nu[[4]]
  Nu5		= Nu[[5]]
  Nu1.trans	= Nu[[7]]
  Nu1.eta	= Nu[[8]]
  p		= Nu1 +1
  p.trans	= matrix(Nu1.trans,ncol=2,nrow=2,byrow=TRUE) +1
  p.eta		= Nu1.eta+1
  a		= Nu5 +1-0.5*(K-1)
  b 		= -Nu2  + 0.25*sum(Nu3*Nu3/Nu4)
  t 		= -2 * Nu4
  m 		= -0.5*Nu3/Nu4

  ## Expectations
  Row.sum.p.trans	= rowSums(p.trans)
  pi.tmp1.trans		= apply(as.matrix(Row.sum.p.trans),1,digamma)
  pi.tmp2.trans		= apply(as.matrix(p.trans),2,digamma)
  tmp3.trans		= pi.tmp2.trans - pi.tmp1.trans
  p.trans 		= .MatToVect.row(p.trans)
  LogPikBarre		= digamma(p)-digamma(sum(p))
  LogPikBarre.eta	= digamma(p.eta)-digamma(sum(p.eta))
  LogPikBarre.trans	= .MatToVect.row(tmp3.trans)
  LambdaBarre    	= (a/b) 
  LambdaMuBarre  	= m*(a/b)
  LambdaMu2Barre 	= (m*m*(a/b))+1/t
  LogLambdaBarre 	= digamma(a)-log(b) 
  PhiBarre		= c(LogPikBarre,LambdaBarre,LambdaMuBarre,LambdaMu2Barre,LogLambdaBarre,LogPikBarre.trans,LogPikBarre.eta)
  
  return(list(LogPikBarre=LogPikBarre,LambdaBarre=LambdaBarre,LambdaMuBarre=LambdaMuBarre,LambdaMu2Barre=LambdaMu2Barre,LogLambdaBarre=LogLambdaBarre,PhiBarre=PhiBarre,
  a=a,b=b,t=t,m=m,p=p,p.trans=p.trans,p.eta=p.eta,LogPikBarre.mat=tmp3.trans,LogPikBarre.eta=LogPikBarre.eta,LogPikBarre.trans=LogPikBarre.trans))
  }



  ## __________________________________________________________
  ##
  ## Name : UBarreVector_hetero
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Calculate the expected value of vector Phi
  ##               in the heterogeneous case
  ##
  ## __________________________________________________________

  .UBarreVectorHMM_hetero<-function(Data,Tau,Mat.transition)
  {
  
  ## Declaration and initialisation of local variables
  Data2 = Data*Data
  nbInd = length(Data)
  K 	= dim(Mat.transition)[2]
  N 	= as.vector(colSums(Tau))[-1]
  S 	= as.vector(t(Tau) %*% Data)[-1]
  V 	= as.vector(t(Tau) %*% Data2)[-1]

  U1       	= c(Tau[1,1],(1-Tau[1,1]))
  U1.trans 	= c((nbInd)*Mat.transition[1,1], (nbInd)*sum(Mat.transition[1,2:K])
		  ,(nbInd)*sum(Mat.transition[(2:K),1]),(nbInd)*sum(Mat.transition[(2:K),(2:K)]))
  U1.eta   	= colSums((nbInd)*Mat.transition)[-1] + Tau[1,-1]
  U2 		= -0.5*V
  U3 		= S
  U4 		= -0.5*N
  U5 		= 0.5*N
  UBarre  	= c(U1,U2,U3,U4,U5,U1.trans,U1.eta)
  UBarregrp 	= UBarre
   
  return(list(U1=U1,U2=U2,U3=U3,U4=U4,U5=U5,UBarre=UBarre,UBarregrp=UBarregrp,U1.trans=U1.trans,U1.eta=U1.eta))
  }




  ## __________________________________________________________
  ##
  ## Name : UBarreVector_homo
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Calculate the expected value of vector Phi
  ##               in the heterogeneous case
  ##
  ## __________________________________________________________

  .UBarreVectorHMM_homo<-function(Data,Tau,Mat.transition)
  {
  
  ## Declaration and initialisation of local variables
  Data2 = Data*Data
  nbInd = length(Data)
  K 	= dim(Mat.transition)[2]
  N 	= as.vector(colSums(Tau))[-1]
  S 	= as.vector(t(Tau) %*% Data)[-1]
  V 	= as.vector(t(Tau) %*% Data2)[-1]
 
  U1		= c(Tau[1,1],(1-Tau[1,1]))
  U1.trans 	= c((nbInd)*Mat.transition[1,1], (nbInd)*sum(Mat.transition[1,2:K])
		  ,(nbInd)*sum(Mat.transition[(2:K),1]),(nbInd)*sum(Mat.transition[(2:K),(2:K)]))
  U1.eta   	= colSums((nbInd)*Mat.transition)[-1] + Tau[1,-1]
  U2 		= -0.5*sum(V)
  U3 		= S
  U4 		= -0.5*N
  U5 		= sum(N)*0.5
  UBarre	= c(U1,U2,U3,U4,U5,U1.trans,U1.eta)
  U2grp		= -0.5*V
  U5grp 	= N*0.5	
  UBarregrp  	= c(U1,U2grp,S,-0.5*N,U5grp,U1.trans,U1.eta)
   
  return(list(U1=U1,U2=U2,U3=U3,U4=U4,U5=U5,UBarre=UBarre,UBarregrp=UBarregrp,U1.trans=U1.trans,U1.eta=U1.eta))
  }


  ## __________________________________________________________
  ##
  ## Name : .UpdateParametersHMM
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Fonction to update the hyperparameters
  ##
  ## __________________________________________________________


  .UpdateParametersHMM<-function(Nu,UBarre)
  {

  ## Update
  Nu1       = Nu[[1]] + UBarre[[1]]
  Nu2       = Nu[[2]] + UBarre[[2]]
  Nu3       = Nu[[3]] + UBarre[[3]]
  Nu4       = Nu[[4]] + UBarre[[4]]
  Nu5       = Nu[[5]] + UBarre[[5]]
  Nu1.trans = Nu[[7]] + UBarre[[8]]
  Nu1.eta   = Nu[[8]] + UBarre[[9]]
  NuTilde   = Nu[[6]] + UBarre[[6]]

  return(list(Nu1=Nu1,Nu2=Nu2,Nu3=Nu3,Nu4=Nu4,Nu5=Nu5,NuTilde=NuTilde,Nu1.trans,Nu1.eta))
  }





  ## __________________________________________________________
  ##
  ## Name : .NonEmpty
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Restrict matrix values
  ##
  ## __________________________________________________________

  .NonEmpty <-function(Tau)
  {

  K = dim(Tau)[2]
  for(k in 1:K)
  {
  Tau[,k] = pmax(Tau[,k],1e-10)
  }
    
  return(Tau)
  }  
  


  ## __________________________________________________________
  ##
  ## Name : .converged
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Test the convergence
  ##
  ## __________________________________________________________


  .converged  <- function(OldParameters,NewParameters,Threshold=1e-6)
  {
  Converged = FALSE
  if(max(abs(OldParameters-NewParameters)) <= Threshold){Converged = TRUE}
  return(Converged)
  }



  ## __________________________________________________________
  ##
  ## Name : .Calcul.constant
  ##
  ## Author : Stevenn Volant
  ##
  ##
  ## Description : Normalising constant of the parameter distribution
  ##
  ## __________________________________________________________


  .Calcul.constant <- function(a,b,m,t,p,p.trans,p.eta,FixedDist=NULL)
  {
  K		= length(p.eta)+1
  p.trans.mat	= matrix(p.trans,nrow=2,ncol=2)
  logB		= sum(lgamma(p))-lgamma(sum(p))
  logB.trans	= apply(lgamma(p.trans.mat),2,sum)-lgamma(colSums(p.trans.mat))
  logB.eta	= sum(lgamma(p.eta))-lgamma(sum(p.eta))
  constante	= -sum(logB)-sum(logB.trans)-sum(logB.eta)-(K-1)*0.5*log(2*pi) +0.5*sum(log(t)) + sum(a*log(b)) - sum(lgamma(a))
  return(constante)
  }




  ## __________________________________________________________
  ##
  ## .ModifTransition
  ##
  ## Author : S. Volant
  ##
  ##	
  ## DESCRIPTION : Transform the transition matrix of Z
  ##
  ## __________________________________________________________
  ##


  .ModifTransition <- function(Omega)
  {

  K 			= dim(Omega)[2]
  Mat.transition 	= matrix(0,ncol=2,nrow=2)
  Omega.tilde    	= matrix(0,ncol=K,nrow=K)

  Mat.transition[1,1] = Omega[1,1]
  Mat.transition[1,2] = 1 - Omega[1,1]
  Mat.transition[2,1] = sum(Omega[2:K,1])/sum(Omega[2:K,])
  Mat.transition[2,2] = 1 - Mat.transition[2,1] 

  eta.tmp = colSums(Omega)[-1]
  eta     = eta.tmp/sum(eta.tmp)

  Omega.tilde[1,1]	= Mat.transition[1,1]
  Omega.tilde[1,2:K]	= Mat.transition[1,2]*eta
  Omega.tilde[2:K,1]	= rep(Mat.transition[2,1],K-1)
  Omega.tilde[2:K,2:K]	= Mat.transition[2,2]*matrix(eta,ncol=K-1,nrow=K-1,byrow=TRUE)

  return(list(Mat.transition=Mat.transition,eta=eta,Omega.tilde=Omega.tilde))
  }

 
 
  ## __________________________________________________________
  ##
  ## .ModifTransition.2clusters
  ##
  ## Author : S. Volant
  ##
  ##	
  ## DESCRIPTION : From S to Z (Transition matrix)
  ##
  ## __________________________________________________________
  ##
 
 
  .ModifTransition.2clusters <- function(Mat.transition,eta)
  {

  K 		= length(eta)+1
  Omega.tilde	= matrix(0,ncol=K,nrow=K)

  Omega.tilde[1,1] 	= Mat.transition[1,1]
  Omega.tilde[1,2:K] 	= Mat.transition[1,2]*eta
  Omega.tilde[2:K,1] 	= rep(Mat.transition[2,1],K-1)
  Omega.tilde[2:K,2:K] 	= Mat.transition[2,2]*matrix(eta,ncol=K-1,nrow=K-1,byrow=TRUE)

  return(list(Omega.tilde=Omega.tilde))
  }

  


  ## __________________________________________________________
  ##
  ## .AverageEstimator
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Compute the averaged estimator of tau.
  ## __________________________________________________________
  ##


  .AverageEstimator<-function(VBEM.object,Weight,kmin,kmax)
  {

  ## Declaration and initialisation of local variables
  nbInd          = length(VBEM.object[[kmin]][[1]][,1])
  TauAverage     = matrix(0,ncol=(kmax-kmin+1),nrow=nbInd)

  ## Gathering the estimation of Tau0 weighting by their posterior weight
  for(k in kmin:kmax)
  {
  TauAverage[,(k-kmin+1)] = as.vector(VBEM.object[[k]][[1]][,1]*Weight[(k-kmin+1)])
  }
  
  TauAverage  = rowSums(TauAverage)
  
  return(TauAverage)
  }



  ## __________________________________________________________
  ##
  ## .MatToVect.row 
  ##
  ## Author : S. Volant
  ##
  ##		
  ## DESCRIPTION : Transform a matrix into a vector
  ## __________________________________________________________
  ##


  .MatToVect.row <- function(Mat)
  {
  
  Mat	= as.matrix(Mat)
  nrow	= dim(Mat)[1]
  ncol	= dim(Mat)[2]
  vect	= Mat[1,]

  if(nrow>=2)
  {
    for(i in 2:nrow)
    {
    vect = c(vect,Mat[i,]) 
    }
  }
  return(vect)
  
  }


  ## __________________________________________________________
  ##
  ## .InitialKmeans 
  ##
  ## Author : S. Volant
  ##
  ##		
  ## DESCRIPTION : Initialise the EM algorithm by the kmeans 
  ## __________________________________________________________
  ##

  .InitialKmeans <- function(data,K,FixedDist.m,FixedDist.var)
  {
  rk	= kmeans(data,K,nstart=25)
  Tau 	= matrix(0,ncol=K,nrow=length(data))

  for(k in 1:K){ Tau[,k] = dnorm(data,rk$centers[k],sd=sd(data[rk$cluster==k]))}

  ind		= which.min(sapply(1:K, FUN=function(k) sum(abs(rk$centers[k] - FixedDist.m))))
  tmp 		= Tau[,ind]
  Tau[,ind]	= Tau[,1]
  Tau[,1] 	= dnorm(data,FixedDist.m,FixedDist.var)
  Tau		= Tau/rowSums(Tau)
  return(Tau)
  }




 ## __________________________________________________________
  ##
  ## .E.step.HMM
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Computes Forward Backward algorithm from C functions
  ## __________________________________________________________
  ##

  .E.step.HMM <- function(Phi,init,Mat.transition,ParamB.vect)
  {

  nbClasse = dim(Phi)[2]
  nbInd    = dim(Phi)[1]

  ## Transform matrix to vectors
  Phi.vect            = .MatToVect(Phi)
  Mat.transition.vect = .MatToVect(Mat.transition)

  ## Forward loop
  resF = .ForwardR(Phi.vect,init,Mat.transition.vect,ParamB.vect)

  ## Backward loop
  resB = .BackwardR(resF$F,Mat.transition.vect,nbClasse)

  ## Turn vector to matrix
  Tau   = matrix(resB$tau,nrow=nbInd,ncol=nbClasse,byrow=FALSE)
  F     = matrix(c(resF$F),nrow=nbInd,ncol=nbClasse)
  G     = matrix(c(resB$G),nrow=nbInd,ncol=nbClasse)

  ## Log-likelihhod
  LambdaParamB 	= resF$LambdaParamB
  loglik 	= -sum(log(LambdaParamB))
  return(list(F=F,G=G,Tau=Tau,loglik=loglik))
  }



  ## __________________________________________________________
  ##
  ## .M.step.HMM  
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Allows one to update the HMM model parameters
  ## __________________________________________________________
  ##

  .M.step.HMM <- function(Tau,F,G,Mat.transition)
  {

  nbInd    = dim(Tau)[1]
  nbClasse = dim(Tau)[2]

  ## Tansition matrix
  Mat.transition.tmp = Mat.transition * (t(F[-nbInd,])%*% (Tau[-1,] / G[-1,]))

  ## Normalise the transition matrix, it must be a stochastic matrix
  Mat.transition = pmax(Mat.transition,0.0001)
  Mat.transition = pmin(Mat.transition,1)
  Mat.transition = Mat.transition.tmp/pmax(rowSums(Mat.transition.tmp),0.0001)
  Mat.transition = t(as.matrix(apply(Mat.transition,1,pmax,0)))
  Mat.transition = t(as.matrix(apply(Mat.transition,1,pmin,1)))

  ## Initial distribution
  val.propre = round(eigen(t(Mat.transition))$values,3)
  pos        = which(val.propre == 1.000) # looking for proper value equals to 1
  init.tmp  = eigen(t(Mat.transition))$vectors[,pos]

  ## Normalise the initial distribution
  init = init.tmp / sum(init.tmp)
  init = as.numeric(init)

  return(list(init=init,Mat.transition=Mat.transition))
  }





  ## __________________________________________________________
  ##
  ## .EM.HMM
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Compute Forward Backward algorithm from C functions and give the update parameters
  ## __________________________________________________________
  ##

  .EM.HMM <- function(Phi,init,Mat.transition,ParamB.vect)
  {

  ## E step
  Estep = .E.step.HMM(Phi,init,Mat.transition,ParamB.vect)

  ## Gathering of E step results
  Tau		= Estep$Tau
  loglik	= Estep$loglik
  Tau		= Tau/rowSums(Tau)
  F		= Estep$F
  G		= Estep$G

  ## M step
  Mstep = .M.step.HMM(Tau,F,G,Mat.transition)

  ## Gathering of M step uptadating
  Mat.transition = Mstep$Mat.transition
  init		 = Mstep$init

  return(list(Tau=Tau,F=F,G=G,init=init,Mat.transition=Mat.transition,loglik=loglik,loglik=loglik))
  }





  ## __________________________________________________________
  ##
  ## .VBHMM_hetero
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Compute the VBEM algorithm with heterogeneous Gaussian
  ## __________________________________________________________
  ##



  .VBHMM_hetero <-function(data,dir=c(1,1),dir.trans=c(1,1),gamma=c(0.00001,0.00001),normal=c(0,0.0001),FixedDist=c(0,1),IterationMax=1000,Threshold=1e-3)
  {

  ## Check
  stopifnot(is.vector(data))
  stopifnot(min(gamma)>0)
  stopifnot(normal[2]>0)
  stopifnot(FixedDist[2]>0)
  stopifnot(length(FixedDist)==2)

  ## Declaration and initialisation of local variables
  K		= length(dir)+1
  nbInd		= length(data)
  data2		= data * data
  Phi		= matrix(0,ncol=K,nrow=nbInd)
  OldParameters = rep(0,((5*K)+1))
     
  ## Initialise the algorithm and define the parameters #####
  NuList 	= .NuVectorHMM_hetero(dir,dir.trans,gamma,normal,K)
  Tau 		= .InitialKmeans(data,K,FixedDist[1],FixedDist[2])
  Init.tau	= .InitialisationHMM(Tau)
  init.omega	= Init.tau$init
  Omega.tmp	= Init.tau$Mat.transition
  tmp		= .ModifTransition(Omega.tmp)

  Mat.transition = tmp$Mat.transition
  Omega		 = tmp$Omega.tilde
  eta		 = tmp$eta


  ######################## EM-algorithm ########################
  for(i in 1:IterationMax)
  {
  UBarre	= .UBarreVectorHMM_hetero(data,Tau,Omega)
  NuTilde	= .UpdateParametersHMM(NuList,UBarre)
  PhiBarreList	= .PhiBarreVectorHMM_hetero(NuTilde,K,FixedDist)
  a		= PhiBarreList$a
  b		= PhiBarreList$b
  m		= PhiBarreList$m
  p		= PhiBarreList$p
  p.trans	= PhiBarreList$p.trans
  p.eta		= PhiBarreList$p.eta
  t		= PhiBarreList$t

  ## Keep the variable in the same order along the loop (order : increasing mean)   
  ordre = order(m)
  m     = m[ordre]
  p.eta = p.eta[ordre]
  t     = t[ordre]
  a     = a[ordre]
  b     = b[ordre]

  ## Gathering of the estimated parameters and the fixed distribution into PhiBarre
  LogPikBarre        	= PhiBarreList$LogPikBarre
  LogPikBarre.trans  	= PhiBarreList$LogPikBarre.trans
  LogPikBarre.eta    	= PhiBarreList$LogPikBarre.eta[ordre]
  LambdaBarre   	= c(1/FixedDist[2],PhiBarreList$LambdaBarre[ordre]) 
  LambdaMuBarre  	= c(FixedDist[1]/FixedDist[2],PhiBarreList$LambdaMuBarre[ordre])
  LambdaMu2Barre 	= c(FixedDist[1]*FixedDist[1]/FixedDist[2],PhiBarreList$LambdaMu2Barre[ordre])
  LogLambdaBarre 	= c(-log(FixedDist[2]),PhiBarreList$LogLambdaBarre[ordre] )

  ## Specific form of PhiBarre for estimating the log likelihood. 
  PhiBarre=c(LogPikBarre[1],LogPikBarre[2]+LogPikBarre.eta,LogPikBarre.trans[1],(LogPikBarre.trans[2]+LogPikBarre.eta),rep(c(LogPikBarre.trans[3],(LogPikBarre.trans[4]		+LogPikBarre.eta)),K-1),LambdaBarre,LambdaMuBarre,LambdaMu2Barre,LogLambdaBarre )

  ## Gathering all parameters into a vector
  NewParameters = c(a,b,m,t,as.vector(p),p.trans,p.eta)

  ## Phi value for each point and for all the estimated distributions
  Phi[,2:K]	= apply(as.matrix(2:K), 1, FUN = function(k) dnorm(data,mean=m[k-1],sd=sqrt(b[k-1]/a[k-1])))
  Phi[,1]	= dnorm(data,mean=FixedDist[1],sd=sqrt(FixedDist[2]))

  ## E step
  EM		= .EM.HMM(Phi,init.omega,Omega,PhiBarre)
  Tau		= EM$Tau
  Tau		= Tau/rowSums(Tau)
  Omega.tmp	= EM$Mat.transition

  ## Shift the transition matrix
  tmp = .ModifTransition(Omega.tmp)

  ## New Transition matrix for S and Z
  Mat.transition = tmp$Mat.transition
  Omega		 = tmp$Omega.tilde
  eta		 = tmp$eta

  ## Initial distribution
  Mat.transition = pmax(Mat.transition,0.0001)
  Mat.transition = pmin(Mat.transition,1)
  eig = eigen(t(Mat.transition))
  val.propre = round(eig$values,3)
  pos        = which(val.propre == 1.000) # looking for proper value equals to 1
  init.tmp  = eig$vectors[,pos]

  ## Normalise the initial distribution
  init = init.tmp / sum(init.tmp)
  init = as.numeric(init)
  init.omega = c(init[1],init[-1]*eta)

  ## About the loglikelihood from the parameters Theta and PhiBarre
  loglik         = EM$loglik

  ## PhiBarre in the right form
  PhiBarre=c(LogPikBarre,LambdaBarre,LambdaMuBarre,LambdaMu2Barre,LogLambdaBarre,LogPikBarre.trans,LogPikBarre.eta )
  
  ### Too small values of Tau are not allowed
  Tau = .NonEmpty(Tau)
  
  ### Testing the algorithm convergence
  Test = .converged(OldParameters,NewParameters,Threshold=Threshold)
  if(Test){break()}

  ## Keep the old parameters	
  OldParameters = NewParameters
  }

  Variance 		= b/a
  parameters 		= list(prop = p/sum(p) , mean = c(FixedDist[1], m), sd = sqrt(c(1/FixedDist[2],Variance)),init=init,eta=eta,init.omega=init.omega)
  hyperparameters	= list(a=a,b=b,t=t,m=m,p=p,p.eta=p.eta,p.trans=p.trans)
  vectors 		= list(NuTilde=NuTilde,PhiBarre=PhiBarre,UBarre=UBarre)
  transition	 	= list(Mat.transition=Mat.transition,Omega=Omega)

  return(list(Tau=Tau,parameters=parameters,hyperparameters=hyperparameters,transition=transition,vectors=vectors,loglik=loglik,Test=Test,Phi=Phi))
  }

 




  ## __________________________________________________________
  ##
  ## .VBHMM_homo
  ##
  ## Author : S. Volant
  ##
  ##
  ## DESCRIPTION : Compute the VBEM algorithm with homogeneous Gaussian
  ## __________________________________________________________
  ##




  .VBHMM_homo <-function(data,dir=c(1,1),dir.trans=c(1,1),gamma=c(0.01,0.01),normal=c(-4,0.01),FixedDist=c(0,1),IterationMax=1000,Threshold=1e-3)
  {

  ## Check
  stopifnot(is.vector(data))
  stopifnot(min(gamma)>0)
  stopifnot(normal[2]>0)
  stopifnot(FixedDist[2]>0)
  stopifnot(length(FixedDist)==2)



  ## Declaration and initialisation of local variables
  K		=length(dir)+1
  nbInd		= length(data)
  data2		= data * data
  Phi		= matrix(0,ncol=K,nrow=nbInd)
  OldParameters = rep(0,((3*K)+5))
     
  ## Initialise the algorithm and define the parameters #####
  NuList 		= .NuVectorHMM_homo(dir,dir.trans,gamma,normal,K)
  Tau 			= .InitialKmeans(data,K,FixedDist[1],FixedDist[2])
  Init.tau 		= .InitialisationHMM(Tau)
  init.omega		= Init.tau$init
  Omega.tmp		= Init.tau$Mat.transition
  tmp			= .ModifTransition(Omega.tmp)
  Mat.transition	= tmp$Mat.transition
  Omega		 	= tmp$Omega.tilde
  eta		 	= tmp$eta


  ######################## EM-algorithm ########################
  for(i in 1:IterationMax)
  {
  UBarre	= .UBarreVectorHMM_homo(data,Tau,Omega)
  NuTilde	= .UpdateParametersHMM(NuList,UBarre)
  PhiBarreList	= .PhiBarreVectorHMM_homo(NuTilde,K,FixedDist)
  a		= PhiBarreList$a
  b		= PhiBarreList$b
  m		= PhiBarreList$m
  p		= PhiBarreList$p
  p.trans	= PhiBarreList$p.trans
  p.eta		= PhiBarreList$p.eta
  t		= PhiBarreList$t

  ## Keep the variable in the same order along the loop (order : increasing mean)   
  ordre = order(m)
  m     = m[ordre]
  p.eta = p.eta[ordre]
  t     = t[ordre]

  ## Gathering of the estimated parameters and the fixed distribution into PhiBarre
  LogPikBarre        	= PhiBarreList$LogPikBarre
  LogPikBarre.trans  	= PhiBarreList$LogPikBarre.trans
  LogPikBarre.eta    	= PhiBarreList$LogPikBarre.eta[ordre]
  LambdaBarre   	= c(1/FixedDist[2],rep(PhiBarreList$LambdaBarre,K-1)) 
  LambdaMuBarre  	= c(FixedDist[1]/FixedDist[2],PhiBarreList$LambdaMuBarre[ordre])
  LambdaMu2Barre 	= c(FixedDist[1]*FixedDist[1]/FixedDist[2],PhiBarreList$LambdaMu2Barre[ordre])
  LogLambdaBarre 	= c(-log(FixedDist[2]),rep(PhiBarreList$LogLambdaBarre,K-1) )

  ## Specific form of PhiBarre for estimating the log likelihood. 
  PhiBarre = c(LogPikBarre[1],LogPikBarre[2]+LogPikBarre.eta,LogPikBarre.trans[1],(LogPikBarre.trans[2]+LogPikBarre.eta),rep(c(LogPikBarre.trans[3],(LogPikBarre.trans[4]	 +LogPikBarre.eta)),K-1),LambdaBarre,LambdaMuBarre,LambdaMu2Barre,LogLambdaBarre )

  ## Gathering all parameters into a vector
  NewParameters = c(a,b,m,t,as.vector(p),p.trans,p.eta)

  ## Phi value for each point and for all the estimated distributions
  for( k in 2:K)
  {
  Phi[,k] = dnorm(data,mean=m[k-1],sd=sqrt(b/a))
  }
  Phi[,1] = dnorm(data,mean=FixedDist[1],sd=sqrt(FixedDist[2]))

  ## E step
  EM		= .EM.HMM(Phi,init.omega,Omega,PhiBarre)
  Tau           = EM$Tau
  Tau           = Tau/rowSums(Tau)
  Omega.tmp	= EM$Mat.transition

  ## Shift the transition matrix
  tmp = .ModifTransition(Omega.tmp)

  ## New Transition matrix for S and Z
  Mat.transition = tmp$Mat.transition
  Omega		 = tmp$Omega.tilde
  eta		 = tmp$eta

  ## Initial distribution
  val.propre = round(eigen(t(Mat.transition))$values,3)
  pos        = which(val.propre == 1.000) # looking for proper value equals to 1
  init.tmp   = eigen(t(Mat.transition))$vectors[,pos]

  ## Normalise the initial distribution
  init		= init.tmp / sum(init.tmp)  
  init		= as.numeric(init)
  init.omega	= c(init[1],init[-1]*eta)

  ## About the loglikelihood from the parameters Theta and PhiBarre
  loglik	= EM$loglik

  ## PhiBarre in the right form
  PhiBarre = c(LogPikBarre,LambdaBarre[c(1,2)],LambdaMuBarre,LambdaMu2Barre,LogLambdaBarre[c(1,2)],LogPikBarre.trans,LogPikBarre.eta )
  
  ## Too small values of Tau is not allowed
  Tau = .NonEmpty(Tau)
 
  ## Testing the algorithm convergence
  Test = .converged(OldParameters,NewParameters,Threshold=Threshold)
  if(Test) break()
 
  ## Keep the old parameters
  OldParameters = NewParameters
 
  }

  Variance 	= b/a
  parameters 	= list(prop = p/sum(p) , mean = c(FixedDist[1], m), sd = sqrt(c(1/FixedDist[2],Variance)),init=init,eta=eta,init.omega=init.omega)
  hyperparameters = list(a=a,b=b,t=t,m=m,p=p,p.eta=p.eta,p.trans=p.trans)
  vectors 	= list(NuTilde=NuTilde,PhiBarre=PhiBarre,UBarre=UBarre,NuList=NuList)
  transition 	= list(Mat.transition=Mat.transition,Omega=Omega)

  return(list(Tau=Tau,parameters=parameters,hyperparameters=hyperparameters,transition=transition,vectors=vectors,loglik=loglik,Test=Test,Phi=Phi))
  }

 

