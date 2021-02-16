#' Mixture likelihood function
#'
#' This function compute the log-likelihood of a mixture model for 
#' a fixed number of Gaussian components.
#' @param y vector of observations.
#' @param w vector of weights; they have to sum to one.
#' @param mu vector of location parameters for the Gaussian components.
#' @param sigma vector of standard deviation parameters for the Gaussian components.
#' @return a scalar representing the log-likelihood of the mixture of Gaussian components.
#' @keywords BayesMIC
#' @export
log_lkd<-function(y,w,mu,sigma) {
  tot=0
  for (i in 1:length(y)) {
    tot=tot+log(sum(w*dnorm(y[i],mean=mu,sd=sigma)))
  }
  return(tot)
}

#' Mixture likelihood function with data augmentation
#'
#' This function compute the log-likelihood of a mixture model for 
#' a fixed number of Gaussian components with data augmentation, i.e.
#' including a latent variable representing the allocation.
#' @param y vector of observations.
#' @param w vector of weights; they have to sum to one.
#' @param mu vector of location parameters for the Gaussian components.
#' @param sigma vector of standard deviation parameters for the Gaussian components.
#' @param S vector of length equal to the vector y; values in 1,...,K, the total number of components. 
#' @return a scalar representing the log-likelihood of the mixture of Gaussian components.
#' @keywords BayesMIC
#' @export
log_lkd_da <- function(y,w,mu,sigma,S)
{
  n <- length(y)
  K <- length(w)
  tot <- 0
  for(i in 1:n){
    for(k in 1:K)
    {
      tot <- tot + (S[i]==k) * ( log(w[k]) + log(dnorm(y[i],mu[k],sigma[k]) ) )
    }		
  }
  return(tot)
}


#' Log-prior density, with a beta-negative binomial prior distribution representing the loss-based prior.
#'
#' This function compute the log-joint prior for the parameters of a mixture model with an  
#' unknown number of Gaussian components.
#' @param w vector of weights; they have to sum to one.
#' @param mu vector of location parameters for the Gaussian components.
#' @param sigma vector of standard deviation parameters for the Gaussian components.
#' @param alpha_bnb shape parameter of the beta-negative-binomial. Default value at 1.
#' @param beta_bnb rate parameter of the beta-negative-binomial. Default value at 1.
#' @param sh shape parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 1.5.
#' @param rt rate parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 0.5.
#' @param alpha_dir parameter of the Dirichlet prior distribution for the weight parameters of the Gaussian components. Default value at 1.
#' @param me_mu mean parameter of the Gaussian prior distribution for the mean parameters of the Gaussian components. Default value at 20.
#' @param sd_mu standard deviation parameter of the Gaussian prior distribution for the standard deviation parameters of the Gaussian components. Default value at 10.
#' @return a scalar representing the log-joint prior density of the parameters of a mixture of Gaussian components.
#' @keywords BayesMIC
#' @importFrom  extraDistr dbnbinom
#' @importFrom mixtools ddirichlet
#' @export
log_prior_bnb<-function(w,mu,sigma,alpha_bnb=1,beta_bnb=1,
                        sh=1.5,rt=0.5,alpha_dir=1,me_mu=20,sd_mu=10) {
  k <- length(w)
  mdp <- dbnbinom(k, 1, alpha = alpha_bnb, beta = beta_bnb, log = TRUE) # beta negative binomial for the number of components
  sp=sum(dgamma(sigma,shape=sh,rate=rt,log=TRUE))   #mean 3, rule out very dense clusters, sd about 2.5
  if(length(w) > 1){
    wp=log(ddirichlet(w,rep(alpha_dir,length(w))))         #cluster weight sum to one else uniform
  } else {
    wp=0}
  mup=sum(dnorm(mu,mean=me_mu,sd=sd_mu,log=TRUE))          #mu is 0-40 at 2 sigma
  return(sum(mup+sp+wp+mdp))
}

#' Log-prior density, with a uniform prior distribution representing the loss-based prior.
#'
#' This function compute the log-joint prior for the parameters of a mixture model with an  
#' unknown number of Gaussian components.
#' @param w vector of weights; they have to sum to one.
#' @param mu vector of location parameters for the Gaussian components.
#' @param sigma vector of standard deviation parameters for the Gaussian components.
#' @param bn_lw lower bound of the uniform prior distribution for the number of Gaussian components. Default value at 0.
#' @param bn_up upper bound of the uniform prior distribution for the number of Gaussian components. Default value at 30.
#' @param sh shape parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 1.5.
#' @param rt rate parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 0.5.
#' @param alpha_dir parameter of the Dirichlet prior distribution for the weight parameters of the Gaussian components. Default value at 1.
#' @param me_mu mean parameter of the Gaussian prior distribution for the mean parameters of the Gaussian components. Default value at 20.
#' @param sd_mu standard deviation parameter of the Gaussian prior distribution for the standard deviation parameters of the Gaussian components. Default value at 10.
#' @return a scalar representing the log-joint prior density of the parameters of a mixture of Gaussian components.
#' @keywords BayesMIC
#' @importFrom mixtools ddirichlet
#' @export
log_prior_unif<-function(w,mu,sigma,bn_lw=0,bn_up=30,sh=1.5,rt=0.5,
                         alpha_dir=1,me_mu=20,sd_mu=10) {
  k <- length(w)
  mdp <- dunif(k, bn_lw,bn_up, log = TRUE) # beta negative binomial for the number of components
  sp=sum(dgamma(sigma,shape=sh,rate=rt,log=TRUE))   #mean 3, rule out very dense clusters, sd about 2.5
  if(length(w) > 1){
    wp=log(ddirichlet(w,rep(alpha_dir,length(w))))         #cluster weight sum to one else uniform
  } else {
    wp=0}
  mup=sum(dnorm(mu,mean=me_mu,sd=sd_mu,log=TRUE))          #mu is 0-40 at 2 sigma
  return(sum(mup+sp+wp+mdp))
}

#' Log-prior density, with a Poisson prior distribution representing the loss-based prior.
#'
#' This function compute the log-joint prior for the parameters of a mixture model with an  
#' unknown number of Gaussian components.
#' @param w vector of weights; they have to sum to one.
#' @param mu vector of location parameters for the Gaussian components.
#' @param sigma vector of standard deviation parameters for the Gaussian components.
#' @param rt_pois rate parameter of the Poisson prior distribution for the number of Gaussian components. Default value at 1.
#' @param sh shape parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 1.5.
#' @param rt rate parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 0.5.
#' @param alpha_dir parameter of the Dirichlet prior distribution for the weight parameters of the Gaussian components. Default value at 1.
#' @param me_mu mean parameter of the Gaussian prior distribution for the mean parameters of the Gaussian components. Default value at 20.
#' @param sd_mu standard deviation parameter of the Gaussian prior distribution for the standard deviation parameters of the Gaussian components. Default value at 10.
#' @return a scalar representing the log-joint prior density of the parameters of a mixture of Gaussian components.
#' @keywords BayesMIC
#' @importFrom mixtools ddirichlet
#' @export
log_prior_pois<-function(w,mu,sigma,rt_pois=1,sh=1.5,rt=0.5,
                         alpha_dir=1,me_mu=20,sd_mu=10) {
  k <- length(w)
  mdp <- dpois(k, rt_pois, log = TRUE) # beta negative binomial for the number of components
  sp=sum(dgamma(sigma,shape=sh,rate=rt,log=TRUE))   #mean 3, rule out very dense clusters, sd about 2.5
  if(length(w) > 1){
    wp=log(ddirichlet(w,rep(alpha_dir,length(w))))         #cluster weight sum to one else uniform
  } else {
    wp=0}
  mup=sum(dnorm(mu,mean=me_mu,sd=sd_mu,log=TRUE))          #mu is 0-40 at 2 sigma
  return(sum(mup+sp+wp+mdp))
}

#' Simulations from the joint posterior distribution of the parameters of a 
#' mixture of Gaussian components with unknown number of components, by using a beta-negative-binomial prior for k. 
#'
#' @param y vector of observations.
#' @param K number of simulations. Default value at 10^5. 
#' @param alpha_bnb shape parameter of the beta-negative-binomial. Default value at 1.
#' @param beta_bnb rate parameter of the beta-negative-binomial. Default value at 1.
#' @param sh shape parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 1.5.
#' @param rt rate parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 0.5.
#' @param alpha_dir parameter of the Dirichlet prior distribution for the weight parameters of the Gaussian components. Default value at 1.
#' @param me_mu mean parameter of the Gaussian prior distribution for the mean parameters of the Gaussian components. Default value at 20.
#' @param sd_mu standard deviation parameter of the Gaussian prior distribution for the standard deviation parameters of the Gaussian components. Default value at 10.
#' @return k, a sample from the posterior distribution of the number of components.
#' @return param, a matrix from the posterior distribution of the weights, the means and the standard deviations of the components.
#' @return olp, a vector a joint prior density of the posterior samples.
#' @return oll, a vector a likelihood values of the posterior samples.
#' @keywords BayesMIC
#' @export
post_bnb <- function(y,K=10^5,alpha_bnb=1,beta_bnb=1,
                     sh=1.5,rt=0.5,alpha_dir=1,me_mu=20,sd_mu=10){

  n=length(y)
  
  #initialise MCMC with one component
  mu=mean(y)
  sigma=sd(y)
  w=1
  
  Sdraw <- rep(1,n)
  G <- 1
  nh <- n
  
  u 		<- c()
  ystar 	<- rep(Inf,n)
  while(sum(abs(ystar)==Inf)!=0)
  {
    for(l in 1:n)
    {
      treshb   =  pnorm(y[l]+1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
      tresha   =  pnorm(y[l]-1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
      
      if(tresha != treshb)
      {
        u[l] = runif(1,tresha , treshb)	
        ystar[l] <- qnorm(u[l],mu[Sdraw[l]], sigma[Sdraw[l]])
      } else {
        ystar[l] <- y[l]
      }
      
    }	# end of the observations loop
  }
  
  oll=log_lkd_da(ystar,w,mu,sigma,Sdraw)
  olp=log_prior_bnb(w,mu,sigma)
  
  #MCMC
  SS=10; Nsamp=K/SS;  
  PL=matrix(NA,Nsamp,3); colnames(PL)<-c("d","olp","oll"); TH=list();
  
  for (k in 1:K) {
    OK=TRUE
    d=length(w)
    if(d==1){
      move=2
      while(move==2 )
        move=sample(1:5,1) #choose a MCMC move (2 birth death, 3 fixed dimension)
    } else {
      move=sample(1:5,1)
    }

    if (move==1) {
      #add a cluster
      i=sample(1:d,1) #pick a component to split - the weight w[i] is shared 
      wn=runif(1,min=0,max=w[i]); wp=c(w,wn); wp[i]=w[i]-wn
      mun=runif(1,min(y)-5,max(y)+5); mup=c(mu,mun)
      sn=rgamma(1,shape=1.5,rate=0.5); sigmap=c(sigma,sn)
      qd=dunif(mun,min(y)-5,max(y)+5,log=TRUE)+dgamma(sn,shape=1.5,rate=0.5,log=TRUE)-log(w[i])
      qn=0
      rhod=-log(d)
      rhon=-log((d+1)*d)
    }
    if (move==2) {
      #kill a cluster
      if (d>1) {
        ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
        #pick component to kill and a component to increment with the weight from 
        #the killed cluster
        wp=w; wp[j]=w[j]+w[i]; wp=wp[-i]
        mup=mu[-i]; sigmap=sigma[-i]
        qd=0
        qn=dunif(mu[i],min(y)-5,max(y)+5,log=TRUE)+dgamma(sigma[i],shape=1.5,rate=0.5,log=TRUE)-log(w[i]+w[j])
        rhod=-log(d*(d-1))
        rhon=-log(d-1)
      } else {
        OK=FALSE
      }
    }
    if (move==3) {
      #fixed dimension mu - simple RW MH
      i=sample(1:d,1); 
      mup=mu; mup[i]=rnorm(1,mean=mu[i],sd=1);
      sigmap=sigma; wp=w;
      qd=qn=rhod=rhon=0      
    }
    if (move==4) {
      #fixed dimension sigma - scaling move
      i=sample(1:d,1); 
      sigmap=sigma; u=runif(1,0.5,2); sigmap[i]=sigma[i]*u
      mup=mu; wp=w;
      qn=-log(u); qd=rhod=rhon=0      
    }
    if (move==5) {
      #fixed dimension w - redistribute the weights between two clusters chosen UAR 
      if (d>1) {
        ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
        wp=w; wp[i]=runif(1,min=0,max=w[i]+w[j]); wp[j]=w[i]+w[j]-wp[i]
        mup=mu; sigmap=sigma;
        qn=qd=rhod=rhon=0
      } else {
        OK=FALSE
      }
    }
    
    # Relabeling
    
    temp_mu <- mup
    mup <- sort(mup)
    wp <- wp[order(temp_mu)]
    sigmap <- sigmap[order(temp_mu)]
    
    do_DA <- FALSE
    if (OK) {
      nll=log_lkd(y,wp,mup,sigmap)
      nlp=log_prior_bnb(wp,mup,sigmap)  
      MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
      if (log(runif(1))<MHR) {
        do_DA <- TRUE
      }
    }
    
    ### Step 2: Multinomial sampler
    
    if(do_DA == TRUE){
      G <- length(wp)
      Sdraw <- c()
      for(l in 1:n)
      {
        probs <- c()
        for(g in 1:G)
        {
          probs[g] <-  dnorm(ystar[l], mup[g], sigmap[g])
        }
        probs[probs<0] <- 0
        if(sum(probs)!=0){
          Sdraw[l] <- sample(1:G,size=1,prob=probs)
        } else {
          Sdraw[l] <- Sdraw[l]
        }
      }
    
      
      ### Step 3: Data augmentation
      
      u 		<- c()
      ystar 	<- c()
      for(l in 1:n)
      {
        treshb   =  pnorm(y[l]+1, mup[Sdraw[l]],sigmap[Sdraw[l]]) 
        tresha   =  pnorm(y[l], mup[Sdraw[l]],sigmap[Sdraw[l]]) 
        u[l] = runif(1,tresha , treshb)
        ystar[l] <- qnorm(u[l],mup[Sdraw[l]], sigmap[Sdraw[l]])
        ystar[l] <- ifelse(abs(ystar[l])==Inf, y[l], ystar[l])
        
      }	# end of the observations loop
      
      # Acceptance / Rejection
      if (OK) {
        nll=log_lkd_da(ystar,wp,mup,sigmap,Sdraw)
        if(is.na(nll)==T){nll = -Inf}
        nlp=log_prior_bnb(wp,mup,sigmap)  
        MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
        if (log(runif(1))<MHR) {
          w=wp; mu=mup; sigma=sigmap
          oll=nll; olp=nlp
        }
      }
    }		
    
    
    # Saving
    if (k%%SS==0) {
      TH[[k/SS]]=list(w=w,mu=mu,sigma=sigma)
      PL[k/SS,]=c(length(w),olp,oll)
    }
    
  }  
  return(list(k=PL[1],param=TH,olp=PL[,2],oll=PL[,3]))
}

#' Simulations from the joint posterior distribution of the parameters of a 
#' mixture of Gaussian components with unknown number of components, by using a uniform prior for k. 
#'
#' @param y vector of observations.
#' @param K number of simulations. Default value at 10^5. 
#' @param bn_lw lower bound of the uniform prior distribution for the number of Gaussian components. Default value at 0.
#' @param bn_up upper bound of the uniform prior distribution for the number of Gaussian components. Default value at 30.
#' @param sh shape parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 1.5.
#' @param rt rate parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 0.5.
#' @param alpha_dir parameter of the Dirichlet prior distribution for the weight parameters of the Gaussian components. Default value at 1.
#' @param me_mu mean parameter of the Gaussian prior distribution for the mean parameters of the Gaussian components. Default value at 20.
#' @param sd_mu standard deviation parameter of the Gaussian prior distribution for the standard deviation parameters of the Gaussian components. Default value at 10.
#' @return k, a sample from the posterior distribution of the number of components.
#' @return param, a matrix from the posterior distribution of the weights, the means and the standard deviations of the components.
#' @return olp, a vector a joint prior density of the posterior samples.
#' @return oll, a vector a likelihood values of the posterior samples.
#' @keywords BayesMIC
#' @export
post_unif <- function(y,K=10^5,bn_lw=0,bn_up=30,sh=1.5,rt=0.5,
                      alpha_dir=1,me_mu=20,sd_mu=10){

  n=length(y)
  
  #initialise MCMC with one component
  mu=mean(y)
  sigma=sd(y)
  w=1
  
  Sdraw <- rep(1,n)
  G <- 1
  nh <- n
  
  u 		<- c()
  ystar 	<- rep(Inf,n)
  while(sum(abs(ystar)==Inf)!=0)
  {
    for(l in 1:n)
    {
      treshb   =  pnorm(y[l]+1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
      tresha   =  pnorm(y[l]-1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
      
      if(tresha != treshb)
      {
        u[l] = runif(1,tresha , treshb)	
        ystar[l] <- qnorm(u[l],mu[Sdraw[l]], sigma[Sdraw[l]])
      } else {
        ystar[l] <- y[l]
      }
      
    }	# end of the observations loop
  }
  
  oll=log_lkd_da(ystar,w,mu,sigma,Sdraw)
  olp=log_prior_unif(w,mu,sigma)
  
  #MCMC
  SS=10; Nsamp=K/SS;  
  PL=matrix(NA,Nsamp,3); colnames(PL)<-c("d","olp","oll"); TH=list();
  
  for (k in 1:K) {
    
    OK=TRUE
    d=length(w)
    if(d==1){
      move=2
      while(move==2 )
        move=sample(1:5,1) #choose a MCMC move (2 birth death, 3 fixed dimension)
    } else {
      move=sample(1:5,1)
    }

    if (move==1) {
      #add a cluster
      i=sample(1:d,1) #pick a component to split - the weight w[i] is shared 
      wn=runif(1,min=0,max=w[i]); wp=c(w,wn); wp[i]=w[i]-wn
      mun=runif(1,min(y)-5,max(y)+5); mup=c(mu,mun)
      sn=rgamma(1,shape=1.5,rate=0.5); sigmap=c(sigma,sn)
      qd=dunif(mun,min(y)-5,max(y)+5,log=TRUE)+dgamma(sn,shape=1.5,rate=0.5,log=TRUE)-log(w[i])
      qn=0
      rhod=-log(d)
      rhon=-log((d+1)*d)
    }
    if (move==2) {
      #kill a cluster
      if (d>1) {
        ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
        #pick component to kill and a component to increment with the weight from 
        #the killed cluster
        wp=w; wp[j]=w[j]+w[i]; wp=wp[-i]
        mup=mu[-i]; sigmap=sigma[-i]
        qd=0
        qn=dunif(mu[i],min(y)-5,max(y)+5,log=TRUE)+dgamma(sigma[i],shape=1.5,rate=0.5,log=TRUE)-log(w[i]+w[j])
        rhod=-log(d*(d-1))
        rhon=-log(d-1)
      } else {
        OK=FALSE
      }
    }
    if (move==3) {
      #fixed dimension mu - simple RW MH
      i=sample(1:d,1); 
      mup=mu; mup[i]=rnorm(1,mean=mu[i],sd=1);
      sigmap=sigma; wp=w;
      qd=qn=rhod=rhon=0      
    }
    if (move==4) {
      #fixed dimension sigma - scaling move
      i=sample(1:d,1); 
      sigmap=sigma; u=runif(1,0.5,2); sigmap[i]=sigma[i]*u
      mup=mu; wp=w;
      qn=-log(u); qd=rhod=rhon=0      
    }
    if (move==5) {
      #fixed dimension w - redistribute the weights between two clusters chosen UAR 
      if (d>1) {
        ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
        wp=w; wp[i]=runif(1,min=0,max=w[i]+w[j]); wp[j]=w[i]+w[j]-wp[i]
        mup=mu; sigmap=sigma;
        qn=qd=rhod=rhon=0
      } else {
        OK=FALSE
      }
    }
    
    # Relabeling
    
    temp_mu <- mup
    mup <- sort(mup)
    wp <- wp[order(temp_mu)]
    sigmap <- sigmap[order(temp_mu)]
    
    do_DA <- FALSE
    if (OK) {
      nll=log_lkd(y,wp,mup,sigmap)
      nlp=log_prior_unif(wp,mup,sigmap)  
      MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
      if (log(runif(1))<MHR) {
        do_DA <- TRUE
      }
    }
    
    ### Step 2: Multinomial sampler
    
    if(do_DA == TRUE){
      G <- length(wp)
      Sdraw <- c()
      for(l in 1:n)
      {
        probs <- c()
        for(g in 1:G)
        {
          probs[g] <-  dnorm(ystar[l], mup[g], sigmap[g])
        }
        probs[probs<0] <- 0
        if(sum(probs)!=0){
          Sdraw[l] <- sample(1:G,size=1,prob=probs)
        } else {
          Sdraw[l] <- Sdraw[l]
        }
      }
      
      ### Step 3: Data augmentation
      
      u 		<- c()
      ystar 	<- c()
      for(l in 1:n)
      {
        treshb   =  pnorm(y[l]+1, mup[Sdraw[l]],sigmap[Sdraw[l]]) 
        tresha   =  pnorm(y[l], mup[Sdraw[l]],sigmap[Sdraw[l]]) 
        u[l] = runif(1,tresha , treshb)
        ystar[l] <- qnorm(u[l],mup[Sdraw[l]], sigmap[Sdraw[l]])
        ystar[l] <- ifelse(abs(ystar[l])==Inf, y[l], ystar[l])
        
      }	# end of the observations loop
      
      # Acceptance / Rejection
      if (OK) {
        nll=log_lkd_da(ystar,wp,mup,sigmap,Sdraw)
        if(is.na(nll)==T){nll = -Inf}
        nlp=log_prior_unif(wp,mup,sigmap)  
        MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
        if (log(runif(1))<MHR) {
          w=wp; mu=mup; sigma=sigmap
          oll=nll; olp=nlp
        }
      }
    }		
    
    
    # Saving
    if (k%%SS==0) {
      TH[[k/SS]]=list(w=w,mu=mu,sigma=sigma)
      PL[k/SS,]=c(length(w),olp,oll)
    }
    
  }  
  return(list(k=PL[1],param=TH,olp=PL[,2],oll=PL[,3]))
}

#' Simulations from the joint posterior distribution of the parameters of a 
#' mixture of Gaussian components with unknown number of components, by using a Poisson prior for k. 
#'
#' @param y vector of observations.
#' @param K number of simulations. Default value at 10^5.
#' @param rt_pois rate parameter of the Poisson prior distribution for the number of Gaussian components. Default value at 1.
#' @param sh shape parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 1.5.
#' @param rt rate parameter of the inverse gamma prior distribution for the standard deviation parameters of the Gaussian components. Default value at 0.5.
#' @param alpha_dir parameter of the Dirichlet prior distribution for the weight parameters of the Gaussian components. Default value at 1.
#' @param me_mu mean parameter of the Gaussian prior distribution for the mean parameters of the Gaussian components. Default value at 20.
#' @param sd_mu standard deviation parameter of the Gaussian prior distribution for the standard deviation parameters of the Gaussian components. Default value at 10.
#' @return k, a sample from the posterior distribution of the number of components.
#' @return param, a matrix from the posterior distribution of the weights, the means and the standard deviations of the components.
#' @return olp, a vector a joint prior density of the posterior samples.
#' @return oll, a vector a likelihood values of the posterior samples.
#' @keywords BayesMIC
#' @export
post_pois <- function(y,K=10^5,rt_pois=1,sh=1.5,rt=0.5,
                      alpha_dir=1,me_mu=20,sd_mu=10){

  n=length(y)
  
  #initialise MCMC with one component
  mu=mean(y)
  sigma=sd(y)
  w=1
  
  Sdraw <- rep(1,n)
  G <- 1
  nh <- n
  
  u 		<- c()
  ystar 	<- rep(Inf,n)
  while(sum(abs(ystar)==Inf)!=0)
  {
    for(l in 1:n)
    {
      treshb   =  pnorm(y[l]+1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
      tresha   =  pnorm(y[l]-1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
      
      if(tresha != treshb)
      {
        u[l] = runif(1,tresha , treshb)	
        ystar[l] <- qnorm(u[l],mu[Sdraw[l]], sigma[Sdraw[l]])
      } else {
        ystar[l] <- y[l]
      }
      
    }	# end of the observations loop
  }
  
  oll=log_lkd_da(ystar,w,mu,sigma,Sdraw)
  olp=log_prior_pois(w,mu,sigma)
  
  #MCMC
  SS=10; Nsamp=K/SS;  
  PL=matrix(NA,Nsamp,3); colnames(PL)<-c("d","olp","oll"); TH=list();
  
  for (k in 1:K) {

    OK=TRUE
    d=length(w)
    if(d==1){
      move=2
      while(move==2 )
        move=sample(1:5,1) #choose a MCMC move (2 birth death, 3 fixed dimension)
    } else {
      move=sample(1:5,1)
    }
    
    if (move==1) {
      #add a cluster
      i=sample(1:d,1) #pick a component to split - the weight w[i] is shared 
      wn=runif(1,min=0,max=w[i]); wp=c(w,wn); wp[i]=w[i]-wn
      mun=runif(1,min(y)-5,max(y)+5); mup=c(mu,mun)
      sn=rgamma(1,shape=1.5,rate=0.5); sigmap=c(sigma,sn)
      qd=dunif(mun,min(y)-5,max(y)+5,log=TRUE)+dgamma(sn,shape=1.5,rate=0.5,log=TRUE)-log(w[i])
      qn=0
      rhod=-log(d)
      rhon=-log((d+1)*d)
    }
    if (move==2) {
      #kill a cluster
      if (d>1) {
        ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
        #pick component to kill and a component to increment with the weight from 
        #the killed cluster
        wp=w; wp[j]=w[j]+w[i]; wp=wp[-i]
        mup=mu[-i]; sigmap=sigma[-i]
        qd=0
        qn=dunif(mu[i],min(y)-5,max(y)+5,log=TRUE)+dgamma(sigma[i],shape=1.5,rate=0.5,log=TRUE)-log(w[i]+w[j])
        rhod=-log(d*(d-1))
        rhon=-log(d-1)
      } else {
        OK=FALSE
      }
    }
    if (move==3) {
      #fixed dimension mu - simple RW MH
      i=sample(1:d,1); 
      mup=mu; mup[i]=rnorm(1,mean=mu[i],sd=1);
      sigmap=sigma; wp=w;
      qd=qn=rhod=rhon=0      
    }
    if (move==4) {
      #fixed dimension sigma - scaling move
      i=sample(1:d,1); 
      sigmap=sigma; u=runif(1,0.5,2); sigmap[i]=sigma[i]*u
      mup=mu; wp=w;
      qn=-log(u); qd=rhod=rhon=0      
    }
    if (move==5) {
      #fixed dimension w - redistribute the weights between two clusters chosen UAR 
      if (d>1) {
        ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
        wp=w; wp[i]=runif(1,min=0,max=w[i]+w[j]); wp[j]=w[i]+w[j]-wp[i]
        mup=mu; sigmap=sigma;
        qn=qd=rhod=rhon=0
      } else {
        OK=FALSE
      }
    }
    
    # Relabeling
    
    temp_mu <- mup
    mup <- sort(mup)
    wp <- wp[order(temp_mu)]
    sigmap <- sigmap[order(temp_mu)]
    
    do_DA <- FALSE
    if (OK) {
      nll=log_lkd(y,wp,mup,sigmap)
      nlp=log_prior_pois(wp,mup,sigmap)  
      MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
      if (log(runif(1))<MHR) {
        do_DA <- TRUE
      }
    }
    
    ### Step 2: Multinomial sampler
    
    if(do_DA == TRUE){
      G <- length(wp)
      Sdraw <- c()
      for(l in 1:n)
      {
        probs <- c()
        for(g in 1:G)
        {
          probs[g] <-  dnorm(ystar[l], mup[g], sigmap[g])
        }
        probs[probs<0] <- 0
        if(sum(probs)!=0){
          Sdraw[l] <- sample(1:G,size=1,prob=probs)
        } else {
          Sdraw[l] <- Sdraw[l]
        }
      }
      
      ### Step 3: Data augmentation
      
      u 		<- c()
      ystar 	<- c()
      for(l in 1:n)
      {
        treshb   =  pnorm(y[l]+1, mup[Sdraw[l]],sigmap[Sdraw[l]]) 
        tresha   =  pnorm(y[l], mup[Sdraw[l]],sigmap[Sdraw[l]]) 
        u[l] = runif(1,tresha , treshb)
        ystar[l] <- qnorm(u[l],mup[Sdraw[l]], sigmap[Sdraw[l]])
        ystar[l] <- ifelse(abs(ystar[l])==Inf, y[l], ystar[l])
        
      }	# end of the observations loop
      
      # Acceptance / Rejection
      if (OK) {
        nll=log_lkd_da(ystar,wp,mup,sigmap,Sdraw)
        if(is.na(nll)==T){nll = -Inf}
        nlp=log_prior_pois(wp,mup,sigmap)  
        MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
        if (log(runif(1))<MHR) {
          w=wp; mu=mup; sigma=sigmap
          oll=nll; olp=nlp
        }
      }
    }		
    
    
    # Saving
    if (k%%SS==0) {
      TH[[k/SS]]=list(w=w,mu=mu,sigma=sigma)
      PL[k/SS,]=c(length(w),olp,oll)
    }
    
  }  
  return(list(k=PL[1],param=TH,olp=PL[,2],oll=PL[,3]))
  
}
