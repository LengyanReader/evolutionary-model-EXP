# Generating P, eigenvector matrix, by col
Gen_P<-function(k){
  I = diag(k)
  m1=matrix(rnorm(k*k),k,k)
  P = I + m1
 # D = diag(egv_H)
 # H = P%*%D%*%solve(P)
  
  return(P)
}

# generating H by eigenvalues and trait dimension k,full rank, not singular
Gen_H<-function(egv_H,k,P){
  #I = diag(k)
 #m1=matrix(rnorm(k*k),k,k)
 # P = I + m1
  D = diag(egv_H)
  H = P%*%D%*%solve(P)
  H<-Re(H)
  
  return(H)
}


# Generating Q, othognal eigenvector matrix, by col
Gen_Q<-function(k){
  A = matrix(rnorm(k*k),k,k)
  Q=qr.Q(qr(A))
  
  return(Q)
}

Gen_H_sym<-function(egv_H,k,Q){ # For H semi def

  H<-Q%*%diag(egv_H)%*%t(Q)

  return(H)
 # return(list(H=B,P=Q))
}


Gen_Sigma2<-function(egv_Sigma2,Q,k){ # for Sigma_x and Sigmax_e, both

  B<-Q%*%diag(egv_Sigma2)%*%t(Q)
  
  return(B )
}

##################################################################################

Gen_Sigma<-function(egv_Sigma,k){ # for sigma_x and sigmax_e, both
  A = matrix(rnorm(k*k),k,k)
  A[lower.tri(A)] <- 0
  diag(A)=egv_Sigma
  
  return(A)
}
#################################################################################


my_chol<-function(M){
  ## function called in estimBM.R evolmodelest.R matrixparametrizations.R
  ## M has to be symmetric--semi--positive--definte
  ## no check done for this done as the function is only called internally
  ## depending on platform two different types of calculations
  ## Debian with pivoting
  ## others without
  ## a different treatment takes place for Debian as the Debian flavor on CRAN
  ## raises an error to calls to chol(), noticed since 2019 XI 27
  ## other platforms do not do this
  ##    mchol<-matrix(NA,nrow=nrow(M),ncol=ncol(M))
  ##    platform<-R.Version()$platform 
  ##    if (grepl("deb", platform, ignore.case = TRUE)){## we are on Debian and some error seems to take place, so we pivot
  ##    	    org_warn<-options("warn")
  ##            options(warn=-1)
  ##            mchol<-chol(M,pivot=TRUE)
  ##            options(warn=org_warn$warn)
  ##    }else{## not Debian
  ##	mchol<-chol(M)
  ##    }
  mchol<-matrix(NA,nrow=nrow(M),ncol=ncol(M)) ## M has to be square
  tryCatch({mchol<-chol(M)},error=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})##,warning=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})
  if (is.na(mchol[1,1])){
    ## an error took place, most probably due to chol() considering M not to be of full rank
    ## therefore setting pivot=TRUE, to deal with such a matrix
    ## chol() throws a warning if M does not have full rank so we suppress it
    org_warn<-options("warn")
    options(warn=-1)
    tryCatch({mchol<-chol(M,pivot=TRUE)},error=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})##,warning=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})
    options(warn=org_warn$warn)    
    if (!is.null(attr(mchol,"pivot"))){attr(mchol,"pivot")<-NULL}
    if (!is.null(attr(mchol,"rank"))){attr(mchol,"rank")<-NULL}
  }
  mchol
}



is_othorgonal<-function(M,k,dif_scale=1e-10){
  dif<-M%*%t(M)%>%as.vector()-diag(1.0,k)%>%as.vector()
  dif<-dif%>%abs()
  
  if(max(dif)<=dif_scale)
    return(TRUE)
  else 
    return(FALSE)
  
}




rotate_cmp<-function(ang,k,k1,k2){
  r_mx<-diag(1,k)
  rt<-matrix(c(cos(ang),sin(ang),-sin(ang),cos(ang)),nrow = 2  ) 
  r_mx[k1,c(k1,k2)]<- rt[1,]
  r_mx[k2,c(k1,k2)]<- rt[2,]
  
  return(r_mx)
}    




Param2HTS<- function(param_HTS,k){
  # param_HTS = param[(k+1):(2*k^2+3*k)]
  
  #H----------------------------------------------
  P<-param_HTS[1:(k*k)]%>%matrix(nrow=k)
  #----------------------------------------------------------------------------------------------------------------------- 
  if(is.singular.matrix(P))return(0)
  egv_H<-param_HTS[(k*k+1):(k*k+k)]
  H<- Gen_H(egv_H,k,P)  #matrix(param_HTS[(k+1-k):(k+k^2-k)],k,k)
  
  #Theta
  Theta = param_HTS[(k*k+k+1):(k*k+2*k)]
  
  # Sigma2_x
  Q_sigma<- param_HTS[(k*k+2*k+1):(2*k*k+2*k)]%>%matrix(nrow=k)    # Q_Sigma2
  if(is.singular.matrix(Q_sigma))return(0)
  if(!is_othorgonal(Q_sigma,k,dif_scale=1e-10))return(0)
  
  egv_Sigma2_x<-param_HTS[(2*k*k+2*k+1):(2*k*k+3*k)]
  Sigma2_x<-Gen_Sigma2(egv_Sigma2_x,Q_sigma,k)
  Sigma_x<-my_chol(Sigma2_x)
  
  # Sigmae_x
  # Sigmae_x =matrix(0,k,k)
  # Sigmae_x[upper.tri(Sigmae_x, diag=TRUE)] = param_HTS[(1.5*k^2+2.5*k+1-k):(2*k^2+3*k-k)]
  #  diag(Sigmae_x) = exp(diag(Sigmae_x)) # to make diagonal non-negative
  #----------------------------------------------------------------------------------------------------------------------- 
  if(!is.singular.matrix(Sigma_x)){ return(list(H=H,Theta=Theta,Sigma_x=Sigma_x#,Sigmae_x=Sigmae_x
  ))}
  else{
    return(0)
  }
  
}


Param2HTS___123<- function(param_HTS,k){
  # param_HTS = param[(k+1):(2*k^2+3*k)]
  
  #H
  H<- matrix(param_HTS[(k+1-k):(k+k^2-k)],k,k)
  diag(H) = exp(diag(H))
  #Theta
  Theta = param_HTS[(k+k^2+1-k):(k^2+2*k-k)]
  # Sigma_x
  Sigma_x =matrix(0,k,k)
  Sigma_x[upper.tri(Sigma_x, diag=TRUE)] = param_HTS[(k^2+2*k+1-k):(1.5*k^2+2.5*k-k)]
  diag(Sigma_x) = exp(diag(Sigma_x)) # to make diagonal non-negative
  
  # Sigmae_x
  Sigmae_x =matrix(0,k,k)
  Sigmae_x[upper.tri(Sigmae_x, diag=TRUE)] = param_HTS[(1.5*k^2+2.5*k+1-k):(2*k^2+3*k-k)]
  diag(Sigmae_x) = exp(diag(Sigmae_x)) # to make diagonal non-negative
  #----------------------------------------------------------------------------------------------------------------------- 
  if(!is.singular.matrix(H) && !is.singular.matrix(Sigma_x) && !is.singular.matrix(Sigmae_x) && is_Cholesky_Factor(Sigma_x) &&  is_Cholesky_Factor(Sigma_x)   ){ return(list(H=H,Theta=Theta,Sigma_x=Sigma_x,Sigmae_x=Sigmae_x))}
  else{
    return(0)
  }
  
}

















#############################################################################################


# A fat function
#Checking if it is a Cholesky Factor , uppertriangle in the context here, consistent with pcmbase
#a) all eigenvalues >0
#b) is upper-triangle?
#c) mx%*%mx' is positive definite
is_Cholesky_Factor<-function(mx){ # suppose the mx is a upper-triangle matrix
  # mx = Sigma_x
  
  cond1 = fastmatrix::is.upper.tri(mx,TRUE)
  
  egv = eigen(mx)$values
  cond2<-length(which(egv>0)) # check if eigenvalues>0
  
  A = mx%*%t(mx) # check if semi-definite matrix
  #is.symmetric.matrix(A)
  cond3 = is.positive.definite(A)
  
  if(cond1>0 & cond2 & cond3){  return(TRUE)
  }else{return(FALSE) }
}


Gen_M_epH<-function(egv_H,k,t){  # desending order
  mx = matrix(0,k,k)
  for( i in c(1:k)){
    for( j in c(1:k)){
      e_sum_t = egv_H[i]+egv_H[j]
      if(e_sum_t==0){
        mx[i,j]=t
      }
      else{
        mx[i,j]=(1 - exp(-e_sum_t*t))/e_sum_t
      }
    }
  }
  return(mx)
}



#----------------------
HTS2VOP<-function(H,Theta,Sigma_x,k,tree,node,tips,pc){
  # VOP: Variance, Omega, Phi
  # node: "1"...
  
  t <- PCMTreeDtNodes(tree)[endNodeLab == node, endTime - startTime]
  #pc <- PCMInfo(X[, tree$tip.label], tree, model.OU.BM)$pc
  
  # active coordinates for tip t and its parent:
  k_c <- pc[, match(node, PCMTreeGetLabels(tree))] #child
  par_node<-tree$edge[,1][tree$edge[,2]==as.integer(node)]%>%as.character()
  k_p <- pc[, match(par_node, PCMTreeGetLabels(tree))]
  # parent
  
  Sigma <- Sigma_x %*% t(Sigma_x)
  res_egv_H <- eigen(H)#$values
  M_epH = Gen_M_epH(res_egv_H$values,k,t)
  M_psp = solve(res_egv_H$vectors)%*%Sigma%*%t(solve(res_egv_H$vectors))
  M_mid = M_epH * M_psp   #------------------------
  
  #-------------------------------------------------------------------------------------------
  V = (res_egv_H$vectors%*%M_mid%*%t(res_egv_H$vectors))[k_c, k_c] 
  V[lower.tri(V)] <- V[upper.tri(V)]
  I_k=diag(1,k)
  Omega = (I_k[k_c,] - expm(-H*t)[k_c,] )%*%Theta#%>%as.matrix() # column-vector
  Phi = expm(-H*t)[k_c, k_p]#%>%as.matrix()
  
  #--------------------
  
  return(list(V=V,Omega=Omega,Phi=Phi))
}