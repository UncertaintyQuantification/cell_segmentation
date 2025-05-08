Matern_5_2_funct<-function(d,beta){
 x=sqrt(5)*beta*d
 (1+x+x^2/3)*exp(-x)
}

Exp_funct<-function(d,beta){
  x=beta*d
  exp(-x)
}

## define loss functions
Neg_log_lik_eigen_with_nugget <- function(param,kernel_type,R01,R02,N, q_X,
                                          X_list,output_mat) {
  output = as.vector(output_mat)
  X=matrix(1,N,1) 
  n1 = nrow(R01)
  n2 = nrow(R02)
  beta=exp(param[1:2])
  nu=exp(param[3])
  if(kernel_type=="matern"){
    R1=Matern_5_2_funct(R01, beta=beta[1])
    R2=Matern_5_2_funct(R02, beta=beta[2])
  }
  if(kernel_type=="exp"){
    R1=Exp_funct(R01, beta=beta[1])
    R2=Exp_funct(R02, beta=beta[2])
  }
  
  eigen_R1=eigen(R1)  
  eigen_R2=eigen(R2)  
  U_x=matrix(NA,N,q_X)
  for(i_q in 1:q_X){
    U_x=as.vector(t(eigen_R1$vectors)%*% X_list[[i_q]]%*%eigen_R2$vectors)
  }
  
  Lambda_tilde_inv=1/(kronecker(eigen_R2$values,eigen_R1$values, FUN = "*")+nu)
  Lambda_tilde_inv_U_x= Lambda_tilde_inv*U_x
  X_R_tilde_inv_X_inv=solve(t(U_x)%*%(Lambda_tilde_inv_U_x))
  output_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat%*%eigen_R2$vectors)
  theta_hat=X_R_tilde_inv_X_inv%*%(t(Lambda_tilde_inv_U_x)%*%output_tilde)
  output_mat_normalized=matrix(output-X%*%theta_hat,n1,n2)
  output_normalize_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat_normalized%*%eigen_R2$vectors)
  S_2=sum((output_normalize_tilde)*Lambda_tilde_inv*output_normalize_tilde  )  
  -(1/2*sum(log(Lambda_tilde_inv))-N/2*log(S_2))
}

lattice_alg <- function(output_mat,input1, input2,kernel_type="matern", 
                        testing_input1=input1, testing_input2=input2,
                        param_ini = c(-2,-2,-3), optim_method='Nelder-Mead'){
  #@output_mat: n1*n2 output matrix
  #@input1:input of the first direction
  #@input2:input of the second direction
  #@kernel_type: default is materm_5_2, can be changed to "exp" (OU).
  #@testing_input1: testing input of the first direction
  #@testing_input2: testing input of the second direction.
  #@param_ini: initial points of optimization. Default is (-2,-2,-3)
  #@optim_method: method used in optim(). Default is "Nelder-Mead" (for large sample size). Can be changed to "L-BFGS-B".
  
  
  
  ## define needed variables
  p = 2 # 2D input
  n1 = length(input1)
  n2 = length(input2)
  N = n1*n2
  output = as.vector(output_mat)
  
  ##constant mean basis
  X=matrix(1,N,1) 
  q_X=dim(X)[2]
  X_list=as.list(1:q_X)
  for(i_q in 1:q_X){
    X_list[[i_q]]=matrix(X[,i_q],n1,n2)
  }
  
  ##input distance
  R01=as.matrix(abs(outer(input1, input1, "-")))
  R02=as.matrix(abs(outer(input2, input2, "-")))
  
  ## optimize parameters
  m_eigen = try(optim(param_ini,Neg_log_lik_eigen_with_nugget,
                      kernel_type=kernel_type,R01=R01,R02=R02,N=N, q_X=q_X,
                      X_list=X_list,output_mat=output_mat,
                      method=optim_method),silent=T)
  while(!is.numeric(m_eigen[[1]])){
    m_eigen=try(optim(param_ini+runif(3), Neg_log_lik_eigen_with_nugget,
                      kernel_type=kernel_type,R01=R01,R02=R02,N=N, q_X=q_X,
                      X_list=X_list,output=output,
                      method=optim_method),silent=T)
  }

  ## build predictive mean
  beta= exp(m_eigen$par[1:p])
  nu=exp(m_eigen$par[p+1])

  if(kernel_type=="matern"){
    R1=Matern_5_2_funct(R01, beta=beta[1])
    R2=Matern_5_2_funct(R02, beta=beta[2])
  }else{
    R1=Exp_funct(R01, beta=beta[1])
    R2=Exp_funct(R02, beta=beta[2])
  }

  eigen_R1=eigen(R1)
  eigen_R2=eigen(R2)

  U_x=matrix(NA,N,q_X)
  for(i_q in 1:q_X){
    U_x=as.vector(t(eigen_R1$vectors)%*% X_list[[i_q]]%*%eigen_R2$vectors)
  }

  Lambda_tilde_inv=1/(kronecker(eigen_R2$values,eigen_R1$values, FUN = "*")+nu)
  Lambda_tilde_inv_U_x= Lambda_tilde_inv*U_x
  X_R_tilde_inv_X_inv=solve(t(U_x)%*%(Lambda_tilde_inv_U_x))
  output_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat%*%eigen_R2$vectors)
  theta_hat=X_R_tilde_inv_X_inv%*%(t(Lambda_tilde_inv_U_x)%*%output_tilde)
  output_mat_normalized=matrix(output-X%*%theta_hat,n1,n2)

  ### test on input,it will be faster if we use filling
  testing_input1=testing_input1
  testing_input2=testing_input2

  X_testing=X

  r01=abs(outer(input1,testing_input1,'-'))
  r02=abs(outer(input2,testing_input2,'-'))
  if(kernel_type=="matern"){
    r1=Matern_5_2_funct(r01, beta[1])
    r2=Matern_5_2_funct(r02, beta[2])
  }else{
    r1 = Exp_funct(r01,beta[1])
    r2 = Exp_funct(r02, beta[2])
  }

  output_normalize_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat_normalized%*%eigen_R2$vectors)
  output_normalized_tilde_lambda_inv_mat=matrix(Lambda_tilde_inv*output_normalize_tilde,n1,n2)
  R_tilde_inv_output_normalize=as.vector((eigen_R1$vectors)%*% output_normalized_tilde_lambda_inv_mat%*%t(eigen_R2$vectors))
  R_tilde_inv_output_normalize_mat=matrix(R_tilde_inv_output_normalize,n1,n2)
  predmean=X_testing%*%theta_hat+as.vector(t(r1)%*%R_tilde_inv_output_normalize_mat%*%r2)
  predmean_mat=matrix(predmean,n1,n2)

  res = list(beta = beta,
             nu=nu,
             pred_mean=predmean_mat)
  return(res)
}
  
  