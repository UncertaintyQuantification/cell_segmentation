###### runnning time of lattice method
library(RobustGaSP)
library(ggplot2)
source('2dim_lattice_func.R')

param_here = c(-2,-2,-3)
beta = exp(param_here[1:2])
nu = exp(param_here[3])

sample_size = seq(10,80,10)
num_iter = 10
fast_comp_time_record_matern <- matrix(NA,num_iter,length(sample_size))
direct_comp_time_record_matern <- matrix(NA,num_iter,length(sample_size))
fast_comp_time_record_exp <- matrix(NA,num_iter,length(sample_size))
direct_comp_time_record_exp <- matrix(NA,num_iter,length(sample_size))

fast_comp_log_lik_record_matern <- matrix(NA,num_iter,length(sample_size))
direct_comp_log_lik_record_matern <- matrix(NA,num_iter,length(sample_size))
fast_comp_log_lik_record_exp <- matrix(NA,num_iter,length(sample_size))
direct_comp_log_lik_record_exp <- matrix(NA,num_iter,length(sample_size))


for(iter in 1:num_iter){
  print(iter)
  for(i in 1:length(sample_size)){
    print(i)
    n1 = sample_size[i]
    n2 = sample_size[i]
    N=n1*n2
    set.seed(iter)
    output = matrix(rnorm(N),n1,n2)
    output_vec = as.vector(output)
    
    input1 = 1:n1
    input2 = 1:n2
    R01=as.matrix(abs(outer(input1, input1, "-")))
    R02=as.matrix(abs(outer(input2, input2, "-")))
    
    
    ### 1. Matern ------------------------------------------
    ####### 1.1 fast computation
    kernel_type_here="matern"
    fast_comp_matern_time <- system.time({
      X=matrix(1,N,1) 
      q_X=dim(X)[2]
      X_list=as.list(1:q_X)
      for(i_q in 1:q_X){
        X_list[[i_q]]=matrix(X[,i_q],n1,n2)
      }
      fast_comp_log_lik_record_matern[iter,i] = Neg_log_lik_eigen_with_nugget(param = param_here,
                                                                         kernel_type=kernel_type_here,R01,R02,
                                                                         N, q_X, X_list,output)
    })
    fast_comp_time_record_matern[iter,i] = fast_comp_matern_time[1]
    
    ####### 1.2 direct computation
    direct_comp_matern_time = system.time({
      R1=Matern_5_2_funct(R01, beta=beta[1])
      R2=Matern_5_2_funct(R02, beta=beta[2])
      R_tilde = kronecker(R1, R2, FUN = "*") + diag(nu, N, N)
      inv_R_tilde = solve(R_tilde)
      one_vec = as.vector(rep(1,N))
      mu_hat = solve(t(one_vec)%*% solve(R_tilde)%*%one_vec)%*%t(one_vec)%*%inv_R_tilde%*%output_vec
      S2 = t(output_vec-c(mu_hat))  %*% inv_R_tilde %*% (output_vec-c(mu_hat))
      direct_comp_log_lik_record_matern[iter,i]=determinant(R_tilde)$modulus[1]/2 + N/2*log(S2)
    })
    direct_comp_time_record_matern[iter,i] = direct_comp_matern_time[1]
    
    ### 2. Exp ------------------------------------------
    ####### 2.1 fast computation
    kernel_type_here="exp"
    fast_comp_exp_time <- system.time({
      X=matrix(1,N,1) 
      q_X=dim(X)[2]
      X_list=as.list(1:q_X)
      for(i_q in 1:q_X){
        X_list[[i_q]]=matrix(X[,i_q],n1,n2)
      }
      fast_comp_log_lik_record_exp[iter,i] = Neg_log_lik_eigen_with_nugget(param = param_here,
                                                                      kernel_type=kernel_type_here,R01,R02,
                                                                      N, q_X, X_list,output)
    })
    fast_comp_time_record_exp[iter,i] = fast_comp_exp_time[1]
    
    ####### 2.2 direct computation
    direct_comp_exp_time = system.time({
      R1=Exp_funct(R01, beta=beta[1])
      R2=Exp_funct(R02, beta=beta[2])
      R_tilde = kronecker(R1, R2, FUN = "*") + diag(nu, N, N)
      inv_R_tilde = solve(R_tilde)
      one_vec = as.vector(rep(1,N))
      mu_hat = solve(t(one_vec)%*% solve(R_tilde)%*%one_vec)%*%t(one_vec)%*%inv_R_tilde%*%output_vec
      S2 = t(output_vec-c(mu_hat))  %*% inv_R_tilde %*% (output_vec-c(mu_hat))
      direct_comp_log_lik_record_exp[iter,i]=determinant(R_tilde)$modulus[1]/2 + N/2*log(S2)
    })
    direct_comp_time_record_exp[iter,i] = direct_comp_exp_time[1]
  }
}


direct_comp_time_record_matern_avg = colMeans(direct_comp_time_record_matern[1:4,])
direct_comp_time_record_exp_avg = colMeans(direct_comp_time_record_exp[1:4,])
fast_comp_time_record_matern_avg = colMeans(fast_comp_time_record_matern[1:4,])
fast_comp_time_record_exp_avg = colMeans(fast_comp_time_record_exp[1:4,])


### label, ticks larger
ggplot_df <- data.frame(
  time=c(direct_comp_time_record_matern_avg,direct_comp_time_record_exp_avg,
         fast_comp_time_record_matern_avg, fast_comp_time_record_exp_avg),
  N = rep(c("10^2","20^2","30^2","40^2","50^2","60^2","70^2","80^2"), 4),
  Methods = rep(c("Direct-Matern","Direct-Exp","Fast-Matern", "Fast-Exp"),each=8)
)
ggplot_df$N = factor(ggplot_df$N,
                     levels = c("10^2","20^2","30^2","40^2","50^2","60^2","70^2","80^2"))


f_1 = ggplot(ggplot_df, aes(x=N, y=time, group=Methods, color=Methods, shape=Methods)) + 
  geom_point(size=1.5) + 
  geom_line() +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Direct-Matern" = "darkgreen",  
      "Direct-Exp" = "orange",    
      "Fast-Matern" = "#578FCA",     
      "Fast-Exp" = "#FFC0CB"      
    )
  ) + 
  scale_x_discrete(
    labels = c(expression(10^{2}), expression(20^{2}), expression(30^{2}),
               expression(40^{2}), expression(50^{2}), expression(60^{2}),
               expression(70^{2}), expression(80^{2}))
  ) +
  labs(y="time (sec)") + 
  theme(legend.position = c(0.25,0.75),
        #legend.justification = c(0,1),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.text = element_text(size = 12),      # size of tick labels
        axis.title = element_text(size = 13))
f_1

f_2 = ggplot(ggplot_df, aes(x=N, y=log10(time), group=Methods, color=Methods,shape=Methods)) + 
  geom_point(size=1.5) + 
  geom_line() +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Direct-Matern" = "darkgreen", 
      "Direct-Exp" = "orange",    
      "Fast-Matern" = "#578FCA",    
      "Fast-Exp" = "#FFC0CB"      
    )
  ) +
  scale_y_continuous(
    breaks = c(-2, 0, 2),
    labels = c(expression(10^{-2}), expression(10^0), expression(10^2))
  ) + 
  scale_x_discrete(
    labels = c(expression(10^{2}), expression(20^{2}), expression(30^{2}),
               expression(40^{2}), expression(50^{2}), expression(60^{2}),
               expression(70^{2}), expression(80^{2}))
  ) + 
  labs(y="time (sec)") + 
  theme(legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        axis.text = element_text(size = 12),      # size of tick labels
        axis.title = element_text(size = 13))
f_2

f_merge <- gridExtra::grid.arrange(f_1,f_2,nrow=1)
ggsave("compare_comp_time.pdf", f_merge, height=2.5, width=8.5)


