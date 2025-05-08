library(lhs)
library(RobustGaSP)
library(plot3D)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(latex2exp)
library(FastGaSP)

source("../src/2dim_lattice_func.R")
source("../src/DMD_alg_Apr20.R")

### 1) lattice_matern has the smallest RMSE across all sigma0
### 2) FMatern performs better than FMOU across all sigma0

#### Generate mean of observations ---------------------------
branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}

n1=100
n2=100
N=n1*n2
p=2 ##2D input
set.seed(1)
LB=c(-5,0)
UB=c(10,15)
range=UB-LB
input1=LB[1]+range[1]*seq(0,1,1/(n1-1)) ##row
input2=LB[2]+range[2]*seq(0,1,1/(n2-1)) ##column
input=cbind(rep(input1,n2),as.vector(t(matrix(input2,n2,n1))))  

output=matrix(NA,N,1)
f=matrix(NA,N,1)
for(i in 1:N){
  f[i]=branin(input[i,])
}

f_mat = matrix(f,n1,n2)
par(mfrow=c(1,1))
image2D(f_mat)


num_repetition = 10
sigma0_list = c(1,5,10)
rmse_lattice_exp <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
rmse_lattice_matern <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
rmse_fmou <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
rmse_pca <- matrix(NA,  nrow=num_repetition, ncol=length(sigma0_list))
rmse_dmd <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))

y_record <- as.list(1:length(sigma0_list))
pred_mean_lattice_exp_record <- as.list(1:length(sigma0_list))
pred_mean_lattice_matern_record <- as.list(1:length(sigma0_list))
pred_mean_fmou_record <- as.list(1:length(sigma0_list))
pred_mean_pca_record <- as.list(1:length(sigma0_list))
pred_mean_dmd_record <- as.list(1:length(sigma0_list))

est_d_record <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))

#### Start 
for(i_sigma0 in 1:length(sigma0_list)){
  sigma_0 = sigma0_list[i_sigma0]
  for(iter in 1:num_repetition){
    set.seed(iter)
    if(iter==1){print(iter)}
    
    #### Generate noisy observations ----------------
    output=f+rnorm(N,mean=0,sd=sigma_0)
    output_mat=matrix(output,n1,n2)
    if(iter==1){
      y_record[[i_sigma0]] = output_mat
    }
    
    #### Lattice - exp kernel -------------------------------------------
    input1_here = input1
    input2_here = input2
    param_ini_here=c(-2,-2,-3)
    kernel_type_here="exp"
    optim_method_here = "L-BFGS-B" # Nelder-Mead
    fit_lattice_exp <- lattice_alg(output_mat, input1_here, input2_here, kernel_type=kernel_type_here,
                                   testing_input1=input1_here, testing_input2=input2_here,
                                   param_ini = param_ini_here, optim_method = optim_method_here)
    pred_mean_lattice_exp = fit_lattice_exp$pred_mean
    rmse_lattice_exp[iter, i_sigma0] = sqrt(mean((f_mat-pred_mean_lattice_exp)^2))
    if(iter==1){
      pred_mean_lattice_exp_record[[i_sigma0]] = pred_mean_lattice_exp
    }
    
    #### Lattice - matern kernel -------------------------------------------
    input1_here = input1
    input2_here = input2
    param_ini_here=c(-2,-2,-3)
    kernel_type_here="matern"
    optim_method_here = "Nelder-Mead" # Nelder-Mead
    fit_lattice_matern <- lattice_alg(output_mat, input1_here, input2_here, kernel_type=kernel_type_here,
                                      testing_input1=input1_here, testing_input2=input2_here,
                                      param_ini = param_ini_here, optim_method = optim_method_here)
    pred_mean_lattice_matern = fit_lattice_matern$pred_mean
    rmse_lattice_matern[iter, i_sigma0] = sqrt(mean((f_mat-pred_mean_lattice_matern)^2))
    if(iter==1){
      pred_mean_lattice_matern_record[[i_sigma0]] = pred_mean_lattice_matern
    }
    
    #### choose number of latent states -------------------------------------------
    svd_output=svd(output_mat)
    U = svd_output$u
    k=n1
    n=n2
    loss_score = NULL
    for(i_d in 1:ceiling(dim(output_mat)[1]*2/3)){
      criteria_val_cur = log(mean((output_mat - U[,1:i_d]%*%t(U[,1:i_d])%*%output_mat)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score = c(loss_score,  criteria_val_cur)
    }
    est_d_here = which.min(loss_score)
    est_d_record[iter, i_sigma0] = est_d_here
    
    #### FMOU -------------------------------------------
    m_fmou <- fmou(output_mat, d=est_d_here)
    fit_fmou <- fit.fmou(m_fmou)
    pred_mean_fmou = fit_fmou$mean_obs
    rmse_fmou[iter, i_sigma0] = sqrt(mean((f_mat-pred_mean_fmou)^2))
    if(iter==1){
      pred_mean_fmou_record[[i_sigma0]] = pred_mean_fmou
    }
    
    
    #### PCA -------------------------------------------
    pred_mean_pca = U[,1:est_d_here] %*% t(U[,1:est_d_here]) %*% output_mat
    rmse_pca[iter, i_sigma0] = sqrt(mean((f_mat-pred_mean_pca)^2))
    if(iter==1){
      pred_mean_pca_record[[i_sigma0]] = pred_mean_pca
    }
    
    #### DMD, fix r -------------------------------------------
    fit_dmd_fix_r <- DMD_alg(output_mat,r=est_d_here, fix_r=T) 
    rmse_dmd[iter, i_sigma0] = sqrt(mean((cbind(output_mat[,1], fit_dmd_fix_r$in_sample_pred) - f_mat)^2))
    if(iter==1){
      pred_mean_dmd_record[[i_sigma0]] = cbind(output_mat[,1], fit_dmd_fix_r$in_sample_pred)
    }
  }
  
}

################################## summary of RMSE 

rmse_summary = rbind(colMeans(rmse_lattice_exp),colMeans(rmse_lattice_matern),
                     colMeans(rmse_fmou), colMeans(rmse_pca), colMeans(rmse_dmd))
rownames(rmse_summary) = c("lattice_exp","lattice_matern","fmou", "PCA", "DMD")
colnames(rmse_summary) = sigma0_list
rmse_summary

# > rmse_summary
#                     1         5        10
# lattice_exp    0.32992148 1.0889037 1.7990646
# lattice_matern 0.07257045 0.2955026 0.5347066
# fmou           0.23933558 1.0913727 2.0494883
# PCA            0.24378955 1.2223282 2.4615920
# DMD            0.26179583 1.3062734 2.5988309

################################## Figure 5(A): violin plot of RMSE 


rmse_boxplot_df <- data.frame(
  RMSE = c(rmse_lattice_exp,rmse_lattice_matern,rmse_fmou,rmse_pca,rmse_dmd),
  Methods = rep(c("Fast-Exp","Fast-Mat","FMOU","PCA","DMD") ,each=30),
  noise_level = rep(rep(sigma0_list, each=10),5)
)
rmse_boxplot_df$noise_level = factor(rmse_boxplot_df$noise_level, 
                                     labels = c(expression(sigma[0] == 1), expression(sigma[0] == 5), expression(sigma[0] == 10)))
rmse_boxplot_df$Methods = factor(rmse_boxplot_df$Methods,
                                 levels = c("Fast-Mat","PCA","FMOU","DMD","Fast-Exp"))



plot_rmse <- ggplot(rmse_boxplot_df, aes(x=Methods, y=RMSE, fill=Methods)) +
  geom_violin() +
  facet_wrap(~noise_level, scales = "free_y",labeller = label_parsed) +
  scale_fill_manual(values = c("Fast-Exp" = "#FFC0CB","Fast-Mat" = "#578FCA",
                               "FMOU" = "#FF8989","PCA" = "#fee440","DMD" = "#66D2CE")) + 
  theme_bw() +
  labs(title = "(A) Branin" ) +
  theme(axis.title.x = element_blank(),
        #axis.text.x   = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12))
plot_rmse


par(mfrow=c(1,3))
image2D(t(f_mat))
image2D(t(y_record[[3]]))
image2D(t(pred_mean_lattice_matern_record[[3]]))

################################## Figure 6(A): obs mean, obs, pred mean
my_palette <- colorRampPalette(c("#8ecae6","#219ebc","#023047","#ffb703","#fb8500"))(100)
mean_df <- expand.grid(input1, input2)
mean_df$z <- f
names(mean_df) = c("input1", "input2","z")

f_1 <- ggplot(mean_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() + 
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(A) Observation mean", x = "", y = "", fill = "") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())

y_image_df <- expand.grid(input1, input2)
y_image_df$z <- as.vector(y_record[[3]])
names(y_image_df) = c("input1", "input2","z")

f_2 <- ggplot(y_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() +  
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(B) Noisy observation", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())
f_2

lattice_matern_image_df <- expand.grid(input1, input2)
lattice_matern_image_df$z <- as.vector(pred_mean_lattice_matern_record[[3]])
names(lattice_matern_image_df) = c("input1", "input2","z")

f_3 <- ggplot(lattice_matern_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() +  # Use geom_tile() or geom_raster()
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(C) Predictive mean", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())
f_3

f_merge <- grid.arrange(f_1,f_2,f_3,nrow=1,
                        top = textGrob("Branin", gp = gpar(fontsize = 12, fontface = "bold"),
                                       vjust=0.2))

#ggsave("signal_obs_pred_branin.pdf",f_merge, width=7.5, height=2.5)

