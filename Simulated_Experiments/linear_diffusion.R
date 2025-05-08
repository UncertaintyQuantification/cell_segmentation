library(ReacTran)
library(deSolve)
library(plot3D)
library(RobustGaSP)
library(FastGaSP)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(latex2exp)
source('../src/2dim_lattice_func.R')
source("../src/DMD_alg_Apr20.R")


#### Generate mean of observations ---------------------------
k = 200 # dimension of observations
n = 200 # time steps

Grid <- setup.grid.1D(N=k, L=1)
pde1D <- function(t, C, parms){
  tran <- tran.1D(C=C, D=D, C.down=Cext, dx=Grid)$dC
  list(tran-Q)
}

D = 1
Q = 0
Cext = 1
times <- seq(0, 0.2, length.out=n)
reality <- ode.1D(y = rep(0, Grid$N), times = times, func = pde1D, 
                  parms = NULL, nspec = 1)
reality <- t(reality[,2:(k+1)]) # mean of obs, each column is associated with a time step


num_repetition = 10
sigma0_list = c(0.05, 0.1, 0.3)
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
  sigma0 = sigma0_list[i_sigma0]
  for(iter in 1:num_repetition){
    set.seed(iter)
    if(iter==5){print(iter)}
    
    #### Generate noisy observations ----------------
    y = reality + matrix(rnorm(n*k,mean=0,sd=sigma0),k,n)
    if(iter==1){
      y_record[[i_sigma0]] = y
    }
    
    #### Lattice - exp kernel -------------------------------------------
    input1_here = 1:k
    input2_here = 1:n
    param_ini_here=c(-2,-2,-3)
    kernel_type_here="exp"
    optim_method_here = "Nelder-Mead" # or "L-BFGS-B"
    
    fit_lattice_exp <- lattice_alg(y, input1_here, input2_here, kernel_type=kernel_type_here,
                                   testing_input1=input1_here, testing_input2=input2_here,
                                   param_ini = param_ini_here, optim_method = optim_method_here)
    
    pred_mean_lattice_exp = fit_lattice_exp$pred_mean
    rmse_lattice_exp[iter, i_sigma0] = sqrt(mean((reality-pred_mean_lattice_exp)^2))
    if(iter==1){
      pred_mean_lattice_exp_record[[i_sigma0]] = pred_mean_lattice_exp
    }
    
    #### Lattice - matern kernel -------------------------------------------
    input1_here = 1:k
    input2_here = 1:n
    param_ini_here=c(-2,-2,-3)
    kernel_type_here="matern"
    optim_method_here = "Nelder-Mead"
    fit_lattice_matern <- lattice_alg(y, input1_here, input2_here, kernel_type=kernel_type_here,
                                        testing_input1=input1_here, testing_input2=input2_here,
                                        param_ini = param_ini_here, optim_method = optim_method_here)
    pred_mean_lattice_matern = fit_lattice_matern$pred_mean
    rmse_lattice_matern[iter, i_sigma0] = sqrt(mean((reality-pred_mean_lattice_matern)^2))
    if(iter==1){
      pred_mean_lattice_matern_record[[i_sigma0]] = pred_mean_lattice_matern
    }
    
    #### choose number of latent states -------------------------------------------
    svd_output=svd(y)
  
    U = svd_output$u
    loss_score = NULL
    for(i_d in 1:ceiling(dim(y)[1]*2/3)){
      criteria_val_cur = log(mean((y - U[,1:i_d]%*%t(U[,1:i_d])%*%y)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score = c(loss_score,  criteria_val_cur)
    }
    est_d = which.min(loss_score)
    est_d_record[iter, i_sigma0] = est_d
    
    #### FMOU -------------------------------------------
    m_fmou <- fmou(y, d=est_d)
    fit_fmou <- fit.fmou(m_fmou, M=100, threshold=10^{-6})
    rmse_fmou[iter, i_sigma0] = sqrt(mean((reality-fit_fmou$mean_obs)^2))
    if(iter==1){
      pred_mean_fmou_record[[i_sigma0]] = fit_fmou$mean_obs
    }
    
    
    #### PCA -------------------------------------------
    pred_mean_pca = U[,1:est_d] %*% t(U[,1:est_d]) %*% y
    rmse_pca[iter, i_sigma0] = sqrt(mean((reality-pred_mean_pca)^2))
    if(iter==1){
      pred_mean_pca_record[[i_sigma0]] = pred_mean_pca
    }
    
    #### DMD -------------------------------------------
    fit_dmd <- DMD_alg(y,r=est_d, fix_r=T) #
    rmse_dmd[iter, i_sigma0] = sqrt(mean((cbind(y[,1], fit_dmd$in_sample_pred) - reality)^2))
    if(iter==1){
      pred_mean_dmd_record[[i_sigma0]] = cbind(y[,1], fit_dmd$in_sample_pred)
    }
  }
} 


### RMSE table --------------------------------------
rmse_summary = rbind(colMeans(rmse_lattice_exp),colMeans(rmse_lattice_matern),
                     colMeans(rmse_fmou), colMeans(rmse_pca), colMeans(rmse_dmd))
rownames(rmse_summary) = c("lattice_exp","lattice_matern","fmou","pca","dmd") 
colnames(rmse_summary) = sigma0_list
rmse_summary



# > rmse_summary
#                    0.05        0.1        0.3
# lattice_exp    0.010072687 0.01375717 0.02222036
# lattice_matern 0.008484715 0.01145377 0.01863110
# fmou           0.010033361 0.01767439 0.04192782
# pca            0.011288589 0.02041124 0.05092774
# dmd            0.015212327 0.02316704 0.05441586

############### Figure 5(B): violin plot of RMSE
rmse_boxplot_df <- data.frame(
  RMSE = c(rmse_lattice_exp,rmse_lattice_matern,rmse_fmou,rmse_pca,rmse_dmd),
  Methods = rep(c("Fast-Exp","Fast-Mat","FMOU","PCA","DMD") ,each=30),
  noise_level = rep(rep(sigma0_list, each=10),5)
) %>%
  dplyr::filter(RMSE!=max(rmse_dmd[,1]))
rmse_boxplot_df$noise_level = factor(rmse_boxplot_df$noise_level, 
                                     labels = c(expression(sigma[0] == 0.05), expression(sigma[0] == 0.1), expression(sigma[0] == 0.3)))
rmse_boxplot_df$Methods = factor(rmse_boxplot_df$Methods,
                                 levels = c("Fast-Mat","PCA","FMOU","DMD","Fast-Exp"))



plot_rmse <- ggplot(rmse_boxplot_df, aes(x=Methods, y=RMSE, fill=Methods)) +
  geom_violin() +
  facet_wrap(~noise_level, scales = "free_y",labeller = label_parsed) +
  scale_fill_manual(values = c("Fast-Exp" = "#FFC0CB","Fast-Mat" = "#578FCA",
                               "FMOU" = "#FF8989","PCA" = "#fee440","DMD" = "#66D2CE")) + 
  theme_bw() +
  labs(title = "(B) Linear diffusion" ) + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12))

(plot_rmse)

############### Figure 6(B): obs mean, obs, pred mean
my_palette <- colorRampPalette(c("#8ecae6","#219ebc","#023047","#ffb703","#fb8500"))(100)
signal_image_df <- expand.grid(1:k, 1:n)
signal_image_df$z <- as.vector(reality)
names(signal_image_df) = c("input1", "input2","z")

f_1 <- ggplot(signal_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() +  # Use geom_tile() or geom_raster()
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(D) Observation mean", x = "", y = "", fill = "") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())

y_image_df <- expand.grid(1:k, 1:n)
y_image_df$z <- as.vector(y_record[[3]])

names(y_image_df) = c("input1", "input2","z")

f_2 <- ggplot(y_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() +  # Use geom_tile() or geom_raster()
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(E) Noisy observation", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())

lattice_matern_image_df <- expand.grid(1:k, 1:n)
lattice_matern_image_df$z <- as.vector(pred_mean_lattice_matern_record[[3]])
names(lattice_matern_image_df) = c("input1", "input2","z")

f_3 <- ggplot(lattice_matern_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() + 
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(F) Predictive mean", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())


f_merge <- grid.arrange(f_1,f_2,f_3,nrow=1,
                        top = textGrob("Linear diffusion", gp = gpar(fontsize = 12, fontface = "bold"),
                                       vjust=0.2))
