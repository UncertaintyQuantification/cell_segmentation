library(png)
library(plot3D)
library(ggplot2)
library(gridExtra)
library(jpeg)
library(FastGaSP)
library(tidyr)
library(latex2exp)
library(grid)
source('../src/2dim_lattice_func.R')
source("../src/DMD_alg_Apr20.R")

# load images
nuclei_figure <- readPNG("../Image_Data/simulation_cells/nuclei_2.png")
nuclei_figure <- nuclei_figure[,,1]
whole_cell_figure <- readJPEG("../Image_Data/simulation_cells/whole_cell_2.jpg")


sigma0_list = c(0.1,0.3,0.5)
num_repetition = 10


######## nuclei_figure ------------------------------------
nuclei_rmse_lattice_exp <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
nuclei_rmse_lattice_matern <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
nuclei_rmse_fmou <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
nuclei_rmse_pca <- matrix(NA,  nrow=num_repetition, ncol=length(sigma0_list))
nuclei_rmse_dmd <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))


nuclei_y_record <- as.list(1:length(sigma0_list))
nuclei_pred_mean_lattice_exp_record <- as.list(1:length(sigma0_list))
nuclei_pred_mean_lattice_matern_record <- as.list(1:length(sigma0_list))
nuclei_pred_mean_fmou_record <- as.list(1:length(sigma0_list))
nuclei_pred_mean_pca_record <- as.list(1:length(sigma0_list))
nuclei_pred_mean_dmd_record <- as.list(1:length(sigma0_list))

nuclei_est_d_record <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))

#### Start 
for(i_sigma0 in 1:length(sigma0_list)){
  sigma_0 = sigma0_list[i_sigma0]
  for(iter in 1:num_repetition){
    set.seed(iter)
    if(iter==1){print(iter)}
    
    #### Generate noisy observations ---------------
    output_mat=nuclei_figure + matrix(rnorm(length(nuclei_figure),sd=sigma_0), nrow(nuclei_figure), ncol(nuclei_figure))
    if(iter==1){
      nuclei_y_record[[i_sigma0]] = output_mat
    }
    
    #### Lattice - exp kernel -------------------------------------------
    input1_here = 1:nrow(nuclei_figure)
    input2_here = 1:ncol(nuclei_figure)
    param_ini_here=c(-2,-2,2)
    kernel_type_here="exp"
    optim_method_here = "Nelder-Mead" # Nelder-Mead
    fit_lattice_exp <- lattice_alg(output_mat, input1_here, input2_here, kernel_type=kernel_type_here,
                                   testing_input1=input1_here, testing_input2=input2_here,
                                   param_ini = param_ini_here, optim_method = optim_method_here)
    pred_mean_lattice_exp = fit_lattice_exp$pred_mean
    nuclei_rmse_lattice_exp[iter, i_sigma0] = sqrt(mean((nuclei_figure-pred_mean_lattice_exp)^2))
    if(iter==1){
      nuclei_pred_mean_lattice_exp_record[[i_sigma0]] = pred_mean_lattice_exp
    }
    
    #### Lattice - matern kernel -------------------------------------------
    input1_here = 1:nrow(nuclei_figure)
    input2_here = 1:ncol(nuclei_figure)
    param_ini_here=c(-2,-2,2)
    kernel_type_here="matern"
    optim_method_here = "Nelder-Mead" # Nelder-Mead
    fit_lattice_matern <- lattice_alg(output_mat, input1_here, input2_here, kernel_type=kernel_type_here,
                                      testing_input1=input1_here, testing_input2=input2_here,
                                      param_ini = param_ini_here, optim_method = optim_method_here)
    pred_mean_lattice_matern = fit_lattice_matern$pred_mean
    nuclei_rmse_lattice_matern[iter, i_sigma0] = sqrt(mean((nuclei_figure-pred_mean_lattice_matern)^2))
    if(iter==1){
      nuclei_pred_mean_lattice_matern_record[[i_sigma0]] = pred_mean_lattice_matern
    }
    
    #### choose number of latent states -------------------------------------------
    svd_output=svd(output_mat)
    U = svd_output$u
    k=nrow(nuclei_figure)
    n=ncol(nuclei_figure)
    loss_score = NULL
    for(i_d in 1:ceiling(dim(output_mat)[1]*2/3)){
      criteria_val_cur = log(mean((output_mat - U[,1:i_d]%*%t(U[,1:i_d])%*%output_mat)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score = c(loss_score,  criteria_val_cur)
    }
    est_d_here = which.min(loss_score)
    nuclei_est_d_record[iter, i_sigma0] = est_d_here
    
    #### FMOU -------------------------------------------
    m_fmou <- fmou(output_mat, d=est_d_here)
    fit_fmou <- fit.fmou(m_fmou)
    pred_mean_fmou = fit_fmou$mean_obs
    nuclei_rmse_fmou[iter, i_sigma0] = sqrt(mean((nuclei_figure-pred_mean_fmou)^2))
    if(iter==1){
      nuclei_pred_mean_fmou_record[[i_sigma0]] = pred_mean_fmou
    }
    
    
    #### PCA -------------------------------------------
    pred_mean_pca = U[,1:est_d_here] %*% t(U[,1:est_d_here]) %*% output_mat
    nuclei_rmse_pca[iter, i_sigma0] = sqrt(mean((nuclei_figure-pred_mean_pca)^2))
    if(iter==1){
      nuclei_pred_mean_pca_record[[i_sigma0]] = pred_mean_pca
    }
    
    #### DMD, fix r -------------------------------------------
    fit_dmd_fix_r <- DMD_alg(output_mat,r=est_d_here, fix_r=T) 
    nuclei_rmse_dmd[iter, i_sigma0] = sqrt(mean((cbind(output_mat[,1], fit_dmd_fix_r$in_sample_pred) - nuclei_figure)^2))
    if(iter==1){
      nuclei_pred_mean_dmd_record[[i_sigma0]] = cbind(output_mat[,1], fit_dmd_fix_r$in_sample_pred)
    }
  }
  
}

nuclei_rmse_summary = rbind(colMeans(nuclei_rmse_lattice_exp),colMeans(nuclei_rmse_lattice_matern),
                            colMeans(nuclei_rmse_fmou), colMeans(nuclei_rmse_pca), colMeans(nuclei_rmse_dmd))
rownames(nuclei_rmse_summary) = c("lattice_exp","lattice_matern","fmou","PCA", "DMD")
colnames(nuclei_rmse_summary) = sigma0_list
nuclei_rmse_summary

# > nuclei_rmse_summary
#                    0.1        0.3        0.5
# lattice_exp    0.02896498 0.04557247 0.05473485
# lattice_matern 0.02688344 0.04404816 0.05392391
# fmou           0.05888946 0.06955734 0.07284038
# PCA            0.05993414 0.07165570 0.07859802
# DMD            0.06034515 0.07326585 0.08255500




nuclei_rmse_boxplot_df <- data.frame(
  RMSE = c(nuclei_rmse_lattice_exp,nuclei_rmse_lattice_matern,nuclei_rmse_fmou,nuclei_rmse_pca,nuclei_rmse_dmd),
  Methods = rep(c("Fast-Exp","Fast-Mat","FMOU","PCA","DMD") ,each=30),
  noise_level = rep(rep(sigma0_list, each=10),5)
)
nuclei_rmse_boxplot_df$noise_level = factor(nuclei_rmse_boxplot_df$noise_level, 
                                            labels = c(expression(sigma[0] == 0.1), expression(sigma[0] == 0.3), expression(sigma[0] == 0.5)))
nuclei_rmse_boxplot_df$Methods = factor(nuclei_rmse_boxplot_df$Methods,
                                        levels = c("Fast-Mat","PCA","FMOU","DMD","Fast-Exp"))

####### Figure 7(A): violin plot of RMSE

nuclei_plot_rmse <- ggplot(nuclei_rmse_boxplot_df, aes(x=Methods, y=RMSE, fill=Methods)) +
  geom_violin() +
  facet_wrap(~noise_level, scales = "free_y",labeller = label_parsed) +
  scale_fill_manual(values = c("Fast-Exp" = "#FFC0CB","Fast-Mat" = "#578FCA",
                               "FMOU" = "#FF8989","PCA" = "#fee440","DMD" = "#66D2CE")) + 
  theme_bw() +
  labs(title = "(A) Cell nuclei" ) +
  theme(axis.title.x = element_blank(),
        #axis.text.x   = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12))
nuclei_plot_rmse

####### Figure 8(A): obs mean, obs, pred mean
input1_nuclei = 1:nrow(nuclei_figure)
input2_nuclei = 1:ncol(nuclei_figure)
my_palette <- colorRampPalette(c("#8ecae6","#219ebc","#023047","#ffb703","#fb8500"))(100)
nuclei_image_nuclei_df <- expand.grid(input1_nuclei, input2_nuclei)
nuclei_image_nuclei_df$z <- as.vector(nuclei_figure)
names(nuclei_image_nuclei_df) = c("input1", "input2","z")

nuclei <- ggplot(nuclei_image_nuclei_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() + 
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(A) Observation mean", x = "", y = "", fill = "") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())

noisy_nuclei_image_nuclei_df <- expand.grid(input1_nuclei, input2_nuclei)
noisy_nuclei_image_nuclei_df$z <- as.vector(nuclei_y_record[[1]])
names(noisy_nuclei_image_nuclei_df) = c("input1", "input2","z")

noisy_nuclei <- ggplot(noisy_nuclei_image_nuclei_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() +  # Use geom_tile() or geom_raster()
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(B) Noisy observation", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())
noisy_nuclei

lattice_matern_nuclei_image_df <- expand.grid(input1_nuclei, input2_nuclei)
lattice_matern_nuclei_image_df$z <- as.vector(nuclei_pred_mean_lattice_matern_record[[1]])
names(lattice_matern_nuclei_image_df) = c("input1", "input2","z")


pred_mean_lattice_matern_nuclei <- ggplot(lattice_matern_nuclei_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() + 
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(C) Predictive mean", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())
pred_mean_lattice_matern_nuclei

nuclei_merge <- grid.arrange(nuclei,noisy_nuclei,pred_mean_lattice_matern_nuclei,nrow=1,
                        top = textGrob("Cell nuclei", gp = gpar(fontsize = 12, fontface = "bold"),
                                       vjust=0.2))
#ggsave("signal_obs_pred_nuclei.pdf",nuclei_merge, width=7.5, height=2.5)

######## whole cell ------------------------------------
wholeCell_rmse_lattice_exp <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
wholeCell_rmse_lattice_matern <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
wholeCell_rmse_fmou <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))
wholeCell_rmse_pca <- matrix(NA,  nrow=num_repetition, ncol=length(sigma0_list))
wholeCell_rmse_dmd <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))


wholeCell_y_record <- as.list(1:length(sigma0_list))
wholeCell_pred_mean_lattice_exp_record <- as.list(1:length(sigma0_list))
wholeCell_pred_mean_lattice_matern_record <- as.list(1:length(sigma0_list))
wholeCell_pred_mean_fmou_record <- as.list(1:length(sigma0_list))
wholeCell_pred_mean_fmatern_record <- as.list(1:length(sigma0_list))
wholeCell_pred_mean_pca_record <- as.list(1:length(sigma0_list))
wholeCell_pred_mean_dmd_record <- as.list(1:length(sigma0_list))


wholeCell_est_d_record <- matrix(NA, nrow=num_repetition, ncol=length(sigma0_list))

#### Start 
for(i_sigma0 in 1:length(sigma0_list)){
  sigma_0 = sigma0_list[i_sigma0]
  for(iter in 1:num_repetition){
    set.seed(iter)
    print(iter)
    
    #### Generate noisy observations ---------------
    output_mat=whole_cell_figure + matrix(rnorm(length(whole_cell_figure),sd=sigma_0), nrow(whole_cell_figure), ncol(whole_cell_figure))
    if(iter==1){
      wholeCell_y_record[[i_sigma0]] = output_mat
    }
    
    #### Lattice - exp kernel -------------------------------------------
    input1_here = 1:nrow(whole_cell_figure)
    input2_here = 1:ncol(whole_cell_figure)
    param_ini_here=c(-2,-2,2)
    kernel_type_here="exp"
    optim_method_here = "Nelder-Mead" # Nelder-Mead
    fit_lattice_exp <- lattice_alg(output_mat, input1_here, input2_here, kernel_type=kernel_type_here,
                                   testing_input1=input1_here, testing_input2=input2_here,
                                   param_ini = param_ini_here, optim_method = optim_method_here)
    pred_mean_lattice_exp = fit_lattice_exp$pred_mean
    wholeCell_rmse_lattice_exp[iter, i_sigma0] = sqrt(mean((whole_cell_figure-pred_mean_lattice_exp)^2))
    if(iter==1){
      wholeCell_pred_mean_lattice_exp_record[[i_sigma0]] = pred_mean_lattice_exp
    }
    
    #### Lattice - matern kernel -------------------------------------------
    input1_here = 1:nrow(whole_cell_figure)
    input2_here = 1:ncol(whole_cell_figure)
    param_ini_here=c(-2,-2,2)
    kernel_type_here="matern"
    optim_method_here = "Nelder-Mead" # Nelder-Mead
    fit_lattice_matern <- lattice_alg(output_mat, input1_here, input2_here, kernel_type=kernel_type_here,
                                      testing_input1=input1_here, testing_input2=input2_here,
                                      param_ini = param_ini_here, optim_method = optim_method_here)
    pred_mean_lattice_matern = fit_lattice_matern$pred_mean
    wholeCell_rmse_lattice_matern[iter, i_sigma0] = sqrt(mean((whole_cell_figure-pred_mean_lattice_matern)^2))
    if(iter==1){
      wholeCell_pred_mean_lattice_matern_record[[i_sigma0]] = pred_mean_lattice_matern
    }
    
    #### choose number of latent states -------------------------------------------
    svd_output=svd(output_mat)
    U = svd_output$u
    k=nrow(whole_cell_figure)
    n=ncol(whole_cell_figure)
    loss_score = NULL
    for(i_d in 1:ceiling(dim(output_mat)[1]*2/3)){
      criteria_val_cur = log(mean((output_mat - U[,1:i_d]%*%t(U[,1:i_d])%*%output_mat)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score = c(loss_score,  criteria_val_cur)
    }
    est_d_here = which.min(loss_score)
    wholeCell_est_d_record[iter, i_sigma0] = est_d_here
    
    #### FMOU -------------------------------------------
    m_fmou <- fmou(output_mat, d=est_d_here)
    fit_fmou <- fit.fmou(m_fmou)
    pred_mean_fmou = fit_fmou$mean_obs
    wholeCell_rmse_fmou[iter, i_sigma0] = sqrt(mean((whole_cell_figure-pred_mean_fmou)^2))
    if(iter==1){
      wholeCell_pred_mean_fmou_record[[i_sigma0]] = pred_mean_fmou
    }
    
    
    #### PCA -------------------------------------------
    pred_mean_pca = U[,1:est_d_here] %*% t(U[,1:est_d_here]) %*% output_mat
    wholeCell_rmse_pca[iter, i_sigma0] = sqrt(mean((whole_cell_figure-pred_mean_pca)^2))
    if(iter==1){
      wholeCell_pred_mean_pca_record[[i_sigma0]] = pred_mean_pca
    }
    
    #### DMD, fix r -------------------------------------------
    fit_dmd_fix_r <- DMD_alg(output_mat,r=est_d_here, fix_r=T) 
    wholeCell_rmse_dmd[iter, i_sigma0] = sqrt(mean((cbind(output_mat[,1], fit_dmd_fix_r$in_sample_pred) - whole_cell_figure)^2))
    if(iter==1){
      wholeCell_pred_mean_dmd_record[[i_sigma0]] = cbind(output_mat[,1], fit_dmd_fix_r$in_sample_pred)
    }
  }
  
}


wholeCell_rmse_summary = rbind(colMeans(wholeCell_rmse_lattice_exp),colMeans(wholeCell_rmse_lattice_matern),
                               colMeans(wholeCell_rmse_fmou), colMeans(wholeCell_rmse_pca), colMeans(wholeCell_rmse_dmd))
rownames(wholeCell_rmse_summary) = c("lattice_exp","lattice_matern","fmou","PCA", "DMD")
colnames(wholeCell_rmse_summary) = sigma0_list
wholeCell_rmse_summary

# > wholeCell_rmse_summary
#                    0.1        0.3        0.5
# lattice_exp    0.03455636 0.06159223 0.07894085
# lattice_matern 0.02956700 0.05534026 0.07374805
# fmou           0.06621916 0.15196315 0.15293420
# PCA            0.06887994 0.15231391 0.15406850
# DMD            0.07162848 0.15287134 0.15539062


######## Figure 7(B): violin plot of RMSE (whole cell)

wholeCell_rmse_boxplot_df <- data.frame(
  RMSE = c(wholeCell_rmse_lattice_exp,wholeCell_rmse_lattice_matern,wholeCell_rmse_fmou,
           wholeCell_rmse_pca,wholeCell_rmse_dmd),
  Methods = rep(c("Fast-Exp","Fast-Mat","FMOU","PCA","DMD") ,each=30),
  noise_level = rep(rep(sigma0_list, each=10),5)
)
wholeCell_rmse_boxplot_df$noise_level = factor(wholeCell_rmse_boxplot_df$noise_level, 
                                               labels = c(expression(sigma[0] == 0.1), expression(sigma[0] == 0.3), expression(sigma[0] == 0.5)))
wholeCell_rmse_boxplot_df$Methods = factor(wholeCell_rmse_boxplot_df$Methods,
                                           levels = c("Fast-Mat","PCA","FMOU","DMD","Fast-Exp"))


wholeCell_plot_rmse <- ggplot(wholeCell_rmse_boxplot_df, aes(x=Methods, y=RMSE, fill=Methods)) +
  geom_violin() +
  facet_wrap(~noise_level, scales = "free_y",labeller = label_parsed) +
  scale_fill_manual(values = c("Fast-Exp" = "#FFC0CB","Fast-Mat" = "#578FCA",
                               "FMOU" = "#FF8989","PCA" = "#fee440","DMD" = "#66D2CE")) + 
  theme_bw() +
  labs(title = "(B) Whole cell" ) +
  theme(axis.title.x = element_blank(),
        #axis.text.x   = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12))
wholeCell_plot_rmse
#ggsave("rmse_wholecell.pdf",wholeCell_plot_rmse ,width=9, height=2)

######## Figure 8(B): whole cell, noisy whole cells, pred mean by lattice matern
input1_wholeCell = 1:nrow(whole_cell_figure)
input2_wholeCell = 1:ncol(whole_cell_figure)

my_palette <- colorRampPalette(c("#8ecae6","#219ebc","#023047","#ffb703","#fb8500"))(100)
wholeCell_image_df <- expand.grid(input1_wholeCell, input2_wholeCell)
wholeCell_image_df$z <- as.vector(whole_cell_figure)
names(wholeCell_image_df) = c("input1", "input2","z")

wholeCell <- ggplot(wholeCell_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() + 
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(D) Observation mean", x = "", y = "", fill = "") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())
wholeCell

noisy_wholeCell_image_df <- expand.grid(input1_wholeCell, input2_wholeCell)
noisy_wholeCell_image_df$z <- as.vector(wholeCell_y_record[[3]])
names(noisy_wholeCell_image_df) = c("input1", "input2","z")

noisy_wholeCell <- ggplot(noisy_wholeCell_image_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() + 
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(E) Noisy observation", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())
noisy_wholeCell

lattice_matern_wholeCell_df <- expand.grid(input1_wholeCell, input2_here)
lattice_matern_wholeCell_df$z <- as.vector(wholeCell_pred_mean_lattice_matern_record[[3]])
names(lattice_matern_wholeCell_df) = c("input1", "input2","z")

pred_mean_lattice_matern_wholeCell <- ggplot(lattice_matern_wholeCell_df, aes(x = input1, y = input2, fill = z)) +
  geom_tile() +  
  scale_fill_gradientn(colors=my_palette) +
  labs(title = "(F) Predictive mean", x = "X", y = "Y", fill = "Z") +
  guides(fill = "none") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank())
pred_mean_lattice_matern_wholeCell

f_merge <- grid.arrange(wholeCell,noisy_wholeCell,pred_mean_lattice_matern_wholeCell,nrow=1,
                        top = textGrob("Whole cell", gp = gpar(fontsize = 12, fontface = "bold"),
                                       vjust=0.2))
#ggsave("signal_obs_pred_wholecell.pdf",f_merge, width=7.5, height=2.5)

