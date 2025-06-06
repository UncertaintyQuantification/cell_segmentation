#Modified Functions - Criterion Smoothing by RobustGaSP Default

separable_GP_param_est <- function(output_mat){
  n1=dim(output_mat)[1]
  n2=dim(output_mat)[2]
  N=n1*n2
  p=2 ##2D input
  set.seed(1)
  
  input1=as.numeric(seq(0,1,1/(n1-1)))
  input2=as.numeric(seq(0,1,1/(n2-1)))
  
  input=cbind(rep(input1,n2),as.vector(t(matrix(input2,n2,n1))))  
  
  gradient_matrix <- gradient(as.matrix(output_mat))
  
  
  ##the range here
  range=c(1,1)
  ##1.numerical gradient estimation
  h1=range[1]/(n1-1)
  h2=range[2]/(n2-1)
  
  output_numer_grad1=matrix(NA,n1,n2)
  output_numer_grad2=matrix(NA,n1,n2)
  
  output_numer_grad1[1,]=(output_mat[2,]-output_mat[1,])/h1
  output_numer_grad1[n1,]=(output_mat[n1,]-output_mat[n1-1,])/h1
  output_numer_grad1[2:(n1-1),]=(output_mat[3:(n1),]-output_mat[1:(n1-2),])/(2*h1)
  
  output_numer_grad2[,1]=(output_mat[,2]-output_mat[,1])/h2
  output_numer_grad2[,n2]=(output_mat[,n2]-output_mat[,n2-1])/h2
  output_numer_grad2[,2:(n2-1)]=(output_mat[,3:n2]-output_mat[,1:(n2-2)])/(2*h2)
  
  output_numer_grad_magnitude=sqrt(output_numer_grad1^2+output_numer_grad2^2)
  
  # image2D(output_numer_grad_magnitude)
  ##2.use numerical grad of post mean of separate model 
  Matern_5_2_funct<-function(d,beta){
    x=sqrt(5)*beta*d
    (1+x+x^2/3)*exp(-x)
  }
  
  Neg_log_lik_eigen_with_nugget <- function(param) {
    ##do not confuse this beta
    beta= exp(param[1:p])
    nu=exp(param[p+1])
    #R1=Exp_funct(R01, beta=beta[1])
    #R2=Exp_funct(R02, beta=beta[2])
    # R1=double_Exp_funct(R01, beta=beta[1])
    # R2=double_Exp_funct(R02, beta=beta[2])
    # R1=Matern_5_2_funct(R01, beta=beta[1])
    # R2=Matern_5_2_funct(R02, beta=beta[2])
    R1=Matern_5_2_funct(R01, beta=beta[1])
    R2=Matern_5_2_funct(R02, beta=beta[2])
    
    eigen_R1=eigen(R1)  
    eigen_R2=eigen(R2)  
    
    U_x=matrix(NA,N,q_X)
    for(i_q in 1:q_X){
      U_x=as.vector(t(eigen_R1$vectors)%*% X_list[[i_q]]%*%eigen_R2$vectors)
    }
    
    Lambda_tilde_inv=1/(kronecker(eigen_R2$values,eigen_R1$values, FUN = "*")+nu)
    #eigen(R_tilde.inv)$values-sort(Lambda_tilde_inv,decreasing=T)
    
    Lambda_tilde_inv_U_x= Lambda_tilde_inv*U_x
    
    X_R_tilde_inv_X_inv=solve(t(U_x)%*%(Lambda_tilde_inv_U_x))
    output_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat%*%eigen_R2$vectors)
    
    theta_hat=X_R_tilde_inv_X_inv%*%(t(Lambda_tilde_inv_U_x)%*%output_tilde)
    
    output_mat_normalized=matrix(as.vector(output_mat)-X%*%theta_hat,n1,n2)
    
    output_normalize_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat_normalized%*%eigen_R2$vectors)
    
    S_2=sum((output_normalize_tilde)*Lambda_tilde_inv*output_normalize_tilde  )  
    
    -(1/2*sum(log(Lambda_tilde_inv))-N/2*log(S_2))
    
  }
  
  
  
  ##mean basis
  X=matrix(1,N,1)  ##use constant mean basis
  q_X=dim(X)[2]
  X_list=as.list(1:q_X)
  for(i_q in 1:q_X){
    X_list[[i_q]]=matrix(X[,i_q],n1,n2)
  }
  ###distance
  R01=as.matrix(abs(outer(input1, input1, "-")))
  R02=as.matrix(abs(outer(input2, input2, "-")))
  
  
  ##initial value 
  param_ini=c(-2,-2,-3)
  Neg_log_lik_eigen_with_nugget(param_ini)
  
  
  
  ##add derivative may make it faster and more stable, but let's use numerical gradient for now
  m_eigen=try(optim(param_ini,method ="L-BFGS-B",#lower = c(-Inf,-Inf,-Inf),
                    Neg_log_lik_eigen_with_nugget),silent=T) ##it will be more accurate to use derivative
  while(!is.numeric(m_eigen[[1]])){
    m_eigen=try(optim(param_ini+runif(3),method ="L-BFGS-B",#lower = c(-Inf,-Inf,-Inf),
                      Neg_log_lik_eigen_with_nugget),silent=T) ##it will be more accurate to use derivative
    
  }
  
  
  
  ##build predictive mean
  
  beta= exp(m_eigen$par[1:p])
  nu=exp(m_eigen$par[p+1])
  
  list(param = c(beta,nu))
}
separable_GP <- function(output_mat, parameters){
  n1=dim(output_mat)[1]
  n2=dim(output_mat)[2]
  N=n1*n2
  p=2 ##2D input
  set.seed(1)
  
  input1=as.numeric(seq(0,1,1/(n1-1)))
  input2=as.numeric(seq(0,1,1/(n2-1)))
  
  input=cbind(rep(input1,n2),as.vector(t(matrix(input2,n2,n1))))  
  
  gradient_matrix <- gradient(as.matrix(output_mat))
  
  
  ##the range here
  range=c(1,1)
  ##1.numerical gradient estimation
  h1=range[1]/(n1-1)
  h2=range[2]/(n2-1)
  
  output_numer_grad1=matrix(NA,n1,n2)
  output_numer_grad2=matrix(NA,n1,n2)
  
  output_numer_grad1[1,]=(output_mat[2,]-output_mat[1,])/h1
  output_numer_grad1[n1,]=(output_mat[n1,]-output_mat[n1-1,])/h1
  output_numer_grad1[2:(n1-1),]=(output_mat[3:(n1),]-output_mat[1:(n1-2),])/(2*h1)
  
  output_numer_grad2[,1]=(output_mat[,2]-output_mat[,1])/h2
  output_numer_grad2[,n2]=(output_mat[,n2]-output_mat[,n2-1])/h2
  output_numer_grad2[,2:(n2-1)]=(output_mat[,3:n2]-output_mat[,1:(n2-2)])/(2*h2)
  
  output_numer_grad_magnitude=sqrt(output_numer_grad1^2+output_numer_grad2^2)
  
  Matern_5_2_funct<-function(d,beta){
    x=sqrt(5)*beta*d
    (1+x+x^2/3)*exp(-x)
  }
  
  Neg_log_lik_eigen_with_nugget <- function(param) {
    ##do not confuse this beta
    beta= exp(param[1:p])
    nu=exp(param[p+1])
    # R1=double_Exp_funct(R01, beta=beta[1])
    # R2=double_Exp_funct(R02, beta=beta[2])
    # R1=Matern_5_2_funct(R01, beta=beta[1])
    # R2=Matern_5_2_funct(R02, beta=beta[2])
    R1=Matern_5_2_funct(R01, beta=beta[1])
    R2=Matern_5_2_funct(R02, beta=beta[2])
    
    eigen_R1=eigen(R1)  
    eigen_R2=eigen(R2)  
    
    U_x=matrix(NA,N,q_X)
    for(i_q in 1:q_X){
      U_x=as.vector(t(eigen_R1$vectors)%*% X_list[[i_q]]%*%eigen_R2$vectors)
    }
    
    Lambda_tilde_inv=1/(kronecker(eigen_R2$values,eigen_R1$values, FUN = "*")+nu)
    #eigen(R_tilde.inv)$values-sort(Lambda_tilde_inv,decreasing=T)
    
    Lambda_tilde_inv_U_x= Lambda_tilde_inv*U_x
    
    X_R_tilde_inv_X_inv=solve(t(U_x)%*%(Lambda_tilde_inv_U_x))
    output_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat%*%eigen_R2$vectors)
    
    theta_hat=X_R_tilde_inv_X_inv%*%(t(Lambda_tilde_inv_U_x)%*%output_tilde)
    
    output_mat_normalized=matrix(as.vector(output_mat)-X%*%theta_hat,n1,n2)
    
    output_normalize_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat_normalized%*%eigen_R2$vectors)
    
    S_2=sum((output_normalize_tilde)*Lambda_tilde_inv*output_normalize_tilde  )  
    
    -(1/2*sum(log(Lambda_tilde_inv))-N/2*log(S_2))
    
  }
  
  
  
  ##mean basis
  X=matrix(1,N,1)  ##use constant mean basis
  q_X=dim(X)[2]
  X_list=as.list(1:q_X)
  for(i_q in 1:q_X){
    X_list[[i_q]]=matrix(X[,i_q],n1,n2)
  }
  ###distance
  R01=as.matrix(abs(outer(input1, input1, "-")))
  R02=as.matrix(abs(outer(input2, input2, "-")))
  
  
  ##initial value 
  param_ini=c(-2,-2,-3)
  Neg_log_lik_eigen_with_nugget(param_ini)
  
  
  beta = parameters[1:2] # adjust this
  nu = parameters[3]
  
  #R1=Exp_funct(R01, beta=beta[1])
  #R2=Exp_funct(R02, beta=beta[2])
  R1=Matern_5_2_funct(R01, beta=beta[1])
  R2=Matern_5_2_funct(R02, beta=beta[2])
  
  eigen_R1=eigen(R1)  
  eigen_R2=eigen(R2)  
  ###get theta
  U_x=matrix(NA,N,q_X)
  for(i_q in 1:q_X){
    U_x=as.vector(t(eigen_R1$vectors)%*% X_list[[i_q]]%*%eigen_R2$vectors)
  }
  
  Lambda_tilde_inv=1/(kronecker(eigen_R2$values,eigen_R1$values, FUN = "*")+nu)
  #eigen(R_tilde.inv)$values-sort(Lambda_tilde_inv,decreasing=T)
  
  Lambda_tilde_inv_U_x= Lambda_tilde_inv*U_x
  
  X_R_tilde_inv_X_inv=solve(t(U_x)%*%(Lambda_tilde_inv_U_x))
  output_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat%*%eigen_R2$vectors)
  
  
  theta_hat=X_R_tilde_inv_X_inv%*%(t(Lambda_tilde_inv_U_x)%*%output_tilde)
  
  
  output_mat_normalized=matrix(as.vector(output_mat)-X%*%theta_hat,n1,n2)
  
  ### test on input,it will be faster if we use filling
  testing_input=input 
  testing_input1=input1
  testing_input2=input2
  
  X_testing=X
  
  r01=abs(outer(input1,testing_input1,'-'))
  r02=abs(outer(input2,testing_input2,'-'))
  r1=Matern_5_2_funct(r01, beta[1])
  r2=Matern_5_2_funct(r02, beta[2])
  
  output_normalize_tilde=as.vector(t(eigen_R1$vectors)%*% output_mat_normalized%*%eigen_R2$vectors)
  output_normalized_tilde_lambda_inv_mat=matrix(Lambda_tilde_inv*output_normalize_tilde,n1,n2)
  R_tilde_inv_output_normalize=as.vector((eigen_R1$vectors)%*% output_normalized_tilde_lambda_inv_mat%*%t(eigen_R2$vectors))
  
  R_tilde_inv_output_normalize_mat=matrix(R_tilde_inv_output_normalize,n1,n2)
  predmean=X_testing%*%theta_hat+as.vector(t(r1)%*%R_tilde_inv_output_normalize_mat%*%r2)
  predmean_mat=matrix(predmean,n1,n2)
  
  
  z_lim=c(min(output_mat,predmean_mat),max(output_mat,predmean_mat) )
  
  #image2D(predmean_mat,x=input1,y=input2,zlim=z_lim,main='GP prediction')
  
  ###build gradient estimate; use KF will save some memory
  ###here for saving the memory we use some same names 
  delta1=10^{-5}
  delta2=10^{-5}
  
  ##plus for 1
  testing_input1_delta=input1+rep(delta1,n1) #t(matrix(c(delta1,0),2,N))
  testing_input2_delta=input2 #t(matrix(c(delta1,0),2,N))
  
  r01_delta=abs(outer(input1,testing_input1_delta,'-'))
  r02_delta=abs(outer(input2,testing_input2_delta,'-'))
  
  r1_delta=Matern_5_2_funct(r01_delta, beta[1])
  r2_delta=Matern_5_2_funct(r02_delta, beta[2])
  
  predmean_plus_1=X_testing%*%theta_hat+as.vector(t(r1_delta)%*%R_tilde_inv_output_normalize_mat%*%r2_delta)
  
  ##minus for 1
  testing_input1_delta=input1-rep(delta1,n1) #t(matrix(c(delta1,0),2,N))
  testing_input2_delta=input2 #t(matrix(c(delta1,0),2,N))
  
  r01_delta=abs(outer(input1,testing_input1_delta,'-'))
  r02_delta=abs(outer(input2,testing_input2_delta,'-'))
  
  r1_delta=Matern_5_2_funct(r01_delta, beta[1])
  r2_delta=Matern_5_2_funct(r02_delta, beta[2])
  
  predmean_minus_1=X_testing%*%theta_hat+as.vector(t(r1_delta)%*%R_tilde_inv_output_normalize_mat%*%r2_delta)
  
  
  predmean_numer_grad1=(predmean_plus_1-predmean_minus_1)/(2*delta1)
  predmean_numer_grad1_mat=matrix(predmean_numer_grad1,n1,n2)
  #image2D(predmean_numer_grad1_mat)
  ##build grad for 2
  ##plus for 2
  testing_input1_delta=input1 #t(matrix(c(delta1,0),2,N))
  testing_input2_delta=input2+rep(delta2,n2) #t(matrix(c(delta1,0),2,N))
  
  r01_delta=abs(outer(input1,testing_input1_delta,'-'))
  r02_delta=abs(outer(input2,testing_input2_delta,'-'))
  
  r1_delta=Matern_5_2_funct(r01_delta, beta[1])
  r2_delta=Matern_5_2_funct(r02_delta, beta[2])
  
  predmean_plus_2=X_testing%*%theta_hat+as.vector(t(r1_delta)%*%R_tilde_inv_output_normalize_mat%*%r2_delta)
  
  
  ##minus for 2
  testing_input1_delta=input1 #t(matrix(c(delta1,0),2,N))
  testing_input2_delta=input2-rep(delta2,n2) #t(matrix(c(delta1,0),2,N))
  
  r01_delta=abs(outer(input1,testing_input1_delta,'-'))
  r02_delta=abs(outer(input2,testing_input2_delta,'-'))
  
  r1_delta=Matern_5_2_funct(r01_delta, beta[1])
  r2_delta=Matern_5_2_funct(r02_delta, beta[2])
  predmean_minus_2=X_testing%*%theta_hat+as.vector(t(r1_delta)%*%R_tilde_inv_output_normalize_mat%*%r2_delta)
  
  
  predmean_numer_grad2=(predmean_plus_2-predmean_minus_2)/(2*delta2)
  predmean_numer_grad2_mat=matrix(predmean_numer_grad2,n1,n2)
  #image2D(predmean_numer_grad2_mat)
  
  predmean_numer_grad_magnitude_mat=sqrt(predmean_numer_grad1_mat^2+predmean_numer_grad2_mat^2)
  #image2D(predmean_numer_grad_magnitude_mat)
  list(
    predmean_mat = predmean_mat,
    grad1 = predmean_numer_grad1_mat,
    grad2 = predmean_numer_grad2_mat,
    grad_magnitude = predmean_numer_grad_magnitude_mat,
    param = c(beta,nu)
  )
}

threshold_image <- function(mat, percentage, count = TRUE) {
  max_value <- max(mat, na.rm = TRUE)
  threshold_value <- percentage * max_value
  thresholded_mat <- ifelse(mat > threshold_value, 1, 0)
  
  if (count) {
    count_mat <- sum(thresholded_mat)
    return(count_mat)
  } else {
    return(thresholded_mat)
  }
}

criterion_1 <- function(predmean_mat, delta = 0.01, nugget = T) {
  
  # Define percentage thresholds
  percentages <- seq(0, 1, by = delta)
  
  # Calculate pixel counts based on thresholding the image
  pixel_counts <- sapply(percentages, function(threshold) {
    threshold_image(mat = predmean_mat, percentage = threshold, count = TRUE)
  })
  
  # Compute the absolute differences in pixel counts
  diff_pixel_counts <- abs(diff(pixel_counts))
  
  
  # Gaussian kernel smoothing 
  #smoothed <- ksmooth(percentages[-1], diff_pixel_counts , kernel = "normal", bandwidth = 2*delta)  # Adjust bandwidth as needed
  
  #Robust GaSP - default
  diff_mod <- rgasp(percentages[-1], diff_pixel_counts, nugget.est = nugget)
  smoothed <- predict(diff_mod, testing_input=as.matrix(percentages[-1]))
  
  #### Using changes in Diff
  diff_pixel_counts<- smoothed$mean
  max_index <- which.max(diff_pixel_counts)
  th<- 0.05*sd(diff_pixel_counts)
  stable_index <- percentages[length(percentages)]
  found_stable <- F
  for (i in (max_index + 1):length(diff_pixel_counts)) {
    # Check if the absolute difference between successive points is less than the threshold
    if (abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1]) < th) {
      stable_index <- i
      found_stable <- T
      break
    }
  }
  
  #default if stable threshold is failed to be found
  #make everything background and have the estimated threshold be the max
  thresholded_image <- matrix(0, nrow=nrow(predmean_mat), ncol=ncol(predmean_mat))
  estimated_percentage <- percentages[length(percentages)]
  
  if (found_stable) {
    estimated_percentage <- percentages[stable_index+1]
    thresholded_image <- threshold_image(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
  }
  
  # Return the thresholded image and pixel counts
  return(list(
    thresholded_image = thresholded_image,
    pixel_counts = pixel_counts,
    diff_pixel_counts = diff_pixel_counts,
    estimated_percentage = estimated_percentage
  ))
}

eliminate_small_areas <- function(GP_masks, size_threshold) {
  unique_labels <- unique(GP_masks[GP_masks > 0])
  filtered_mask <- GP_masks
  nrow_mask <- nrow(GP_masks)
  ncol_mask <- ncol(GP_masks)
  
  for (label in unique_labels) {
    label_mask <- (GP_masks == label)
    area <- sum(label_mask)
    on_boundary <- any(label_mask[1, ]) || any(label_mask[nrow_mask, ]) || 
      any(label_mask[, 1]) || any(label_mask[, ncol_mask])
    
    if (area < size_threshold && !on_boundary) {
      # Remove the label by setting it to 0 (background)
      filtered_mask[label_mask] <- 0
    }
    
    if (on_boundary && area < size_threshold/5) {
      filtered_mask[label_mask] <- 0
    }
  }
  
  return(filtered_mask)
}
generate_GP_Masks <- function(file_path, delta = 0.01, remove_size_threshold = 50, 
                              nugget = T) {
  
  img <- image_read(file_path)
  img_info <- image_info(img)
  img_width <- img_info$width
  img_height <- img_info$height
  
  # Dynamically determine row and column proportions
  row_proportion <- get_proportion(img_height)
  col_proportion <- get_proportion(img_width)
  
  # Calculate the size of each cropped piece based on the determined proportions
  crop_width <- img_width * col_proportion
  crop_height <- img_height * row_proportion
  
  # Determine the number of pieces in each dimension
  num_pieces_x <- ceiling(1 / col_proportion)
  num_pieces_y <- ceiling(1 / row_proportion)
  
  ori_images <- list()
  processed_images <- list()
  thresholded1_images <- list()
  crit_1_opt_thresholds <- list()
  connected_parts_count <- list()
  
  parameters <- NULL  # Parameters estimated for the first cropped piece
  count <- 1
  
  # Process each piece individually without combining them
  for (i in 1:num_pieces_x) {
    for (j in 1:num_pieces_y) {
      x_offset <- (i - 1) * crop_width
      y_offset <- (j - 1) * crop_height
      
      # Crop the image using the calculated dimensions and offsets
      cropped_img <- image_crop(img, geometry_area(crop_width, crop_height, x_offset, y_offset))
      img_matrix <- as.numeric(cropped_img[[1]])[,,1]
      
      # If this is the first piece, estimate GP parameters
      if (i == 1 && j == 1) {
        parameters <- separable_GP_param_est(img_matrix)
      }
      
      # Apply the GP function on the matrix using the estimated parameters
      separable_GP_info <- separable_GP(img_matrix, parameters$param)
      predmean_mat <- separable_GP_info$predmean_mat
      
      # Store processed images and original images
      processed_images[[count]] <- predmean_mat
      ori_images[[count]] <- img_matrix
      
      # Apply Criterion 1 thresholding to the predicted mean matrix
      criterion_1_info <- criterion_1(predmean_mat, delta, nugget)
      thresholded1_img <- criterion_1_info$thresholded_image
      
      thresholded1_images[[count]] <- thresholded1_img
      
      # Count the number of connected components in the thresholded image
      connected_parts_count[[count]] <- length(unique(as.vector(bwlabel(thresholded1_img))))
      
      # Store the optimal threshold for this piece
      crit_1_opt_thresholds[[count]] <- criterion_1_info$estimated_percentage
      
      count <- count + 1
    }
  }
  
  # Outlier Detection: Identify pieces with anomalous connected parts count
  # Need to be changed?
  mean_connected_parts <- mean(unlist(connected_parts_count))
  sd_connected_parts <- sd(unlist(connected_parts_count))
  outlier_threshold <- 2
  outliers <- which(abs(unlist(connected_parts_count) - mean_connected_parts) > outlier_threshold * sd_connected_parts)
  
  # Reprocess outliers by applying the mean threshold
  for (outlier_index in outliers) {
    rethresholded_image <- threshold_image(mat = processed_images[[outlier_index]], 
                                           percentage = mean(as.numeric(crit_1_opt_thresholds)[-outliers]), 
                                           count = FALSE)
    
    # Replace the original thresholded image with the rethresholded image
    thresholded1_images[[outlier_index]] <- rethresholded_image
    
    #Update outlier threshold to average
    crit_1_opt_thresholds[outlier_index] <- mean(as.numeric(crit_1_opt_thresholds)[-outliers])
  }
  
  # Initialize combined matrices for final image
  combined_predmean <- matrix(0, nrow = img_height, ncol = img_width)
  combined_thresholded1 <- matrix(0, nrow = img_height, ncol = img_width)
  
  # Combine all processed images into the full matrix
  count <- 1
  for (i in 1:num_pieces_x) {
    for (j in 1:num_pieces_y) {
      x_offset <- (i - 1) * crop_width
      y_offset <- (j - 1) * crop_height
      
      # Retrieve the processed images and thresholded images after outlier handling
      predmean_mat <- processed_images[[count]]
      thresholded1_img <- thresholded1_images[[count]]
      
      # Place each piece in the combined matrices
      piece_height <- nrow(predmean_mat)
      piece_width <- ncol(predmean_mat)
      combined_predmean[y_offset + 1:piece_height, x_offset + 1:piece_width] <- predmean_mat
      combined_thresholded1[y_offset + 1:piece_height, x_offset + 1:piece_width] <- thresholded1_img
      
      count <- count + 1
    }
  }
  
  # Create distance map and apply watershed segmentation
  dist_map <- distmap(as.Image(combined_thresholded1))
  segmented_image <- watershed(dist_map)
  GP_masks_raw <- segmented_image@.Data
  GP_masks <- eliminate_small_areas(GP_masks_raw, remove_size_threshold)
  
  # Return both the individual processed pieces and the combined result
  return(list(
    ori_images = ori_images,
    processed_images = processed_images,
    crit_1_opt_thresholds = crit_1_opt_thresholds,
    connected_parts_count = connected_parts_count,
    outliers = outliers,
    combined_predmean = combined_predmean,
    combined_thresholded1 = combined_thresholded1,
    GP_masks = GP_masks
  ))
}


find_boundaries <- function(mask) {
  boundary <- matrix(0, nrow = nrow(mask), ncol = ncol(mask))
  
  for (i in 2:(nrow(mask) - 1)) {
    for (j in 2:(ncol(mask) - 1)) {
      if (mask[i, j] > 0) {
        # Check only the four direct neighbors (up, down, left, right)
        if (mask[i - 1, j] != mask[i, j] || mask[i + 1, j] != mask[i, j] ||
            mask[i, j - 1] != mask[i, j] || mask[i, j + 1] != mask[i, j]) {
          boundary[i, j] <- 1
        }
      }
    }
  }
  
  return(boundary)
}

# Automatically determine row and column proportions based on target piece size (200-400 pixels)
get_proportion <- function(size, target_min = 200, target_max = 400) {
  divisors <- c(1/4, 1/3, 1/2, 1)  
  # Some Common divisors(Assume there is no common figure larger than 1600)
  for (div in divisors) {
    piece_size <- size * div
    if (piece_size >= target_min && piece_size <= target_max) {
      return(div)
    }
  }
  # Default to 0.25 if no suitable divisor is found
  return(1/4)
}

generate_GP_Masks_test <- function(file_path, delta = 0.01, remove_size_threshold = 50,
                                   nugget = T) {
  
  img <- image_read(file_path)
  img_info <- image_info(img)
  img_width <- img_info$width
  img_height <- img_info$height
  
  # Dynamically determine row and column proportions
  row_proportion <- get_proportion(img_height)
  col_proportion <- get_proportion(img_width)
  
  # Calculate integer crop sizes for each piece based on proportions
  crop_width <- as.integer(img_width * col_proportion)
  crop_height <- as.integer(img_height * row_proportion)
  
  # Determine the number of pieces in each dimension
  num_pieces_x <- floor(img_width / crop_width)
  num_pieces_y <- floor(img_height / crop_height)
  
  # Adjust crop dimensions if they don’t add up to the original image dimensions
  crop_width <- img_width %/% num_pieces_x
  crop_height <- img_height %/% num_pieces_y
  
  # Initialize combined matrices for final image
  combined_predmean <- matrix(0, nrow = img_height, ncol = img_width)
  combined_thresholded1 <- matrix(0, nrow = img_height, ncol = img_width)
  
  ori_images <- list()
  processed_images <- list()
  thresholded1_images <- list()
  crit_1_opt_thresholds <- list()
  connected_parts_count <- list()
  
  parameters <- NULL  # Parameters estimated for the first cropped piece
  count <- 1
  
  # Process each piece individually and place directly in the combined matrices
  for (i in 1:num_pieces_x) {
    for (j in 1:num_pieces_y) {
      # Set offsets for cropping
      x_offset <- (i - 1) * crop_width
      y_offset <- (j - 1) * crop_height
      
      # Adjust dimensions for the last piece to fill the remaining pixels
      piece_width <- ifelse(i == num_pieces_x, img_width - x_offset, crop_width)
      piece_height <- ifelse(j == num_pieces_y, img_height - y_offset, crop_height)
      
      # Crop the image using the calculated dimensions and offsets
      cropped_img <- image_crop(img, geometry_area(crop_width, crop_height, x_offset, y_offset))
      img_matrix <- as.numeric(cropped_img[[1]])[,,1]
      
      # If this is the first piece, estimate GP parameters
      if (i == 1 && j == 1) {
        parameters <- separable_GP_param_est(img_matrix)
      }
      
      # Apply the GP function on the matrix using the estimated parameters
      separable_GP_info <- separable_GP(img_matrix, parameters$param)
      predmean_mat <- separable_GP_info$predmean_mat
      
      # Store processed images and original images
      processed_images[[count]] <- predmean_mat
      ori_images[[count]] <- img_matrix
      
      # Apply Criterion 1 thresholding to the predicted mean matrix
      criterion_1_info <- criterion_1(predmean_mat, delta, nugget)
      thresholded1_img <- criterion_1_info$thresholded_image
      thresholded1_images[[count]] <- thresholded1_img
      
      # Count the number of connected components in the thresholded image
      connected_parts_count[[count]] <- length(unique(as.vector(bwlabel(thresholded1_img))))
      
      # Store the optimal threshold for this piece
      crit_1_opt_thresholds[[count]] <- criterion_1_info$estimated_percentage
      
      count <- count + 1
    }
  }
  
  # Outlier Detection: Identify pieces with anomalous connected parts count
  # Need to be changed?
  mean_connected_parts <- mean(unlist(connected_parts_count))
  sd_connected_parts <- sd(unlist(connected_parts_count))
  outlier_threshold <- 2
  outliers <- which(abs(unlist(connected_parts_count) - mean_connected_parts) > outlier_threshold * sd_connected_parts)
  
  # Reprocess outliers by applying the mean threshold
  for (outlier_index in outliers) {
    rethresholded_image <- threshold_image(mat = processed_images[[outlier_index]], 
                                           percentage = mean(as.numeric(crit_1_opt_thresholds)[-outliers]), 
                                           count = FALSE)
    
    #update threshold
    crit_1_opt_thresholds[outlier_index] <- mean(as.numeric(crit_1_opt_thresholds)[-outliers])
    #If vast majority pixels are foreground, revert to background
    if(sum(as.vector(rethresholded_image)) > 0.99*(nrow(rethresholded_image)*ncol(rethresholded_image))) {
      rethresholded_image <- matrix(0, nrow=nrow(rethresholded_image), ncol=ncol(rethresholded_image))
      crit_1_opt_thresholds[outlier_index] <- 1
      
    }
    
    thresholded1_images[[outlier_index]] <- rethresholded_image
  }
  
  # Initialize combined matrices for final image
  combined_predmean <- matrix(0, nrow = img_height, ncol = img_width)
  combined_thresholded1 <- matrix(0, nrow = img_height, ncol = img_width)
  
  # Combine all processed images into the full matrix
  count <- 1
  for (i in 1:num_pieces_x) {
    for (j in 1:num_pieces_y) {
      x_offset <- (i - 1) * crop_width
      y_offset <- (j - 1) * crop_height
      
      # Retrieve the processed images and thresholded images after outlier handling
      predmean_mat <- processed_images[[count]]
      thresholded1_img <- thresholded1_images[[count]]
      
      # Place each piece in the combined matrices
      piece_height <- nrow(predmean_mat)
      piece_width <- ncol(predmean_mat)
      combined_predmean[y_offset + 1:piece_height, x_offset + 1:piece_width] <- predmean_mat
      combined_thresholded1[y_offset + 1:piece_height, x_offset + 1:piece_width] <- thresholded1_img
      
      count <- count + 1
    }
  }
  
  # Create distance map and apply watershed segmentation
  dist_map <- distmap(as.Image(combined_thresholded1))
  segmented_image <- watershed(dist_map)
  GP_masks_raw <- segmented_image@.Data
  GP_masks <- eliminate_small_areas(GP_masks_raw, remove_size_threshold)
  
  # Return both the individual processed pieces and the combined result
  return(list(
    ori_images = ori_images,
    processed_images = processed_images,
    crit_1_opt_thresholds = crit_1_opt_thresholds,
    connected_parts_count = connected_parts_count,
    outliers = outliers,
    combined_predmean = combined_predmean,
    combined_thresholded1 = combined_thresholded1,
    GP_masks = GP_masks
  ))
}