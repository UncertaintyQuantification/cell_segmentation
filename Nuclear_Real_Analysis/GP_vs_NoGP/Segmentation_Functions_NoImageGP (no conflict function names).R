#Segmentation with no GP

threshold_image2 <- function(mat, percentage, count = TRUE) {
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

criterion_1_2 <- function(predmean_mat, delta = 0.01, nugget = T) {
  
  # Define percentage thresholds
  percentages <- seq(0, 1, by = delta)
  
  # Calculate pixel counts based on thresholding the image
  pixel_counts <- sapply(percentages, function(threshold) {
    threshold_image2(mat = predmean_mat, percentage = threshold, count = TRUE)
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
    thresholded_image <- threshold_image2(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
  }
  
  # Return the thresholded image and pixel counts
  return(list(
    thresholded_image = thresholded_image,
    pixel_counts = pixel_counts,
    diff_pixel_counts = diff_pixel_counts,
    estimated_percentage = estimated_percentage
  ))
}

eliminate_small_areas2 <- function(GP_masks) {
  unique_labels <- unique(GP_masks[GP_masks > 0])
  label_counts <- table(as.vector(GP_masks))[names(table(as.vector(GP_masks))) != 0]
  filtered_mask <- GP_masks
  nrow_mask <- nrow(GP_masks)
  ncol_mask <- ncol(GP_masks)
  
  mean_obj_size <- mean(label_counts)
  
  for (label in unique_labels) {
    label_mask <- (GP_masks == label)
    area <- sum(label_mask)
    on_boundary <- any(label_mask[1, ]) || any(label_mask[nrow_mask, ]) || 
      any(label_mask[, 1]) || any(label_mask[, ncol_mask])
    
    if (area < mean_obj_size*0.15 && !on_boundary) {
      # Remove the label by setting it to 0 (background)
      filtered_mask[label_mask] <- 0
    }
    
    if (on_boundary && area < mean_obj_size*0.05) {
      filtered_mask[label_mask] <- 0
    }
  }
  
  return(filtered_mask)
}
generate_GP_Masks2 <- function(file_path, delta = 0.01, 
                              nugget = T) {
  
  img <- image_read(file_path)
  img_info <- image_info(img)
  img_width <- img_info$width
  img_height <- img_info$height
  
  # Dynamically determine row and column proportions
  row_proportion <- get_proportion2(img_height)
  col_proportion <- get_proportion2(img_width)
  
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
      
      predmean_mat <- img_matrix
      
      # Store processed images and original images
      processed_images[[count]] <- predmean_mat
      ori_images[[count]] <- img_matrix
      
      # Apply Criterion 1 thresholding to the predicted mean matrix
      criterion_1_info <- criterion_1_2(predmean_mat, delta, nugget)
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
    rethresholded_image <- threshold_image2(mat = processed_images[[outlier_index]], 
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
  segmented_image <- EBImage::watershed(dist_map)
  GP_masks_raw <- segmented_image@.Data
  GP_masks <- eliminate_small_areas2(GP_masks_raw)
  
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


find_boundaries2 <- function(mask) {
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
get_proportion2 <- function(size, target_min = 200, target_max = 400) {
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

generate_GP_Masks_test2 <- function(file_path, delta = 0.01,
                                   nugget = T) {
  
  img <- image_read(file_path)
  img_info <- image_info(img)
  img_width <- img_info$width
  img_height <- img_info$height
  
  # Dynamically determine row and column proportions
  row_proportion <- get_proportion2(img_height)
  col_proportion <- get_proportion2(img_width)
  
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
      
      
      predmean_mat <- img_matrix
      
      # Store processed images and original images
      processed_images[[count]] <- predmean_mat
      ori_images[[count]] <- img_matrix
      
      # Apply Criterion 1 thresholding to the predicted mean matrix
      criterion_1_info <- criterion_1_2(predmean_mat, delta, nugget)
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
    rethresholded_image <- threshold_image2(mat = processed_images[[outlier_index]], 
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
  segmented_image <- EBImage::watershed(dist_map)
  GP_masks_raw <- segmented_image@.Data
  GP_masks <- eliminate_small_areas2(GP_masks_raw)
  
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