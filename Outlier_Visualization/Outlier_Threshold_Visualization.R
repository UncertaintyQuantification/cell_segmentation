####Outlier Threshold Visualization ####

#### Import Libraries and Functions ####
library(magick)
library(pracma)
library(lhs)
library(RobustGaSP)
library(plot3D)
library(dplyr)
library(ggplot2)
library(EBImage)
#Ensure directory is set to file folder
source("../src/Modified_Functions_RGasp.R")

#Image 5 has outlier panels
file_path5 <- "../Image_Data/cropped_pieces/cropped_img_5.jpg"
img5 <- image_read(file_path5)
ori_img_matrix5 <- as.numeric(img5[[1]])[,,1]
image2D(ori_img_matrix5, zlim=c(0,1))

#### Use Separable Kernel GP and Cropped Images to generate GP Masks ####
gp_masks_result <- generate_GP_Masks_test(file_path5)
GP_masks <- gp_masks_result$GP_masks
image2D(GP_masks, main = "Cropped Image 5")
gp_masks_result$outliers

#### Plotting Threshold for Outlier Sub-Images ####
test <- gp_masks_result$combined_predmean

delta = 0.01 #change
percentages <- seq(0, 1, by = delta) #percentage thresholds

#Full Predictive Mean
image2D(test, zlim=c(0,max(as.vector(test))), xlab="", ylab="",
        axes=F)

#Select outlier sub-image
empty_indx <- gp_masks_result$outliers

#First predictive mean
predmean_mat <- gp_masks_result$processed_images[[empty_indx]]
image2D(predmean_mat, zlim=c(0,1.1), xlab="", ylab="", axes=F)

#Plot thresholds for empty sub-images; FIGURE S1
for (k in 1:length(empty_indx)) {
  predmean_mat <- gp_masks_result$processed_images[[empty_indx[k]]]
  
  #Get pixel counts
  pixel_counts <- sapply(percentages, function(threshold) {
    threshold_image(mat = predmean_mat, percentage = threshold, count = TRUE)
  })
  
  # Compute the absolute differences in pixel counts
  diff_pixel_counts <- abs(diff(pixel_counts))
  
  #Robust GaSP - default
  diff_mod <- rgasp(percentages[-1], diff_pixel_counts, nugget.est = T)
  smoothed <- predict(diff_mod, testing_input=as.matrix(percentages[-1]))
  
  diff_pixel_counts<- smoothed$mean
  second_diff <- c()
  abs_second_diff <- c()
  max_index <- which.max(diff_pixel_counts)
  th<- 0.05*sd(diff_pixel_counts)
  for (i in 2:length(diff_pixel_counts)) {
    crit_bool <- T
    second_diff[i] <- diff_pixel_counts[i] - diff_pixel_counts[i - 1]
    abs_second_diff[i] <- abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1])
  }
  
  found_stable <- F
  for (i in (max_index + 1):length(diff_pixel_counts)) {
    # Check if the absolute difference between successive points is less than the threshold
    if (abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1]) < th) {
      stable_index <- i
      found_stable <- T
      break
    }
  }
  
  image2D(gp_masks_result$processed_images[[empty_indx[k]]], zlim=c(0,max(as.vector(test))), main=paste0(empty_indx[k]))
  plot(percentages[2:length(percentages)], abs_second_diff,
       xlab = "Thresholds", ylab="Absolute Second Difference",
       pch=16, cex = 0.5, cex.lab=1.2)
  lines(percentages[2:length(percentages)], abs_second_diff)
  if (found_stable) {
    abline(v = percentages[stable_index+1], lty = "dotted", lwd = 2, col = "red")
    estimated_percentage <- percentages[stable_index+1]
    thresholded_image <- threshold_image(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
    image2D(thresholded_image, col=c("white", "grey"), xlab="", ylab="", axes=F, colkey=F)
  }
  
  if (!found_stable) {
    thresholded_image <- threshold_image(mat = predmean_mat, percentage = 1, count = FALSE)
    image2D(thresholded_image, col=c("white", "grey"), xlab="", ylab="", axes=F, colkey=F)
  }
}

estimated_percentage #get first threshold

####Flood fill to detect number of objects in outlier####
image2D(bwlabel(thresholded_image), col = c("white", hcl.colors(600)), xlab="", ylab="", colkey=F)
#object count
outlier_count <- length(unique(as.vector(bwlabel(threshold_image(mat = gp_masks_result$processed_images[[gp_masks_result$outliers]], 
                                percentage = estimated_percentage, count = FALSE)))))
outlier_count-1 #537 objects total - 1 for background ID = 536

####Re-thresholding based off mean threshold####
rethresholded_image <- threshold_image(mat = gp_masks_result$processed_images[[gp_masks_result$outliers]], 
                                     percentage = gp_masks_result$crit_1_opt_thresholds[[gp_masks_result$outliers]], 
                                     count = FALSE)
image2D(rethresholded_image, col=c("white", "grey"), xlab="", ylab="", colkey=F, axes=F)
#note that the average threshold is...
gp_masks_result$crit_1_opt_thresholds[[gp_masks_result$outliers]]

#Flood fill for rethresholded image
image2D(bwlabel(rethresholded_image),col = c("white", hcl.colors(40)), axes=F, xlab="", ylab="", colkey=F)
#object count
outlier_recount <- length(unique(as.vector(bwlabel(rethresholded_image))))
outlier_recount - 1 #22 objects total - 1 for background ID = 21

#Flood fill for binary matrix before flood fill
prefinal_ff <- bwlabel(rethresholded_image)
remove_ids <- as.numeric(names(table(as.vector(prefinal_ff)))[table(as.vector(prefinal_ff)) < 50])
prefinal_ff[which(prefinal_ff %in% remove_ids, arr.ind=T)] <- 0

#"Pre"-Final result (before watershed segmentation)
image2D(prefinal_ff,col = c("white", hcl.colors(40)), axes=F, xlab="", ylab="", colkey=F)

outlier_prefinal_count <- length(unique(as.vector(prefinal_ff)))
outlier_prefinal_count - 1 #12 objects total - 1 for background ID = 11

#Objects in final result (after watershed)
#The final result is not split into sub-images, so the indices need to be retrieved
#to output the corresponding sub-matrix within the final result
img <- image_read(file_path5)
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
    
# Crop the final result using the calculated dimensions and offsets
x_split <- seq(1, img_width-1, length.out = num_pieces_x + 1) #subtract 1 from image width to match img_crop behavior
y_split <- seq(1, img_height, length.out = num_pieces_y + 1)

# Used to locate cropped area for final masks
#for (i in 1:(length(x_split)-1)) {
#  for (j in 1:(length(y_split)-1)) {
#    image2D(GP_masks[y_split[j]:y_split[j+1], x_split[i]:x_split[i+1]],
#            main=paste(i,j))
#  }
#}

# Manually determined from code above
x_indx <- 2
y_indx <- 2

#Final result
final_crop_mask <- GP_masks[y_split[y_indx]:y_split[y_indx+1],
                            x_split[x_indx]:x_split[x_indx+1]]
image2D(final_crop_mask, col = c("white", hcl.colors(40)), axes=F, xlab="", 
        ylab="", colkey=F)

#retrieve final object count
outlier_final_count <- length(unique(as.vector(final_crop_mask)))
outlier_final_count - 1 #13 objects total - 1 for background = 12 cell objects
