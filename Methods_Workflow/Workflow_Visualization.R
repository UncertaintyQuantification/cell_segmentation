####Workflow Visualization####

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

#### Display the original Image ####
file_path <- "../Image_Data/cropped_pieces/cropped_img_1.jpg"
img <- image_read(file_path)
print(img)
ori_img_matrix <- as.numeric(img[[1]])[,,1]

#### Use Separable Kernel GP and Cropped Images to generate GP Masks ####
gp_masks_result <- generate_GP_Masks_test(file_path) # solves the numeric flow issue, but still has learning part
GP_masks <- gp_masks_result$GP_masks
image2D(GP_masks, axes=F, xlab="", ylab="")

####Workflow Figure Images####
test <- gp_masks_result$combined_predmean

#Original Image
image2D(ori_img_matrix, zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F)

#Original Cropped into 4 (step 1)
image2D(ori_img_matrix[1:nrow(ori_img_matrix)/2, 1:ncol(ori_img_matrix)/2], zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F, colkey=F)
image2D(ori_img_matrix[1:(nrow(ori_img_matrix)/2), (ncol(ori_img_matrix)/2):ncol(ori_img_matrix)], zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F, colkey=F)
image2D(ori_img_matrix[(nrow(ori_img_matrix)/2):nrow(ori_img_matrix), 1:ncol(ori_img_matrix)/2], zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F, colkey=F)
image2D(ori_img_matrix[(nrow(ori_img_matrix)/2):nrow(ori_img_matrix), (ncol(ori_img_matrix)/2):ncol(ori_img_matrix)],zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F, colkey=F)

##Predmean calculated for cropped pieces (step 2)
image2D(test[1:nrow(ori_img_matrix)/2, 1:ncol(ori_img_matrix)/2], zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F)
image2D(test[1:(nrow(ori_img_matrix)/2), (ncol(ori_img_matrix)/2):ncol(ori_img_matrix)], zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F)
image2D(test[(nrow(ori_img_matrix)/2):nrow(ori_img_matrix), 1:ncol(ori_img_matrix)/2], zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F)
image2D(test[(nrow(ori_img_matrix)/2):nrow(ori_img_matrix), (ncol(ori_img_matrix)/2):ncol(ori_img_matrix)],zlim=c(0,max(as.vector(test))), xlab="", ylab="", axes=F)

#Threshold set for cropped predmean (step 3)
bin <- gp_masks_result$combined_thresholded1
image2D(bin[1:nrow(ori_img_matrix)/2, 1:ncol(ori_img_matrix)/2], zlim=c(0,max(as.vector(test))), col=rev(gray(seq(0, 1, length.out = 100))), xlab="", ylab="", axes=F, colkey=F)
image2D(bin[1:(nrow(ori_img_matrix)/2), (ncol(ori_img_matrix)/2):ncol(ori_img_matrix)], zlim=c(0,max(as.vector(test))), col=rev(gray(seq(0, 1, length.out = 100))), xlab="", ylab="", axes=F, colkey=F)
image2D(bin[(nrow(ori_img_matrix)/2):nrow(ori_img_matrix), 1:ncol(ori_img_matrix)/2], zlim=c(0,max(as.vector(test))), col=rev(gray(seq(0, 1, length.out = 100))), xlab="", ylab="", axes=F, colkey=F)
image2D(bin[(nrow(ori_img_matrix)/2):nrow(ori_img_matrix), (ncol(ori_img_matrix)/2):ncol(ori_img_matrix)],zlim=c(0,max(as.vector(test))), col=rev(gray(seq(0, 1, length.out = 100))), xlab="", ylab="", axes=F, colkey=F)

#Binary image stitched back together (step 4)
image2D(bin,zlim=c(0,max(as.vector(test))), col=rev(gray(seq(0, 1, length.out = 100))), xlab="", ylab="", axes=F, colkey=F)

#Distance matrix (step 5) - not displayed in workflow figure but used for watershed
image2D(distmap(bin), xlab="", ylab="", axes=F)

#GP Masks by watershed (final)
library(randomcoloR)
image2D(GP_masks, col=c("white",distinctColorPalette(length(unique(as.vector(GP_masks)))-1)), xlab="", ylab="", axes=F, colkey=F)
