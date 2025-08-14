#### Criterion 1 Visualization ####

#### Import Libraries and Functions ####
library(magick)
library(pracma)
library(lhs)
library(RobustGaSP)
library(plot3D)
library(dplyr)
library(ggplot2)
library(EBImage)
library(plotly)
#Ensure directory is set to file folder
source("../src/Modified_Functions_RGasp.R")

#### Display the original Image ####
file_path <- "../Image_Data/cropped_pieces/cropped_img_4.jpg"
img <- image_read(file_path)
ori_img_matrix <- as.numeric(img[[1]])[,,1]

#### Use Separable Kernel GP and Cropped Images to generate GP Masks ####
gp_masks_result <- generate_GP_Masks_test(file_path)
GP_masks <- gp_masks_result$GP_masks
image2D(GP_masks)

#For displaying criterion 1, only a portion of the original predictive mean needs to be utilized
predmean_mat <- gp_masks_result$combined_predmean[1:200, 1:300]
image2D(predmean_mat)

#### Criterion Examples - 0.05*SD ####
delta = 0.01

# Define percentage thresholds
percentages <- seq(0, 1, by = delta)

# Calculate pixel counts based on thresholding the image
pixel_counts <- sapply(percentages, function(threshold) {
  threshold_image(mat = predmean_mat, percentage = threshold, count = TRUE)
})

# Compute the absolute differences in pixel counts
diff_pixel_counts <- abs(diff(pixel_counts))

diff_mod <- rgasp(percentages[-1], diff_pixel_counts, nugget.est=T)
smoothed <- predict(diff_mod, testing_input=as.matrix(percentages[-1]))

#### Using changes in Diff
diff_pixel_counts<- smoothed$mean

second_diff <- c()
abs_second_diff <- c()
max_index <- which.max(diff_pixel_counts)
th<- 0.05*sd(diff_pixel_counts)
for (i in 2:length(diff_pixel_counts)) {
  second_diff[i] <- diff_pixel_counts[i] - diff_pixel_counts[i - 1]
  abs_second_diff[i] <- abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1])
}
for (i in (max_index + 1):length(diff_pixel_counts)) {
  # Check if the absolute difference between successive points is less than the threshold
  if (abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1]) < th) {
    stable_index <- i
    break
  }
}

#Optimal threshold - set by criterion 1
estimated_percentage <- percentages[stable_index+1]
thresholded_image <- threshold_image(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
image2D(thresholded_image, col=c("white", "grey"))

#Get height matrix
height_matrix <- predmean_mat

# Set threshold
threshold <- matrix(percentages[stable_index+1]*max(as.vector(predmean_mat)))

# Create interactive 3D surface plot with threshold plane
surface_plot <- plot_ly() %>%
  # Add the main surface
  add_surface(
    z = height_matrix,
    name = "Height Surface"
  ) %>%
  # Add the threshold plane
  add_surface(
    z = matrix(threshold, 
               nrow = nrow(height_matrix), 
               ncol = ncol(height_matrix)),  # Make plane full size
    colorscale = list(c(0,1), c("red","red")),
    opacity = 0.2,  # Adjust transparency
    name = "Optimal Threshold",
    showscale = FALSE,
    hoverinfo = "none"
  ) %>%
  # Add text annotation for threshold
  add_annotations(
    x = dim(predmean_mat)[2],  # Adjusted for your matrix width
    y = dim(predmean_mat)[1],  # Adjusted for your matrix height
    z = threshold,
    text = "Optimal Threshold",
    showarrow = TRUE,
    arrowhead = 2,
    arrowsize = 1,
    arrowcolor = "red"
  ) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      yaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      zaxis = list(title = "Predictive Mean of Intensity")
    )
  )

# Display the plot
print(surface_plot)

#### Criterion Examples - 0.1*SD ####
delta = 0.01

# Define percentage thresholds
percentages <- seq(0, 1, by = delta)

# Calculate pixel counts based on thresholding the image
pixel_counts <- sapply(percentages, function(threshold) {
  threshold_image(mat = predmean_mat, percentage = threshold, count = TRUE)
})

# Compute the absolute differences in pixel counts
diff_pixel_counts <- abs(diff(pixel_counts))

diff_mod <- rgasp(percentages[-1], diff_pixel_counts, nugget.est=T)
smoothed <- predict(diff_mod, testing_input=as.matrix(percentages[-1]))

#### Using changes in Diff
diff_pixel_counts<- smoothed$mean

second_diff <- c()
abs_second_diff <- c()
max_index <- which.max(diff_pixel_counts)
th<- 0.1*sd(diff_pixel_counts)
for (i in 2:length(diff_pixel_counts)) {
  second_diff[i] <- diff_pixel_counts[i] - diff_pixel_counts[i - 1]
  abs_second_diff[i] <- abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1])
}
for (i in (max_index + 1):length(diff_pixel_counts)) {
  # Check if the absolute difference between successive points is less than the threshold
  if (abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1]) < th) {
    stable_index <- i
    break
  }
}

#Optimal threshold - set by criterion 1
estimated_percentage <- percentages[stable_index+1]
thresholded_image <- threshold_image(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
image2D(thresholded_image, col=c("white", "grey"))

#Get height matrix
height_matrix <- predmean_mat

# Set threshold
threshold <- matrix(percentages[stable_index+1]*max(as.vector(predmean_mat)))

# Create interactive 3D surface plot with threshold plane
surface_plot <- plot_ly() %>%
  # Add the main surface
  add_surface(
    z = height_matrix,
    name = "Height Surface"
  ) %>%
  # Add the threshold plane
  add_surface(
    z = matrix(threshold, 
               nrow = nrow(height_matrix), 
               ncol = ncol(height_matrix)),  # Make plane full size
    colorscale = list(c(0,1), c("red","red")),
    opacity = 0.2,  # Adjust transparency
    name = "Optimal Threshold",
    showscale = FALSE,
    hoverinfo = "none"
  ) %>%
  # Add text annotation for threshold
  add_annotations(
    x = dim(predmean_mat)[2],  # Adjusted for your matrix width
    y = dim(predmean_mat)[1],  # Adjusted for your matrix height
    z = threshold,
    text = "Optimal Threshold",
    showarrow = TRUE,
    arrowhead = 2,
    arrowsize = 1,
    arrowcolor = "red"
  ) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      yaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      zaxis = list(title = "Predictive Mean of Intensity")
    )
  )

# Display the plot
print(surface_plot)


#### Criterion Examples - 2.5*SD ####
delta = 0.01

# Define percentage thresholds
percentages <- seq(0, 1, by = delta)

# Calculate pixel counts based on thresholding the image
pixel_counts <- sapply(percentages, function(threshold) {
  threshold_image(mat = predmean_mat, percentage = threshold, count = TRUE)
})

# Compute the absolute differences in pixel counts
diff_pixel_counts <- abs(diff(pixel_counts))

diff_mod <- rgasp(percentages[-1], diff_pixel_counts, nugget.est=T)
smoothed <- predict(diff_mod, testing_input=as.matrix(percentages[-1]))

#### Using changes in Diff
diff_pixel_counts<- smoothed$mean

second_diff <- c()
abs_second_diff <- c()
max_index <- which.max(diff_pixel_counts)
th<- 2.5*sd(diff_pixel_counts)
for (i in 2:length(diff_pixel_counts)) {
  second_diff[i] <- diff_pixel_counts[i] - diff_pixel_counts[i - 1]
  abs_second_diff[i] <- abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1])
}
for (i in (max_index + 1):length(diff_pixel_counts)) {
  # Check if the absolute difference between successive points is less than the threshold
  if (abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1]) < th) {
    stable_index <- i
    break
  }
}

#Optimal threshold - set by criterion 1
estimated_percentage <- percentages[stable_index+1]
thresholded_image <- threshold_image(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
image2D(thresholded_image, col=c("white", "grey"))

#Get height matrix
height_matrix <- predmean_mat

# Set threshold
threshold <- matrix(percentages[stable_index+1]*max(as.vector(predmean_mat)))

# Create interactive 3D surface plot with threshold plane
surface_plot <- plot_ly() %>%
  # Add the main surface
  add_surface(
    z = height_matrix,
    name = "Height Surface"
  ) %>%
  # Add the threshold plane
  add_surface(
    z = matrix(threshold, 
               nrow = nrow(height_matrix), 
               ncol = ncol(height_matrix)),  # Make plane full size
    colorscale = list(c(0,1), c("red","red")),
    opacity = 0.2,  # Adjust transparency
    name = "Optimal Threshold",
    showscale = FALSE,
    hoverinfo = "none"
  ) %>%
  # Add text annotation for threshold
  add_annotations(
    x = dim(predmean_mat)[2],  # Adjusted for your matrix width
    y = dim(predmean_mat)[1],  # Adjusted for your matrix height
    z = threshold,
    text = "Optimal Threshold",
    showarrow = TRUE,
    arrowhead = 2,
    arrowsize = 1,
    arrowcolor = "red"
  ) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      yaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      zaxis = list(title = "Predictive Mean of Intensity")
    )
  )

# Display the plot
print(surface_plot)

#### Criterion Examples - 0.01*SD ####
delta = 0.01

# Define percentage thresholds
percentages <- seq(0, 1, by = delta)

# Calculate pixel counts based on thresholding the image
pixel_counts <- sapply(percentages, function(threshold) {
  threshold_image(mat = predmean_mat, percentage = threshold, count = TRUE)
})

# Compute the absolute differences in pixel counts
diff_pixel_counts <- abs(diff(pixel_counts))

diff_mod <- rgasp(percentages[-1], diff_pixel_counts, nugget.est=T)
smoothed <- predict(diff_mod, testing_input=as.matrix(percentages[-1]))

#### Using changes in Diff
diff_pixel_counts<- smoothed$mean

second_diff <- c()
abs_second_diff <- c()
max_index <- which.max(diff_pixel_counts)
th<- 0.01*sd(diff_pixel_counts)
for (i in 2:length(diff_pixel_counts)) {
  second_diff[i] <- diff_pixel_counts[i] - diff_pixel_counts[i - 1]
  abs_second_diff[i] <- abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1])
}
for (i in (max_index + 1):length(diff_pixel_counts)) {
  # Check if the absolute difference between successive points is less than the threshold
  if (abs(diff_pixel_counts[i] - diff_pixel_counts[i - 1]) < th) {
    stable_index <- i
    break
  }
}

#Optimal threshold - set by criterion 1
estimated_percentage <- percentages[stable_index+1]
thresholded_image <- threshold_image(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
image2D(thresholded_image, col=c("white", "grey"))

#Get height matrix
height_matrix <- predmean_mat

# Set threshold
threshold <- matrix(percentages[stable_index+1]*max(as.vector(predmean_mat)))

# Create interactive 3D surface plot with threshold plane
surface_plot <- plot_ly() %>%
  # Add the main surface
  add_surface(
    z = height_matrix,
    name = "Height Surface"
  ) %>%
  # Add the threshold plane
  add_surface(
    z = matrix(threshold, 
               nrow = nrow(height_matrix), 
               ncol = ncol(height_matrix)),  # Make plane full size
    colorscale = list(c(0,1), c("red","red")),
    opacity = 0.2,  # Adjust transparency
    name = "Optimal Threshold",
    showscale = FALSE,
    hoverinfo = "none"
  ) %>%
  # Add text annotation for threshold
  add_annotations(
    x = dim(predmean_mat)[2],  # Adjusted for your matrix width
    y = dim(predmean_mat)[1],  # Adjusted for your matrix height
    z = threshold,
    text = "Optimal Threshold",
    showarrow = TRUE,
    arrowhead = 2,
    arrowsize = 1,
    arrowcolor = "red"
  ) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      yaxis = list(
        showticklabels = FALSE,
        showgrid = TRUE,
        zeroline = TRUE,
        showline = TRUE,
        title = ""
      ),
      zaxis = list(title = "Predictive Mean of Intensity")
    )
  )

# Display the plot
print(surface_plot)

