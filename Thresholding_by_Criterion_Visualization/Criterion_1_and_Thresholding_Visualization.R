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

#### Print Criterion Examples ####
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

#select points for display
indx_3 <- 20
indx_4 <- 80
selected_pts <- c(min(second_diff, na.rm=T), max(second_diff, na.rm=T), second_diff[indx_3], second_diff[indx_4])

#first pt - minimum of second difference
indx_1 <- which(second_diff == selected_pts[1])
pct_1 <- percentages[2:length(percentages)][indx_1]
thresholded_image1 <- threshold_image(mat = predmean_mat, percentage = pct_1, count = FALSE)
image2D(thresholded_image1, col=c("white", "grey"))

#second pt - maximum of second difference
indx_2 <- which(second_diff == selected_pts[2])
pct_2 <- percentages[2:length(percentages)][indx_2]
thresholded_image2 <- threshold_image(mat = predmean_mat, percentage = pct_2, count = FALSE)
image2D(thresholded_image2, col=c("white", "grey"))

#third pt - second difference at index=20
pct_3 <- percentages[2:length(percentages)][indx_3]
thresholded_image3 <- threshold_image(mat = predmean_mat, percentage = pct_3, count = FALSE)
image2D(thresholded_image3, col=c("white", "grey"))

#fourth pt - second difference at index=80
pct_4 <- percentages[2:length(percentages)][indx_4]
thresholded_image4 <- threshold_image(mat = predmean_mat, percentage = pct_4, count = FALSE)
image2D(thresholded_image4, col=c("white", "grey"))

#Optimal threshold - set by criterion 1
estimated_percentage <- percentages[stable_index+1]
thresholded_image <- threshold_image(mat = predmean_mat, percentage = estimated_percentage, count = FALSE)
image2D(thresholded_image, col=c("white", "grey"))


### PLOT SECOND ABSOLUTE DIFFERENCE

par(mgp = c(2.1, 0.3, 0))  # Adjusts the axis spacing: titles, labels, and lines
##### Plot with formulas in axes
plot(percentages[2:length(percentages)], abs_second_diff, 
     xlab = expression(~alpha[m]), 
     ylab = expression("|" ~ Delta * c[k]^"*" ~ (alpha[m]) - Delta * c[k]^"*" ~ (alpha[m-1]) ~ "|"), 
     pch = 16, 
     cex = 0.5, 
     font.lab = 2,    
     bty = "n",       
     ylim = c(0, 4000), 
     axes = FALSE)    

#Add x-axis 
axis(1, at = seq(0, 1, by = 0.25), labels = TRUE, lwd.ticks = 0, lwd = 1, cex.axis = 0.8)
#Add y-axis 
axis(2, at = seq(0, 4000, by = 1000), labels = TRUE, lwd.ticks = 0, lwd = 1, las = 1, cex.axis = 0.8)

#denote locations of percentages of interest
lines(percentages[2:length(percentages)], abs_second_diff, lwd=2)
 points(pct_4, abs_second_diff[indx_4], col="orange", pch = 16)
 points(pct_3, abs_second_diff[indx_3], col="yellow", pch = 16)
 points(pct_2, abs_second_diff[indx_2], col="green", pch = 16)
 points(pct_1, abs_second_diff[indx_1], col="blue", pch = 16)
 abline(v = percentages[stable_index+1], lty = "dotted", lwd = 2, col = "red")
 abline(h=0.05*sd(diff_pixel_counts), col="purple", lty = "dotted", lwd = 2)
 
 text(x = 0.94, y = 0.05 * sd(diff_pixel_counts), 
      labels = expression(0.05 ~ tau[k]), 
      col = "purple", cex = 1.2, pos = 3)

### Plot with pixel intensity histogram
 
df <- as.data.frame(as.vector(predmean_mat))

ggplot(data = df, aes(x=as.vector(predmean_mat))) + geom_histogram(color="black", fill="lightblue") +
   geom_vline(xintercept=percentages[stable_index+1]*max(as.vector(predmean_mat)), col="red") +
   annotate("text", x=percentages[stable_index+1]*max(as.vector(predmean_mat))+0.015, y=9000, label="Optimal Threshold", angle=90, col="red", cex=5) +
   geom_vline(xintercept=pct_1*max(as.vector(predmean_mat)), col="blue", size = 0.75, linetype = "dashed") +
   geom_vline(xintercept=pct_2*max(as.vector(predmean_mat)), col="green", size = 0.75, linetype = "dashed") +
   geom_vline(xintercept=pct_3*max(as.vector(predmean_mat)), col="yellow", size = 0.75, linetype = "dashed") +
   geom_vline(xintercept=pct_4*max(as.vector(predmean_mat)), col="orange", size = 0.75, linetype = "dashed") +
   xlab("Pixel Intensity Value") +
   ylab("Frequency") +
   theme_minimal() + theme(text = element_text(size = 15))

### Plot optimal threshold against pixel intensity heights

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
      zaxis = list(title = "Height")
    )
  )

# Display the plot
print(surface_plot)
