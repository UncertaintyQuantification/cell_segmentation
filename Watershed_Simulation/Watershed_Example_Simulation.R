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
library(reshape2)
#Ensure directory is set to file folder
source("../src/Modified_Functions_RGasp.R")

#### Display the original Image ####
file_path <- "../Image_Data/cropped_pieces/cropped_img_2.jpg"
img <- image_read(file_path)
ori_img_matrix <- as.numeric(img[[1]])[,,1]

#### Use Separable Kernel GP and Cropped Images to generate GP Masks ####
gp_masks_result <- generate_GP_Masks_test(file_path)
GP_masks <- gp_masks_result$GP_masks
image2D(GP_masks)

#Crop to retrieve two cells
gp_mask_select <- GP_masks[344:412, 694:735]
image2D(gp_mask_select) #Result from using watershed() function

#Write to CSV - only run when diagonal map is not provided and needs to be prepped
#write.csv(gp_mask_select, "different_size_diag_cropping_2.csv")

#Get binary mask
bin_two <- matrix(0, nrow=nrow(gp_mask_select), ncol=ncol(gp_mask_select))
bin_two[which(gp_mask_select > 0, arr.ind=T)] <- 1
image2D(bin_two, col=c("white", "grey"), xlab="", ylab="", axes=F, colkey=F)

#Get distance matrix
dist_map <- distmap(bin_two)
image2D(dist_map, xlab="", ylab="", axes=F)

#Get negative distance matrix
neg_dist <- -dist_map
image2D(neg_dist, xlab="", ylab="", axes=F, cex=2)

#observe diagnonal - prepared in Excel and provided
modified_neg_dist <- as.matrix(read.csv("different_size_diag_cropping_2.csv"))[,-1]
image2D(modified_neg_dist) # the "diagonal" elements are 1

#View diagonal heights
heights <- neg_dist[which(modified_neg_dist==1, arr.ind=T)]
plot(heights, xlab="Position", ylab="Height")

#View diagonal on distance matrix
neg_dist_diag <- neg_dist
neg_dist_diag[which(modified_neg_dist==1, arr.ind=T)] <- 1
image2D(neg_dist_diag)

# Function to simulate water filling up to a certain level
add_water <- function(df, water_level) {
  df$water <- pmax(water_level, df$height)
  return(df)
}

#Generate dataframe with height and position of diagonal
indices <- which(modified_neg_dist==1, arr.ind=T)
landscape_df <- data.frame(
  x = c(1:length(neg_dist[indices])),
  height = neg_dist[indices]
)

generate_watershed_plot <- function(water_level) {
  filled_landscape <- add_water(landscape_df, water_level)
  
  ggplot(filled_landscape, aes(x = x)) +
    # Plot the water level first (below the landscape line)
    geom_ribbon(aes(ymin = height, ymax = water), 
                fill = "lightblue", alpha = 0.7) +
    # Plot the landscape on top
    geom_line(aes(y = height), color = "brown", size = 1) +
    # Plot the current water level line
    geom_hline(yintercept = water_level, 
               color = "blue", linetype = "dashed") +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 20, color = "black")) +
    coord_cartesian(ylim = c(-16, 1)) #+
}

#Initial plot - labeled minimums, no water level
#see where the IDs switch based on negative distance map
switch_vis <- cbind(neg_dist[indices], gp_mask_select[indices])
df_switch <- as.data.frame(switch_vis)
colnames(df_switch) <- c("height", "label")
df_switch <- df_switch %>% group_by(label) %>% 
  summarize(min_height = min(height))

min_heights <- rbind(switch_vis[switch_vis[,1] == df_switch$min_height[2] & switch_vis[,2] == df_switch$label[2],1],
                     switch_vis[switch_vis[,1] == df_switch$min_height[3] & switch_vis[,2] == df_switch$label[3],1])
min_positions <- rbind(which(switch_vis[,1] == min_heights[1]), which(switch_vis[,1] == min_heights[2]))

local_min <- as.data.frame(cbind(min_positions, min_heights))
colnames(local_min) <- c("x", "y")

filled_landscape <- add_water(landscape_df, -15)

ggplot(filled_landscape, aes(x = x)) +
  # Plot the landscape on top
  geom_line(aes(y = height), color = "brown", size = 1) +
  
  # Customize the plot
  theme_minimal() +
  coord_cartesian(ylim = c(-16, 1)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 20, color = "black")) +
  geom_point(data=local_min, aes(x=x, y=y), col=c("purple", "orange"), cex=2)

#Watershed plots - labeled water level before watershed line is formed
water_levels_demo <- c(-13)

print(generate_watershed_plot(water_levels_demo))

#Watershed line plot
#see where the IDs switch based on negative distance map
switch_vis <- cbind(neg_dist[indices], gp_mask_select[indices])
switch_vis
#switches at approx height = -10.198039 and position 14

watershed_height <- -10.198039
watershed_position <- 14

#Euclidean distances between the two minimums:
#purple
purple_distance <- sqrt((local_min[1,1]-watershed_position)^2 + (local_min[1,2] - watershed_height)^2)
purple_distance #purple is closer

#orange
orange_distance <- sqrt((local_min[2,1]-watershed_position)^2 + (local_min[2,2] - watershed_height)^2)
orange_distance

final_plot <- generate_watershed_plot(round(watershed_height,2)) +
  geom_vline(xintercept = watershed_position, 
             color = "red", 
             linetype = "dashed") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 20, color = "black"))

print(final_plot)

#Full plot - water level is at top
final_plot <- generate_watershed_plot(0) +
  geom_vline(xintercept = 14, 
             color = "red", 
             linetype = "dashed") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 20, color = "black"))

print(final_plot)

#Watershed result
watershed_result <- watershed(dist_map)
image2D(watershed_result, col = hcl.colors(100, "plasma"), xlab="", ylab="",
        axes=F, colkey=F)

#### Watershed Overlay Plots ####
ori_with_diagonal <- ori_img_matrix[344:412, 694:735]

#with diagonal
image2D(ori_with_diagonal, col=grey.colors(length(table(as.vector(ori_with_diagonal)))),
        colkey=F, axes=F, xlab="", ylab="")
lines(x=indices[,1]/nrow(ori_with_diagonal), y=indices[,2]/ncol(ori_with_diagonal),
      col="red", lwd=2)

#without diagonal
image2D(ori_with_diagonal, col=c(grey.colors(length(table(as.vector(ori_with_diagonal))))),
        colkey=F, axes=F, xlab="", ylab="")

#with water levels
all_levels <- c(-13, watershed_height, 0)
for (level in all_levels) {
  image2D(ori_with_diagonal, col=grey.colors(length(table(as.vector(ori_with_diagonal)))),
          colkey=F, axes=F, xlab="", ylab="")
  
  #purple cell indx
  cell_purple_indx <- which(neg_dist <= level & watershed_result == 1, arr.ind=T)
  points(cell_purple_indx[,1]/nrow(ori_with_diagonal),
         cell_purple_indx[,2]/ncol(ori_with_diagonal),
         col = rgb(0.5, 0, 0.5, alpha = 0.3), pch = 16, cex = 2
  )
  
  #orange cell indx
  cell_orange_indx <- which(neg_dist < level & watershed_result == 2, arr.ind=T)
  points(cell_orange_indx[,1]/nrow(ori_with_diagonal),
         cell_orange_indx[,2]/ncol(ori_with_diagonal),
         col = rgb(1, 0.647, 0, alpha = 0.3), pch = 16, cex = 2
  )
}
