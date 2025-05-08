#### Import Libraries and Functions ####
library(magick)
library(pracma)
library(RobustGaSP)
library(plot3D)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EBImage)
#Ensure directory is set to file folder
source("../src/Modified_Functions_RGasp.R")

#### Define Folder Paths ####
base_dir <- getwd()  # Uses the current working directory

folders <- paste0("../Image_Data/nuclear_test_images/nuclei_figure_", 1:5)  # Automatically generate folder names

# This is to solve the python and R image visualization issue
process_image_mask <- function(image_path) {
  imageJ_mask <- readImage(image_path)
  imageJ_mask <- imageJ_mask@.Data
  
  # Extract unique values and map them to consecutive integers
  unique_values <- unique(as.vector(imageJ_mask))
  value_to_integer <- setNames(seq(0, length(unique_values) - 1), unique_values)
  
  # Convert image to integer mask
  integer_mask <- matrix(value_to_integer[as.character(imageJ_mask)], nrow = nrow(imageJ_mask), ncol = ncol(imageJ_mask))
  
  # Rotate the mask
  rotated_mask <- t(integer_mask)[, nrow(integer_mask):1]
  
  # Mirror the rotated mask
  mirrored_mask <- rotated_mask[, ncol(rotated_mask):1]
  
  return(mirrored_mask)
}

#### Use Average Precision under certain IoU as metric 
compute_ious <- function(true_mask, pred_mask) {
  # Get unique labels, excluding background
  unique_true_labels <- unique(true_mask[true_mask > 0])
  unique_pred_labels <- unique(pred_mask[pred_mask > 0])
  
  # Create binary masks for each label
  true_regions <- lapply(unique_true_labels, function(label) true_mask == label)
  pred_regions <- lapply(unique_pred_labels, function(label) pred_mask == label)
  
  # Precompute areas of each region 
  true_region_sums <- sapply(true_regions, sum)
  pred_region_sums <- sapply(pred_regions, sum)
  
  # Initialize IoU matrix
  ious <- matrix(0, nrow = length(unique_true_labels), ncol = length(unique_pred_labels))
  rownames(ious) <- unique_true_labels
  colnames(ious) <- unique_pred_labels
  
  for (i in seq_along(true_regions)) {
    ious[i, ] <- sapply(seq_along(pred_regions), function(j) {
      intersection <- sum(true_regions[[i]] & pred_regions[[j]])
      union <- true_region_sums[i] + pred_region_sums[j] - intersection
      if (union > 0) intersection / union else 0
    })
  }
  
  return(ious)
}
compute_ap_from_ious <- function(ious, threshold = 0.5) {
  tp <- 0
  fp <- 0
  fn <- 0
  
  # Track which rows (true labels) and columns (predicted labels) have been matched
  true_matched <- rep(FALSE, nrow(ious))
  pred_matched <- rep(FALSE, ncol(ious))
  
  # Iterate through the IoU matrix to count TP, FP, and FN (We assume the max iou is the best fit)
  for (i in seq_len(nrow(ious))) {
    row_max_iou <- max(ious[i, ])
    max_j <- which.max(ious[i, ])
    
    # True Positive if IoU meets or exceeds the threshold
    if (row_max_iou >= threshold) {
      tp <- tp + 1
      true_matched[i] <- TRUE
      pred_matched[max_j] <- TRUE
    }
  }
  
  # Count FN and FP based on unmatched rows and columns
  fn <- sum(!true_matched)  # False Negatives (unmatched true labels)
  fp <- sum(!pred_matched)  # False Positives (unmatched predicted labels)
  
  if (tp + fp > 0) {
    precision <- tp / (tp + fp + fn)
  } else {
    precision <- 0
  }
  
  return(list(
    precision = precision,
    tp = tp,
    fp = fp,
    fn = fn
  ))
}
compute_ap_table <- function(ious_gp, ious_imagej) {
  # Define IoU thresholds for AP calculation
  thresholds <- seq(0.5, 0.80, by = 0.05)
  
  # Initialize vectors to store precision values
  gp_precisions <- numeric(length(thresholds))
  imagej_precisions <- numeric(length(thresholds))
  is_gp_better <- logical(length(thresholds))  
  
  # Loop through each threshold to compute Average Precision (AP)
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    
    # Compute AP for GP-based method
    ap_result_gp <- compute_ap_from_ious(ious_gp, threshold = threshold)
    gp_precisions[i] <- ap_result_gp$precision
    
    # Compute AP for ImageJ-based method
    ap_result_imagej <- compute_ap_from_ious(ious_imagej, threshold = threshold)
    imagej_precisions[i] <- ap_result_imagej$precision
    
    # Compare methods: TRUE if GP has higher precision than ImageJ
    is_gp_better[i] <- gp_precisions[i] > imagej_precisions[i]
  }
  
  ap_table <- data.frame(
    Threshold = thresholds,
    GP_Method_AP = gp_precisions,  
    ImageJ_Method_AP = imagej_precisions,  
    GP_Better = is_gp_better  
  )
  
  return(ap_table)  
}

#### Save Mask Figures ####
save_boundary_figures <- function(folder_path, ori_img_matrix, GP_masks, ImageJ_Masks, True_Mask) {
  # Find Boundaries
  GP_boundaries <- find_boundaries(GP_masks)
  ImageJ_boundaries <- find_boundaries(ImageJ_Masks)
  True_boundaries <- find_boundaries(True_Mask)
  
  #### Save GP Masks Figure
  png(file.path(folder_path, "GP_boundaries.png"), width = 800, height = 800)
  par(mar = c(0, 0, 0, 0))
  image2D(ori_img_matrix, colkey = FALSE, axes = FALSE)
  boundary_coords <- which(GP_boundaries == 1, arr.ind = TRUE)
  points(boundary_coords[, 1] / nrow(ori_img_matrix),
         boundary_coords[, 2] / ncol(ori_img_matrix),
         col = "black", pch = 19, cex = 0.5)
  dev.off()
  
  #### Save ImageJ Masks Figure
  png(file.path(folder_path, "ImageJ_boundaries.png"), width = 800, height = 800)
  par(mar = c(0, 0, 0, 0))
  image2D(ori_img_matrix, colkey = FALSE, axes = FALSE)
  ImageJ_boundary_coords <- which(ImageJ_boundaries == 1, arr.ind = TRUE)
  points(ImageJ_boundary_coords[, 1] / nrow(ori_img_matrix), 
         ImageJ_boundary_coords[, 2] / ncol(ori_img_matrix), 
         col = "black", pch = 19, cex = 0.5)
  dev.off()
  
  #### Save True Masks Figure
  png(file.path(folder_path, "True_boundaries.png"), width = 800, height = 800)
  par(mar = c(0, 0, 0, 0))
  image2D(ori_img_matrix, colkey = FALSE, axes = FALSE)
  True_boundary_coords <- which(True_boundaries == 1, arr.ind = TRUE)
  points(True_boundary_coords[, 1] / nrow(ori_img_matrix), 
         True_boundary_coords[, 2] / ncol(ori_img_matrix), 
         col = "black", pch = 19, cex = 0.5)
  dev.off()
}

#### Compute and Save IoUs ####
save_ious <- function(folder_path, True_Mask, GP_masks, ImageJ_Masks) {
  # Compute IoUs
  ious_gp <- compute_ious(True_Mask, GP_masks)
  ious_imagej <- compute_ious(True_Mask, ImageJ_Masks)
  
  # Save IoU Data
  write.csv(ious_gp, file = file.path(folder_path, "ious_gp.csv"), row.names = TRUE)
  write.csv(ious_imagej, file = file.path(folder_path, "ious_imagej.csv"), row.names = TRUE)
  
  return(compute_ap_table(ious_gp, ious_imagej))
}

#### Process Each Image ####
process_image <- function(folder_name, base_dir) {
  folder_path <- file.path(base_dir, folder_name)
  
  file_path <- file.path(folder_path, "original_fig.png")
  ImageJ_Masks_path <- file.path(folder_path, "original_ImageJ_masks.tif")
  True_Masks_path <- file.path(folder_path, "original_true_masks.png")
  
  img <- image_read(file_path)
  ori_img_matrix <- as.numeric(img[[1]])[, , 1]
  
  # Generate GP Masks
  gp_masks_result <- generate_GP_Masks_test(file_path, nugget = T) # Default is Nugget = F
  GP_masks <- gp_masks_result$GP_masks 
  
  # Process ImageJ and True Masks
  ImageJ_Masks <- process_image_mask(ImageJ_Masks_path)
  True_Mask <- process_image_mask(True_Masks_path)
  
  # Save Boundary Figures
  save_boundary_figures(folder_path, ori_img_matrix, GP_masks, ImageJ_Masks, True_Mask)
  
  # Compute and Save IoUs
  ap_table <- save_ious(folder_path, True_Mask, GP_masks, ImageJ_Masks)
  
  # Add folder name as label
  ap_table$Pair <- folder_name  
  return(ap_table)
}

#### Process All Images ####
all_ap_tables <- lapply(folders, process_image, base_dir = base_dir)
combined_table <- bind_rows(all_ap_tables)

#### Read IoU CSV files to skip the codes above and save into a single table ####
skip = T # Change this
if (skip){
  iou_data <- lapply(folders, function(folder_name) {
    folder_path <- file.path(base_dir, folder_name)
    
    # Read GP and ImageJ IoU CSVs from the folder
    ious_gp <- read.csv(file.path(folder_path, "ious_gp.csv"), row.names = 1)
    ious_imagej <- read.csv(file.path(folder_path, "ious_imagej.csv"), row.names = 1)
    
    # Compute AP Table
    ap_table <- compute_ap_table(ious_gp, ious_imagej)
    
    # Add folder name as label
    ap_table$Pair <- folder_name  
    
    return(ap_table)
  })
  combined_table <- bind_rows(iou_data)
}

#### Plot AP ####
plot_data <- combined_table %>%
  pivot_longer(cols = c(GP_Method_AP, ImageJ_Method_AP), 
               names_to = "Method", values_to = "AP") %>%
  mutate(Method = recode(Method, "GP_Method_AP" = "GP", "ImageJ_Method_AP" = "ImageJ"),
         Pair = recode(Pair,
                       "nuclei_figure_1" = "Image 1",
                       "nuclei_figure_2" = "Image 2",
                       "nuclei_figure_3" = "Image 3",
                       "nuclei_figure_4" = "Image 4",
                       "nuclei_figure_5" = "Image 5"))

# Create and save the boxplot comparing GP vs. ImageJ
plot_box <- ggplot(plot_data, aes(x = as.factor(Threshold), y = AP, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2, aes(color = Method)) +
  labs(x = "Threshold", 
       y = "Average Precision (AP)", 
       fill = "Method",
       color = "Method") +
  theme_classic() +
  scale_fill_manual(values = c("GP" = "blue", "ImageJ" = "red")) +
  scale_color_manual(values = c("GP" = "blue", "ImageJ" = "red")) +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20, color = "black"),
    legend.position = "top",
    legend.justification = "right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20, color = "black"),
    plot.margin = margin(30, 30, 30, 30)
  )

print(plot_box)
ggsave("mean_ap_boxplot_comparison.png", plot = plot_box, width = 10, height = 8, dpi = 300)

# Create and save the line plot showing individual test image performance
plot_line <- ggplot(plot_data, aes(
  x = Threshold, 
  y = AP, 
  color = Pair, 
  linetype = Method, 
  group = interaction(Pair, Method)
)) +
  geom_line(size = 1.5) +  
  geom_point(size = 3) +  
  labs(
    x = "Threshold",
    y = "Average Precision (AP)",
    color = "Images", 
    linetype = "Method"
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "Image 1" = "blue",
    "Image 2" = "green",
    "Image 3" = "purple",
    "Image 4" = "orange",
    "Image 5" = "brown"
  )) +
  scale_linetype_manual(values = c("GP" = "solid", "ImageJ" = "dashed")) +
  guides(
    linetype = guide_legend(
      override.aes = list(size = 1.5, linetype = c("solid", "dashed"), color = c("black", "black"))
    )
  ) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "right",
    legend.spacing.x = unit(0.4, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(plot_line)
ggsave("ap_comparison_plot.png", plot = plot_line, width = 10, height = 8, dpi = 300)
