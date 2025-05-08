# Unsupervised Cellular Boundary Detection by Fast Gaussian Processes
This repository contains the code for the paper *Unsupervised Cellular Boundary Detection by Fast Gaussian Processes* by Laura Baracaldo, Blythe King, Haoran Yan, Yizi Lin, Nina Miolane, and Mengyang Gu.

This folder contains data and code to reproduce the results from the paper.

There are 9 folders.

1. src -- This folder contains the needed R functions to perform the cell image segmentation by fast Gaussian processes. In R, we use generate_GP_Masks_test() to call these functions. This folder also contains the functions for comparing smoothing results for simulated experiments.
2. Image_Data -- Contains the raw images used for the demonstration of methods and real data analysis. The images in cropped_pieces are smaller images used for in-text figures and examples. The images in nuclear_test_figures are utilized in the real data analysis for the nuclear channel. The images in whole_cell_test_images are utilized in the real data analysis for the whole cell channel. Additionally, the test_images folders contain the generated true masks and ImageJ masks used for results comparison. Finally, simulated_cells contains the cell images used in Section 4.
3. Methods_Workflow -- Contains the code needed to reproduce the methods workflow given in Figure 1 in Section 3 and provides image output at each major step of the proposed segmentation method.
4. Thresholding_by_Criterion_Visualization -- Contains the code needed to reproduce Figures 2 and 3 in Section 3.2 and provides a visualization of how the proposed criterion can be used to set a threshold for an image.
5. Watershed_Simulation -- Contains a 1-D and 2-D simulation of how watershed segments two connecting cells in our analysis. This code can be used to reproduce Figure 4 in Section 3.3.
6. Simulated_Experiments -- Contains code to retrieve numerical results comparing accuracy and computational speed of different image smoothing algorithms. This code can be used to reproduce the figures in Section 4.
7. Nuclear_Real_Analysis -- Contains the code for reproducing the numerical results as well as the test image segmentation results for the nuclear cell channel. This code can be used to reproduce Figure 10 in Section 5.3.
8. Whole_Cell_Real_Analysis -- Contains the code for reproducing the numerical results as well as the test image segmentation results for the whole cell channel. This code can be used to reproduce Figure 11 in Section 5.4.
9. Outlier_Visualization -- Contains a demonstration of how the segmentation workflow handles sub-images with large amounts of noise affecting object detection. The code can be used to reproduce Figure S1 in the supplementary materials.

