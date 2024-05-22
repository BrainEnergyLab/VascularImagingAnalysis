# VascularImagingAnalysis
This repository contains the code used to extract the in vivo vascular responses recorded using two-photon microscopy or optical imaging spectroscopy for the following papers:
1. Bonnar, O., Shaw, K., Anderle, S., Grijseels, D. M., Clarke, D., Bell, L., King, S. L., & Hall, C. N. (2023). APOE4 expression confers a mild, persistent reduction in neurovascular function in the visual cortex and hippocampus of awake mice. Journal of cerebral blood flow and metabolism : official journal of the International Society of Cerebral Blood Flow and Metabolism, 43(11), 1826–1841. https://doi.org/10.1177/0271678X231172842
2. Anderle, A., Bonnar, O., Henderson, J., Shaw, K., Chagas, A., McMullan, L., & Hall, C. N. (In Preparation). APOE4 and sedentary lifestyle synergistically impair neurovascular function in visual cortex of awake mice.
The author of these MATLAB codes is Dr Kira Shaw (who was a Postdoctoral Researcher in Catherine Hall's Brain Energy Lab). 

# Data extraction
The following MATLAB scripts were run to extract the data from 2P generated TIF (image) files into continuous traces to represent individual vessel diameter, velocity or haematocrit, or regional haemodynamic traces (e.g. SO2, Hbt, Hbr, HbO, flux, speed):

# 1. xyFWHM
  This MATLAB function loads in a TIF file XY recording of single vessel branches (labelled with a fluorescent dye) and extracts the diameter of each branch using a full width half maximum calculation. To do this a skeleton is created along the centre of the vessel, and a perpendicular line generated for each pixel of the skeleton. The intensity plot is calculated for each skeleton pixel's perpendicular line, which distinguishes between the vessel (bright pixels) and background (dark pixels). The full width half maximum calculation is applied to each intensity plot to generate the vessel diameter. A 2D variable 'cont_diam' is created (e.g. 25 x 1000) with the first dimension being the diameter for each skeleton point along the branch, and the last dimension the diameter for each data frame. Averaging across the first dimension gives you a vessel diameter reading per time frame (which can be linked to the timing onset of external stimulations or locomotion where applicable). This continuous trace is suitable for further analysis to assess vascular diameter changes over time (e.g. spontaneous vasomotion, dilation propogations along the vessel branch, stimulus- or neuronally- induced dilations, etc). 
  
# 2. linescanDiamAnalysis
  This MATLAB function loads a TIF file linescan recording taken from a single capillary (labelled with a fluorescent dye) and extracts the diameter using a FWHM calculation for each row of the linescan. Linescan images should be cropped to centre around the portion of the linescan which crosses the vessel, so that the fluorescence intensity along each row represents the change between vessel edges and background. The tif file should be labelled 'diam.TIF'. A continuous diameter trace over time is generated. 
   
# 3. linescanVelocityAnalysis
  This MATLAB function loads a TIF file linescan recording taken from a single capillary (labelled with a fluorescent dye) and extracts the red blood cell velocity. The TIF file should be cropped around the portion of the linescan which travels down the centre of the vessel, and stored as 'RBCV.tif'. The code builds on that freely shared by Patrick Drew's lab, and uses a Radon transform (across a 40ms time window) to calculate the RBCV. A continuous RBCV trace over time is generated. 
   
# 4. extractHaematocrit
  This MATLAB function loads the same TIF file as the linescanVelocityAnalysis ('RBCV.tif'), but simply calculates the percentage of pixels which are fluorescent (vascular lumen) vs not fluorescent (RBCs). The haematocrit value is a measure of the ratio of RBCs to the total volume of blood. A continuous haematocrit trace over time is generated. 
   
# 5. extractProbeData
  This MATLAB function extracts the regional haemodynamic traces including Hbr, HbO, flux, speed, SO2 recorded by the Moors OXY-VMS probe (and loaded into Microsoft Excel), and additionally calculates Hbt and CMRO2 from these traces. A continuous trace over time is generated for each of these measures. 
   
