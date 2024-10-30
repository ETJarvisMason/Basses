# Basses

This is data and code that contributed to the following paper:

Environment-driven trends in larval abundance predict fishery recruitment in two saltwater basses (in Review, ICES J of Mar Sci) 

*Erica T. Jarvis Mason<sup>1,2</sup>, William Watson<sup>2</sup>, Eric J. Ward <sup>3</sup>, Andrew R. Thompson<sup>2</sup>, Brice X. Semmens<sup>1</sup>*

<sup>1</sup>Scripps Institution of Oceanography, University of California San Diego, California, 92037 USA  
<sup>2</sup>Southwest Fisheries Science Center, NOAA Fisheries, California, 92037 USA  
<sup>3</sup>Northwest Fisheries Science Center, NOAA Fisheries, California, 98112 USA  

citation:  

**R-script files:**

1. BSB_sdmTMB_IndicesOfAbundance_Masonetal.R (sdmTMB model to generate standardized indices of larval abundance for Barred Sand Bass) Includes the following sections:
   - Define a geostatistical mesh
   - Test models accounting for various factors that can affect catchability, but not abundance
   - Check model parsimony (AIC) and predictive power (8-fold cross validation)
   - Create a grid for making predictions
   - Calculate index of abundance (with and without years with no July survey)
   - Calculate predictions for plotting spatial distribution of larval densities
   - Plot predicted larval densities through time and space
   - Calculate a center of gravity (COG) through time
   - Contains script for Figures 2, 3, and 5 (relies on loading Kelp Bass results from "01_tidy_data/" once generated with KB_sdmTMB_IndicesOfAbundance_Masonetal.R), as well as Supplemental Figures S8, S9 (map of spatiotemporal random effects), and S11A (temporal trends of standardized index of abundance containing estimates for years with no July survey)  

2. KB_sdmTMB_IndicesOfAbundance_Masonetal.R (sdmTMB model to generate standardized indices of larval abundance for Kelp Bass) Includes the following sections:
   - Define a geostatistical mesh
   - Test models accounting for various factors that can affect catchability, but not abundance
   - Check model parsimony and predictive power
   - Create a grid for making predictions
   - Calculate index of abundance (with and without years with no July survey)
   - Calculate predictions for plotting spatial distribution of larval densities
   - Plot predicted larval densities through time and space
   - Calculate a center of gravity (COG) through time
   - Contains script for Figures 2C,D, 4, 5C,D, and Supplemental Figures S10 (map of spatiotemporal random effects) and S11B (temporal trends of standardized index of abundance containing estimates for years with no July survey)  

3. CrossCorrelation_Masonetal.R (runs cross-correlation analyses between the larval indices of abundance generated above and 2 different sources of catch data (total catch estimates from all fishing modes and harvested catch reported on CPFV logbooks)
   - Cross-correlation for Barred Sand Bass
   - Cross-correlation for Kelp Bass
   - Contains script for Figure 4 and Supplemental Figure S3 (catch trends)
  
4. Basses_sdmTMB_EnvironmentalInfluence_Masonetal.R (sdmTMB model to test the effects of a variety of environmental covariates on each species larval density)
   - Test a suite of environmental covariates over the entire study period (1963-2016) - For Barred Sand Bass and for Kelp Bass
   - Check model parsimony (AIC) and predictive power (8-fold cross validation) for Time Period 1 - For Barred Sand Bass and for Kelp Bass
   - Test the same suite of environmental covariates, adding two more, over a shorter time period (1984-2016) - For Barred Sand Bass and for Kelp Bass
   - Check model parsimony (AIC) and predictive power (8-fold cross validation) for Time Period 2 - For Barred Sand Bass and for Kelp Bass
   - Compare covariate effects across species and time periods
   - Contains script for Figure 7 and for Supplemental Figures S12,13,14,15 (conditional effects plots)


