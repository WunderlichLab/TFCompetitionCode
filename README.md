# TFCompetitionCode
Code for analyzing reporter expression, modeling TF binding probability, and generating figures in Waymack, et al., 2021. 
Matlab analysis of transcriptional reporter activity
Before running any of the reporter analysis scripts, you will first need to track individual transcription spots using the LivemRNA package available on Github at  https://github.com/GarciaLab/mRNADynamics with instructions for these scripts at Image processing Nikon instructions (https://docs.google.com/document/d/1t3_gyt9EjffZVtf5o6vI93UdrGLoCPojmTlul4FEL0g/edit?usp=sharing). Fill out the DataStatus excel sheet for each movie that you produce and want analyzed. Create a separate sheet for each construct, naming that sheet the name of the construct. Each column from B on should be filled in for a separate movie. The Prefix row and CompileParticles rows must be filled in for that movie’s data to be used in any of the following scripts. Fill in Prefix row with ‘Prefix=’ followed by the name you have given that movie. Ex: for the movie named ‘2018-03-5-Kr1_Kr1’ ,the Prefix row would be filled out as: Prefix=’2018-03-05-Kr1_Kr1’ To have an embryo’s data used in analysis, the CompiledParticles row must be filled out “ApproveAll’
Running "RW" scripts for compiling and plotting data as in Waymack, et al, 2021:
1.	Run “RunComparingSpotCorrAdjRW” to create SpotDiff.m files for all approved movies – The most important thing this function does is use the ElapsedTime of your movie, the nuclear tracking information (‘XXX_lin.mat’ and APDetection.mat) , and the spot tracking information (CompiledParticles) to determine for each time frame whether a given nucleus exists in the movie or not and if so, whether or not it has one or two transcriptional spots (and the corresponding fluorescence values for those spots at that time point). a. You will first need to edit the ConstructList at the beginning of the code to list the names of your constructs as you have them in DataStatus.xls b. You can also perform the underlying function ComparingSpotsCorrRWAdj for each individual movie- to do so you need to load the CompiledParticles.m, APDetection.m, and ‘XXX_lin.mat’ files from that movie’s folder, and manually run the save function at the end of the code.
2.	Run “RunBurstPropertiesSlope.m” to create BurstPropertiesSlope.mat files for all approved movies. This is the code that goes through all of the fluorescence traces to determine where transcriptional bursts are happening and records the different properties of each burst (i.e. burst size, frequency, duration, etc.) a. Here you must also edit the ConstructList at the beginning b. If you are using a different microscope, you will also want to edit the ON and OFF thresholds in SlopeBurstCalling.m (lines 26-27) to correspond to the FRNAP (or negative FRNAP for OFF threshold) calculated for your microscope as done in Lammers, et al., 2018

Once the two above files have been created (SpotDiff.m and BurstPropertiesSlope.mat) for all movies you wish to analyze, you can run the following scripts:

Comparing total mRNA produced by different constructs - TotalExpressionLevelsRW

1.	This groups all of the integrated fluorescence values associated with each allele (recorded in CompiledParticles) by AP bin, embryo, and construct. The output is a structure, AvgProdAllAP, that contains a few different things, the most useful being an array of all integrated fluorescence values of a given construct (ie across all the movies you have of that construct) organized by AP bin (each column is a different AP bin). It also uses this data to calculate the 95% CI (as well as underlying SD and SE) for each construct at each AP bin.
2.	The very end of the calculations section removes data points from the final structure where there was only 1 allele at that AP bin for the entire construct. This was done to avoid the issue of trying to calculate error bars with only one data point. If you’re not planning on imaging many embryos (>= 3) you might want to remove this section.

Compare expression boundaries between reporter configurations – APExpressionBoundsRW
1.	You will first need to run AvgmRNAProdbyEmbryoRW to generate the data structure (AvgProdAllAP.m) used for this analysis. 
2.	This script will find the positional boundaries where reporter expression is greater than or equal to half maximum expression of the homozygous configuration of the reporter construct. It uses bootstrapping to estimate error bars on these boundaries. 

Compare levels of competition between reporter configurations – CompetitionLevels_BootStrap_RW
1.	You will first need to run AvgmRNAProdbyEmbryoRW to generate the data structure (AvgProdAllAP.m) used for this analysis. 
2.	This calculates the percent higher expression in hemizygote embryos compared to homozygous embryos for the different reporters and insertion sites used. It also estimates errorbars using bootstrapping. 

Calculating fraction of nuclei with negative covariance as fx of egg length – NegCovarRatesRW
1.	This script is a subset of the TotalNoiseCalcRW script used in Waymack,et al., 2020 to analyze the total expression noise associated with different enhancer constructs. This script performs the same calculations but just retains the plotting code for looking at negative covariance rates across the AP axis for the different enhancer reporters. 

Calculating reporter distance – ReporterDistDistribution_RW
1.	Calculates the distance between MS2 reporters in a nucleus across the time of nc14. 
2.	Plots this information as in Supplemental Figure 6A (all average distances for each reporter construct) and 6B (Average reporter distances for each reporter construct broken down into nuclei that did or did not have negative covariance), and Supplemental Figure 8 (Distribution of distances between reporters for all time points for each reporter construct). 


Python TF binding probability model - Clean pBound Analysis

Modeling of TF binding probability was done using Python code, found in the TF_pBound Jupyter Notebook, run on a high performance computing cluster. The bash code used to call the modeling code is also in this notebook, along with the code to generate model-related figures. 
The Bcd gradient used in this code is from Fowlkes, et al., 2008, but different TF gradients could be used. You first need to generate the input parameters for the model (i.e. the range of Es, C, and T you want to use). The resulting file is then used as input for the modeling code, which again, is called using the bash script.   
