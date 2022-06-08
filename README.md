# Flux_Plasticity
# Plant growth and metabolic flux plasticity
Predicting plasticity of rosette growth and metabolic fluxes in Arabidopsis thaliana

## 1. System Requirements

1.1 Hardware requirements
All codes require only a standard computer with enough RAM to support the in-memory operations.

1.2 Software requirements

1.2.1 Operation System (OS) Requirements
All codes can be run under Linux or Windows OS. The code here has been tested on Ubuntu 14.04 LTS. 

1.2.2 Programming software
R language and Matlab are needed to run all codes. The code here has been tested on R (version 3.4.4) and Matlab (version 8.5 (R2015a)). 

## 2. Installation Guide

All installation guide present below is only for Linux OS. For Windows OS, it might be slightly different. 

2.1 R language
R language can be installed automatically under Linux OS using the demand below: 

sudo apt-get install r-base

In order to install R, you need the password of sudo account. The installation could take few minutes on a normal desktop computer. 

After successful install R, you could simply type 'R' in the terminal to open a workspace under R, then run the code written in R language. 

The netGS depends on two R packages: rrBLUP and R.matlab, which can be installed using the commands below in R workspace:

install.packages(rrBLUP)
install.packages(R.matlab)

2.2 Matlab

Matlab can be installed follow the instruction in graphical user interface (GUI) after download Matlab installation compress file. A license is needed to install Matlab. The installation could take up to half an hour on a normal desktop computer. 

After successful install Matlab, you could simply type 'matlab' in the terminal to open a workspace under Matlab, then run the code written in Matlab language. 

The netGS depends on Matlab toolbox cobra, solver glpk and cvx. 

To install cobra toolbox, you can clone the repository in the terminal using:

git clone --depth=1 https://github.com/opencobra/cobratoolbox.git cobratoolbox

then change to the folder cobratoolbox/ and run the command below in Matlab:

initCobraToolbox

To install the solver glpk in Matlab, simply download glpkcc.mex* files and unzip it, then add the path where you put glpkmex to Matlab. 

To install the solver cvx in Matlab, similarly, download the package and unzip it, then add the path where you put cvx to Matlab, and run the below command in Matlab:

cvx_setup

## 3. Folder instructions

### 3.1 FW (Fresh weight and its plasticity)
•	GxE analysis in FW (FW-h2.csv)
•	BLUP of FW and its fold change (FW-BLUP.csv, FCFW-BLUP.csv)

•	Predictability of fresh weight and its plasticity (rrBLUP/)
•	Genetic correlation between fresh weight and its plasticity (gcorr/)

More information on how to run the Matlab code of netGS, can be found https://github.com/Hao-Tong/netGS.

### 3.2 netGS across environments

All codes and data are under the folder named ‘netGS_env’. 
Please change the path name in the first section of the code to the path on your computer. If you want to use a new folder, please make sure the results of rrBLUP in previous section (fluxpredict_*.csv and fluxprediction_cor.csv) are copied here or indicate the path to access these results. 

•	FluxDist.m
This code should be run in Matlab. 
This code is used to estimate the reference flux distribution of Col-0 in two environments (optimal and low N) based on FBA and minimization of quadratic program. 

The output including: the reference flux distribution in optimal N condition both in .mat format (fluxcol0_optN.mat) and .csv format (fluxcol0_optN.csv), the reference flux distribution in low N condition both in .mat format (fluxcol0_lowN.mat) and .csv format (fluxcol0_lowN.csv), and the index for nonzero flux in reference flux distribution in optimal N (nonzeroid.csv). 

The code could take few seconds for the example on a normal desktop computer. 

•	Biomass.m 
This code should be run in Matlab. 
This code is used to estimate the genotype flux distribution in steady-state in low N condition by minimization of quadratic program. The biomass flux included in this flux distribution is used as the final biomass prediction in low N condition. 

The output including: the final flux distribution for each replicate and each fold in cross validations (biomasspredict_lowN_r*_f*.csv).

The code could take up to an hour for the example on a normal desktop computer.

•	Correlation.R
This code should be run in R. 
This code is used to check the correlation coefficient between predicted biomass in low N condition in netGS and measured biomass in low N condition as the netGS prediction accuracy across environments. 

The output including: the correlation coefficient for each replicate and each fold in cross validations (biomcorr_lowN.csv).

The code could take few seconds for the example on a normal desktop computer. 

### 3.3 netGS robustness

All codes and data are under the folder named ‘netGS_robust’. 
Please change the path name in the first section of the code to the path on your computer. 

•	FluxDist.m
This code should be run in Matlab. 
This code is used to test the robustness of flux distribution of Col-0 based on FBA and minimization of quadratic program. 

The output including: the reference flux distribution both in .mat format (fluxcol0.mat) and .csv format (fluxcol0.csv), the random sampled reference flux distribution both in .mat format (fluxcol0_sample.mat) and .csv format (fluxcol0_sample.csv), the random sampled reference flux distribution in steady-state as robustness in both in .mat format (fluxcol0_robust.mat) and .csv format (fluxcol0_robust.csv), and the index for nonzero flux in reference flux distribution (nonzeroid.csv). 

The code could take few minutes for the example on a normal desktop computer. 

## 4. Reference

Predicting plasticity of rosette growth and metabolic fluxes in Arabidopsis thaliana, under review.

Please see the Methods section in this paper for the model details in mathematical equations. 
Any further questions: tong@mpimp-golm.mpg.de
