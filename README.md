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

## 3. Instructions

### 3.1 FW (Fresh weight and its plasticity)
•	GxE analysis in FW (FW-h2.csv)

•	BLUP of FW and its fold change (FW-BLUP.csv, FCFW-BLUP.csv)

•	Predictability of fresh weight and its plasticity (rrBLUP/)

•	Genetic correlation between fresh weight and its plasticity (gcorr/)

### 3.2 netGS across environments


More information on how to run the Matlab code of netGS, can be found https://github.com/Hao-Tong/netGS.

## 4. Reference

Predicting plasticity of rosette growth and metabolic fluxes in Arabidopsis thaliana, under review.

Please see the Methods section in this paper for the model details in mathematical equations. 
Any further questions: tong@mpimp-golm.mpg.de
