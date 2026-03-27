# Codes for the Article "NPMAG"

## Overview
This repository contains the **MATLAB programs/codes** required for the manuscript: **"NPMAG:Normalized Polarization Magnetic Amplitude Gradient Operator for High-Resolution Induced Polarization Imaging Overcoming High-Resistivity Layer"** . 

To ensure the reproducibility of our research, the corresponding datasets (including forward modeling simulation data and field comparative experimental data) have been permanently hosted on **Zenodo**. 

*Note: The field survey was conducted in an area with a concealed sulfide deposit in Jilin Province, China. As the area belongs to a commercial mining right, true geographic coordinates of the survey lines are withheld due to confidentiality requirements.*

## 💻 Repository Content (Codes)
This GitHub repository provides the MATLAB scripts implementing the operator/algorithm proposed in the article.


## 📊 Dataset Access (Zenodo)
The datasets required to run these codes are hosted on Zenodo and are saved in MATLAB `.mat` format. 
**Download the data here: https://doi.org/10.5281/zenodo.19249784**

The Zenodo dataset includes:
1. **Forward Modeling Data (`Forward_Modeling_Data.mat`)**: Synthetic data generated for the simulation and numerical verification of the operator.
2. **Field Electrical Data (`Electrical_Data.mat`)**: Acquired using the **ELOG-4 High-Power 3D Induced Polarization (IP) System** (Guoke Instruments Co., Ltd.).
3. **Field Magnetic Data (`Magnetic_Data.mat`)**: Three-component magnetic field data acquired using a **self-developed unshielded SERF (Spin-Exchange Relaxation-Free) magnetometer system**.

*Usage Tip: Please download the `.mat` files from Zenodo and place them in the `Data/` folder of this repository before running the scripts.*

## Requirements
* MATLAB (Tested on R2022a or newer)


## Contact
For any academic inquiries, please contact:Cheng Linhan at chenglh22@mails.jlu.edu.cn.
