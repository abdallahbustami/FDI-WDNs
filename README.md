# Water Distribution Network Stealthy Attack Simulation
This repository contains the implementation of Full-Stealth False Data Injection Attacks (FS-FDI) for water distribution networks as detailed in our paper _"Breaking the Flow and the Bank: Stealthy Cyberattacks on Water Network Hydraulics"_.

## Overview
This code simulates stealthy cyberattacks against water distribution networks, focusing on hydraulic operations and intrusion detection evasion. It demonstrates how carefully crafted sensor manipulations can impact system operations while remaining undetected.

## Main Files
- **FS_FDI_Net1.m**: Main script implementing the FS-FDI attack on the Net1 benchmark network.
- **Net1.inp**: EPANET input file for the Net1 benchmark network.
- **ExtractEPANET.m**: Extracts hydraulic parameters from EPANET files.
- **PipeHeadLoss_Linearization.m**: Implements piecewise linearization of the pipe head-loss equations.
- **PumpCurveApproximate.m**: Approximates pump characteristic curves for optimization.
- **PumpPowerConsumptionApprox.m**: Calculates power consumption for different pump operating points.
- **ConfigurationConstants.m**: Contains configuration parameters and unit conversion factors.
- **AdjustingBMatrixWithFlowJunctions.m**: Adjusts the incidence matrix for proper flow balance.
- **analyzeNetwork.m**: Provides a basic network analysis (e.g., connectivity and pump relationships).
- **checkNetworkBalance.m**: Validates the mass and energy balances of the network.

## Usage
1. **Setup:**
   - Ensure MATLAB (with the Optimization Toolbox) and the EPANET-MATLAB toolkit are installed.
   - Place all the files in one directory along with the EPANET input file (`Net1.inp`).
     
2. **Run the Simulation:**
   - Open MATLAB and navigate to the repository directory.
   - Execute the main script: `FS_FDI_Net1`
   - The script will simulate the network operation, perform state estimation, inject the FS-FDI attack, and run intrusion detection.

3. **Customization:**
   - Modify attack parameters (e.g., `attack_magnitude`, `att_start`, `att_end`) in `FS_FDI_Net1.m` to test different scenarios.
   - To simulate other attack types (HA-FDI, HU-FDI, R-FDI), adjust the corresponding sections in the code.

## Notes
- The current implementation is focused on the FS-FDI attack for research purposes. Future updates will include a more generalized implementation covering all attack types and larger network scenarios.
- While this code demonstrates FS-FDI on a benchmark network (Net1), it can be adapted for vulnerability assessment of real water distribution systems.
- For additional guidance on using the code and for detailed examples, please refer to the documentation within the source code.

## Contact
For questions or collaboration, please contact: abdallah.b.alalem.albustami@vanderbilt.edu
