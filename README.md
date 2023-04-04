# SoOp-R Retrieval Sensitivity Analysis for Subsurface Soil Moisture and Vegetation Water Content Using SOCRATES-Retrieval

## Overview
SOCRATES-Retrieval is a retrieval module of SOCRATES (Signals of Opportunity Constellation and Remote sensing Analysis Tool for Earth Science) to simulate the retrieval of subsurface soil moisture (RZSM) and vegetation water content (VWC) using synthetic observations of multi-frequency/polarimetric signals of opportunity (SoOp). This simulation tool integrates the SCoBi model and the Principle of Maximum Entropy (POME) model to provide a bistatic SoOp-R environment and resulting measurements in a comprehensive manner. The current release of SOCRATES-Retrieval boasts the following capabilities:

  - Soil moisutre profile and VWC retrieval using multi-frequency/polarimetric SoOp.
  - Sensitivity analysis of retrieval accuracy to SoOp system parameters.
  - Antenna property realizations including gain pattern.
  - Geophysical data realizations including soil moisture and vegetation.

Use SOCRATES for more functionality such as generic SoOp coverage analysis and trade study of SoOp constellation design.

## Installation

SOCRATES-Retrieval supports Linux. See the manual, "SOCRATES-Retrieval\_user\_manual.pdf" in the ./doc folder for details.

To install from the Git repository:

  1. Clone the code from Git: `git clone git@github.itap.purdue.edu:RadioNavigationLab/SOCRATES-Retrieval.git`.
  2. Open a terminal on the ./code directory and run `make` to compile the simulator as a default mode.
  3. An executable file "SOCRATES-Retrieval" should be created in the top directory.

## Getting Started
See the manual, "SOCRATES-Retrieval\_user\_manual.pdf" in the ./doc folder.

## Update History

  - v1.0.0 - (Mar 09, 2023) First commit.

## Citation
S. Kim, J. L. Garrison and M. Kurum, "Retrieval of Subsurface Soil Moisture and Vegetation Water Content From Multifrequency SoOp Reflectometry: Sensitivity Analysis," in IEEE Transactions on Geoscience and Remote Sensing, vol. 61, pp. 1-18, 2023, Art no. 4502818, doi: 10.1109/TGRS.2023.3284800.

Seho Kim, James L. Garrison, and Mehmet Kurum, "SoOp-R Retrieval Sensitivity Analysis for Subsurface Soil Moisture and Vegetation Water Content Using SOCRATES-Retrieval" [Source Code], Code Ocean, DOI:https://doi.org/10.24433/CO.5405959.v1
