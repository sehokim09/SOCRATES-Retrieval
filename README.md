# SOCRATES-Retrieval
SoOp Constellation Remote sensing Analysis Tool for Earth Science - Retrieval module

SOCRATES-Retrieval is a retrieval module of SOCRATES to simulate the retrieval of subsurface soil moisture (RZSM) and vegetation water content (VWC) using synthetic observations of multi-frequency/polarimetric signals of opportunity (SoOp). This simulation tool integrates the SCoBi model and the Principle of Maximum Entropy (POME) model to provide a bistatic SoOp-R environment and resulting measurements in a comprehensive manner. The current release of SOCRATES-Retrieval boasts the following capabilities:

  - Soil moisutre profile and VWC retrieval using multi-frequency/polarimetric SoOp.
  - Sensitivity analysis of retrieval accuracy to SoOp system parameters.
  - Antenna property realizations including gain pattern.
  - Geophysical data realizations including soil moisture and vegetation.

Use SOCRATES for more functionality such as generic SoOp coverage analysis and trade study of SoOp constellation design.

## Installation

SOCRATES-Retrieval supports Linux.

To install from the Git repository:

  1. Clone the code from Git: `git clone git@github.itap.purdue.edu:RadioNavigationLab/SOCRATES-Retrieval.git`.
  2. Open a terminal on the top-level directory and run `make` to compile the simulator as a default mode. See the manual, "SOCRATES-Retrieval\_user\_manual.pdf" in the ./doc folder for details.

## Getting Started
See the manual, "SOCRATES-Retrieval\_user\_manual.pdf" in the ./doc folder.

## Update History

  - v1.0.0 - (Mar 09, 2023) First commit.

## Citation

Seho Kim, James L. Garrison, and Mehmet Kurum, “Retrieval of Subsurface Soil Moisture and Vegetation Water Content from Multi-frequency SoOp-Reflectometry: Sensitivity Analysis,” IEEE Transactions of Geoscience and Remote Sensing, DOI: https://doi.org/10.36227/techrxiv.22114433.v1, Manuscript submitted.

