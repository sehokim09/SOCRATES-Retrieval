 /* Copyright (C) 2023 Seho Kim

	This file is part of SOCRATES-Retrieval.

    SOCRATES-Retrieval is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SOCRATES-Retrieval is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>. */
	
#ifndef __DEFINES_H__
#define __DEFINES_H__

#define DEBUG 0
#define SA_TUNING 0

#ifndef TRUE
	#define TRUE 1
#endif

#ifndef FALSE
	#define FALSE 0
#endif

#ifndef POSITIVE
	#define POSITIVE 1
#endif

#ifndef NEGATIVE
	#define NEGATIVE 0
#endif

#define EPS 1E-15
#define EL 0.5772156649015329
#define POS_ZERO 0

#define COMPLETE 110

//#define LIGHTSPEED 3e8 // speed of light [m/s]
#define LIGHTSPEED 299792458 // speed of light [m/s]
#define BOLTZMANN 1.380649E-23 // Boltzmann's constant [J/K]

#define FILLVALUE -9999
#define POME_FAILURE -1
#define POME_DRY 1
#define POME_WET 2
#define POME_DYNAMIC 3
#define SOILMOISTMAX_EFF 1.5
#define SOILMOISTMIN_EFF -0.15
#define VWCMAX 3.0
#define VWCMIN 0.0

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define ERRNO(e) {printf("Error: %s\n", strerror(errno)); exit(ERRCODE);}

#define MIN(X,Y) ((X)<(Y) ? (X):(Y))
#define MAX(X,Y) ((X)>(Y) ? (X):(Y))

// NetCDf
#define SINGLE_DIM_NANE "single"
#define SOILSEPARATES_DIM_NAME "soilSeparates"
#define SAMPLES_DIM_NAME "samples"
#define LAYERS_DIM_NAME "layers"
#define SUBLAYERS_DIM_NAME "sublayers"
#define SOOP_DIM_NAME "SoOps"
#define RETRIEVAL_DIM_NAME "retrievals"
#define CHAR_DIM_NAME "chars"
#define VEGS_DIM_NAME "vegs"
#define VEGVARS_DIM_NAME "vegVars"

#define LAT_VAR_NAME "lat"
#define LON_VAR_NAME "lon"
#define UTCDATE_VAR_NAME "utcDate"
#define UTCTIME_VAR_NAME "utcTime"
#define IND_PRODUCTTIME_VAR_NAME "productTimeIndex"
#define IND_OBSTIME_VAR_NAME "obsTimeIndex"
#define SOILDEPTH_VAR_NAME "soilDepth"
#define SOILDEPTHPOME_VAR_NAME "soilDepthPOME"
#define SOILTYPE_VAR_NAME "soilType"
#define SOILMOIST_VAR_NAME "soilMoist"
#define SOILTEMP_VAR_NAME "soilTemp"
#define AIRTEMPMIN_VAR_NAME "airTempMin"
#define PREC_VAR_NAME "prec"
#define VEGTYPE_VAR_NAME "vegType"
#define VWCREF_VAR_NAME "VWCref"
#define VEGVAR_VAR_NAME "vegVar"

#define FLAGQCFORWARD_VAR_NAME "flagQcForward"
#define FLAGQCRET_VAR_NAME "flagQcRet"
#define PENDEPTHH_VAR_NAME "penDepthH"
#define PENDEPTHV_VAR_NAME "penDepthV"
#define PENDEPTHHRET_VAR_NAME "penDepthHRet"
#define PENDEPTHVRET_VAR_NAME "penDepthVRet"
#define REFLCOEFCOR_R_B_VAR_NAME "reflCoefCoR_R_B"
#define REFLCOEFCOI_R_B_VAR_NAME "reflCoefCoI_R_B"
#define REFLCOEFXR_R_B_VAR_NAME "reflCoefXR_R_B"
#define REFLCOEFXI_R_B_VAR_NAME "reflCoefXI_R_B"
#define REFLCOEFCOR_X_B_VAR_NAME "reflCoefCoR_X_B"
#define REFLCOEFCOI_X_B_VAR_NAME "reflCoefCoI_X_B"
#define REFLCOEFXR_X_B_VAR_NAME "reflCoefXR_X_B"
#define REFLCOEFXI_X_B_VAR_NAME "reflCoefXI_X_B"
#define REFLCOEFCOR_R_V_VAR_NAME "reflCoefCoR_R_V"
#define REFLCOEFCOI_R_V_VAR_NAME "reflCoefCoI_R_V"
#define REFLCOEFXR_R_V_VAR_NAME "reflCoefXR_R_V"
#define REFLCOEFXI_R_V_VAR_NAME "reflCoefXI_R_V"
#define REFLCOEFCOR_X_V_VAR_NAME "reflCoefCoR_X_V"
#define REFLCOEFCOI_X_V_VAR_NAME "reflCoefCoI_X_V"
#define REFLCOEFXR_X_V_VAR_NAME "reflCoefXR_X_V"
#define REFLCOEFXI_X_V_VAR_NAME "reflCoefXI_X_V"

#define REFLCO_R_B_VAR_NAME "reflect_CO_R_B"
#define REFLX_R_B_VAR_NAME "reflect_X_R_B"
#define REFLCO_X_B_VAR_NAME "reflect_CO_X_B"
#define REFLX_X_B_VAR_NAME "reflect_X_X_B"

#define REFLCO_R_V_VAR_NAME "reflect_CO_R_V"
#define REFLX_R_V_VAR_NAME "reflect_X_R_V"
#define REFLCO_X_V_VAR_NAME "reflect_CO_X_V"
#define REFLX_X_V_VAR_NAME "reflect_X_X_V"

#define REFLX_R_B_REF_VAR_NAME "reflect_X_R_B_ref"
#define REFLX_R_B_EST_VAR_NAME "reflect_X_R_B_est"
#define REFLX_R_B_ERR_VAR_NAME "reflect_X_R_B_err"
#define REFLX_X_B_REF_VAR_NAME "reflect_X_X_B_ref"
#define REFLX_X_B_EST_VAR_NAME "reflect_X_X_B_est"
#define REFLX_X_B_ERR_VAR_NAME "reflect_X_X_B_err"
#define REFLCO_R_B_REF_VAR_NAME "reflect_CO_R_B_ref"
#define REFLCO_R_B_EST_VAR_NAME "reflect_CO_R_B_est"
#define REFLCO_R_B_ERR_VAR_NAME "reflect_CO_R_B_err"
#define REFLCO_X_B_REF_VAR_NAME "reflect_CO_X_B_ref"
#define REFLCO_X_B_EST_VAR_NAME "reflect_CO_X_B_est"
#define REFLCO_X_B_ERR_VAR_NAME "reflect_CO_X_B_err"

#define REFLX_R_V_REF_VAR_NAME "reflect_X_R_V_ref"
#define REFLX_R_V_EST_VAR_NAME "reflect_X_R_V_est"
#define REFLX_R_V_ERR_VAR_NAME "reflect_X_R_V_err"
#define REFLX_X_V_REF_VAR_NAME "reflect_X_X_V_ref"
#define REFLX_X_V_EST_VAR_NAME "reflect_X_X_V_est"
#define REFLX_X_V_ERR_VAR_NAME "reflect_X_X_V_err"
#define REFLCO_R_V_REF_VAR_NAME "reflect_CO_R_V_ref"
#define REFLCO_R_V_EST_VAR_NAME "reflect_CO_R_V_est"
#define REFLCO_R_V_ERR_VAR_NAME "reflect_CO_R_V_err"
#define REFLCO_X_V_REF_VAR_NAME "reflect_CO_X_V_ref"
#define REFLCO_X_V_EST_VAR_NAME "reflect_CO_X_V_est"
#define REFLCO_X_V_ERR_VAR_NAME "reflect_CO_X_V_err"

#define REFLERR_VAR_NAME "reflErr"
#define TYPEPOME_VAR_NAME "typePOME"
#define SURPOME_VAR_NAME "surPOME"
#define BOTTPOME_VAR_NAME "bottPOME"
#define AVGPOME_VAR_NAME "avgPOME"
#define VODH_VAR_NAME "VODH"
#define VODV_VAR_NAME "VODV"
#define VWCREFMODEL_VAR_NAME "VWCrefmodel"
#define ARGHR_VAR_NAME "argHr"
#define ARGHI_VAR_NAME "argHi"
#define ARGVR_VAR_NAME "argVr"
#define ARGVI_VAR_NAME "argVi"
#define E_C_CAN_R_VAR_NAME "e_c_can_r"
#define E_C_CAN_I_VAR_NAME "e_c_can_i"

#define TYPEPOMERET_VAR_NAME "typePOMERet"
#define SURPOMERET_VAR_NAME "surPOMERet"
#define BOTTPOMERET_VAR_NAME "bottPOMERet"
#define AVGPOMERET_VAR_NAME "avgPOMERet"
#define VWCINIT_VAR_NAME "VWCinit"
#define VWCRET_VAR_NAME "VWCRet"
#define SOILMOISTPOME_VAR_NAME "soilMoistPOME"
#define SOILMOISTRET_VAR_NAME "soilMoistRet"
#define SA_TIME_VAR_NAME "sa_time"
#define SA_ITER_VAR_NAME "sa_iter"
#define SA_JUMP_VAR_NAME "sa_jump"
#define SA_COSTMIN_VAR_NAME "sa_costMin"

#define INC_MUOS_VAR_NAME "incMUOS"
#define INC_ORBCOMM_VAR_NAME "incORBCOMM"
#define INC_GPS_VAR_NAME "incGPS"

#define RUNTIME "runtime"

#define QC_GOOD 0b0
#define QC_MISSING 0b1
#define QC_FROZEN 0b10

// MPI
#define MASTER 0

// Simulation mode
#define FORWARD 0
#define INVERSE 1

// USDA Soil Texture
#define SAND 0
#define LOAMY_SAND 1
#define SANDY_LOAM 2
#define SILT_LOAM 3
#define SILT 4
#define LOAM 5
#define SANDY_CLAY_LOAM 6
#define SILTY_CLAY_LOAM 7
#define CLAY_LOAM 8
#define SANDY_CLAY 9
#define SILTY_CLAY 10
#define CLAY 11

// Number of Vegetation Types
#define NUM_VEGTYPE 5

// Vegetation Scatterer Type
#define DISK 0
#define CYLINDER 1
#define SPHERE 2

// Antenna Pattern
#define IDEAL 1
#ifndef USER_DEFINED
	#define USER_DEFINED 2
#endif
#define GAUSSIAN 3

#define ANT_PAT_TH_RANGE_DEG 180
#define ANT_PAT_PH_RANGE_DEG 360

// Air Dielectric Constant
#define EPS_DIEL_AIR 1.0

#endif /* __DEFINES_H__ */
