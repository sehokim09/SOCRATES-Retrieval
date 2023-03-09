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
	
#ifndef __TYPES_H__
#define __TYPES_H__

#include <stdint.h>
#include <complex.h>
#include <stddef.h>

#undef I
#undef I1
#define I1 _Complex_I

struct MultiLayerType{
	// Parameters for EXPO SMP model (used in SCoBi)
	double zA; // Air Layer Thickness [m]
	double zB; // The Bottom-most Layer Thickness [m]
	double zS; // Total Layer Depth [m]
	long FitFcn; // Fitting Function for Dielectric Profile

	// Parameters for POME SMP model
	double sur; // surface effective SM for boundary condition
	double bott; // bottom-most effective SM for boundary condition
	double avg; // mean effective SM
	double D; // total depth [m]
	double z_infl; // inflection depth (cm)
	double infl_val; // inflection value
	double *SMP; // SM profile

	// Common Parameters
	double delZ; // Layer Discretization [m]
	double *z; // Layer Profile [m]
	long NSublayer; // Number of Discretized Layers (=zS/delZ)
	long NSublayer_extend;
	long NSublayer_POME;
};

struct NcDataType{
	char filename[30];
	long startYear;
	long startMonth;
	long endYear;
	long endMonth;

	// USCRN input
	int utcDateID;
	int utcTimeID;
	int soilDepthID;
	int soilTypeID;
	int soilMoistID;
	int soilTempID;
	int airTempMinID;
	int precID;
	int vegTypeID;
	int incMUOSID;
	int incORBCOMMID;
	int incGPSID;
	int VWCrefID;
	int vegVarID;
	size_t soilDepthDim;
	size_t soilTypeDim;
	size_t soilMoistDim;
	size_t soilSeparateDim;
	size_t sampleDim;
	size_t layerDim;
	size_t charDim;
	size_t vegDim;
	size_t vegVarDim;
	long *utcDate;
	long *utcTime;
	float *soilDepth;
	float soilType[3];
	float *soilMoist;
	float *soilTemp;
	float *airTempMin;
	float *prec;
	char *vegType;
	double *incMUOS;
	double *incORBCOMM;
	double *incGPS;
	double *VWCref;
	double *vegVar;

	// Forward Output
	int ncForwardID;
	int flagQcForwardID;
	int soilDepthPOMEID;
	int soilMoistPOMEID;
	int penDepthHID;
	int penDepthVID;
	int reflCoefCoR_R_B_ID;
	int reflCoefCoI_R_B_ID;
	int reflCoefXR_R_B_ID;
	int reflCoefXI_R_B_ID;
	int reflCoefCoR_X_B_ID;
	int reflCoefCoI_X_B_ID;
	int reflCoefXR_X_B_ID;
	int reflCoefXI_X_B_ID;
	int reflCoefCoR_R_V_ID;
	int reflCoefCoI_R_V_ID;
	int reflCoefXR_R_V_ID;
	int reflCoefXI_R_V_ID;
	int reflCoefCoR_X_V_ID;
	int reflCoefCoI_X_V_ID;
	int reflCoefXR_X_V_ID;
	int reflCoefXI_X_V_ID;
	int reflCo_R_B_ID;
	int reflX_R_B_ID;
	int reflCo_X_B_ID;
	int reflX_X_B_ID;
	int reflCo_R_V_ID;
	int reflX_R_V_ID;
	int reflCo_X_V_ID;
	int reflX_X_V_ID;
	int typePOMEID;
	int surPOMEID;
	int bottPOMEID;
	int avgPOMEID;
	int VODHID;
	int VODVID;
	int VWCrefmodelID;
	int argHrID;
	int argHiID;
	int argVrID;
	int argViID;
	size_t soopDim;
	size_t sublayerDim;
	int *flagQcForward;
	float *soilDepthPOME;
	double *soilMoistPOME;
	double *penDepthH;
	double *penDepthV;
	double *reflCoefCoR_R_B;
	double *reflCoefCoI_R_B;
	double *reflCoefXR_R_B;
	double *reflCoefXI_R_B;
	double *reflCoefCoR_X_B;
	double *reflCoefCoI_X_B;
	double *reflCoefXR_X_B;
	double *reflCoefXI_X_B;
	double *reflCoefCoR_R_V;
	double *reflCoefCoI_R_V;
	double *reflCoefXR_R_V;
	double *reflCoefXI_R_V;
	double *reflCoefCoR_X_V;
	double *reflCoefCoI_X_V;
	double *reflCoefXR_X_V;
	double *reflCoefXI_X_V;
	double *reflCo_R_B;
	double *reflX_R_B;
	double *reflCo_X_B;
	double *reflX_X_B;
	double *reflCo_R_V;
	double *reflX_R_V;
	double *reflCo_X_V;
	double *reflX_X_V;
	int *typePOME;
	double *surPOME;
	double *bottPOME;
	double *avgPOME;
	double *VODH;
	double *VODV;
	double *VWCrefmodel;
	double *argHr;
	double *argHi;
	double *argVr;
	double *argVi;	

	// Retrieval Output
	int flagQcRetID;
	int ncRetID;
	int iProductTimeID;
	int iObsTimeID;
	int reflX_R_B_RefID;
	int reflX_R_B_EstID;
	int reflX_R_B_ErrID;
	int reflX_X_B_RefID;
	int reflX_X_B_EstID;
	int reflX_X_B_ErrID;
	int reflCo_R_B_RefID;
	int reflCo_R_B_EstID;
	int reflCo_R_B_ErrID;
	int reflCo_X_B_RefID;
	int reflCo_X_B_EstID;
	int reflCo_X_B_ErrID;
	int reflX_R_V_RefID;
	int reflX_R_V_EstID;
	int reflX_R_V_ErrID;
	int reflX_X_V_RefID;
	int reflX_X_V_EstID;
	int reflX_X_V_ErrID;
	int reflCo_R_V_RefID;
	int reflCo_R_V_EstID;
	int reflCo_R_V_ErrID;
	int reflCo_X_V_RefID;
	int reflCo_X_V_EstID;
	int reflCo_X_V_ErrID;
	int reflErrID;
	int soilMoistRetID;
	int penDepthHRetID;
	int penDepthVRetID;
	int typePOMERetID;
	int surPOMERetID;
	int bottPOMERetID;
	int avgPOMERetID;
	int VWCRetID;
	int saTimeID;
	int saIterID;
	int saJumpID;
	int saCostMinID;
	int runtimeID;
	size_t retDim;
	int *flagQcRet;
	int *iProductTime;
	int *iObsTime;
	double *soilMoistRet;
	double *penDepthHRet;
	double *penDepthVRet;
	double *reflX_R_B_Ref;
	double *reflX_R_B_Est;
	double *reflX_R_B_Err;
	double *reflX_X_B_Ref;
	double *reflX_X_B_Est;
	double *reflX_X_B_Err;
	double *reflCo_R_B_Ref;
	double *reflCo_R_B_Est;
	double *reflCo_R_B_Err;
	double *reflCo_X_B_Ref;
	double *reflCo_X_B_Est;
	double *reflCo_X_B_Err;
	double *reflX_R_V_Ref;
	double *reflX_R_V_Est;
	double *reflX_R_V_Err;
	double *reflX_X_V_Ref;
	double *reflX_X_V_Est;
	double *reflX_X_V_Err;
	double *reflCo_R_V_Ref;
	double *reflCo_R_V_Est;
	double *reflCo_R_V_Err;
	double *reflCo_X_V_Ref;
	double *reflCo_X_V_Est;
	double *reflCo_X_V_Err;
	double reflErr;
	int *typePOMERet;
	double *surPOMERet;
	double *bottPOMERet;
	double *avgPOMERet;
	double *VWCRet;
	size_t nCycle;
	double *saTime;
	long *saIter;
	long *saJump;
	double *saCostMin;
	double runtime;
};

struct VegDataType{
	char label[20];
	long nLayer;
	double *thickness;
	double *typkndThick;
	double depth;
	long *nKind;
	char ***kind;
	int **typknd;
	int **typknd_new;
	long *part_type_idx;
	double ***dsty;
	double ***dim1;
	double ***dim2;
	double ***dim3;
	double complex ****e_c;
	double ***beginAng;
	double ***endAng;
	long ***shape;
	long nKindMax;
	long nTypeMax;
	double complex ****dKz; // Vegetation Effect
	double complex **dKzNc; // Vegetation Effect for NetCDF
	double ***VWC; // Vegetation Water Content of each type and kind, VWC = VWC[kg/m3]*depth in kg/m2
	double VWC_total; // Total vegetation water content, sum of all type and kind's
};

struct VegKindType{
	char ID[2];
	double density; // [particles/m^3]
	double dim[3]; // [m]
	double mg; // Gravimetric Water Content
	double mv; // Volumetric Water Content
	double complex *e_c;
	double beginAng; // [deg]
	double endAng; // [deg]
	double VWC; // Vegetation water content, VWC = density*volume*water content (*depth) in kg/m3
};

struct VegModelType{
	long shape;
	double theta1;
	double theta2;
	double tin;
	double pin;
	double phi;
	double freq;
	double complex epsilon;
	// Disk
	double semiMajorx;
	double semiMinorx;
	double thickness;
	// Cylinder
	double radius;
	double length;
};

struct GndDataType{
	char filename[30];
	double depth; // Layer Depth [m]
	double layer_bottom; // Layer Bottom Depth [m]
	double layer_thickness; // Layer Thickness [m]
};

struct SttDataType{
	long ExistVeg;
};

struct LocalGndType{
	// Static Data
	float sand; // sand textural component of a soil by weight (0<=S<=1)
	float clay; // clay textural component of a soil by weight (0<=C<=1)
	float bulk; // bulk density of soil sample in [g/cc] (grams/cubic centimeter)
	float porosity; // soil porosity (GLDAS 0.25 deg resolution) (NOT USED)
	uint8_t landMask;
	uint8_t landType;

	// Look-up data
	float fc; // soil field capacity (m^3/m^3)
	float poro; // soil porosity (m^3/m^3)
	float resid; // residual soil moisture (m^3/m^3)
	float wp; // wilting point (m^3/m^3)
	uint8_t soilType;

	// Dynamic Data
	uint8_t snowCover;
	uint8_t freezeThaw;
	float soilTemp1; // temperature in [C]

	// Ground Data
	float *VSM; // volumetric water content in g*cm^-3 (NLayer)
	double VWC; // vegetation water content in kg/m2

	// Dielectric Profile Realization Method
	long model;

	// Dielectric Constant
	double complex *e_c; // Complex dielectric constant (NLayer)
	double complex *e_c_profile; // Dielectric Constant Profile Generated by Fitting Function (MultiLayer->NSublayer)

	// Soil Moisture Profile
	double *SMP;

	// Reflection
	double h; // Effective Roughness Parameter h = (2*RMSH*k0)^2, k0 = TwoPi/Wavelength[cm] in [cm]
	double complex FVV; // Standard Fresnel coefficient
	double complex FHH;
	double complex FLR;
	double complex FRR;
	double gammaV; // reflectivity (Vertical polarization)
	double gammaH; // reflectivity (Horizontal polarization)
	double complex RcV; // reflection coefficient
	double complex RcH;

	// Penetration depth
	double pdv;
	double pdh;
};

struct BistaticType{
	double r_tr; // Slant Range from Transmitter to Receiver [m]
	double r_ts; // Slant Range from Transmitter to Specular Point [m]
	double r_sr; // Slant Range from Specular Point to Receiver [m]
	double r_tsr; // r_ts + r_sr
	double idn[3]; // propagation vector (i_d^-)
	double isn[3]; // propagation vector (i_s^-)
	double osp[3]; // propagation vector (o_s^+)
	double Tgs[3][3]; // Transformation matrix for transforming a vector from the ground frame to local (specular) ground system
	double Tgr[3][3]; // Transformation matrix for transforming a vector from the ground frame to receiver system
	double Tgt[3][3]; // Transformation matrix for transforming a vector from the ground frame to transmitters system
	double complex u_tr[2][2]; // Rotation Matrix (Transmitter to Receiver)
	double complex u_ts[2][2]; // Rotation Matrix (Transmitter to Specular Point)
	double complex u_sr[2][2]; // Rotation Matrix (Specular Point to Receiver)
	double complex g_rt[2][2]; // Normalized Voltage Antenna Pattern of Receiver in the transmitter direction
	double complex g_rs[2][2]; // Normalized Voltage Antenna Pattern of Receiver in the specular direction
	char polRx; // Receive Antenna Polarization
	char polTx; // Transmit Antenna Polarization
	double AngT2S_sf_th0;
	double complex ArgH;
	double complex ArgV;

	double Freq; // [Hz]
	double Wavelength; // [m]
	double Wavenumber;
	double th; // Incidence Angle [rad]
	double EIRP_Tx; // [dB]
	double EIRP_Tx_Lin; // Decibel to Natural
	double G_Rx; // [dB]
	double G_Rx_Lin; // Decibel to Natural
	double t_coh; // Coherent Integration Time [sec]
	double Bandwidth; // Noise bandwidth of the front-end [Hz]
	double Tn; // Noise temerature [K]
	long N_inc; // Number of Independent Samples
};

struct FixedObsType{
	long ID;

	double Freq; // [Hz]
	double Wavelength; // [m]
	double EIRP_Tx; // [dB]
	char polRx; // Receive Antenna Polarization
	char polTx; // Transmit Antenna Polarization
	double G_Rx; // [dB]
	double h; // Effective Roughness Parameter h = (2*RMSH*k0)^2, k0 = TwoPi/Wavelength[cm] in [cm]
	double t_coh; // Coherent Integration Time [sec]
	double Bandwidth; // Noise bandwidth of the front-end [Hz]
	double Tn; // Noise temerature [K]
	double sres_req; // Spatial Resolution Requirement [m]
	double t_inc; // Max Incoherent Integration Time [sec] - set by sres_req
	long N_inc; // Number of Independent Samples - floor(t_inc/t_coh)

	// Bistatic Geometry
	double rTx; // Transmitter range from Earth's center [m]
	double hTx; // Transmitter altitude [m]
	double elTx; // Tx Elevation = 90 - Incidence angle [rad]
	double thTx; // Tx Incidence angle [rad]
	double phTx; // Tx Azimuth angle [rad]
	double hRx; // Receiver altitude [m]
	double elRx; // Rx Elevation = 90 - Incidence angle [rad]
	double thRx; // Rx Antenna looking angle (angle of incidence) [rad]
	double phRx; // Rx Azimuth angle of receiver position [rad]
	double vcircRx; // Rx Circular Velocity [m/s] in ECEF frame
	double S0x; // Distance of specular point away from the receiver ground projection [m]
	double x1; // Distance to the center of the ellipse [m]
	long Nfz; // Number of Fresnel zones
	double ax1; // semi-major axis of Fresnel zone [m]
	double by1; // semi-minor axis of Fresnel zone [m]
	long AntTag; // Antenna Type - USER_DEFINED, IDEAL
	char AntFileName[4][30]; // Antenna Pattern File Name
	double ***AntPattern; // Antenna Pattern
	double AntPatRes; // Antenna Pattern Resolution [deg]
	double hpbw; // Receive Antenna Half-Power Bandwidth [deg]
	double SLL; // Side-lobe Level of Antenna Pattern [dB]
	double XPL; // Cross-Polarization Level of Antenna Pattern [dB]
	double AngT2R_rf[2];
	double AngS2R_rf[2];
	double AngT2S_sf[2];
	char orbitName[5];
};

struct InputType{
	struct FixedObsType *F;
	double RMSH; // Root-mean-square height of roughness of the soil surface [cm]
};

struct OutputType{
	// Factor K
	double complex K;
	double complex Kd;
	double complex Kc;

	// Field
	// Direct
	double complex b_d1[2];
	double complex b_d2[2];
	// Specular
	double complex b_coh1[2];
	double complex b_coh2[2];
	double complex b_coh1b[2];
	double complex b_coh2b[2];
	double complex b_coh1v[2];
	double complex b_coh2v[2];

	// Power - Reflectivity
	// Direct
	double P_d1[4];
	double P_d2[4];
	// Specular
	double P_coh1[4];
	double P_coh2[4];
	double P_coh1_err[4]; // error included 
	double P_coh2_err[4]; // error included
	double P_coh1b[4];
	double P_coh2b[4];
	double P_coh1v[4];
	double P_coh2v[4];
	double Gamma_Std_mod_lin[2];
	double Gamma_Std_mod_fused[2];
	double Gamma_Std_mod_fused_avg[2];

	// Power
	double PWR_d1[2];
	double PWR_d2[2];
	double PWR_coh1[2];
	double PWR_coh2[2];

	// SNR
	double SNR0[2]; // for detectability
	double SNR_d1[2];
	double SNR_d2[2];
	double SNR_coh1[2];
	double SNR_coh2[2];

	// Delay and Doppler
	double delayD; // delay for direct signal [sec]
	double delayR; // delay for reflected signal [sec]
	double dopplerD; // doppler for direct signal [Hz]
	double dopplerR; // doppler for reflected signal [Hz]

	// Differential Phase
	double phaseRD;

	// Optical depth
	double tauH;
	double tauV;

	// Range Corrected Gain
	double RCG_R; // in dB
	double RCG_D; // in dB
	double Gamma_eff;
};

struct SimAnnealType{
	// Inputs
	long nDim;
	long nIter;
	long nTemp;
	long nTrap;
	long nGen; // for ASA
	long nAccepted; // for ASA
	double *c; // for ASA
	double c_cost; // for ASA
	double f_min;
	double f_err;
	double *lb;
	double *ub;
	double T0; // Initial temperature
	double k_accept;
	double *x_opt;
	double f;
	double f_init;
	double f_opt;
	double *T_gen;
	double *s;
	// Outputs
	long iIter;
	long cnt_trapped;
};

struct RetrievalType{
	double *VSM_ref;
	double *VSM_est;
	double *VSM_err;
	double *R_CO_ref;
	double *R_CO_est;
	double *R_X_ref;
	double *R_X_est;
	double stddev;
	long period;
	long window;
	double complex *ArgH;
	double complex *ArgV;
	double VWC_ref;
	// POME model
	double sur_ref;
	double avg_ref;
	double bott_ref;
	double sur_est;
	double avg_est;
	double bott_est;
	double sur_err;
	double avg_err;
	double bott_err;
};

/* Types from 42 */
struct RandomProcessType {
   /* For UniformRandom */
   long Index;
   long PreviousValue;
   long LookupTable[32];

   /* For GaussianRandom */
   long HaveSavedValue;
   double SavedValue;
};

#endif
