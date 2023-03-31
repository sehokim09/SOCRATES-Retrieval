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

#ifndef __SMATRET_H__
#define __SMATRET_H__

/* Disable extern keyword to declare globals */
#ifdef DECLARE_GLOBALS
   #define EXTERN
#else
   #define EXTERN extern
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "defines.h"
#include "types.h"
#include "util.h"
#include "mpi.h"
#include "netcdf.h"
#include "sys/stat.h"

EXTERN char InPath[30];
EXTERN char ForwardPath[40];
EXTERN char InversePath[40];
EXTERN char InverseSubPath[50];
EXTERN char AntPath[30];

EXTERN char SimMode;
EXTERN char EndFlag;
EXTERN struct InputType *InputSoOp;
EXTERN struct FixedObsType *Fixed;
EXTERN struct BistaticType *Bistatic;
EXTERN struct SttDataType *SttData;
EXTERN struct GndDataType *GndData;
EXTERN struct VegDataType *VegData;
EXTERN struct VegKindType *VegKind;
EXTERN struct VegModelType *VegModel;
EXTERN struct MultiLayerType *MultiLayer;
EXTERN struct LocalGndType *Gnd;
EXTERN struct OutputType *Out;
EXTERN long NData;
EXTERN long NSoOp, NPol, flagMixPol;
EXTERN char RxPol[2];
EXTERN long NLayer; // Number of Soil Layers
EXTERN long NVeg, NKind; // Number of Vegetation Types, Kinds

EXTERN struct SimAnnealType *SimAnneal;
EXTERN struct RetrievalType *Ret;
EXTERN struct NcDataType *Nc;
EXTERN double **b_LUT; // Look-up table for VOD = b * VWC

/* Math Basics */
EXTERN double Pi, TwoPi, HalfPi, D2R, R2D;

/* Physical Constants */
EXTERN double EarthRad, EarthMu, EarthW;

/* Variables from 42 */
// Master Random Process
EXTERN struct RandomProcessType *RNG;

/* Functions for Signals of Opportunity */
// obs.c
void updateProgress(long cnt);

// geometry.c
void updateBistaticFixed(struct FixedObsType *F);
void calcTxGeometry(struct FixedObsType *F);
void calcGeometryFixed(struct FixedObsType *F);

// antenna.c
void calcAntPolRot(void);
void calcMueller(double complex u[2][2], double complex U[4][4]);
void LoadAntPattern_Fixed(struct FixedObsType *F);
void LoadAntPattern_Fixed_GG(struct FixedObsType *F);
void LookupAntPattern(double ***antPat, double res, double th, double ph, double complex g[2][2]);

// ground.c
void calcRcMulti(void);
double complex calcDielDebye(double freq);
void calcDielMironov(void);
double complex calcDielUlabyElRayesMv(double freq, double mv);
void calcMLDielProfile(void);
int SMP_POME_main(double sur, double bott, double avg);
void LookupSoiltype(void);
void qsimp_ThAvg(double a, double b, double eps, struct VegModelType *V, double complex s[4]);
void qsimp_pxf(double phi, double tol, struct VegModelType *V, double complex s[4]);
void calcPropagationNetCDF(long iSoOp, long iVeg, long iSample);

// product.c
void calcSpecularTermForward(void);
void calcSpecularTermSimple(struct FixedObsType *F, double VWC);

// init.c
void Welcome(void);
void InitSOCRATES_Ret(int argc,char **argv);

// retrieval.c
void writeNcRet(void);
void writeNcForward(void);
void createNcRetPeriod(long year, char *stationName, long period, long window);
void createNcForward(long year, char *stationName);
void updateNcVeg(long iSample);
void readNcForward(long year, char *stationName);
void readNcUSCRN(long year, char *stationName);
void ObsFixedNetCDF(void);
void InitRetrieval(void);
void InitSA(void);
void RetrievalNetCDF(void);

// 42kit.c
FILE *FileOpen(const char *Path, const char *File, const char *CtrlCode);
struct RandomProcessType *CreateRandomProcess(long Seed);
void DestroyRandomProcess(struct RandomProcessType *RP);
double UniformRandom(struct RandomProcessType *RP);
double GaussianRandom(struct RandomProcessType *RP);
void MxM(double A[3][3], double B[3][3], double C[3][3]);
void MxV (double M[3][3], double V[3], double W[3]);
double CopyUnitV(double V[3], double W[3]);
double signum(double x);

/////////////////////////////////////////////

#endif /* __SMATRET_H__ */
