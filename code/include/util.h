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
	
#ifndef __UTIL_H__
#define __UTIL_H__

long decodeStr(char *s);
void readLineParsor(FILE *infile, double *data, const char *delimiter);
FILE *FileOpen(const char *Path, const char *File, const char *CtrlCode);
struct RandomProcessType *CreateRandomProcess(long Seed);
void DestroyRandomProcess(struct RandomProcessType *RP);
double UniformRandom(struct RandomProcessType *RP);
double GaussianRandom(struct RandomProcessType *RP);
long UniformRandomDiscrete(struct RandomProcessType *RP, int lb, int ub);
int msta1(double x,int mp);
int msta2(double x,int n,int mp);
int cbessjy01(double complex z, double complex *cj0, double complex *cj1, double complex *cy0,double complex *cy1,double complex *cj0p,
    double complex *cj1p,double complex *cy0p,double complex *cy1p);
int cbessjyna(int n, double complex z, int *nm, double complex *cj, double complex *cy, double complex *cjp, double complex *cyp);
int cbessh(int n, double complex z, int *nm, double complex *ch, double complex *chp);
double unifpdf(double x, double a, double b);
double roundzero(double x);
double newatan2(double y, double x);
double atan2TwoPi(double y, double x);
double complex sqrte(double complex z);
void sph2car(double theta, double phi, double AR, double AT, double AP, double axis[3]);
double complex VoV_complex(double complex A[3], double complex B[3]);
void M4xV4_complex(double complex M[4][4], double complex V[4], double complex W[4]);
void M2xV2_complex(double complex M[2][2], double complex V[2], double complex W[2]);
void MTxV_complex (double M[3][3], double V[3], double complex W[3]);
void M2xM2_complex(double complex A[2][2], double complex B[2][2], double complex C[2][2]);
void MxM(double A[3][3], double B[3][3], double C[3][3]);
void MxV (double M[3][3], double V[3], double W[3]);
double CopyUnitV(double V[3], double W[3]);
double signum(double x);
void RotX(double th, double R[3][3]);
void RotY(double th, double R[3][3]);
void RotZ(double th, double R[3][3]);

#endif
