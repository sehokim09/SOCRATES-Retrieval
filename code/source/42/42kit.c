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

#include "socratesret.h"


/*    This file is distributed with 42,                               */
/*    the (mostly harmless) spacecraft dynamics simulation            */
/*    created by Eric Stoneking of NASA Goddard Space Flight Center   */

/*    Copyright 2010 United States Government                         */
/*    as represented by the Administrator                             */
/*    of the National Aeronautics and Space Administration.           */

/*    No copyright is claimed in the United States                    */
/*    under Title 17, U.S. Code.                                      */

/*    All Other Rights Reserved.                                      */

///////////////////////////////////////////////////// From iokit.c of 42
FILE *FileOpen(const char *Path, const char *File, const char *CtrlCode)
{
      FILE *FilePtr;
      char FileName[1024];

      strcpy(FileName,Path);
      strcat(FileName,File);
      FilePtr=fopen(FileName,CtrlCode);
      if(FilePtr == NULL) {
         printf(">> Error opening %s: %s\n",FileName, strerror(errno));
         exit(1);
      }
      return(FilePtr);
}

////////////////////////////////////////////////////// From sigkit.c of 42
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-(1.2E-7))

/**********************************************************************/
struct RandomProcessType *CreateRandomProcess(long Seed)
{
      long j;
      long k;
      long Value;
      struct RandomProcessType *RP;

      RP = (struct RandomProcessType *)
         calloc(1,sizeof(struct RandomProcessType));

      Value = Seed;

      /* Make sure Value is in valid range */
      if (Value <= 0) {
         if (-Value < 1) Value = 1;
         else Value = -Value;
      }

      /* Do eight warmup cycles */
      for(j=0;j<8;j++) {
         /* Multiplicative Congruential Algorithm */
         k = Value/IQ;
         Value = IA*(Value-k*IQ)-IR*k;
         if (Value < 0) Value += IM;
      }
      /* Fill Lookup Table */
      for(j=0;j<NTAB;j++) {
         /* Multiplicative Congruential Algorithm */
         k = Value/IQ;
         Value = IA*(Value-k*IQ)-IR*k;
         if (Value < 0) Value += IM;

         /* Place in Table */
         RP->LookupTable[j] = Value;
      }
      RP->Index = (long) (RP->LookupTable[0]/NDIV);
      RP->PreviousValue = Value;

      return(RP);
}
/**********************************************************************/
void DestroyRandomProcess(struct RandomProcessType *RP)
{
      free(RP);
}
/**********************************************************************/
/* This portable random number generator is taken from Numerical      */
/* Recipes in C, where it is called ran1.  It has been rewritten for  */
/* readability.                                                       */
/* It returns a uniform random deviate between 0.0 and 1.0.           */
double UniformRandom(struct RandomProcessType *RP)
{
      long k;
      long Value;
      long ReadValue;
      double Output;

      Value = RP->PreviousValue;

      /* Make sure Value is in valid range */
      if (Value <= 0) {
         if (-Value < 1) Value = 1;
         else Value = -Value;
      }

      /* Multiplicative Congruential Algorithm */
      k = Value/IQ;
      Value = IA*(Value-k*IQ)-IR*k;
      if (Value < 0) Value += IM;

      /* Read from Table */
      ReadValue = RP->LookupTable[RP->Index];

      /* Remember value for next time */
      RP->PreviousValue = Value;

      /* Replace read-out value with fresh one */
      RP->LookupTable[RP->Index] = Value;
      /* Compute index for next table lookup */
      RP->Index = (long) (ReadValue/NDIV);

      /* Scale output  to (0:1) */
      Output = (double) (ReadValue*AM);
      if (Output > RNMX) Output = RNMX;
      return(Output);

}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef RNMX

/**********************************************************************/
/*  This function returns a Gaussian random deviate with zero mean,   */
/*  unit variance, using the "UniformRandom" function above to obtain */
/*  uniform deviates.  Algorithm based on Numerical Recipes Sec. 7-2  */
double GaussianRandom(struct RandomProcessType *RP)
{
      double x,y,r2,a;

      if (RP->HaveSavedValue) {
         RP->HaveSavedValue = 0;
         return(RP->SavedValue);
      }
      else {
         /* Throw (uniform) darts until hit in unit circle */
         do {
            x = 2.0*UniformRandom(RP)-1.0;
            y = 2.0*UniformRandom(RP)-1.0;
            r2 = x*x+y*y;
         } while (r2 >= 1.0 || r2 == 0.0);
         a = sqrt(-2.0*log(r2)/r2);
         /* Use one now, save one for next time */
         RP->HaveSavedValue = 1;
         RP->SavedValue=x*a;
         return(y*a);
      }
}

/////////////////////////////////////////////////////////// From mathkit.c of 42
/**********************************************************************/
/*   3x3 Matrix Product                                               */
void MxM(double A[3][3], double B[3][3], double C[3][3])
{

      C[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0];
      C[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1];
      C[0][2]=A[0][0]*B[0][2]+A[0][1]*B[1][2]+A[0][2]*B[2][2];
      C[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0];
      C[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1];
      C[1][2]=A[1][0]*B[0][2]+A[1][1]*B[1][2]+A[1][2]*B[2][2];
      C[2][0]=A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0];
      C[2][1]=A[2][0]*B[0][1]+A[2][1]*B[1][1]+A[2][2]*B[2][1];
      C[2][2]=A[2][0]*B[0][2]+A[2][1]*B[1][2]+A[2][2]*B[2][2];
}
/**********************************************************************/
/*  3x3 Matrix times 3x1 Vector                                       */
void MxV (double M[3][3], double V[3], double W[3])
{
      W[0]=V[0]*M[0][0]+V[1]*M[0][1]+V[2]*M[0][2];
      W[1]=V[0]*M[1][0]+V[1]*M[1][1]+V[2]*M[1][2];
      W[2]=V[0]*M[2][0]+V[1]*M[2][1]+V[2]*M[2][2];
}
/**********************************************************************/
/*  Copy and normalize a 3-vector.  Return its magnitude              */
double CopyUnitV(double V[3], double W[3])
{
      double A;

      A=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
      if (A > 0.0) {
         W[0] = V[0]/A;
         W[1] = V[1]/A;
         W[2] = V[2]/A;
      }
      else {
         printf("Attempted divide by zero in COPYUNITV (Line %d of mathkit.c)\n",__LINE__);
         W[0] = 0.0;
         W[1] = 0.0;
         W[2] = 0.0;
      }
      return(A);
}
/**********************************************************************/
double signum(double x)
{
      return(x>=0 ? 1.0 : -1.0);
}
////////////////////////////////////////////////////////////////////////////////////////////