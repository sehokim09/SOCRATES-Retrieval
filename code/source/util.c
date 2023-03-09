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

long decodeStr(char *s)
{
    if (!strcmp(s,"FALSE")) return FALSE;
    else if (!strcmp(s,"TRUE")) return TRUE;

    else if (!strcmp(s,"USER_DEFINED")) return USER_DEFINED;
    else if (!strcmp(s,"IDEAL")) return IDEAL;
    else if (!strcmp(s,"GAUSSIAN")) return GAUSSIAN;
    else {
        printf(">> Bogus input %s in decodeStr (init.c:%d)\n",s,__LINE__);
        exit(1);
    }
}

void readLineParsor(FILE *infile, double *data, const char *delimiter)
{
    char *token;
    char line[3600];
    int i;

    fgets(line, sizeof line, infile);
    token = strtok(line, delimiter);
    i = 0;
    while(token!=NULL){	
        data[i++] = atof(token);
        token = strtok(NULL, delimiter);
    }
}


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

//////////////////////////////////////////////////////////////////////////

long UniformRandomDiscrete(struct RandomProcessType *RP, int lb, int ub)
{
    return (long) floor(lb+(ub-lb+1)*UniformRandom(RP));
}

int msta1(double x,int mp)
{
    double a0,f0,f1,f;
    int i,n0,n1,nn;

    a0 = fabs(x);
    n0 = (int)(1.1*a0)+1;
    f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-mp;
    n1 = n0+5;
    f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-mp;
    for (i=0;i<20;i++) {
        nn = n1-(n1-n0)/(1.0-f0/f1);
        f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-mp;
        if (abs(nn-n1) < 1) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}

int msta2(double x,int n,int mp)
{
    double a0,ejn,hmp,f0,f1,f,obj;
    int i,n0,n1,nn;

    a0 = fabs(x);
    hmp = 0.5*mp;
    ejn = 0.5*log10(6.28*n)-n*log10(1.36*a0/n);
    if (ejn <= hmp) {
        obj = mp;
        n0 = (int)(1.1*a0);
        if (n0 < 1) n0 = 1;
    }
    else {
        obj = hmp+ejn;
        n0 = n;
    }
    f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-obj;
    n1 = n0+5;
    f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-obj;
    for (i=0;i<20;i++) {
        nn = n1-(n1-n0)/(1.0-f0/f1);
        f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-obj;
        if (abs(nn-n1) < 1) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn+10;
}

int cbessjy01(double complex z, double complex *cj0, double complex *cj1, double complex *cy0,double complex *cy1,double complex *cj0p,
    double complex *cj1p,double complex *cy0p,double complex *cy1p)
{
    double complex z1,z2,cr,cp,cs,cp0,cq0,cp1,cq1,ct1,ct2,cu;
    double a0,w0,w1;
    int k,kz;

    static double a[] = {
        -7.03125e-2,
         0.112152099609375,
        -0.5725014209747314,
         6.074042001273483,
        -1.100171402692467e2,
         3.038090510922384e3,
        -1.188384262567832e5,
         6.252951493434797e6,
        -4.259392165047669e8,
         3.646840080706556e10,
        -3.833534661393944e12,
         4.854014686852901e14,
        -7.286857349377656e16,
         1.279721941975975e19};
    static double b[] = {
         7.32421875e-2,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1,
         5.513358961220206e2,
        -1.825775547429318e4,
         8.328593040162893e5,
        -5.006958953198893e7,
         3.836255180230433e9,
        -3.649010818849833e11,
         4.218971570284096e13,
        -5.827244631566907e15,
         9.476288099260110e17,
        -1.792162323051699e20};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947,
         1.215978918765359e2,
        -3.302272294480852e3,
         1.276412726461746e5,
        -6.656367718817688e6,
         4.502786003050393e8,
        -3.833857520742790e10,
         4.011838599133198e12,
        -5.060568503314727e14,
         7.572616461117958e16,
        -1.326257285320556e19};
    static double b1[] = {
        -0.1025390625,
         0.2775764465332031,
        -1.993531733751297,
         2.724882731126854e1,
        -6.038440767050702e2,
         1.971837591223663e4,
        -8.902978767070678e5,
         5.310411010968522e7,
        -4.043620325107754e9,
         3.827011346598605e11,
        -4.406481417852278e13,
         6.065091351222699e15,
        -9.833883876590679e17,
         1.855045211579828e20};

    a0 = cabs(z);
    z2 = z*z;
    z1 = z;
    if (a0 == 0.0) {
        *cj0 = 1;
        *cj1 = 0;
        *cy0 = -1e308;
        *cy1 = -1e308;
        *cj0p = 0;
        *cj1p = 0.5;
        *cy0p = 1e308;
        *cy1p = 1e308;
        return 0;
    }
    if (creal(z) < 0.0) z1 = -z;
    if (a0 <= 12.0) {
        *cj0 = 1;
        cr = 1;
        for (k=1;k<=40;k++) {
            cr *= -0.25*z2/(double)(k*k);
            *cj0 += cr;
            if (cabs(cr) < cabs(*cj0)*EPS) break;
        }
        *cj1 = 1;
        cr = 1;
        for (k=1;k<=40;k++) {
            cr *= -0.25*z2/(k*(k+1.0));
            *cj1 += cr;
            if (cabs(cr) < cabs(*cj1)*EPS) break;
        }
        *cj1 *= 0.5*z1;
        w0 = 0.0;
        cr = 1;
        cs = 0;
        for (k=1;k<=40;k++) {
            w0 += 1.0/k;
            cr *= -0.25*z2/(double)(k*k);
            cp = cr*w0;
            cs += cp;
            if (cabs(cp) < cabs(cs)*EPS) break;
        }
        *cy0 = M_2_PI*((clog(0.5*z1)+EL)* (*cj0)-cs);
        w1 = 0.0;
        cr = 1;
        cs = 1;
        for (k=1;k<=40;k++) {
            w1 += 1.0/k;
            cr *= -0.25*z2/(k*(k+1.0));
            cp = cr*(2.0*w1+1.0/(k+1.0));
            cs += cp;
            if (cabs(cp) < cabs(cs)*EPS) break;
        }
        *cy1 = M_2_PI*((clog(0.5*z1)+EL)*(*cj1)-1.0/z1-0.25*z1*cs);
    }
    else {
        if (a0 >= 50.0) kz = 8;         // can be changed to 10
        else if (a0 >= 35.0) kz = 10;   //   "      "     "  12
        else kz = 12;                   //   "      "     "  14
        ct1 = z1 - M_PI_4;
        cp0 = 1;
        for (k=0;k<kz;k++) {
            cp0 += a[k]*cpow(z1,-2.0*k-2.0);
        }
        cq0 = -0.125/z1;
        for (k=0;k<kz;k++) {
            cq0 += b[k]*pow(z1,-2.0*k-3.0);
        }
        cu = csqrt(M_2_PI/z1);
        *cj0 = cu*(cp0*ccos(ct1)-cq0*csin(ct1));
        *cy0 = cu*(cp0*csin(ct1)+cq0*ccos(ct1));
        ct2 = z1 - 0.75*M_PI;
        cp1 = 1;
        for (k=0;k<kz;k++) {
            cp1 += a1[k]*cpow(z1,-2.0*k-2.0);
        }
        cq1 = 0.375/z1;
        for (k=0;k<kz;k++) {
            cq1 += b1[k]*cpow(z1,-2.0*k-3.0);
        }
        *cj1 = cu*(cp1*ccos(ct2)-cq1*csin(ct2));
        *cy1 = cu*(cp1*csin(ct2)+cq1*ccos(ct2));
    }
    if (creal(z) < 0.0) {
        if (cimag(z) < 0.0) {
            *cy0 -= 2.0*I1*(*cj0);
            *cy1 = -(*cy1-2.0*I1*(*cj1));
        }
        else if (cimag(z) > 0.0) {
            *cy0 += 2.0*I1*(*cj0);
            *cy1 = -(*cy1+2.0*I1*(*cj1));
        }
        *cj1 = -*cj1;
    }
    *cj0p = -*cj1;
    *cj1p = *cj0-*cj1/z;
    *cy0p = -*cy1;
    *cy1p = *cy0-*cy1/z;
    return 0;
}

int cbessjyna(int n, double complex z, int *nm, double complex *cj, double complex *cy, double complex *cjp, double complex *cyp)
{
    double complex cbj0,cbj1,cby0,cby1,cj0,cjk,cj1,cf,cf1,cf2;
    double complex cs,cg0,cg1,cyk,cyl1,cyl2,cylk,cp11,cp12,cp21,cp22;
    double complex ch0,ch1,ch2;
    double complex cbjp, cbyp;
    double a0,yak,ya1,ya0,wa;
    int m,k,lb,lb0;
    
    if (n < 0) return 1;
    a0 = cabs(z);
    *nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            cj[k] = 0;
            cy[k] = -1e308;
            cjp[k] = 0;
            cyp[k] = 1e308;
        }
        cj[0] = 1;
        if(n>0)
            cjp[1] = 0.5;
        return 0;
    }
    if(n<1){
        cbessjy01(z,&cj[0],&cbj1,&cy[0],&cby1,&cjp[0],&cbjp,&cyp[0],&cbyp);
        cbj0 = cj[0];
        cby0 = cy[0];
    }
    else{
        cbessjy01(z,&cj[0],&cj[1],&cy[0],&cy[1],&cjp[0],&cjp[1],&cyp[0],&cyp[1]);
        cbj0 = cj[0];
        cbj1 = cj[1];
        cby0 = cy[0];
        cby1 = cy[1];
    }
    if (n <= 1) return 0;
    if (n < (int)0.25*a0) {
        cj0 = cbj0;
        cj1 = cbj1;
        for (k=2;k<=n;k++) {
            cjk = 2.0*(k-1.0)*cj1/z-cj0;
            cj[k] = cjk;
            cj0 = cj1;
            cj1 = cjk;
        }
    }
    else {
        m = msta1(a0,200);
        if (m < n) *nm = m;
        else m = msta2(a0,n,15);
        cf2 = 0;
        cf1 = 1.0e-100;
        for (k=m;k>=0;k--) {
            cf = 2.0*(k+1.0)*cf1/z-cf2;
            if (k <=*nm) cj[k] = cf;
            cf2 = cf1;
            cf1 = cf;
        }
        if (cabs(cbj0) > cabs(cbj1)) cs = cbj0/cf;
        else cs = cbj1/cf2;
        for (k=0;k<=*nm;k++) {
            cj[k] *= cs;
        }
    }
    for (k=2;k<=*nm;k++) {
        cjp[k] = cj[k-1]-(double)k*cj[k]/z;
    }
    ya0 = cabs(cby0);
    lb = 0;
    cg0 = cby0;
    cg1 = cby1;
    for (k=2;k<=*nm;k++) {
        cyk = 2.0*(k-1.0)*cg1/z-cg0;
        yak = cabs(cyk);
        ya1 = cabs(cg0);
        if ((yak < ya0) && (yak < ya1)) lb = k;
        cy[k] = cyk;
        cg0 = cg1;
        cg1 = cyk;
    }
    lb0 = 0;
    if ((lb > 4) && (cimag(z) != 0.0)) {
        while (lb != lb0) {
            ch2 = 1;
            ch1 = 0;
            lb0 = lb;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*k*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp12 = ch0;
            cp22 = ch2;
            ch2 = 0;
            ch1 = 1;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*k*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp11 = ch0;
            cp21 = ch2;
            if (lb == *nm)
                cj[lb+1] = 2.0*lb*cj[lb]/z-cj[lb-1];
            if (cabs(cj[0]) > cabs(cj[1])) {
                cy[lb+1] = (cj[lb+1]*cby0-2.0*cp11/(M_PI*z))/cj[0];
                cy[lb] = (cj[lb]*cby0+2.0*cp12/(M_PI*z))/cj[0];
            }
            else {
                cy[lb+1] = (cj[lb+1]*cby1-2.0*cp21/(M_PI*z))/cj[1];
                cy[lb] = (cj[lb]*cby1+2.0*cp22/(M_PI*z))/cj[1];
            }
            cyl2 = cy[lb+1];
            cyl1 = cy[lb];
            for (k=lb-1;k>=0;k--) {
                cylk = 2.0*(k+1.0)*cyl1/z-cyl2;
                cy[k] = cylk;
                cyl2 = cyl1;
                cyl1 = cylk;
            }
            cyl1 = cy[lb];
            cyl2 = cy[lb+1];
            for (k=lb+1;k<n;k++) {
                cylk = 2.0*k*cyl2/z-cyl1;
                cy[k+1] = cylk;
                cyl1 = cyl2;
                cyl2 = cylk;
            }
            for (k=2;k<=*nm;k++) {
                wa = cabs(cy[k]);
                if (wa < cabs(cy[k-1])) lb = k;
            }
        }
    }
    for (k=2;k<=*nm;k++) {
        cyp[k] = cy[k-1]-(double)k*cy[k]/z;
    }
    return 0;
}

int cbessh(int n, double complex z, int *nm, double complex *ch, double complex *chp)
{
    int err, k, nm1;
    int n1 = n+1;
    int n2 = n+2;
    double complex cj[n2], cy[n2], cjp[n2], cyp[n2], ch1;

    err = cbessjyna(n1, z, &nm1, &cj[0], &cy[0], &cjp[0], &cyp[0]);
    *nm = nm1 - 1;
    for(k=0; k<=*nm; k++){
        ch[k] = cj[k] + I1 * cy[k];
    }
    ch1 = cj[nm1] + I1 * cy[nm1];

    for(k=1; k<=*nm; k++){
        chp[k-1] = (double)(k-1)*ch[k-1]/z - ch[k];
    }
    chp[*nm] = (double)(*nm)*ch[*nm]/z - ch1;
    return err;
}

/* unifpdf (from Matlab statistics toolbox)
  UNIFPDF Uniform (continuous) probability density function (pdf).
   Y = UNIFPDF(X,A,B) returns the continuous uniform pdf on the
   interval [A,B] at the values in X. By default A = 0 and B = 1.
   The size of Y is the common size of the input arguments. A scalar input
   functions as a constant matrix of the same size as the other inputs.

   Reference:
      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
      Functions", Government Printing Office, 1964, 26.1.34.
   Copyright 1993-2000 The MathWorks, Inc.
   $Revision: 2.9 $  $Date: 2000/05/26 18:53:54 $
*/
double unifpdf(double x, double a, double b)
{
    double y;
    if(x>=a && x<= b && a < b)
        y = 1 / (b - a);
    else
        y = NAN;
    return y;    
}

double trapzxy(double x1, double x2, double y1, double y2)
{
    return 0.5 * (y2 + y1) * (x2 - x1);
}

#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, long n)
{
    double x,tnm,sum,del;
    static double s;
    long it,j;

    if (n == 1) {
            return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } else {
            for (it=1,j=1;j<n-1;j++) it <<= 1;
            tnm=it;
            del=(b-a)/tnm;              /*This is the spacing of points to be
                                          added. */
            x=a+0.5*del;
            for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
            s=0.5*(s+(b-a)*sum/tnm);  /*This replaces s by its refined value.*/
            return s;
    }
}
#undef FUNC

double qtrap(double (*func)(double), double a, double b, double eps, double imax)
{
    double trapzd(double (*func)(double), double a, double b, long n);
    long i;
    double s,olds;

    olds = -1.0e30;
    for (i=1;i<=imax;i++) {
        s=trapzd(func,a,b,i);
        if (i > 5)
            if (fabs(s-olds) < eps*fabs(olds) || (s == 0.0 && olds == 0.0))
                return s;
        olds=s;
    }
    printf(">> Too many steps in routine qtrap\n");
    return 0.0;
}

/* Returns the integral of the function func from a to b. The parameters eps can be set to the
desired fractional accuracy and imax so that 2 to the power imax-1 is the maximum allowed
number of steps. Integration is performed by Simpson's rule. W.H.Press et al, Numerical Recipes in C */
double qsimp(double (*func)(double), double a, double b, double eps, double imax)
{
    double trapzd(double (*func)(double), double a, double b, long n);
    long i;
    double s,st,ost=0.0,os=0.0;
    for (i=1; i<=imax;i++) {
        st = trapzd(func,a,b,i);
        s = (4.0*st - ost)/3.0;
        if (i > 5) // Avoid spurious early convergence.
            if (fabs(s-os) < eps*fabs(os) || (s == 0.0 && os == 0.0))
                return s;
        os = s;
        ost = st;
    }
    printf(">> Too many steps in routine qsimp\n");
    return 0.0;
}

double roundzero(double x)
{
	return (fabs(x) < EPS ? 0 : x);
}

double newatan2(double y, double x)
{
    if(fabs(x) == 0 && fabs(y) == 0){
        if(signbit(x) != POS_ZERO && signbit(y) == POS_ZERO)
            return (atan2(y,x) - Pi);
        else if(signbit(x) != POS_ZERO && signbit(y) != POS_ZERO)
            return (atan2(y,x) + Pi);
        else
            return atan2(y,x);
    }
    else
        return atan2(y,x);
}

double atan2TwoPi(double y, double x)
{
    double res; 
    res = atan2(y,x);
    if(res<0)
        res += TwoPi;
    return res;
}

/* Return evanescent SQRT for waves problems */
double complex sqrte(double complex z)
{
	double complex y;

	if(cimag(z)==0 && creal(z)<0)
		y = -I1 * csqrt(abs(z));
	else
		y = csqrt(z);

	return y;
}

/* Conversion from Spherical to Cartesian coordinate system */
void sph2car(double theta, double phi, double AR, double AT, double AP, double axis[3])
{
	double st, ct, sp, cp;

	st = sin(theta);
	ct = cos(theta);
	sp = sin(phi);
	cp = cos(phi);

	axis[0] = AR*st*cp + AT*ct*cp - AP*sp;
	axis[1] = AR*st*sp + AT*ct*sp + AP*cp;
	axis[2] = AR*ct - AT*st;
}

/* Complex Vector Dot Product */
double complex VoV_complex(double complex A[3], double complex B[3])
{
    return(conj(A[0])*conj(B[0])+conj(A[1])*conj(B[1])+conj(A[2])*conj(B[2]));
//      return(A[0]*conj(B[0])+A[1]*conj(B[1])+A[2]*conj(B[2]));
}

/*  4x4 Matrix times 4x1 Vector */
void M4xV4_complex(double complex M[4][4], double complex V[4], double complex W[4])
{
      W[0] = V[0]*M[0][0] + V[1]*M[0][1] + V[2]*M[0][2] + V[3]*M[0][3];
      W[1] = V[0]*M[1][0] + V[1]*M[1][1] + V[2]*M[1][2] + V[3]*M[1][3];
      W[2] = V[0]*M[2][0] + V[1]*M[2][1] + V[2]*M[2][2] + V[3]*M[2][3];
      W[3] = V[0]*M[3][0] + V[1]*M[3][1] + V[2]*M[3][2] + V[3]*M[3][3];
}

/*  2x2 Matrix times 2x1 Vector */
void M2xV2_complex(double complex M[2][2], double complex V[2], double complex W[2])
{
      W[0] = V[0]*M[0][0] + V[1]*M[0][1];
      W[1] = V[0]*M[1][0] + V[1]*M[1][1];
}

/*  Transpose of 3x3 Matrix times 3x1 Complex Vector */
void MTxV_complex(double M[3][3], double V[3], double complex W[3])
{
      W[0] = M[0][0]*V[0] + M[1][0]*V[1] + M[2][0]*V[2];
      W[1] = M[0][1]*V[0] + M[1][1]*V[1] + M[2][1]*V[2];
      W[2] = M[0][2]*V[0] + M[1][2]*V[1] + M[2][2]*V[2];
}

/*  2x2 Complex Matrix Product */
void M2xM2_complex(double complex A[2][2], double complex B[2][2], double complex C[2][2])
{
      C[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0];
      C[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1];
      C[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0];
      C[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1];
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
void RotX(double th, double R[3][3])
{
	R[0][0] = 1;
	R[0][1] = 0;
	R[0][2] = 0;
	R[1][0] = 0;
	R[1][1] = cos(th);
	R[1][2] = -sin(th);
	R[2][0] = 0;
	R[2][1] = sin(th);
	R[2][2] = cos(th);
}

void RotY(double th, double R[3][3])
{
	R[0][0] = cos(th);
	R[0][1] = 0;
	R[0][2] = sin(th);
	R[1][0] = 0;
	R[1][1] = 1;
	R[1][2] = 0;
	R[2][0] = -sin(th);
	R[2][1] = 0;
	R[2][2] = cos(th);
}

void RotZ(double th, double R[3][3])
{
	R[0][0] = cos(th);
	R[0][1] = -sin(th);
	R[0][2] = 0;
	R[1][0] = sin(th);
	R[1][1] = cos(th);
	R[1][2] = 0;
	R[2][0] = 0;
	R[2][1] = 0;
	R[2][2] = 1;
}