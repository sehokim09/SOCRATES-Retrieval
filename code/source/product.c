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

void calcSpecularTermForward(void)
{
	double complex g_t0[2][2] = {{1, 0}, {0, 1}}, G_t0[4][4]; // Ideal transmitter antenna pattern
	double complex G_r[4][4]; // receiver antenna pattern
	double complex U_ts[4][4], U_sr[4][4];
	double complex e_t10[2] = {1, 0}, e_t20[2] = {0, 1}; 
	double complex E_t10[4] = {1, 0, 0, 0}, E_t20[4] = {0, 0, 0, 1}; // Transmitter Pol State
	double complex E_t1[4], E_t2[4];
	double complex Q[4][4] = {{1, 0, 0, 0}, {0, 0, 0, 1}, {0, 1, 1, 0}, {0, -I1, -I1, 0}};
	double complex r_sb[2][2] = {{1, 0}, {0, 1}}, R_sb[4][4];
	double complex temp1_V2[2], temp2_V2[2], temp1_V4[4], temp2_V4[4];

	M4xV4_complex(Q, E_t10, E_t1);
	M4xV4_complex(Q, E_t20, E_t2);
	
	calcMueller(Bistatic->g_rs, G_r);
	calcMueller(g_t0, G_t0);
	calcMueller(Bistatic->u_ts, U_ts);
	calcMueller(Bistatic->u_sr, U_sr);

	// Reflection Coefficient
	calcRcMulti();
	r_sb[0][0] = Gnd->RcV;
	r_sb[1][1] = Gnd->RcH;
	calcMueller(r_sb, R_sb);
	
	// Field - Bare soil
	// b_coh1b = g_r * u_sr * r_sb * u_ts * g_t * e_t1 ;
	// b_coh2b = g_r * u_sr * r_sb * u_ts * g_t * e_t2 ;

	M2xV2_complex(g_t0, e_t10, temp1_V2);
	M2xV2_complex(Bistatic->u_ts, temp1_V2, temp2_V2);
	M2xV2_complex(r_sb, temp2_V2, temp1_V2);
	M2xV2_complex(Bistatic->u_sr, temp1_V2, temp2_V2);
	M2xV2_complex(Bistatic->g_rs, temp2_V2, Out->b_coh1b);

	M2xV2_complex(g_t0, e_t20, temp1_V2);
	M2xV2_complex(Bistatic->u_ts, temp1_V2, temp2_V2);
	M2xV2_complex(r_sb, temp2_V2, temp1_V2);
	M2xV2_complex(Bistatic->u_sr, temp1_V2, temp2_V2);
	M2xV2_complex(Bistatic->g_rs, temp2_V2, Out->b_coh2b);

	// Power - Bare soil
	// P_coh1b = G_r * U_sr * R_sb * U_ts * G_t * E_t1 ;
	// P_coh2b = G_r * U_sr * R_sb * U_ts * G_t * E_t2 ;
	
	// P_coh1b = G_r0 * U_sr * R_sb * U_ts * G_t0 * E_t1
	M4xV4_complex(G_t0, E_t1, temp1_V4);
	M4xV4_complex(U_ts, temp1_V4, temp2_V4);
	M4xV4_complex(R_sb, temp2_V4, temp1_V4);
	M4xV4_complex(U_sr, temp1_V4, temp2_V4);
	M4xV4_complex(G_r, temp2_V4, temp1_V4);
	Out->P_coh1b[0] = creal(temp1_V4[0]);
	Out->P_coh1b[1] = creal(temp1_V4[1]);
	Out->P_coh1b[2] = creal(temp1_V4[2]);
	Out->P_coh1b[3] = creal(temp1_V4[3]);
	
	// P_coh2b = G_r0 * U_sr * R_sb * U_ts * G_t0 * E_t2
	M4xV4_complex(G_t0, E_t2, temp1_V4);
	M4xV4_complex(U_ts, temp1_V4, temp2_V4);
	M4xV4_complex(R_sb, temp2_V4, temp1_V4);
	M4xV4_complex(U_sr, temp1_V4, temp2_V4);
	M4xV4_complex(G_r, temp2_V4, temp1_V4);
	Out->P_coh2b[0] = creal(temp1_V4[0]);
	Out->P_coh2b[1] = creal(temp1_V4[1]);
	Out->P_coh2b[2] = creal(temp1_V4[2]);
	Out->P_coh2b[3] = creal(temp1_V4[3]);

	Out->P_coh1[0] = Out->P_coh1b[0];
	Out->P_coh1[1] = Out->P_coh1b[1];
	Out->P_coh1[2] = Out->P_coh1b[2];
	Out->P_coh1[3] = Out->P_coh1b[3];

	Out->P_coh2[0] = Out->P_coh2b[0];
	Out->P_coh2[1] = Out->P_coh2b[1];
	Out->P_coh2[2] = Out->P_coh2b[2];
	Out->P_coh2[3] = Out->P_coh2b[3];

	// Vegetation, if any
	if(SttData->ExistVeg){
		double complex t_sv[2][2], temp1_M22[2][2], r_sv[2][2], R_sv[4][4];
		
		// Optical depth
		Out->tauH = 2 * cos(Bistatic->th) * cimag(Bistatic->ArgH);
		Out->tauV = 2 * cos(Bistatic->th) * cimag(Bistatic->ArgV);

		t_sv[0][0] = cexp(I1*Bistatic->ArgV);
		t_sv[0][1] = 0;
		t_sv[1][0] = 0;
		t_sv[1][1] = cexp(I1*Bistatic->ArgH);
		M2xM2_complex(t_sv, r_sb, temp1_M22);
		M2xM2_complex(temp1_M22, t_sv, r_sv);
		calcMueller(r_sv, R_sv);

		// Field - Vegetation
		// b_coh1v = g_r * u_sr * r_sv * u_ts * g_t * e_t1 ;
		// b_coh2v = g_r * u_sr * r_sv * u_ts * g_t * e_t2 ;

		M2xV2_complex(g_t0, e_t10, temp1_V2);
		M2xV2_complex(Bistatic->u_ts, temp1_V2, temp2_V2);
		M2xV2_complex(r_sv, temp2_V2, temp1_V2);
		M2xV2_complex(Bistatic->u_sr, temp1_V2, temp2_V2);
		M2xV2_complex(Bistatic->g_rs, temp2_V2, Out->b_coh1v);

		M2xV2_complex(g_t0, e_t20, temp1_V2);
		M2xV2_complex(Bistatic->u_ts, temp1_V2, temp2_V2);
		M2xV2_complex(r_sv, temp2_V2, temp1_V2);
		M2xV2_complex(Bistatic->u_sr, temp1_V2, temp2_V2);
		M2xV2_complex(Bistatic->g_rs, temp2_V2, Out->b_coh2v);

		// Power - Vegetation
		// P_coh1v = G_r * U_sr * R_sv * U_ts * G_t * E_t1 ;
		// P_coh2v = G_r * U_sr * R_sv * U_ts * G_t * E_t2 ;
		
		M4xV4_complex(G_t0, E_t1, temp1_V4);
		M4xV4_complex(U_ts, temp1_V4, temp2_V4);
		M4xV4_complex(R_sv, temp2_V4, temp1_V4);
		M4xV4_complex(U_sr, temp1_V4, temp2_V4);
		M4xV4_complex(G_r, temp2_V4, temp1_V4);

		Out->P_coh1v[0] = creal(temp1_V4[0]);
		Out->P_coh1v[1] = creal(temp1_V4[1]);
		Out->P_coh1v[2] = creal(temp1_V4[2]);
		Out->P_coh1v[3] = creal(temp1_V4[3]);

		M4xV4_complex(G_t0, E_t2, temp1_V4);
		M4xV4_complex(U_ts, temp1_V4, temp2_V4);
		M4xV4_complex(R_sv, temp2_V4, temp1_V4);
		M4xV4_complex(U_sr, temp1_V4, temp2_V4);
		M4xV4_complex(G_r, temp2_V4, temp1_V4);

		Out->P_coh2v[0] = creal(temp1_V4[0]);
		Out->P_coh2v[1] = creal(temp1_V4[1]);
		Out->P_coh2v[2] = creal(temp1_V4[2]);
		Out->P_coh2v[3] = creal(temp1_V4[3]);

		Out->P_coh1[0] = Out->P_coh1v[0];
		Out->P_coh1[1] = Out->P_coh1v[1];
		Out->P_coh1[2] = Out->P_coh1v[2];
		Out->P_coh1[3] = Out->P_coh1v[3];

		Out->P_coh2[0] = Out->P_coh2v[0];
		Out->P_coh2[1] = Out->P_coh2v[1];
		Out->P_coh2[2] = Out->P_coh2v[2];
		Out->P_coh2[3] = Out->P_coh2v[3];
	}
}

void calcSpecularTermSimple(struct FixedObsType *F, double VWC)
{
	double complex g_t0[2][2] = {{1, 0}, {0, 1}}, G_t0[4][4]; // Ideal transmitter antenna pattern
	double complex G_r[4][4]; // receiver antenna pattern
	double complex U_ts[4][4], U_sr[4][4];
	double complex E_t10[4] = {1, 0, 0, 0}; // Transmitter Pol State
	double complex E_t1[4];
	double complex Q[4][4] = {{1, 0, 0, 0}, {0, 0, 0, 1}, {0, 1, 1, 0}, {0, -I1, -I1, 0}};
	double complex r_sb[2][2] = {{1, 0}, {0, 1}}, R_sb[4][4];
	double complex temp1_V4[4], temp2_V4[4];

	M4xV4_complex(Q, E_t10, E_t1);
	
	calcMueller(Bistatic->g_rs, G_r);
	calcMueller(g_t0, G_t0);
	calcMueller(Bistatic->u_ts, U_ts);
	calcMueller(Bistatic->u_sr, U_sr);

	// Reflection Coefficient
	calcRcMulti();
	r_sb[0][0] = Gnd->RcV;
	r_sb[1][1] = Gnd->RcH;

	calcMueller(r_sb, R_sb);

	// Power - Bare soil
	// P_coh1b = G_r * U_sr * R_sb * U_ts * G_t * E_t1 ;
	// P_coh2b = G_r * U_sr * R_sb * U_ts * G_t * E_t2 ;
	// P_coh1b = G_r0 * U_sr * R_sb * U_ts * G_t0 * E_t1
	M4xV4_complex(G_t0, E_t1, temp1_V4);
	M4xV4_complex(U_ts, temp1_V4, temp2_V4);
	M4xV4_complex(R_sb, temp2_V4, temp1_V4);
	M4xV4_complex(U_sr, temp1_V4, temp2_V4);
	M4xV4_complex(G_r, temp2_V4, temp1_V4);
	Out->P_coh1[0] = creal(temp1_V4[0]);
	Out->P_coh1[1] = creal(temp1_V4[1]);

	// Vegetation, if any
	if(SttData->ExistVeg){
		double VOD_Co, VOD_X;
								
		if(F->polRx == 'R'){
			VOD_Co = b_LUT[F->ID][0] * VWC; // RR-pol
			VOD_X = b_LUT[F->ID][1] * VWC; // RL-pol
			Out->P_coh1[0] *= exp(-2*VOD_Co/cos(Bistatic->th));
			Out->P_coh1[1] *= exp(-2*VOD_X/cos(Bistatic->th));
		}
		else if(F->polRx == 'X'){
			VOD_Co = b_LUT[F->ID][2] * VWC; // RX-pol
			VOD_X = b_LUT[F->ID][3] * VWC; // RY-pol
			Out->P_coh1[0] *= exp(-2*VOD_Co/cos(Bistatic->th));
			Out->P_coh1[1] *= exp(-2*VOD_X/cos(Bistatic->th));
		}
		else{
			printf(">> Wrong Rx polarization.\n");
			exit(1);
		}
	}
}