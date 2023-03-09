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

void updateBistaticFixed(struct FixedObsType *F)
{
	Bistatic->Freq = F->Freq; // Transmitter Frequency [Hz]
	Bistatic->polTx = F->polTx; // Transmit Antenna Polarization
	Bistatic->Wavelength = F->Wavelength; // [m]
	Bistatic->Wavenumber = TwoPi / Bistatic->Wavelength;
	Bistatic->th = F->thTx;
	Bistatic->EIRP_Tx = F->EIRP_Tx; // Transmitter EIRP [dB]
	Bistatic->EIRP_Tx_Lin = pow(10, F->EIRP_Tx/10);
	Bistatic->G_Rx = F->G_Rx; // Receive Antenna Gain [dB]
	Bistatic->G_Rx_Lin= pow(10, F->G_Rx/10);
	Bistatic->t_coh = F->t_coh;
	Bistatic->Bandwidth = F->Bandwidth;
	Bistatic->Tn = F->Tn;
	Bistatic->N_inc = F->N_inc;
}

/* Calculates the Fresnel zones for the specular reflection */
void calcFresnelZones(struct FixedObsType *F)
{
	double dist, hd, hs, rd2, rs, del0, sinTh, cosTh, secTh, del, a, b, c, nump, denp, p;
	long nn;

	// Transmitter/Receiver ground range
	dist = F->hTx * tan(F->thTx);
	// Distance of specular point away from the receiver ground projection
	F->S0x = dist / (1 + F->hTx/F->hRx);

	// Parameters for Fresnel zone calculation
	hd = F->hTx - F->hRx;
	hs = F->hTx + F->hRx;
	rd2 = sqrt(hd*hd + dist*dist);
	rs = sqrt(hs*hs + dist*dist);
	del0 = rs - rd2; // the shortest distance after the direct path

	sinTh = hd / rd2;
	cosTh = dist / rd2;
	secTh = 1/ cosTh;

	// Iterate for the decided number of Fresnel zones
	// This iteration may be meaningful, when incoherent contribution is considered. Nfz is always 1 for the current version.
	for(nn = 0; nn < F->Nfz; nn++){
		// the path for the n-th Fresnel zone
		del = del0 + nn * Bistatic->Wavelength / 2;

		a = (dist * secTh + del) / 2;
		b = sqrt(del*del + 2*dist*del*secTh) / 2;
		c = (F->hTx + F->hRx) / 2;

		nump = c * (a*a - b*b) * sinTh * cosTh;
		denp = (b*b*cosTh*cosTh + a*a*sinTh*sinTh);
		p = -nump / denp;
		// distance to the center of the ellipse
		F->x1 = dist/2 + p;
		// semi-minor axis
		F->by1 = b * sqrt(1 - c*c/denp);
		// semi-major axis
		F->ax1 = F->by1 * a / sqrt(denp);
	}
}

/* Calculates the transmitter geometry-related transformation matrices */
void calcTxGeometry(struct FixedObsType *F)
{
	double ux[3] = {1, 0, 0}, uy[3] = {0, 1, 0}, uz[3] = {0, 0, 1}; // Ground Reference Coordinate System (ENU)
	double phInTx, AntRotZ_Tx[3][3], AntRotY_Tx[3][3], AntRot_Tx[3][3];

	/* .. Incoming Signal .. */
	// Azimuth angle (in standard spherical coords) of incoming signal
	phInTx = HalfPi - F->phTx + Pi;
	
	/* .. Tx Antenna Rotation Matrices .. */
	RotZ(phInTx, AntRotZ_Tx);
	RotY(Pi - F->thTx, AntRotY_Tx);
	MxM(AntRotZ_Tx, AntRotY_Tx, AntRot_Tx);

	/* .. Tx Antenna Coordinate Systems .. */
	double uxt[3], uyt[3], uzt[3];
	MxV(AntRot_Tx, ux, uxt);
	MxV(AntRot_Tx, uy, uyt);
	MxV(AntRot_Tx, uz, uzt);

	/* .. Transformation .. */
	// Transformation matrix for transforming a vector from the ground frame to transmitter system
	Bistatic->Tgt[0][0] = uxt[0]; Bistatic->Tgt[0][1] = uxt[1]; Bistatic->Tgt[0][2] = uxt[2];
	Bistatic->Tgt[1][0] = uyt[0]; Bistatic->Tgt[1][1] = uyt[1]; Bistatic->Tgt[1][2] = uyt[2];
	Bistatic->Tgt[2][0] = uzt[0]; Bistatic->Tgt[2][1] = uzt[1]; Bistatic->Tgt[2][2] = uzt[2];
}

/* Calculates the bistatic geometry-related transformation matrices, direction vectors, and rotation matrices for fixed observation */
void calcGeometryFixed(struct FixedObsType *F)
{
	double ux[3] = {1, 0, 0}, uy[3] = {0, 1, 0}, uz[3] = {0, 0, 1}; // Ground Reference Coordinate System (ENU)
	//double posGnd[3] = {0, 0, 0}; // Center of reference coordinate system
	double AntRotZ_Rx[3][3], AntRotY_Rx[3][3], AntRot_Rx[3][3], AntRotZ_Tx[3][3];
	double rd, posTx[3], posSp[3], posRx[3] = {0, 0, 0};//, ht;
	double tempV3[3] = {0, 0, 0};
	double phRx, phTx;
	double osp_rf[3], idn_rf[3], isn_sf[3];
	double th0, ph0;

	phRx = HalfPi - F->phRx;
	phTx = HalfPi - F->phTx;

	// Antenna Rotation Matrix for Receiver
	RotZ(phRx, AntRotZ_Rx);
	RotY(Pi - F->thRx, AntRotY_Rx);
	MxM(AntRotZ_Rx, AntRotY_Rx, AntRot_Rx);
	// Antenna Rotation Matrix for Transmitter
	RotZ(phTx, AntRotZ_Tx);

	/* .. Transmitter Position and Related Locations .. */
	// Slant range (Approximated)
	rd = sqrt(F->rTx*F->rTx - EarthRad*EarthRad*cos(F->elTx)*cos(F->elTx)) - EarthRad*sin(F->elTx);
	// Transmitter Antenna position
	posTx[0] = rd * sin(F->thTx) * cos(phTx);
	posTx[1] = rd * sin(F->thTx) * sin(phTx);
	posTx[2] = rd * cos(F->thTx);
	// Transmitter altitude
	F->hTx = posTx[2];
	// Specular reflection point and 1st Fresnel zone ellipse
	F->Nfz = 1; // number of fresnel zone is 1 since only the specular contribution is considered
	calcFresnelZones(F);
	// Specular reflection point position
	tempV3[0] = F->S0x;
	MxV(AntRotZ_Tx, tempV3, posSp);
	
	/* .. Receiver Position and Pointing Locations .. */
	// Receiver Antenna position
	posRx[2] = F->hRx;

	/* .. Propagation Vectors .. */
	double RT[3], ST[3], RS[3];
	long i;
	for(i=0; i<3; i++){
		// Transmitter to Receiver
		RT[i] = posRx[i] - posTx[i];
		// Tansmitter to Specular point
		ST[i] = posSp[i] - posTx[i];
		// Specular point to Receiver
		RS[i] = posRx[i] - posSp[i];
	}
	Bistatic->r_tr = CopyUnitV(RT, Bistatic->idn);
	Bistatic->r_ts = CopyUnitV(ST, Bistatic->isn);
	Bistatic->r_sr = CopyUnitV(RS, Bistatic->osp);
	Bistatic->r_tsr = Bistatic->r_ts + Bistatic->r_sr;

	/* .. Coordinate Systems .. */
	double uxr[3], uyr[3], uzr[3], uxs[3], uys[3], uzs[3];
	// Receiver Antenna Coordinate System
	MxV(AntRot_Rx, ux, uxr);
	MxV(AntRot_Rx, uy, uyr);
	MxV(AntRot_Rx, uz, uzr);
	// Specular Point Coordinate System
	MxV(AntRotZ_Tx, ux, uxs);
	MxV(AntRotZ_Tx, uy, uys);
	MxV(AntRotZ_Tx, uz, uzs);

	/* .. Transformations .. */
	// Transformation matrix for transforming a vector from the ground frame to local (specular) ground system
	Bistatic->Tgs[0][0] = uxs[0]; Bistatic->Tgs[0][1] = uxs[1]; Bistatic->Tgs[0][2] = uxs[2];
	Bistatic->Tgs[1][0] = uys[0]; Bistatic->Tgs[1][1] = uys[1]; Bistatic->Tgs[1][2] = uys[2];
	Bistatic->Tgs[2][0] = uzs[0]; Bistatic->Tgs[2][1] = uzs[1]; Bistatic->Tgs[2][2] = uzs[2];
	
	// Transformation matrix for transforming a vector from the ground frame to receiver system
	Bistatic->Tgr[0][0] = uxr[0]; Bistatic->Tgr[0][1] = uxr[1]; Bistatic->Tgr[0][2] = uxr[2];
	Bistatic->Tgr[1][0] = uyr[0]; Bistatic->Tgr[1][1] = uyr[1]; Bistatic->Tgr[1][2] = uyr[2];
	Bistatic->Tgr[2][0] = uzr[0]; Bistatic->Tgr[2][1] = uzr[1]; Bistatic->Tgr[2][2] = uzr[2];

	/* .. Angles .. */
	// The incidence angle on the receiver in receiver coordinates
	// propagation vector from Transmitter to reciever in receiver antenna system
	MxV(Bistatic->Tgr, Bistatic->idn, idn_rf);

	// Transmitter to Receiver
	// off-axis angle of zr towards Transmitter
	th0 = acos(-idn_rf[2]);
	// orientation - azimuth
	ph0 = atan2TwoPi(-idn_rf[1], -idn_rf[0]);

	F->AngT2R_rf[0] = th0;
	F->AngT2R_rf[1] = ph0;
	//printf("T->R (ph,th) = %lf %lf\n", ph0*R2D, th0*R2D);
	// Propagation vector from specular point to receiver in receiver antenna system
	MxV(Bistatic->Tgr, Bistatic->osp, osp_rf);
	// Specular Point to Receiver
	// off-axis angle of zr towards specular point
	th0 = acos(-osp_rf[2]);
	// orientation - azimuth
	ph0 = atan2TwoPi(-osp_rf[1], -osp_rf[0]);

	F->AngS2R_rf[0] = th0;
	F->AngS2R_rf[1] = ph0;
	//printf("S->R (ph,th) = %lf %lf\n", ph0*R2D, th0*R2D);

	// The incidence angle in specular frame
	// propagation vector from Transmitter to ground in local (specular) ground system
	MxV(Bistatic->Tgs, Bistatic->isn, isn_sf);
	// Transmitter to Specular Point
	// off-axis angle of zs towards Transmitter
	th0 = acos(-isn_sf[2]);
	// orientation - azimuth
	ph0 = atan2TwoPi(-isn_sf[1], -isn_sf[0]);

	F->AngT2S_sf[0] = th0;
	F->AngT2S_sf[1] = ph0;

	Bistatic->AngT2S_sf_th0 = F->AngT2S_sf[0];
}