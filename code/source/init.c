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

void Welcome(void)
{
	int retval;

	printf("/************************************************************/\n");
	printf("/*                                                          */\n");
	printf("/*                 |  SOCRATES-retrieval  |                 */\n");
	printf("/*                                                          */\n");
	printf("/************************************************************/\n");
	// Make an output directories if not exist
	retval = mkdir("./results/", 0777);
	if(retval != 0 && errno != EEXIST){
		printf(">> mkdir fail: %s\n", "./results/");
		exit(ERRCODE);            
	}
	if(SimMode == FORWARD){
		printf(">> Simulation Mode: Forward\n");
		// Make an output directories if not exist
		retval = mkdir(ForwardPath, 0777);
		if(retval != 0 && errno != EEXIST){
			printf(">> mkdir fail: %s\n", ForwardPath);
			exit(ERRCODE);            
		}
		else{
			printf(">> Output dir: %s\n", ForwardPath);
		}
	}
	else if(SimMode == INVERSE){
		printf(">> Simulation Mode: Inverse\n");
		// Make an output directories if not exist
		retval = mkdir(InversePath, 0777);
		if(retval != 0 && errno != EEXIST){
			printf(">> mkdir fail: %s\n", InversePath);
			exit(ERRCODE);            
		}
		else{
			printf(">> Output dir: %s\n", InversePath);
		}
	}
	printf(">> Number of Input Soil Layers is %ld\n", NLayer);
	printf(">> Ground Multi-layered Structure: POME model, # sublayers = %ld + %ld\n", MultiLayer->NSublayer_POME, MultiLayer->NSublayer_extend);
	printf(">> Initiating Fixed Geometry Observation ...\n");
}

void InitConstants(void)
{
	// Physical constants
	EarthRad = 6.378145E6;
	EarthMu = 3.986004E14;
	EarthW = 7.292115E-5;

	// Math
	Pi = 4.0*atan(1.0);
	TwoPi = 2.0*Pi;
	HalfPi = 0.5*Pi;
	R2D = 180.0/Pi;
	D2R = Pi/180.0;
}

void freeInputSoOp(void)
{
	free(InputSoOp->F);	
	free(InputSoOp);
}

void InitVegetation(void)
{
	long iVeg, iKind, jKind, iLayer, iType;
	long sTypknd[NUM_VEGTYPE], num_types, num_kinds;
	long part_type, part_type_idx[NUM_VEGTYPE] = {0}, part_kind;
	long iSoOp, i;
	
	for(iVeg=0; iVeg<NVeg; iVeg++){
		// Reset parameters
		for(iType=0; iType<NUM_VEGTYPE; iType++)
			sTypknd[iType] = 0;
		num_types = 0;
		num_kinds = 0;
		// Initialization
		VegData[iVeg].typknd = (int **) calloc(VegData[iVeg].nLayer, sizeof(int *));
		for(iLayer=0; iLayer<VegData[iVeg].nLayer; iLayer++){
			VegData[iVeg].typknd[iLayer] = (int *) calloc(NUM_VEGTYPE, sizeof(int));
			for(iKind=0; iKind<VegData[iVeg].nKind[iLayer]; iKind++){
				switch(VegData[iVeg].kind[iLayer][iKind][0]){
					case 'L':
						VegData[iVeg].typknd[iLayer][0]++;
						break;
					case 'B':
						VegData[iVeg].typknd[iLayer][1]++;
						break;
					case 'T':
						VegData[iVeg].typknd[iLayer][2]++;
						break;
					case 'N':
						VegData[iVeg].typknd[iLayer][3]++;
						break;
					case 'W':
						VegData[iVeg].typknd[iLayer][4]++;
						break;
					default:
						printf(">> Wrong kind of vegetation...\n"); exit(1);	
				}
			}
			for(iType=0; iType<NUM_VEGTYPE; iType++){
				sTypknd[iType] += VegData[iVeg].typknd[iLayer][iType];
				if(num_kinds < VegData[iVeg].typknd[iLayer][iType])
					num_kinds = VegData[iVeg].typknd[iLayer][iType];
			}
			VegData[iVeg].depth += VegData[iVeg].thickness[iLayer];
		}

		VegData[iVeg].part_type_idx = (long *) calloc(NUM_VEGTYPE, sizeof(long));
		for(iType=0; iType<NUM_VEGTYPE; iType++){
			if(sTypknd[iType] > 0){
				part_type_idx[iType] = num_types;
				VegData[iVeg].part_type_idx[iType] = num_types;
				num_types++;
			}
		}
		VegData[iVeg].nKindMax = num_kinds;
		VegData[iVeg].nTypeMax = num_types;
		VegData[iVeg].typknd_new = (int **) calloc(VegData[iVeg].nLayer, sizeof(int *));
		for(iLayer=0; iLayer<VegData[iVeg].nLayer; iLayer++)
			VegData[iVeg].typknd_new[iLayer] = (int *) calloc(VegData[iVeg].nTypeMax, sizeof(int));
		VegData[iVeg].dsty = (double ***) calloc(num_kinds, sizeof(double **));
		VegData[iVeg].dim1 = (double ***) calloc(num_kinds, sizeof(double **));
		VegData[iVeg].dim2 = (double ***) calloc(num_kinds, sizeof(double **));
		VegData[iVeg].dim3 = (double ***) calloc(num_kinds, sizeof(double **));
		VegData[iVeg].e_c = (double complex ****) calloc(num_kinds, sizeof(double complex ***));
		VegData[iVeg].beginAng = (double ***) calloc(num_kinds, sizeof(double **));
		VegData[iVeg].endAng = (double ***) calloc(num_kinds, sizeof(double **));
		VegData[iVeg].shape = (long ***) calloc(num_kinds, sizeof(long **));
		VegData[iVeg].VWC = (double ***) calloc(num_kinds, sizeof(double **));
		for(iKind=0; iKind<num_kinds; iKind++){
			VegData[iVeg].dsty[iKind] = (double **) calloc(num_types, sizeof(double *));
			VegData[iVeg].dim1[iKind] = (double **) calloc(num_types, sizeof(double *));
			VegData[iVeg].dim2[iKind] = (double **) calloc(num_types, sizeof(double *));
			VegData[iVeg].dim3[iKind] = (double **) calloc(num_types, sizeof(double *));
			VegData[iVeg].e_c[iKind] = (double complex ***) calloc(num_types, sizeof(double complex **));
			VegData[iVeg].beginAng[iKind] = (double **) calloc(num_types, sizeof(double *));
			VegData[iVeg].endAng[iKind] = (double **) calloc(num_types, sizeof(double *));
			VegData[iVeg].shape[iKind] = (long **) calloc(num_types, sizeof(long *));
			VegData[iVeg].VWC[iKind] = (double **) calloc(num_types, sizeof(double *));
			for(iType=0; iType<num_types; iType++){
				VegData[iVeg].dsty[iKind][iType] = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
				VegData[iVeg].dim1[iKind][iType] = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
				VegData[iVeg].dim2[iKind][iType] = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
				VegData[iVeg].dim3[iKind][iType] = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
				VegData[iVeg].e_c[iKind][iType] = (double complex **) calloc(VegData[iVeg].nLayer, sizeof(double complex *));
				VegData[iVeg].beginAng[iKind][iType] = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
				VegData[iVeg].endAng[iKind][iType] = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
				VegData[iVeg].shape[iKind][iType] = (long *) calloc(VegData[iVeg].nLayer, sizeof(long));
				VegData[iVeg].VWC[iKind][iType] = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
				for(iLayer=0; iLayer<VegData[iVeg].nLayer; iLayer++){
					VegData[iVeg].e_c[iKind][iType][iLayer] = (double complex *) calloc(NSoOp, sizeof(double complex));
				}
			}
		}

		// Find layer indices for each kind
		VegData[iVeg].typkndThick = (double *) calloc(NKind, sizeof(double));
		for(iKind=0; iKind<NKind; iKind++){
			for(iLayer=0; iLayer<VegData[iVeg].nLayer; iLayer++){
				for(jKind=0; jKind<VegData[iVeg].nKind[iLayer]; jKind++){
					if(!strcmp(VegKind[iKind].ID, VegData[iVeg].kind[iLayer][jKind])){
						VegData[iVeg].typkndThick[iKind] += VegData[iVeg].thickness[iLayer];
					}
				}
			}
		}

		for(iLayer=0; iLayer<VegData[iVeg].nLayer; iLayer++){
			part_kind = 0;
			for(jKind=0; jKind<NKind; jKind++){
				if(jKind>0 && (VegKind[jKind].ID[0] != VegKind[jKind-1].ID[0]))
					part_kind = 0;
				for(iKind=0; iKind<VegData[iVeg].nKind[iLayer]; iKind++){
					if(!strcmp(VegKind[jKind].ID, VegData[iVeg].kind[iLayer][iKind])){
						switch(VegData[iVeg].kind[iLayer][iKind][0]){
							case 'L':
								part_type = part_type_idx[0];
								VegData[iVeg].shape[part_kind][part_type][iLayer] = DISK;
								VegData[iVeg].typknd_new[iLayer][part_type] = VegData[iVeg].typknd[iLayer][0];
								break;
							case 'B':
								part_type = part_type_idx[1];
								VegData[iVeg].shape[part_kind][part_type][iLayer] = CYLINDER;
								VegData[iVeg].typknd_new[iLayer][part_type] = VegData[iVeg].typknd[iLayer][1];
								break;
							case 'T':
								part_type = part_type_idx[2];
								VegData[iVeg].shape[part_kind][part_type][iLayer] = CYLINDER;
								VegData[iVeg].typknd_new[iLayer][part_type] = VegData[iVeg].typknd[iLayer][2];
								break;
							case 'N':
								part_type = part_type_idx[3];
								VegData[iVeg].shape[part_kind][part_type][iLayer] = CYLINDER;
								VegData[iVeg].typknd_new[iLayer][part_type] = VegData[iVeg].typknd[iLayer][3];
								break;
							case 'W':
								part_type = part_type_idx[4];
								VegData[iVeg].shape[part_kind][part_type][iLayer] = DISK;
								VegData[iVeg].typknd_new[iLayer][part_type] = VegData[iVeg].typknd[iLayer][4];
								break;
						}						
						VegData[iVeg].dsty[part_kind][part_type][iLayer] = VegKind[jKind].density * VegData[iVeg].depth / VegData[iVeg].typkndThick[jKind];	
						VegData[iVeg].dim1[part_kind][part_type][iLayer] = VegKind[jKind].dim[0];	
						VegData[iVeg].dim2[part_kind][part_type][iLayer] = VegKind[jKind].dim[1];	
						VegData[iVeg].dim3[part_kind][part_type][iLayer] = VegKind[jKind].dim[2];
						VegData[iVeg].beginAng[part_kind][part_type][iLayer] = VegKind[jKind].beginAng;
						VegData[iVeg].endAng[part_kind][part_type][iLayer] = VegKind[jKind].endAng;
						VegData[iVeg].VWC[part_kind][part_type][iLayer] = VegKind[jKind].VWC * VegData[iVeg].typkndThick[jKind];
						for(iSoOp=0; iSoOp<NSoOp; iSoOp++)
							VegData[iVeg].e_c[part_kind][part_type][iLayer][iSoOp] = VegKind[jKind].e_c[iSoOp];
						part_kind++;
						break;
					}
				}	
			}
		}
/*
		VegData[iVeg].dKz = (double complex ****)calloc(NSoOp, sizeof(double complex ***));
		for(iSoOp=0; iSoOp<NSoOp; iSoOp++){
			VegData[iVeg].dKz[iSoOp] = (double complex ***)calloc(2, sizeof(double complex **));
			for(i=0; i<2; i++){
				VegData[iVeg].dKz[iSoOp][i] = (double complex **)calloc(90, sizeof(double complex *));
				for(iAng=0; iAng<90; iAng++)
					VegData[iVeg].dKz[iSoOp][i][iAng] = (double complex *)calloc(VegData[iVeg].nLayer, sizeof(double complex));
			}
		}
*/
		VegData[iVeg].dKzNc = (double complex **)calloc(2, sizeof(double complex *));
		for(i=0; i<2; i++){
			VegData[iVeg].dKzNc[i] = (double complex *)calloc(VegData[iVeg].nLayer, sizeof(double complex));
		}
		
		VegModel = (struct VegModelType *) calloc(1, sizeof(struct VegModelType));
		if (VegModel == NULL) {
			printf("VegModel calloc returned null pointer.  Bailing out!\n");
			exit(1);
		}
	}
}

void InitFixedObs(void)
{
	long ISoOp;

	Fixed = (struct FixedObsType *) calloc(NSoOp, sizeof(struct FixedObsType));
	if (Fixed == NULL) {
		printf("Fixed Observation calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	for(ISoOp=0; ISoOp<NSoOp; ISoOp++){
		Fixed[ISoOp].ID = ISoOp;
		Fixed[ISoOp].Freq = InputSoOp->F[ISoOp].Freq; // [Hz]
		Fixed[ISoOp].Wavelength = LIGHTSPEED / Fixed[ISoOp].Freq;
		Fixed[ISoOp].polRx = InputSoOp->F[ISoOp].polRx; // Receive Antenna Polarization
		Fixed[ISoOp].polTx = InputSoOp->F[ISoOp].polTx; // Transmit Antenna Polarization
		Fixed[ISoOp].hpbw = InputSoOp->F[ISoOp].hpbw; // Half-power beanwidth of the receive antenna pattern [deg]
		Fixed[ISoOp].SLL = InputSoOp->F[ISoOp].SLL; // First Sidelobe level of the receive antenna pattern [dB]
		Fixed[ISoOp].XPL = InputSoOp->F[ISoOp].XPL; // Cross-polarization level of the receive antenna pattern [dB]
		Fixed[ISoOp].AntPatRes = InputSoOp->F[ISoOp].AntPatRes; // Receive antenna pattern resolution [deg]
		Fixed[ISoOp].G_Rx = InputSoOp->F[ISoOp].G_Rx; // [dB]
		Fixed[ISoOp].rTx = InputSoOp->F[ISoOp].rTx; // Transmitter range from Earth's center [m]
		Fixed[ISoOp].hTx = InputSoOp->F[ISoOp].hTx; // Transmitter altitude [m]
		Fixed[ISoOp].phTx = 0; // Tx Azimuth angle [rad]
		Fixed[ISoOp].hRx = InputSoOp->F[ISoOp].hRx; // Receiver altitude [m]
		Fixed[ISoOp].phRx = 0; // Rx Azimuth angle of receiver position [rad]

		Fixed[ISoOp].AntTag = InputSoOp->F[ISoOp].AntTag;
		if(Fixed[ISoOp].AntTag == USER_DEFINED){
			strcpy(Fixed[ISoOp].AntFileName[0], InputSoOp->F[ISoOp].AntFileName[0]);
			strcpy(Fixed[ISoOp].AntFileName[1], InputSoOp->F[ISoOp].AntFileName[1]);
			strcpy(Fixed[ISoOp].AntFileName[2], InputSoOp->F[ISoOp].AntFileName[2]);
			strcpy(Fixed[ISoOp].AntFileName[3], InputSoOp->F[ISoOp].AntFileName[3]);
			LoadAntPattern_Fixed(&Fixed[ISoOp]);
		}
		else if(Fixed[ISoOp].AntTag == GAUSSIAN){
			LoadAntPattern_Fixed_GG(&Fixed[ISoOp]);
		}
		Fixed[ISoOp].h = ( 2*InputSoOp->RMSH*TwoPi / (Fixed[ISoOp].Wavelength *100) ) * ( 2*InputSoOp->RMSH*TwoPi / (Fixed[ISoOp].Wavelength *100) );
		strcpy(Fixed[ISoOp].orbitName, InputSoOp->F[ISoOp].orbitName);
	}

	// Random Process
	RNG = CreateRandomProcess((long)time(NULL));
}

void LoadVegetation(void)
{
	FILE *infile;
    char junk[120], newline, response[120], kind[2];
	long iVeg, iLayer, iKind, iSoOp;
	double freq;// temp1, temp2;

	/* .. Read from file Inp_Fixed.txt */
    infile=FileOpen(InPath,"Inp_Veg.txt","r");

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

	// Number of Vegetation Types
	fscanf(infile, "%ld %[^\n] %[\n]", &NVeg, junk, &newline);

	VegData = (struct VegDataType *) calloc(NVeg, sizeof(struct VegDataType));
	if (VegData == NULL) {
		printf("VegData calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}

	for(iVeg=0; iVeg<NVeg; iVeg++){
		fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
		// Vegetation Type (Label)
		fscanf(infile,"%s %[^\n] %[\n]",VegData[iVeg].label,junk,&newline);
		// Number of Layers
		fscanf(infile, "%ld %[^\n] %[\n]", &VegData[iVeg].nLayer, junk, &newline);
		VegData[iVeg].thickness = (double *) calloc(VegData[iVeg].nLayer, sizeof(double));
		if (VegData[iVeg].thickness == NULL) {
			printf("VegData[iVeg].thickness calloc returned null pointer.  Bailing out!\n");
			exit(1);
		}
		VegData[iVeg].nKind = (long *) calloc(VegData[iVeg].nLayer, sizeof(long));
		if (VegData[iVeg].nKind == NULL) {
			printf("VegData[iVeg].nKind calloc returned null pointer.  Bailing out!\n");
			exit(1);
		}
		VegData[iVeg].kind = (char ***) calloc(VegData[iVeg].nLayer, sizeof(char**));

		for(iLayer=0; iLayer<VegData[iVeg].nLayer; iLayer++){
			fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
			// Thickness [m]
			fscanf(infile, "%lf %[^\n] %[\n]", &VegData[iVeg].thickness[iLayer], junk, &newline);
			// Number of Included Kinds
			fscanf(infile, "%ld ", &VegData[iVeg].nKind[iLayer]);
			
			VegData[iVeg].kind[iLayer] = (char **) calloc(VegData[iVeg].nKind[iLayer], sizeof(char*));
			if (VegData[iVeg].kind[iLayer] == NULL) {
				printf("VegData[iVeg].kind calloc returned null pointer.  Bailing out!\n");
				exit(1);
			}
			// Particle ID
			for(iKind=0; iKind<VegData[iVeg].nKind[iLayer]; iKind++){
				VegData[iVeg].kind[iLayer][iKind] = (char *) calloc(2, sizeof(char));
			
				fscanf(infile,"%s ", kind);
				strcpy(VegData[iVeg].kind[iLayer][iKind], kind);
			}
			fscanf(infile,"%[^\n] %[\n]", junk, &newline);
		}
	}

	fscanf(infile,"%[^\n] %[\n]",junk,&newline);
	fscanf(infile,"%[^\n] %[\n]",junk,&newline);
	fscanf(infile,"%[^\n] %[\n]",junk,&newline);
	
	// Number of Kinds
	fscanf(infile, "%ld %[^\n] %[\n]", &NKind, junk, &newline);

	VegKind = (struct VegKindType *) calloc(NKind, sizeof(struct VegKindType));
	if (VegKind == NULL) {
		printf("VegKind calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}

	for(iKind=0; iKind<NKind; iKind++){
		fscanf(infile,"%[^\n] %[\n]",junk,&newline);
		// Particle ID
		fscanf(infile,"%s %[^\n] %[\n]", VegKind[iKind].ID, junk,&newline);
		// Density [particles/m^3]
		fscanf(infile, "%lf %[^\n] %[\n]", &VegKind[iKind].density, junk, &newline);
		// Dimensions [m, m, m]
		fscanf(infile, "%lf %lf %lf %[^\n] %[\n]", &VegKind[iKind].dim[0], &VegKind[iKind].dim[1], &VegKind[iKind].dim[2], junk, &newline);
		// Dielectric Constant (real, imag)
		//fscanf(infile, "%lf %lf %[^\n] %[\n]", &temp1, &temp2, junk, &newline);
		//VegKind[iKind].e_c = temp1 + I1*temp2;
		// Gravimetric Water Content
		fscanf(infile, "%lf %[^\n] %[\n]", &VegKind[iKind].mv, junk, &newline);
		// Interval Angle of Orientation [deg] (min, max)
		fscanf(infile, "%lf %lf %[^\n] %[\n]", &VegKind[iKind].beginAng, &VegKind[iKind].endAng, junk, &newline);
		// degee to radian
		VegKind[iKind].beginAng *= D2R;
		VegKind[iKind].endAng *= D2R;

		// Vegetation Water Content
		VegKind[iKind].VWC = VegKind[iKind].density * VegKind[iKind].mv * 1000 * Pi * VegKind[iKind].dim[0] * VegKind[iKind].dim[1] * VegKind[iKind].dim[2];

		// Dielectric Constant
		VegKind[iKind].e_c = (double complex *) calloc(NSoOp, sizeof(double complex));

		for(iSoOp=0; iSoOp<NSoOp; iSoOp++){
			freq = Fixed[iSoOp].Freq;
			switch(VegKind[iKind].ID[0]){
				case 'L':
					VegKind[iKind].e_c[iSoOp] = calcDielUlabyElRayesMv(freq, VegKind[iKind].mv);
					break;
				case 'B':
					VegKind[iKind].e_c[iSoOp] = calcDielUlabyElRayesMv(freq, VegKind[iKind].mv);
					break;
				case 'T':
					VegKind[iKind].e_c[iSoOp] = calcDielUlabyElRayesMv(freq, VegKind[iKind].mv);
					break;
				case 'N':
					VegKind[iKind].e_c[iSoOp] = calcDielUlabyElRayesMv(freq, VegKind[iKind].mv);
					break;
				case 'W':
					VegKind[iKind].e_c[iSoOp] = calcDielDebye(freq);
					break;
				default:
					printf(">> Wrong kind of vegetation...\n"); exit(1);	
			}
		}
	}
}

void LoadFixed(void)
{
	FILE *infile;
    char junk[120], newline, response[120];
	char AntFileName[4][30], orbitName[5];
	long ISoOp;
	double hRx, G_Rx, hpbw, SLL, XPL, antPatRes;
	long AntTag;

	/* .. Read from file Inp_Fixed.txt */
    infile=FileOpen(InPath,"Inp_Fixed.txt","r");

	fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
	fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
	// Receiver Altitude [km]
	fscanf(infile, "%lf %[^\n] %[\n]", &hRx, junk, &newline);
	hRx *= 1E3; // [km] to [m]
	// Receive Antenna Gain [dB] and Polarization (R,L,X,Y)
	fscanf(infile, "%lf %s %[^\n] %[\n]", &G_Rx, RxPol, junk, &newline);
	if(!strcmp(RxPol, "RX")){
		NPol = 2;
	}
	else{
		NPol = 1;
	}
	// Receive Antenna Pattern
	fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
	AntTag = decodeStr(response);
	// (GAUSSIAN) HPBW[deg], Sidelobe[dB], X-pol[dB], Res[deg]
	fscanf(infile, "%lf %lf %lf %lf %[^\n] %[\n]", &hpbw, &SLL, &XPL, &antPatRes, junk, &newline);
	// (USER_DEFINED) File Name for Antenna Pattern
	fscanf(infile,"%s %[^\n] %[\n]",AntFileName[0],junk,&newline);
	fscanf(infile,"%s %[^\n] %[\n]",AntFileName[1],junk,&newline);
	fscanf(infile,"%s %[^\n] %[\n]",AntFileName[2],junk,&newline);
	fscanf(infile,"%s %[^\n] %[\n]",AntFileName[3],junk,&newline);
	// (For retrieval only) Orbit type
	fscanf(infile,"%s %[^\n] %[\n]", orbitName,junk,&newline);

	fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
	// Number of Transmitters
	fscanf(infile, "%ld %[^\n] %[\n]", &NSoOp, junk, &newline);

	InputSoOp->F = (struct FixedObsType *) calloc(NSoOp, sizeof(struct FixedObsType));
	if (InputSoOp->F == NULL) {
		printf("InputSoOp->F calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}

	for(ISoOp=0; ISoOp<NSoOp; ISoOp++){
		fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
		// Transmitter Frequency [Hz]
		fscanf(infile, "%lf %[^\n] %[\n]", &InputSoOp->F[ISoOp].Freq, junk, &newline);
		InputSoOp->F[ISoOp].Freq *= 1E6; // [MHz] to [Hz]
		// Transmitter Altitude [km]
		fscanf(infile, "%lf %[^\n] %[\n]", &InputSoOp->F[ISoOp].hTx, junk, &newline);
		InputSoOp->F[ISoOp].hTx *= 1E3; // [km] to [m]
		// Transmitter Range to Earth Center [m]
		InputSoOp->F[ISoOp].rTx = InputSoOp->F[ISoOp].hTx + EarthRad;
		// Transmitter Polarization (R,L,X,Y)
		fscanf(infile, "%c %[^\n] %[\n]", &InputSoOp->F[ISoOp].polTx, junk, &newline);
	}

	for(ISoOp=0; ISoOp<NSoOp; ISoOp++){
		// Receiver Altitude [km]
		InputSoOp->F[ISoOp].hRx = hRx;
		// Receive Antenna Gain [dB]
		InputSoOp->F[ISoOp].G_Rx = G_Rx;
		// Receive Antenna Pattern
		InputSoOp->F[ISoOp].AntTag = AntTag;
		// (GAUSSIAN) HPBW[deg], Sidelobe[dB], X-pol[dB], Res[deg]
		InputSoOp->F[ISoOp].hpbw = hpbw;
		InputSoOp->F[ISoOp].SLL = SLL;
		InputSoOp->F[ISoOp].XPL = XPL;
		InputSoOp->F[ISoOp].AntPatRes = antPatRes;
		// (USER_DEFINED) File Name for Antenna Pattern
		strcpy(InputSoOp->F[ISoOp].AntFileName[0], AntFileName[0]);
		strcpy(InputSoOp->F[ISoOp].AntFileName[1], AntFileName[1]);
		strcpy(InputSoOp->F[ISoOp].AntFileName[2], AntFileName[2]);
		strcpy(InputSoOp->F[ISoOp].AntFileName[3], AntFileName[3]);
		// (For retrieval only) Orbit Type
		strcpy(InputSoOp->F[ISoOp].orbitName, orbitName);
	}
}

void LoadMain(void)
{
	FILE *infile;
    char junk[120], newline, response[120];
	long iLayer;
	double RMSH;

    infile=FileOpen(InPath,"Inp_Main.txt","r");

	fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
	    
	/* .. Number of Input Soil Layers .. */
	NLayer = 5;
	GndData = (struct GndDataType *) calloc(NLayer,sizeof(struct GndDataType));
	if (GndData == NULL) {
		printf("GndData calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	/* .. Multilayered Ground Structure .. */
	MultiLayer = (struct MultiLayerType *) calloc(1,sizeof(struct MultiLayerType));
	if (MultiLayer == NULL) {
		printf("MultiLayer calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	// Total depths [cm]
	fscanf(infile,"%lf %[^\n] %[\n]", &MultiLayer->D,junk,&newline);
	MultiLayer->D *= 1E-2; // [cm] -> [m]
	MultiLayer->D -= 0.05;
	// Layer discretization [mm]
	fscanf(infile,"%lf %[^\n] %[\n]", &MultiLayer->delZ,junk,&newline);
	MultiLayer->delZ *= 1E-3; // [mm] -> [m]
	// Initial inflection point depth [cm]
	fscanf(infile,"%lf %[^\n] %[\n]", &MultiLayer->z_infl,junk,&newline);
	MultiLayer->z_infl *= 1E-2; // [cm] -> [m]
	MultiLayer->z_infl -= 0.05;

	MultiLayer->NSublayer_POME = (long) round(MultiLayer->D/MultiLayer->delZ) + 1;
	MultiLayer->NSublayer_extend = (long) round(0.05/MultiLayer->delZ);
	MultiLayer->NSublayer = MultiLayer->NSublayer_POME + MultiLayer->NSublayer_extend;
	
	MultiLayer->z = (double *) calloc(MultiLayer->NSublayer, sizeof(double));
	if (MultiLayer->z == NULL) {
		printf("MultiLayer->z calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	MultiLayer->SMP = (double *) calloc(MultiLayer->NSublayer_POME, sizeof(double));
	if (MultiLayer->SMP == NULL) {
		printf("MultiLayer->SMP calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	for(iLayer=0; iLayer<MultiLayer->NSublayer; iLayer++){
		MultiLayer->z[iLayer] = iLayer*MultiLayer->delZ;
	}

	/* .. Data Source - USCRN .. */
	Nc = (struct NcDataType *) calloc(1, sizeof(struct NcDataType));
	if (Nc == NULL) {
		printf("Nc calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	fscanf(infile,"%[^\n] %[\n]",junk,&newline);
	// .. USCRN Station List File Name .. //
	fscanf(infile,"%s %[^\n] %[\n]",Nc->filename,junk,&newline);
	// .. Start Year and End Year .. //
	fscanf(infile,"%ld %ld %[^\n] %[\n]",&Nc->startYear, &Nc->endYear, junk,&newline);
	// .. Start Month and End Month .. //
	fscanf(infile,"%ld %ld %[^\n] %[\n]",&Nc->startMonth, &Nc->endMonth, junk,&newline);

    /* .. Static Data .. */
	SttData = (struct SttDataType *) calloc(1,sizeof(struct SttDataType));
	if (SttData == NULL) {
		printf("SttData calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
    // Data On/off
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
	fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
	SttData->ExistVeg = decodeStr(response);
	// Surface Roughness RMSH
	fscanf(infile,"%lf %[^\n] %[\n]",&RMSH,junk,&newline);
	InputSoOp->RMSH = RMSH;

	fclose(infile);
}

void InitSOCRATES_Ret(int argc,char **argv)
{
	int iLayer;
	/* .. Set directories .. */
	sprintf(InPath, "../code/inputs/");
	sprintf(ForwardPath, "../results/forward/");
	sprintf(InversePath, "../results/inverse/");
	strcpy(AntPath, InPath);
        
	/* .. Set simulation Mode .. */
	if (argc > 1){
		if(!strcmp(argv[1], "f")){
			SimMode = FORWARD;
		}
		else if(!strcmp(argv[1], "i")){
			SimMode = INVERSE;
			if (argc > 2) sprintf(InverseSubPath, "%s%s/", InversePath, argv[2]);
			else printf(">> Error: Missing second argument. Specify your output directory.\n");
		}
		else{
			printf(">> Error: Wrong first argument. Specify your simulation mode.\n");
		}
	}
	else{
		printf(">> Error: Arguments missing. Specify your simulation mode.\n");
		exit(1);
	}

	/* .. Initialize Math Constants .. */
	InitConstants();

	/* .. Initialize Input/Output Parameters Structure .. */
	InputSoOp = (struct InputType *) calloc(1, sizeof(struct InputType));
	if (InputSoOp == NULL) {
		printf("InputSoOp calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	Out = (struct OutputType *) calloc(1, sizeof(struct OutputType));
	if (Out == NULL) {
		printf("Output calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	
	/* .. Read Inp_Main.txt */
	LoadMain();

	/* .. Read Inp_Fixed.txt .. */
	LoadFixed();

	/* .. Initialize Local Ground Parameters .. */
	Gnd = (struct LocalGndType *) calloc(1, sizeof(struct LocalGndType));
	if (Gnd == NULL) {
		printf("Local Gnd parameter calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	Gnd->VSM = (float *)calloc(NLayer, sizeof(float));
	if (Gnd->VSM == NULL) {
		printf("Gnd->VSM calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	for(iLayer=0; iLayer<NLayer; iLayer++)
		Gnd->VSM[iLayer] = NAN;
	Gnd->e_c = (double complex *)calloc(NLayer, sizeof(double complex));
	if (Gnd->e_c == NULL) {
		printf("Gnd->e_c calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	Gnd->e_c_profile = (double complex *)calloc(MultiLayer->NSublayer, sizeof(double complex));
	if (Gnd->e_c_profile == NULL) {
		printf("Gnd->e_c_profile calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}
	Gnd->SMP = (double *)calloc(MultiLayer->NSublayer, sizeof(double));
	if (Gnd->SMP == NULL) {
		printf("Gnd->SMP calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}

	/* .. Initialize Bistatic Parameters .. */
	Bistatic = (struct BistaticType *)calloc(1, sizeof(struct BistaticType));
	if (Bistatic == NULL) {
		printf("Bistatic parameter calloc returned null pointer.  Bailing out!\n");
		exit(1);
	}

	/* .. Initialize Fixed Observation Parameters .. */
	InitFixedObs();
	
	/* .. Initialize Vegetation if exist .. */
	if(SttData->ExistVeg){
		LoadVegetation();
		InitVegetation();
	}

	/* .. Free Input variables .. */
	freeInputSoOp();
}