/***************************************************************************/
/*                                                                         */
/*  SoOp Constellation and Remote sensing Analysis Tool for Earth Science  */
/*                                                                         */
/*                          < Retrieval Module >                           */
/*                                                                         */
/***************************************************************************/
/*                                                                         */
/*   Written by Seho Kim @ Purdue University                               */
/*   Integrating "SCoBi" by M.Kurum @ Mississippi State University         */
/*                                                                         */
/***************************************************************************/
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

#define DECLARE_GLOBALS
#include "socratesret.h"
#undef DECLARE_GLOBALS

int main(int argc,char **argv)
{
	int provided, taskid;

	/* .. Initialize MPI .. */
	MPI_Init_thread(&argc,&argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	
	/* .. Initialize Simulator .. */
	InitSMAT_ret(argc, argv);
	if(taskid == MASTER)
		Welcome();

	/* .. Run forward/inverse model .. */
	switch(SimMode){
		case FORWARD:
			ObsFixedNetCDF();
			break;
		case INVERSE:
			RetrievalNetCDF();
			break;
		default:
			break;
	}

	/* .. Finalize Simulator .. */
	if(taskid == MASTER)
		printf(">> Simulation End.\n");

	/* .. Finalize MPI .. */
	MPI_Finalize();
	
	return 0;
}

