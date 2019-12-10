//this code deiscretizes a heat/diffusion equation with defined ICs and BCs. (for the cylinderical case
//numerical and exact solutions for some specified BC and IC have been compared
//The current 2-D code has been parallelized in one dimension of z. 
//Author: Mohammad Reza Ghanbari


// z.main.c

#include "mpi.h"
#include <stdio.h>
#include<malloc.h>
#include <stdlib.h>
#include <math.h>

#include </home/mghanbar/Math578/lab9a/setup.h>
#include </home/mghanbar/Math578/lab9a/mainMR.h>
#include </home/mghanbar/Math578/lab9a/mainWR.h>


int main(int argc, char **argv)

{

	int nPROC; //nWRs+1;
	int myID,ierr;
	double  tt0, tt1;

	MPI_Status stat;

	ierr = MPI_Init(&argc,&argv);
	ierr = MPI_Comm_size(MPI_COMM_WORLD,&nPROC);   //..returns nPROC
	ierr = MPI_Comm_rank(MPI_COMM_WORLD,&myID);    //...assigns my ID
	
    //-----nPROC is specified at miprun; see makefile-----

	int mster = 0;     //---RANK = 0
	int nWRs = nPROC -1;    //---NUMBER OF WORKERS

	//----------------start 0,...,nWRs tasks----------------

	printf("number of workers is = %d", nWRs);

	if(myID == mster)

	{
		tt0 = MPI_Wtime();    //...start CPU time on MR

        MASTER(nWRs,myID);

		tt1 = MPI_Wtime();    //...end timer

		printf("\n\n\nMR <<<<<timing>>>>> on nWRs %d is = <<<<<%lf>>>>>\n\n\n", nWRs, tt1-tt0);

	}

		else
	{

		WORKER(nWRs,myID);
		printf("\n\n\nworker number %d is done\n\n\n", myID);

		if(ierr!=0)
		{
			printf("worker number %d was ended", myID);
		}

		
	}

	ierr = MPI_Finalize();

	return 0;
}


