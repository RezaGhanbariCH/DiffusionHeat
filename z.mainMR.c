// MASTER

#include "mpi.h"
#include <stdio.h>
#include<malloc.h>
#include <stdlib.h>
#include <math.h>



#include </home/mghanbar/Math578/lab9a/io.h>
#include </home/mghanbar/Math578/lab9a/messaging.h>

double COMPARISON(double *,double *, double **, int, int,double, double);

void MASTER(int nWRs, int myID)

{
	int i,Mp1,Mp2,MR,MZt,M,maxSteps, nsteps;
	double D,Rin,Rout,z0,zl,dtfactor;
	double dtout,tend;
	double dr,dz;	
	double dtEXPL,dt;
    double time1;
    double tout;
	double Mz,Mr;

	double **U,**Fr,**Fz, Ue;

    double *R, *Z;//, Fr[n][n],Fz[n][n];

    int *iparms,Niparms,Nparms;
	double *parms;
	int ierr;
	
    Nparms = 10;
	//Niparms = 1;

	parms = (double*)malloc(sizeof(double)*(Nparms)); 
	//iparms = (int*)malloc(sizeof(int)*(Niparms)); 

//----------reading variables from a file---------


INPUT(&Rin, &Rout, &z0, &zl, &D, &dtfactor, &dtout, &tend, &Mr, &Mz);	  //---read runtime parameters---


//	Mr = (int)Mr1;

    MR = Mr*(Rout-Rin); // =(Rout-Rin)/dr
	dr = 1.0/Mr;

	MZt = Mz*(zl-z0); // =(zl-z0)/dz
    M=MZt/nWRs ;    //---mesh size for each node
	dz = 1.0/Mz;
	
	Mp1 = MR+1;
	
	Mp2 = MZt+1;

	dtEXPL=((dr*dr*dz*dz)/(dr*dr+dz*dz))*(1/(2*D));
	dt = dtfactor*dtEXPL;

    printf("\n\n\n\n -----M is------ %d", M);

	R = (double*)malloc(sizeof(double)*(Mp1+1));
	Z = (double*)malloc(sizeof(double)*(Mp2+1));

	U = malloc((Mp2+1)*sizeof(double*));
    double *t = malloc((Mp2+1)*(Mp1+1)*sizeof(double));
    for(i = 0; i <= Mp2; i++)
    U[i] = t + i * (Mp1+1);
	
	Fr = malloc((Mp2+1)*sizeof(double*));
    double *t1 = malloc((Mp2+1)*(Mp1+1)*sizeof(double));
    for(i = 0; i <=Mp2; i++)
    Fr[i] = t1 + i * (Mp1+1);

	Fz = malloc((Mp2+1)*sizeof(double*));
    double *t2 = malloc((Mp2+1)*(Mp1+1)*sizeof(double));
    for(i = 0; i <=Mp2; i++)
    Fz[i] = t2 + i * (Mp1+1);   

	MESH(R,Z,MR,MZt,dr,dz,Rin,Rout,z0,zl,nWRs,myID);

    parms[8] = Mr;
	parms[0] = Rin; parms[1] = Rout; parms[2] = z0; parms[3] = zl; parms[4] = D; parms[5] = dtfactor; parms[6] = dtout; parms[7] = tend; parms[9] = Mz;
  
//	ierr = MPI_Bcast(iparms,Niparms,MPI_INT,1,MPI_COMM_WORLD);
	ierr = MPI_Bcast(parms,Nparms,MPI_DOUBLE,0,MPI_COMM_WORLD); 

	nsteps = 0; 
	time1 = 0;
    tout = dtout; 
	maxSteps = tend/dt + 1;

	RECV_output_MPI(nWRs,MR,M,U); //---receive output at time 0 from workers*/
		
	OUTPUT(R,Z,U,MR,MZt,time1);

	for(nsteps=1;nsteps<=maxSteps;nsteps++)
	{
 
//------------synchornize everyone--------------

		ierr = MPI_Barrier(MPI_COMM_WORLD);//---synchrnozing all of the workers

		time1+=dt;

		if(time1>= tout)
		{    

//--------------------receive output from Workers;----------------

            RECV_output_MPI(nWRs,MR,M,U);
           // output(X,U,Mt,time1); //----profile at time tout-----
			tout+=dtout;
		}
		  
	}

	ierr = MPI_Barrier(MPI_COMM_WORLD);
	RECV_output_MPI(nWRs,MR,M,U);//---final receive from all workers---
 //  OUTPUT(R,Z,U,MR,MZt,time1);
	
	printf("\n\n\n\n\n <<<<<maximum error>>>>> is <<<<<%e>>>>>\n\n\n\n\n",COMPARISON(R,Z,U,MR,MZt,time1,D));//-----error would be printed out------
}

double COMPARISON(double *R, double *Z, double **U, int MR, int MZt, double time1, double D)	
{
int i,j;
int Mp1,Mp2;

Mp1 = MR+1;

Mp2 = MZt+1;
double error;
FILE *fouti1;
double max;
double Ue; 

fouti1 = fopen("values","w");

 for (i = 0; i <= Mp1; i++) {
				for (j = 0; j <= Mp2; j++) 
		{

		Ue = exp(-time1)*log10(R[i])*sin(Z[j]);;

		error =  fabs(U[j][i]-Ue);

		if(error>max)

			max = error;

  fprintf(fouti1, "%le\t%le\t%le\t%le\t%le\n", R[i],Z[j], U[j][i], Ue, error);
		}
	}

	
		return max;
}

