//WORKERS

#include <stdio.h>
#include<malloc.h>
#include <stdlib.h>
#include <math.h>

#include </home/mghanbar/Math578/lab9a/update.h>



void WORKER(int nWRs,int myID)

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

	double *R, *Z;

    int *iparms,Niparms,Nparms;
	double *parms;
	int ierr;
	int NodeRT, NodeLT;
	double lz0,lzl;

    Nparms = 10;
//	Niparms = 1;
	parms = (double*)malloc(sizeof(double)*(Nparms)); 
//	iparms = (int*)malloc(sizeof(int)*(Niparms)); 


//----------reading variables from a file---------


//---------unpacking initial iparms and parms-----------
//    ierr = MPI_Bcast(iparms,Niparms,MPI_INT,1,MPI_COMM_WORLD);   
	ierr = MPI_Bcast(parms,Nparms,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   
    Mr= parms[8];
	Rin = parms[0];  Rout= parms[1]; z0 = parms[2] ; zl = parms[3]; D = parms[4] ; dtfactor = parms[5] ; dtout = parms[6]; tend = parms[7]; Mz = parms[9];

	    NodeRT = myID+1;//---deteriming the right node of the current node---
		NodeLT = myID-1;//---deteriming the left node of the current node---

	
    MR = Mr*(Rout-Rin); // =(Rout-Rin)/dr
	dr = 1.0/(double)Mr; 

	MZt = Mz*(zl-z0); // (zl-z0)/dz
	M = MZt/nWRs;
	dz = 1.0/(double)Mz;
	
	Mp1 = MR+1;
	
	Mp2 = M+1;

	dtEXPL=((dr*dr*dz*dz)/(dr*dr+dz*dz))*(1/(2*D));
	dt = dtfactor*dtEXPL;
	R = (double*)malloc(sizeof(double)*(Mp1+1));
	Z = (double*)malloc(sizeof(double)*(Mp2+1));
	
	U = malloc((Mp2+1)*sizeof(double*));
    double *t = malloc((Mp2+1)*(Mp1+1)*sizeof(double));
    for(i = 0; i <= Mp2; i++)
    U[i] = t + i * (Mp1+1);
	
	Fr = malloc((Mp2+1)*sizeof(double*));
    double *t1 = malloc((Mp2+1)*(Mp1+1)*sizeof(double));
    for(i = 0; i <= Mp2; i++)
    Fr[i] = t1 + i * (Mp1+1);

	Fz = malloc((Mp2+1)*sizeof(double*));
    double *t2 = malloc((Mp2+1)*(Mp1+1)*sizeof(double));
    for(i = 0; i <= Mp2; i++)
    Fz[i] = t2 + i * (Mp1+1);

	lz0 = z0 + (myID-1)*M*dz;

	lzl = z0 + myID*M*dz;

	if(myID==nWRs)
		lzl = zl;


   MESH(R,Z,MR,M,dr,dz,Rin,Rout,lz0,lzl,nWRs,myID); 


//-----------initialize-------

	nsteps = 0; 
	time1 = 0; 
	tout = dtout; 
	maxSteps = tend/dt + 1;

	INIT(R,Z,U,MR,M,time1,D);

	SEND_output_MPI(myID,MR,M,U);//---send initial U to Master----




// --------begin timestepping----------

	for(nsteps=1;nsteps<=maxSteps;nsteps++)
	{
 
	
//------------synchornize everyone--------------

		ierr = MPI_Barrier(MPI_COMM_WORLD);//---synchronization of all of the workers---

//------------Exchange "boundary" values---------

        EXCHANGE_bry_MPI(nWRs,myID,NodeRT,NodeLT,MR,M,U);//---boundary values are exchnaged---
    
		FLUX(U,Z,MR,M,dr,dz,D,Fr,Fz,time1,nWRs,myID); 
		
        PDE(Fr,Fz,MR,M,dr,dz,dt,U,R,D,time1);

		time1+=dt;

		if(time1>= tout)
		{    

//--------------------send output to the MASTER----------------

           SEND_output_MPI(myID,MR,M,U);
			tout+=dtout;
		}
		  
	}

	EXCHANGE_bry_MPI(nWRs,myID,NodeRT,NodeLT,MR,M,U);//---final boundary values are exchanged---
	ierr = MPI_Barrier(MPI_COMM_WORLD);//---synchronization of all of the workers---

	SEND_output_MPI(myID,MR,M,U);//---final Us are sent to Master---

//------------end if timestepping-----------------

}


