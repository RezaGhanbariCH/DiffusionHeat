#include "mpi.h"


void RECV_output_MPI(int nWRs,int MR, int M,double **U)
{

    int i,I2,J2;
	int msgtag,Ime,ierr;

//-------------only MR does this---------------

	I2 = MR + 2;    //size of U(i)
	J2 = M + 2;
    MPI_Status stat;
	MPI_Datatype mytype;
    ierr = MPI_Type_contiguous(I2,MPI_DOUBLE,&mytype);
	ierr = MPI_Type_commit(&mytype);

//-------------------receive values from everybody for output--------------------

	for(i=1;i<=nWRs;i++)
		{

		Ime=(i-1)*M;
	  	msgtag = 1000 + i;

				MPI_Recv(&U[Ime][0],J2, mytype, i, msgtag, MPI_COMM_WORLD, &stat);//---receive local U from the Workers by Master---
	}
	ierr = MPI_Type_free(&mytype);
}

void SEND_output_MPI(int myID,int MR, int M,double **U)
{

    int I2,J2;
	int msgtag, mster,ierr;

//-------------every WR does this---------------

	I2 = MR + 2;    //size of U(i)
    
    MPI_Datatype mytype;
    ierr = MPI_Type_contiguous(I2,MPI_DOUBLE,&mytype);
	ierr = MPI_Type_commit(&mytype);

//-------------------everybody sends values to the Master for the output--------------------

    mster = 0;
	J2 = M + 2;
	msgtag = 1000 + myID;
 
				MPI_Send(&U[0][0],J2, mytype, mster, msgtag, MPI_COMM_WORLD);//---send local U to the Master from Workers---
	
	ierr = MPI_Type_free(&mytype);
}


void EXCHANGE_bry_MPI(int nWRs,int myID,int NodeRT,int NodeLT,int MR,int M,double **U)
{

//---------Exchange "boundary" values btn neighbors---------

//--------------every WR does this----------------

    int I2,Ib,Ib1,ierr;
     
    int msgtag;

    int msgrt,msglt;

	msgrt = 10;
	msglt = 20;

    I2 =MR+2;

	Ib = M;

	Ib1 = Ib +1;

    MPI_Status stat;
	MPI_Datatype mytype;
    ierr = MPI_Type_contiguous(I2,MPI_DOUBLE,&mytype);
	ierr = MPI_Type_commit(&mytype);


//-------sending data from nodes to their left ones-------*except for node = 1, all would do this------

    if(myID != 1)
    {
	msgtag = msglt;
    ierr =  MPI_Send(&U[1][0],1, mytype, NodeLT, msgtag, MPI_COMM_WORLD); //---sending bv to the left node----
	}
    
	if(myID != nWRs)
    {
    msgtag = msglt;
	ierr =  MPI_Recv(&U[Ib1][0],1, mytype, NodeRT, msgtag, MPI_COMM_WORLD, &stat);//---receiving bv from the right node----
	}

    if(myID != nWRs)
    {
	msgtag = msgrt;
    ierr =  MPI_Send(&U[Ib][0],1, mytype, NodeRT, msgtag, MPI_COMM_WORLD); //---sending bv to the right node----
	}

	if(myID != 1)
    {
    msgtag = msgrt;
	ierr =  MPI_Recv(&U[0][0],1, mytype, NodeLT, msgtag, MPI_COMM_WORLD, &stat);//---receiving bv from the left node----
	}

ierr = MPI_Type_free (&mytype);
}






