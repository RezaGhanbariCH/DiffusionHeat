#include <math.h>


void FLUX(double **U, double *Z, int MR,int MZ, double dr, double dz, double D, double **Fr, double **Fz, double time1,int nWRs, int myID)
{

int i,j;
int Mp1,Mp2;

Mp1 = MR+1;

Mp2 = MZ+1;

	for(i=0;i<=Mp1;i++){

		if(myID==1)
		{
		U[0][i] = 0;
	    Fz[1][i] = -2*D*(U[1][i]-U[0][i])*(1/dz);
		}

		else
		{
			Fz[1][i] = -D*(U[1][i]-U[0][i])*(1/dz);
		}

		if(myID==nWRs)
		{
			U[Mp2][i] = 0;
			Fz[Mp2][i] = -2*D*(U[Mp2][i]-U[Mp2-1][i])*(1/dz);
		}

		else
		{
			Fz[Mp2][i] = -D*(U[Mp2][i]-U[Mp2-1][i])*(1/dz);
		}
	}


	for(j=0;j<=Mp2;j++){
		U[j][0] = 0;
		U[j][Mp1]=exp(-time1)*log10(2.0)*sin(Z[j]);}

	for(j=0;j<=Mp2;j++){
		Fr[j][1]=-2*D*(U[j][1]-U[j][0])*(1/dr);
		Fr[j][Mp1]=-2*D*(U[j][Mp1]-U[j][Mp1-1])*(1/dr);}

	for(j=0;j<=Mp2;j++){
	    for(i=2;i<Mp1;i++){
            Fr[j][i]=-D*(U[j][i]-U[j][i-1])*(1/dr);}
	}
			
	for(i=0;i<=Mp1;i++){
        for(j=2;j<Mp2;j++){
            Fz[j][i]=-D*(U[j][i]-U[j-1][i])*(1/dz);}
	}
	}

void PDE(double **Fr,double **Fz, int MR, int MZ, double dr,double dz, double dt, double **U, double *R, double D, double time1)
{
	
int i,j;
int Mp1,Mp2;

Mp1 = MR+1;

Mp2 = MZ+1;

	for(i=1;i<Mp1;i++){
			for(j=1;j<Mp2;j++){
		U[j][i]+=dt*((1/dr)*((R[i]-(dr/2))/R[i])*Fr[j][i]-(1/dr)*((R[i]+(dr/2))/R[i])*Fr[j][i+1]+(1/dz)*Fz[j][i]-(1/dz)*Fz[j+1][i]);}

	}
	}
