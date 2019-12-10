void MESH(double *R, double *Z, int MR, int MZ, double dr,double dz,double Rin, double Rout, double z0, double zl,int nWRs, int myID)
{

int i,j;
int Mp1,Mp2;
	
Mp1 = MR+1;

Mp2 = MZ+1;

R[0] = Rin;
R[1] = Rin+dr/2;
R[Mp1] = Rout;

for(i=2;i<Mp1;i++)
	R[i] =R[i-1]+dr;


	if(myID == 0 || myID == 1)
		{

          Z[0] = z0;
	    }
	else
	    {
		  Z[0] = z0-dz/2;
	    }

	    Z[1]=z0+dz/2;

		if(myID == nWRs || myID == 0)
	    {
          Z[Mp2] = zl;
			
	    }
        else 
	    {
		  Z[Mp2] = zl+dz/2;
	    }
          for(i=2;i<=MZ;i++){

		     	Z[i] = Z[i-1]+dz;}

	}
void INIT(double *R, double *Z, double **U, int MR, int MZ, double time1, double D)
{

int i,j;
int Mp1,Mp2;

Mp1 = MR+1;

Mp2 = MZ+1;

    for(i=0;i<=Mp1;i++){
       for(j=0;j<=Mp2;j++){
        U[j][i]= log10(R[i])*sin(Z[j]);}
	}
	}
