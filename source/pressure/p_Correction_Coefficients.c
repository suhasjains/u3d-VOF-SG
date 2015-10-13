#include"pressure.h"

void pressure_correction_east()
{
	double d;	/* rho*dx */
	double rhoe;
	
  
  	for (i=2;i<=nxm;i++)
	  {
		ieast = i+1;
		iwest = i-1;
		
	  	for (j=2;j<=nym;j++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	{
	      
		      		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
		      		
				rhoe = 0.5 * rho[i][j][k] +   0.5 * rho[ieast][j][k];
				
				/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
	
	      			d = rhoe*s;
	      			
/*	      			if(rhoe!=1)	      		printf("rhoe = %e\n",rhoe);*/
	      			
	      			Fe[i][j][k] = rhoe*u[i][j][k][l]*s; 
	      
	      			ae[i][j][k]  = d*s*apu[i][j][k];
	      			aw[ieast][j][k]= ae[i][j][k];
              		}
	      	  }
	     }
}

void pressure_correction_north()
{
	double d;	/* rho*dx */
	double rhon;
	
  
  for (j=2;j<=nym;j++)
  {
	jnorth = j+1;
	jsouth = j-1;

  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	  		
	  		s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
	  		
	  		rhon = 0.5 * rho[i][j][k] +   0.5 * rho[i][jnorth][k];
			
			/* Harmonic Interpolation of density*/
/*			rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/

			d = rhon * s;
			
			Fn[i][j][k] = rhon*v[i][j][k][l]*s;
      	
	      		an[i][j][k]  = d*apv[i][j][k]*s;
      			as[i][jnorth][k] = an[i][j][k];
      			
     		}
    	}
  }
}

void pressure_correction_top()
{
	double d;	/* rho*dx */
	double rhot;
	
  
  for (k=2;k<=nzm;k++)
  {
	ktop = k+1;
	kbottom = k-1;
	
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
	  		
	  		s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
			
			rhot = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][ktop];
			
			/* Harmonic Interpolation of density*/
/*			rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/

			d = rhot * s;
			
			Ft[i][j][k] = rhot*w[i][j][k][l]*s;
		
	      		at[i][j][k]  = d*apw[i][j][k]*s;
	      		ab[i][j][ktop] = at[i][j][k];
	      		
     		}
    	}
  }
}

void pressure_correction_source()
{
	maxMassResidual=0; /* To check for continuity */
  
  	for (i=2;i<=nxm;i++)
	  {
	  	for (j=2;j<=nym;j++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	 {
		  	 	ieast 	= i + 1;
      				iwest 	= i - 1;
      				jnorth	= j + 1;
      				jsouth	= j - 1;
      				ktop  	= k + 1;
      				kbottom	= k - 1;
            
      				su[i][j][k] = Fe[iwest][j][k] - Fe[i][j][k] + Fn[i][jsouth][k] - Fn[i][j][k] + Ft[i][j][kbottom] - Ft[i][j][k] ;
      				ap[i][j][k] = ae[i][j][k] + aw[i][j][k] + an[i][j][k] + as[i][j][k] + at[i][j][k] + ab[i][j][k];
      				maxMassResidual	= maxMassResidual + fabs(su[i][j][k]);	/* Checking continuity */
      				pp[i][j][k][l] = 0;
    		 }
    	}
  	 }
}

