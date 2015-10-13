#include"update.h"

//Updates the values of the pressure
void update_Pressure_Velocity()
{
	double ppe,ppw,ppn,pps,ppt,ppb; /* Linear interpolation of pprime */
	double rhoe,rhon,rhot;
	
	
  
  for(i=2;i<=nxm;i++)
  {

    	for(j=2;j<=nym;j++)
	  {
	
	  	for(k=2;k<=nzm;k++)
		  {
			
      		
			p[i][j][k][l] = p[i][j][k][l] + pUnderRelaxCoeff*(pp[i][j][k][l] - ppo);
      
		  } 
      	}
  }
  
  
	for  (k=2; k<=nzm; k++)
	{
	 	for (i=2; i<=nxm-1;i++)
	 	{
		  	for(j=2; j<=nym; j++)
		  	{
			  	s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
			  	
			  	rhoe = 0.5 * rho[i][j][k] +   0.5 * rho[ieast][j][k];
	  	
			  	/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
						
        			u[i][j][k][l] = u[i][j][k][l] + apu[i][j][k] * (pp[i][j][k][l] - pp[i+1][j][k][l]) * s;
        			
				Fe[i][j][k] = rhoe*u[i][j][k][l]*s; 
	 	 	}
		}
	}


	//updates v velocity
	for  (i=2; i<=nxm; i++)
	{
	 	for (j=2; j<=nym-1;   j++)
	  	{
	  	  	for(k=2; k<=nzm; k++)
          		{
          			s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
          			
          			rhon = 0.5 * rho[i][j][k] +   0.5 * rho[i][jnorth][k];
          	
          			/* Harmonic Interpolation of density*/
/*				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/
				
          			v[i][j][k][l] = v[i][j][k][l] + apv[i][j][k] * (pp[i][j][k][l] - pp[i][j+1][k][l]) * s;
          	
          			Fn[i][j][k] = rhon*v[i][j][k][l]*s;
	  		}
	  	}
	  }


	//updates w velocity
	for  (j=2; j<=nym; j++)
	{
	 	for (k=2; k<=nzm-1;   k++)
	  	{
	  		for(i=2; i<=nxm; i++)
          		{
          			s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
          			
          			rhot = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][ktop];
          	
          			/* Harmonic Interpolation of density*/
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/
		
          			w[i][j][k][l] = w[i][j][k][l] + apw[i][j][k] * (pp[i][j][k][l] - pp[i][j][k+1][l]) * s;
		
				Ft[i][j][k] = rhot*w[i][j][k][l]*s;
  	  		}
  	  	}
  	 }
}


