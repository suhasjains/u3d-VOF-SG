void TemperatureBoundaryConditions()
{

	/* Bottom boundary conditions */
  	j = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	  	
	  		/*Inverse of Prandtl number */
			Prr = Kcond[i][j-1][k]/(Cp[i][j-1][k]*mu[i][j-1][k]);
	
			diff = mu[i][j-1][k]*Prr*(xf[i]-xf[i-1])*(zf[k]-zf[k-1])/(yc[j]-yc[j-1]);
			
			as[i][j][k] = 0;
		
		
			//Wall 
/*    			apT[i][j][k] = apT[i][j][k] + diff;*/
/*	    		sT[i][j][k]  = sT[i][j][k]  + diff*T[i][j-1][k][l];*/
/*			*/
			//zero Gradient
			T[i][j-1][k][l] = T[i][j][k][l];


	  	}
  	  }
  	
	    
  	  /* Top boundary conditions */
  	j = nym;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	  	
	  		/*Inverse of Prandtl number */
			Prr = Kcond[i][j+1][k]/(Cp[i][j+1][k]*mu[i][j+1][k]);

			diff = mu[i][j+1][k]*Prr*(xf[i]-xf[i-1])*(zf[k]-zf[k-1])/(yc[j+1]-yc[j]);
			
			an[i][j][k] = 0;
		
			//Wall
/*  			apT[i][j][k] = apT[i][j][k] + diff;*/
/*	    		sT[i][j][k]  = sT[i][j][k]  + diff*T[i][j+1][k][l];*/


			//zero Gradient
			T[i][j+1][k][l] = T[i][j][k][l];

			
	  	}
	  }
  
  	/* West boundary conditions */
	i = 2;
  	for (j=2;j<=nym;j++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	  	
	  		/*Inverse of Prandtl number */
			Prr = Kcond[i-1][j][k]/(Cp[i-1][j][k]*mu[i-1][j][k]);

			diff = mu[i-1][j][k]*Prr*(yf[j]-yf[j-1])*(zf[k]-zf[k-1])/(xc[i]-xc[i-1]);
			
			aw[i][j][k] = 0;
		
			//wall
    			apT[i][j][k] = apT[i][j][k] + diff;
    			sT[i][j][k]  = sT[i][j][k]  + diff*T[i-1][j][k][l];

			//zero Gradient
/*			T[i-1][j][k][l] = T[i][j][k][l];*/


		}
    }
          
  	/* East boundary conditions */
 	i = nxm;
  	for (j=2;j<=nym;j++)
  	{
  		for(k=2;k<=nzm;k++)
  		{
  		
  			/*Inverse of Prandtl number */
			Prr = Kcond[i+1][j][k]/(Cp[i+1][j][k]*mu[i+1][j][k]);
  		
			diff = mu[i+1][j][k]*Prr*(yf[j]-yf[j-1])*(zf[k]-zf[k-1])/(xc[i+1]-xc[i]);
			
			ae[i][j][k] = 0;
		
			//Wall
    			apT[i][j][k] = apT[i][j][k] + diff;
    			sT[i][j][k]  = sT[i][j][k]  + diff*T[i+1][j][k][l];
	
			//zero Gradient
/*			T[i+1][j][k][l] = T[i][j][k][l];*/

  		}
  	}  
  
 	/* Back bounday conditions */
	k = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
	  	
	  		/*Inverse of Prandtl number */
			Prr = Kcond[i][j][k-1]/(Cp[i][j][k-1]*mu[i][j][k-1]);
	  		
			diff = mu[i][j][k-1]*Prr*(yf[j]-yf[j-1])*(xf[i]-xf[i-1])/(zc[k]-zc[k-1]);
			
			ab[i][j][k] = 0;
		
			//Wall
/*    			apT[i][j][k] = apT[i][j][k] + diff;*/
/*    			sT[i][j][k]  = sT[i][j][k]  + diff*T[i][j][k-1][l];*/

			//zero Gradient
			T[i][j][k-1][l] = T[i][j][k][l];

		}
    	}
          
  	/* Front boundary conditions */
  	k = nzm;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
	  	
	  		/*Inverse of Prandtl number */
			Prr = Kcond[i][j][k+1]/(Cp[i][j][k+1]*mu[i][j][k+1]);

			diff = mu[i][j][k+1]*Prr*(yf[j]-yf[j-1])*(xf[i]-xf[i-1])/(zc[k+1]-zc[k]);
			
			at[i][j][k] = 0;
			
			//Wall
/*	    		apT[i][j][k] = apT[i][j][k] + diff;*/
/*    			sT[i][j][k]  = sT[i][j][k]  + diff*T[i][j][k+1][l];*/
		
			//Zero Gradient
			T[i][j][k+1][l] = T[i][j][k][l];
		  		
			
		}
    	}

}

void setting_Temperature()
{
  j=ny;						// Top boundary 
  for (i=2;i<=nxm;i++)
   for(k=2;k<=nzm;k++)
    T[i][j][k][2] = TopBoundTemp;

   
   j=1;						// Bottom boundary 
  for (i=2;i<=nxm;i++)
   for(k=2;k<=nzm;k++)
    T[i][j][k][2] = BottomBoundTemp;

   
   k=nz;					// Front boundary 
  for (i=2;i<=nxm;i++)
   for(j=2;j<=nym;j++)
    T[i][j][k][2] = FrontBoundTemp;
   
   k=1;						// Back boundary 
  for (i=2;i<=nxm;i++)
   for(j=2;j<=nym;j++)
    T[i][j][k][2] = BackBoundTemp;
   
   
   i=nx;					// Right boundary 
  for (j=2;j<=nym;j++)
   for(k=2;k<=nzm;k++)
    T[i][j][k][2] = RightBoundTemp;
   
   i=1;						// Left boundary 
  for (j=2;j<=nym;j++)
   for(k=2;k<=nzm;k++)
    T[i][j][k][2] = LeftBoundTemp;
   
   
}

