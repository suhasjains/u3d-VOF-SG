void PsiEqnCoeff()
{

	
	
	
	//Constructing convective fluxes along east face
	for  (i=wcID[e][f][g]-1; i<=ecID[e][f][g]; i++)
 	{
		ieast = i+1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;
		fxp = 1.0 - fxe;
		
/*		if(i>=nxm)	continue;*/
/*		if(i<2)		continue;*/

	 	for (j=scID[e][f][g]; j<=ncID[e][f][g]; j++)
   		{
	  		for(k=bcID[e][f][g]; k<=tcID[e][f][g]; k++)
			{
		
		
				/* Evaluation of cell face area */
				s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
				
		
				upwind_implicit(Fe);

			      		
		      		ae[i][j][k]     =  convf/2.0;
		      		aw[ieast][j][k] = -convp/2.0;
		      		
		      		if(Fe[i][j][k]>0)       ePsi[i][j][k] = Psi[i][j][k][l-1];
		      		else                    ePsi[i][j][k] = Psi[ieast][j][k][l-1];
		      		
		      		
/*		      		if(ePsi[i][j][k]>1)        ePsi[i][j][k] = 1.0;*/
/*		      		if(ePsi[i][j][k]<0)        ePsi[i][j][k] = 0.0;    */
		      		
		      		wPsi[ieast][j][k] = ePsi[i][j][k]; 
		      		
/*		      		printf("ePsi = %e \n",ePsi[i][j][k]);*/
		      		
				
			}
  		}
	}
	
	
	
	//Constructing convective fluxes along north face
	for(j=scID[e][f][g]-1;j<=ncID[e][f][g];j++)
	{
		jnorth = j+1;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;
	    	fyp = 1.0 - fyn;
	    	
/*	    	if(j>=nym)	continue;*/
/*	    	if(j<2)		continue;*/

		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				
				s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);

      				upwind_implicit(Fn);
      			
      					
      				an[i][j][k]	 =  convf/2.0;
      				as[i][jnorth][k] = -convp/2.0;
      				
      				if(Fn[i][j][k]>0)       nPsi[i][j][k] = Psi[i][j][k][l-1];
		      		else                    nPsi[i][j][k] = Psi[i][jnorth][k][l-1];
		      		
/*		      		if(nPsi[i][j][k]>1)        nPsi[i][j][k] = 1.0;*/
/*		      		if(nPsi[i][j][k]<0)        nPsi[i][j][k] = 0.0;*/
		      		
		      		sPsi[i][jnorth][k] = nPsi[i][j][k];

				
      
			}
		}
	}
	
	
	//Constructing convective fluxes along top face
	for(k=bcID[e][f][g]-1;k<=tcID[e][f][g];k++)
	{
		ktop = k+1;
		dzpt = zc[ktop]-zc[k];
		fzt = (zf[k] - zc[k])/dzpt;
		fzp = 1.0 - fzt;
		
/*		if(k>=nzm)	continue; */
/*		if(k<2)		continue;*/
		
		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
				s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
      
      				upwind_implicit(Ft);
      			
      				at[i][j][k]   	=  convf/2.0;
      				ab[i][j][ktop]	= -convp/2.0;
      				
      				if(Ft[i][j][k]>0)       tPsi[i][j][k] = Psi[i][j][k][l-1];
		      		else                    tPsi[i][j][k] = Psi[i][j][ktop][l-1];
		      		
/*		      		if(tPsi[i][j][k]>1)        tPsi[i][j][k] = 1.0;*/
/*		      		if(tPsi[i][j][k]<0)        tPsi[i][j][k] = 0.0;*/
		      				      		
		      		bPsi[i][j][ktop] = tPsi[i][j][k];

      				
      							
				
			}
		}
	}
	
	//Adding unsteady contribution to source terms
	for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
			
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				apt	= vol/dt;
				
				if(simulationTime==dt)
				{
					sP[i][j][k]	= sP[i][j][k] + apt*Psi[i][j][k][l-1];
					apP[i][j][k]	= apP[i][j][k] + apt;  
				}
				else
				{
					sP[i][j][k]	= sP[i][j][k] + (1.0 + Gamma)*apt*Psi[i][j][k][l-1] - 0.5*Gamma*apt*Psi[i][j][k][l-2];
					apP[i][j][k]	= apP[i][j][k] + (1.0 + 0.5*Gamma)*apt;  
				}			
			}
		}
	}
	
	
	
	
	


}
