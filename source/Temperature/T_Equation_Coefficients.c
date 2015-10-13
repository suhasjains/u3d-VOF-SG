void TEqnCoeff()
{

	
	
	
	//Constructing convective fluxes along east face 
	for  (i=wcID[e][f][g]-1; i<=ecID[e][f][g]; i++)
 	{
		ieast = i+1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;
		fxp = 1.0 - fxe;
		
		if(i>=nxm)	continue;
		if(i<2)		continue;

	 	for (j=scID[e][f][g]; j<=ncID[e][f][g]; j++)
   		{
	  		for(k=bcID[e][f][g]; k<=tcID[e][f][g]; k++)
			{
		
		
				/* Evaluation of cell face area */
				s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
				
		
				
				/*Interpolating viscosity field*/	
				mue = fxp*mu[i][j][k] + fxe*mu[ieast][j][k];
				
				/*Interpolating conductivity field*/
				Kconde = fxp*Kcond[i][j][k] + fxe*Kcond[ieast][j][k];
				
				/*Interpolating conductivity field*/
				Cpe = fxp*Cp[i][j][k] + fxe*Cp[ieast][j][k];
				
/*				printf("Kconde = %e \n",Kconde);*/
				
				/*Inverse of Prandtl number */
				Prr = Kconde/(Cpe*mue);	
				
/*				printf("Prr = %e \n",Prr);			*/
				
				/* Evaluation of diffusive term */
				diff = mue*Prr*s/dxpe;
		
/*				printf("diff = %e \n",diff);*/
		
				upwind_implicit(Fe);

			      		
		      		ae[i][j][k]     =  convf - diff;
		      		aw[ieast][j][k] = -convp - diff;
		      		
/*		      		printf("Kconde = %e \n",convp);*/
		      		
		      		
				fTuds = convp*T[i][j][k][l] + convf*T[ieast][j][k][l];
				fTcds = Fe[i][j][k]*(T[ieast][j][k][l]*fxe + T[i][j][k][l]*fxp);
				
				sT[i][j][k]     = sT[i][j][k]     + DCvalue*(fTuds - fTcds);
      				sT[ieast][j][k] = sT[ieast][j][k] - DCvalue*(fTuds - fTcds);
      				
/*      				printf("Kconde = %e \n",fTuds);*/
				
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
	    	
	    	if(j>=nym)	continue;
	    	if(j<2)		continue;

		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				
				s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
      
				
				/*Interpolating viscosity field*/
				mun = fyp*mu[i][j][k] + fyn*mu[i][jnorth][k];
				
				/*Interpolating conductivity field*/
				Kcondn = fyp*Kcond[i][j][k] + fyn*Kcond[i][jnorth][k];
				
				/*Interpolating conductivity field*/
				Cpn = fyp*Cp[i][j][k] + fyn*Cp[i][jnorth][k];

				/*Inverse of Prandtl number */
				Prr = Kcondn/(Cpn*mun);
								
				/* Evaluation of diffusive term */
		      		diff = mun*Prr*s/dypn;
      	
      				upwind_implicit(Fn);
      			
      					
      				an[i][j][k]	 =  convf - diff;
      				as[i][jnorth][k] = -convp - diff;

				
      				fTuds = convp*T[i][j][k][l] + convf*T[i][jnorth][k][l];
				fTcds = Fn[i][j][k]*(T[i][jnorth][k][l]*fyn + T[i][j][k][l]*fyp);

				sT[i][j][k]      = sT[i][j][k]      + DCvalue*(fTuds - fTcds);
      				sT[i][jnorth][k] = sT[i][jnorth][k] - DCvalue*(fTuds - fTcds);
      				
/*      				printf("Kconde = %e \n",diff);*/
      
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
		
		if(k>=nzm)	continue; 
		if(k<2)		continue;
		
		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
				s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
      
      
				/*Interpolating viscosity field*/				
				mut = fzp*mu[i][j][k] + fzt*mu[i][j][ktop];
				
				/*Interpolating conductivity field*/
				Kcondt = fzp*Kcond[i][j][k] + fzt*Kcond[i][j][ktop];
				
				/*Interpolating conductivity field*/
				Cpt = fzp*Cp[i][j][k] + fzt*Cp[i][j][ktop];
				
				/*Inverse of Prandtl number */
				Prr = Kcondt/(Cpt*mut);
								
				/* Evaluation of diffusive term */
			      	diff = mut*Prr*s/dzpt;
      
      				upwind_implicit(Ft);
      			
      				at[i][j][k]   	=  convf - diff;
      				ab[i][j][ktop]	= -convp - diff;


      				fTuds = convp*T[i][j][k][l] + convf*T[i][j][ktop][l];
				fTcds = Ft[i][j][k]*(T[i][j][ktop][l]*fzt + T[i][j][k][l]*fzp);
      							
      				sT[i][j][k]      = sT[i][j][k]      + DCvalue*(fTuds - fTcds);
      				sT[i][j][ktop] 	 = sT[i][j][ktop]   - DCvalue*(fTuds - fTcds);
      				
/*      				printf("Kconde = %e \n",diff);*/
      				
      							
				
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
				
				apt	= rho[i][j][k]*vol/dt;
				
				if(simulationTime==dt)
				{
					sT[i][j][k]	= sT[i][j][k] + apt*T[i][j][k][l-1];
					apT[i][j][k]	= apT[i][j][k] + apt;  
				}
				else
				{
					sT[i][j][k]	= sT[i][j][k] + (1.0 + Gamma)*apt*T[i][j][k][l-1] - 0.5*Gamma*apt*T[i][j][k][l-2];
					apT[i][j][k]	= apT[i][j][k] + (1.0 + 0.5*Gamma)*apt;  
				}			
			}
		}
	}
	
	
	
	
	


}
