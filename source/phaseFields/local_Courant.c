void localCourantNumber()
{


	double rhoe,rhon,rhot;	

	for(i=2;i<=nxm-1;i++)	//Calculation of local Courant number only at internal faces
  	{
		ieast 	= i + 1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;

    		for(j=2;j<=nym;j++)
		{
			jnorth	= j + 1;
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				dzpt = zc[ktop]-zc[k];
				fzt = (zf[k] - zc[k])/dzpt;
				fzp = 1.0 - fzt;

			
				
				//Calculating density at cell interface
				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];
				
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/

				s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
	      			vole = dxpe * s;


 				//Sum of courant numbers of outflow faces of donor cell
				Ce[i][j][k] = fabs(Fe[i][j][k]/(rhoe*vole))*dt;

/*				printf("Ce=%e\n",Ce[i][j][k]);*/
			}
		}
	}
	
	
		for(i=2;i<=nxm;i++)	//Calculation of local Courant number only at internal faces
  	{
		ieast 	= i + 1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;

    		for(j=2;j<=nym-1;j++)
		{
			jnorth	= j + 1;
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				dzpt = zc[ktop]-zc[k];
				fzt = (zf[k] - zc[k])/dzpt;
				fzp = 1.0 - fzt;

			
				
				//Calculating density at cell interface
				rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];
				
/*				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/


				s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
				voln = s * dypn;


 				//Sum of courant numbers of outflow faces of donor cell
				Cn[i][j][k] = fabs(Fn[i][j][k]/(rhon*voln))*dt;
				
/*				printf("Ce=%e\n",Ce[i][j][k]);*/
			}
		}
	}
	
	for(i=2;i<=nxm;i++)	//Calculation of local Courant number only at internal faces
  	{
		ieast 	= i + 1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;

    		for(j=2;j<=nym;j++)
		{
			jnorth	= j + 1;
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
	
	  		for(k=2;k<=nzm-1;k++)
			{
				ktop  	= k + 1;
				dzpt = zc[ktop]-zc[k];
				fzt = (zf[k] - zc[k])/dzpt;
				fzp = 1.0 - fzt;

			
				
				//Calculating density at cell interface
				rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];
				
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/


				s = (xf[i]-xf[i-1])*(yf[j]-yf[j-1]);
				volt = s * dzpt;


 				//Sum of courant numbers of outflow faces of donor cell
				Ct[i][j][k] = fabs(Ft[i][j][k]/(rhot*volt))*dt;

				
				
/*				printf("Ce=%e\n",Ce[i][j][k]);*/
			}
		}
	}
	
	
	for(i=2;i<=nxm;i++)	//Calculation of local Courant number only at internal faces
  	{

    		for(j=2;j<=nym;j++)
		{
	
	  		for(k=2;k<=nzm;k++)
			{
			
				COutD[i][j][k] = Ce[i][j][k] + Cn[i][j][k] + Ct[i][j][k];
/*				printf("COutD=%lf\n",COutD[i][j][k]);*/
/*				printf("Ce=%e\n",Ce[i][j][k]);*/
/*				printf("Cn=%e\n",Cn[i][j][k]);*/
/*				printf("Ct=%e\n",Ct[i][j][k]);*/
			}
		}
	}
	
	
} 
				
				
