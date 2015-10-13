void surfaceTension()
{


	double rhoe,rhon,rhot;

/*	for(i=2;i<=nxm-1;i++)	//Calculation of components of Surface Tension*/
/*	{*/
/*		ieast 	= i + 1;*/
/*		iwest 	= i-1;*/
/*		dx = xf[i]-xf[iwest];*/

/*    		for(j=2;j<=nym;j++)*/
/*		{*/
/*	  		for(k=2;k<=nzm;k++)*/
/*			{*/
/*				*/
/*				if(curvature==2)	xST[i][j][k] = -6.0*sigma*Kappa[i][j][k]*xGradAlphaS[i][j][k][l]*alphaS[i][j][k][l]*(1.0 - alphaS[i][j][k][l]);*/
/*				else if(curvature==1)	xST[i][j][k] = -sigma*Kappa[i][j][k]*(eAlpha[i][j][k] - wAlpha[i][j][k])/dx;*/
/*				*/
/*				else if(curvature==3)	xST[i][j][k] = -sigma*KappaHF[i][j][k]*(eAlpha[i][j][k] - wAlpha[i][j][k])/dx;*/
/*				*/
/*				xST[i][j][k] = (xSTe[i][j][k] + xSTe[iwest][j][k])/2.0;  */
/*				*/
/*				*/
/*			}*/
/*		}*/
/*	}*/
/*	*/
/*	for(i=2;i<=nxm;i++)	//Calculation of components of Surface Tension*/
/*	{*/

/*    		for(j=2;j<=nym-1;j++)*/
/*		{*/
/*			jnorth	= j + 1;*/
/*			jsouth = j-1;*/
/*			dy = yf[j]-yf[jsouth];*/
/*	*/
/*	  		for(k=2;k<=nzm;k++)*/
/*			{*/
/*			*/
/*				if(curvature==2)	yST[i][j][k] = -6.0*sigma*Kappa[i][j][k]*yGradAlphaS[i][j][k][l]*alphaS[i][j][k][l]*(1.0 - alphaS[i][j][k][l]);*/
/*				else if(curvature==1)	yST[i][j][k] = -sigma*Kappa[i][j][k]*(nAlpha[i][j][k] - sAlpha[i][j][k])/dy;*/
/*				*/
/*				else if(curvature==3)	yST[i][j][k] = -sigma*KappaHF[i][j][k]*(nAlpha[i][j][k] - sAlpha[i][j][k])/dy;*/

/*				yST[i][j][k] = (ySTn[i][j][k] + ySTn[i][jsouth][k])/2.0; */

/*			}*/
/*		}*/
/*	}*/
/*	*/
/*	*/
/*	for(i=2;i<=nxm;i++)	//Calculation of components of Surface Tension*/
/*	{*/

/*    		for(j=2;j<=nym;j++)*/
/*		{*/
/*	*/
/*	  		for(k=2;k<=nzm-1;k++)*/
/*			{*/
/*				ktop  	= k + 1;*/
/*				kbottom = k-1;*/
/*				dz = zf[k]-zf[kbottom];*/
/*				*/
/*				if(curvature==2)	zST[i][j][k] = -6.0*sigma*Kappa[i][j][k]*zGradAlphaS[i][j][k][l]*alphaS[i][j][k][l]*(1.0 - alphaS[i][j][k][l]);*/
/*				else if(curvature==1)	zST[i][j][k] = -sigma*Kappa[i][j][k]*(tAlpha[i][j][k] - bAlpha[i][j][k])/dz;*/
/*				*/
/*				else if(curvature==3)	zST[i][j][k] = -sigma*KappaHF[i][j][k]*(tAlpha[i][j][k] - bAlpha[i][j][k])/dz;*/

/*				zST[i][j][k] = (zSTt[i][j][k] + zSTt[i][j][kbottom])/2.0; */

/*			}*/
/*		}*/
/*	}		*/
	

	for(i=2;i<=nxm-1;i++)	//Calculation of components of Surface Tension
	{
		ieast 	= i + 1;
		iwest 	= i-1;
		dx = xf[i]-xf[iwest];
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;

    		for(j=2;j<=nym;j++)
		{
	  		for(k=2;k<=nzm;k++)
			{
				/* Arithmetic Interpolation of density*/
				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];
			
	
				/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
	
	
/*				xSTe[i][j][k] = (xST[i][j][k] + xST[ieast][j][k])/2.0;*/
				xSTe[i][j][k] = -sigma*eKappaHF[i][j][k]*(alphaS[ieast][j][k][l]-alphaS[i][j][k][l])/(dxpe); 
/*				xSTe[i][j][k] = -sigma*0.5*(alpha[ieast][j][k][l]-alpha[i][j][k][l])/(dxpe); */

			}
		}
	}
	
	for(i=2;i<=nxm;i++)	//Calculation of components of Surface Tension
	{

    		for(j=2;j<=nym-1;j++)
		{
			jnorth	= j + 1;
			jsouth = j-1;
			dy = yf[j]-yf[jsouth];
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
	
	  		for(k=2;k<=nzm;k++)
			{
				/* Arithmetic Interpolation of density*/
				rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];

				/* Harmonic Interpolation of density*/
/*				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/
			
	
/*				ySTn[i][j][k] = (yST[i][j][k] + yST[i][jnorth][k])/2.0;*/
				ySTn[i][j][k] = -sigma*nKappaHF[i][j][k]*(alphaS[i][jnorth][k][l] - alphaS[i][j][k][l])/(dypn); 
/*				ySTn[i][j][k] = -sigma*0.5*(alpha[i][jnorth][k][l] - alpha[i][j][k][l])/(dypn); */
			}
		}
	}
	
	
	for(i=2;i<=nxm;i++)	//Calculation of components of Surface Tension
	{

    		for(j=2;j<=nym;j++)
		{
	
	  		for(k=2;k<=nzm-1;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;
				dz = zf[k]-zf[kbottom];
				dzpt = zc[ktop]-zc[k];
				fzt = (zf[k] - zc[k])/dzpt;
				fzp = 1.0 - fzt; 
				
				/* Arithmetic Interpolation of density*/
				rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];
				
				/* Harmonic Interpolation of density*/
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/
				
	
/*				zSTt[i][j][k] = (zST[i][j][k] + zST[i][j][ktop])/2.0;*/
				zSTt[i][j][k] = -sigma*tKappaHF[i][j][k]*(alphaS[i][j][ktop][l] - alphaS[i][j][k][l])/(dzpt); 
/*				zSTt[i][j][k] = -sigma*0.5*(alpha[i][j][ktop][l] - alpha[i][j][k][l])/(dzpt); 		*/
			}
		}
	}
	
		


}
