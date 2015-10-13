void alphaEqnCorrCoeff()
{

	double alphaDc,alphaAc,alphaD,alphaA;
	double negError,posError;
	double delAlpha;
	double Betac;
	
	
	double Ape,Apw,Apn,Aps,Apt,Apb;
	double Ae,Aw,An,As,At,Ab;

	double rhoe,rhon,rhot;



	//Constructing alpha fluxes for faces along x direction
	for  (i=2; i<=nxm-1; i++)
 	{
		ieast 	= i+1;
		iwest 	= i-1;
		ieasteast = i+2;
		dx = xf[i]-xf[iwest];
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;

	 	for (j=2; j<=nym; j++)
   		{
	  		for(k=2; k<=nzm; k++)
			{


				/* Arithmetic Interpolation of density*/
				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];
				
				
				/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
				
				//Calculation of present and previous time step Donor cell alpha
				if(Fe[i][j][k]>=0)
				{
					alphaDc = tauAlpha[i][j][k][l];
					alphaAc = tauAlpha[ieast][j][k][l];
					alphaD = tauAlpha[i][j][k][l-1];
					alphaA = tauAlpha[ieast][j][k][l-1];
				}
				else			
				{
					alphaDc = tauAlpha[ieast][j][k][l];
					alphaAc = tauAlpha[i][j][k][l];
					alphaD = tauAlpha[ieast][j][k][l-1];
					alphaA = tauAlpha[i][j][k][l-1];
				}
			
				//Calculation of Error in each cell 	
				negError = std::max(-alphaDc,0.0);
				posError = std::max(alphaDc-1.0,0.0);
				
/*				if(negError!=0)	printf("negError = %e \n",negError);*/
/*				if(posError!=0) printf("posError = %e \n",posError);*/
				
				//if(maxnegError<negError)	maxnegError = negError;

				
				delAlpha = (alphaA + alphaAc)/2.0 - (alphaD + alphaDc)/2.0; 
				
				
/*				printf("delAlpha=%e\n",delAlpha);*/

				//Calculation of error in Beta
				if(negError!=0.0)
				{
					if(delAlpha>negError)		Betac = std::min(negError*(2.0 + Ce[i][j][k] - 2.0*Ce[i][j][k]*eBeta[i][j][k])/(2.0*Ce[i][j][k]*(delAlpha-negError)),eBeta[i][j][k]);
					else if(delAlpha<=negError)	Betac = 0.0;
				}
			
				else if(posError!=0.0)
				{
					if(delAlpha<-posError)		Betac = std::min(posError*(2.0 + Ce[i][j][k] - 2.0*Ce[i][j][k]*eBeta[i][j][k])/(2.0*Ce[i][j][k]*(-delAlpha-posError)),eBeta[i][j][k]);
					else if(delAlpha>=-posError)	Betac = 0.0;

					else 				Betac = std::min(posError*(2.0 + Ce[i][j][k] - 2.0*Ce[i][j][k]*eBeta[i][j][k])/(2.0*Ce[i][j][k]*(-delAlpha-posError)),eBeta[i][j][k]);
				}

				else	continue;
				
				
				
/*				printf("Betac = %e delAlpha = %e posError = %e\n",Betac,delAlpha,posError);*/

				printf("Ce=%e\n",Ce[i][j][k]);

				//Calculation of new Beta 
				eBeta[i][j][k] = eBeta[i][j][k] - Betac;
				
				
				
				
				//Correction of Alpha at Cell face
				if(Fe[i][j][k]>=0)	eAlpha[i][j][k] = (1.0 - eBeta[i][j][k])*tauAlpha[i][j][k][l-1] + eBeta[i][j][k]*tauAlpha[ieast][j][k][l-1];
				else			eAlpha[i][j][k] = (1.0 - eBeta[i][j][k])*tauAlpha[ieast][j][k][l-1] + eBeta[i][j][k]*tauAlpha[i][j][k][l-1];
		

				wAlpha[ieast][j][k] = eAlpha[i][j][k];


				//Calculation of coefficients
				if(Fe[i][j][k]>=0)
				{
					Ape = 1.0 - eBeta[i][j][k]; 	  
					Apw = eBeta[i][j][k];
					Ae  = eBeta[i][j][k];
					Aw  = 1.0 - eBeta[i][j][k];
				}

				else
				{
					Ape = eBeta[i][j][k];
					Apw = 1.0 - eBeta[i][j][k];
					Ae  = 1.0 - eBeta[i][j][k];
					Aw  = eBeta[i][j][k];
				}				
				
				aE[i][j][k][l] 	= Ae*Fe[i][j][k]/(2*rhoe);
				aW[ieast][j][k][l] = -Aw*Fe[i][j][k]/(2*rhoe);

				aPE[i][j][k][l]	= Ape*Fe[i][j][k]/(2*rhoe);
				aPW[ieast][j][k][l]= -Apw*Fe[i][j][k]/(2*rhoe);

			}
		}
	}
			

	//Constructing alpha fluxes for faces along y direction
	for(j=2;j<=nym-1;j++)
	{
		jnorth = j+1;
		jsouth = j-1;
		jnorthnorth = j+2;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;	
    		fyp = 1.0 - fyn;

		for(i=2;i<=nxm;i++)
		{
			for(k=2;k<=nzm;k++)
			{

				/* Arithmetic Interpolation of density*/
				rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];


				/* Harmonic Interpolation of density*/
/*				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/
				
				//Calculation of present and previous Donor cell alpha
				if(Fn[i][j][k]>=0)
				{
					alphaDc = tauAlpha[i][j][k][l];
					alphaAc = tauAlpha[i][jnorth][k][l];
					alphaD = tauAlpha[i][j][k][l-1];
					alphaA = tauAlpha[i][jnorth][k][l-1];
				}
				else			
				{
					alphaDc = tauAlpha[i][jnorth][k][l];
					alphaAc = tauAlpha[i][j][k][l];
					alphaD = tauAlpha[i][jnorth][k][l-1];
					alphaA = tauAlpha[i][j][k][l-1];
				}


				//Calculation of Error in each cell 	
				negError = std::max(-alphaDc,0.0);
				posError = std::max(alphaDc-1.0,0.0);


				delAlpha = (alphaA + alphaAc)/2.0 - (alphaD + alphaDc)/2.0; 

				//Calculation of error in Beta
				if(negError!=0.0)
				{
					if(delAlpha>negError)		Betac = std::min(negError*(2.0 + Cn[i][j][k] - 2.0*Cn[i][j][k]*nBeta[i][j][k])/(2.0*Cn[i][j][k]*(delAlpha-negError)),nBeta[i][j][k]);
					else if(delAlpha<=negError)	Betac = 0.0;
				}
			
				else if(posError!=0.0)
				{
					if(delAlpha<-posError)		Betac = std::min(posError*(2.0 + Cn[i][j][k] - 2.0*Cn[i][j][k]*nBeta[i][j][k])/(2.0*Cn[i][j][k]*(-delAlpha-posError)),nBeta[i][j][k]);
					else if(delAlpha>=-posError)	Betac = 0.0;
					else 				Betac = std::min(posError*(2.0 + Cn[i][j][k] - 2.0*Cn[i][j][k]*nBeta[i][j][k])/(2.0*Cn[i][j][k]*(-delAlpha-posError)),nBeta[i][j][k]);
				}

				else	continue;

				//Calculation of new Beta 
				nBeta[i][j][k] = nBeta[i][j][k] - Betac;
				
/*				printf("Cn=%e\n",Cn[i][j][k]);*/
				
/*				printf("Betac = %e delAlpha = %e posError = %e\n",Betac,delAlpha,posError);*/
				
/*				printf("Betac = %e\n",Betac);*/
				
				//Correction of Alpha at Cell face
				if(Fn[i][j][k]>=0)	nAlpha[i][j][k] = (1.0 - nBeta[i][j][k])*tauAlpha[i][j][k][l-1] + nBeta[i][j][k]*tauAlpha[i][jnorth][k][l-1];
				else			nAlpha[i][j][k] = (1.0 - nBeta[i][j][k])*tauAlpha[i][jnorth][k][l-1] + nBeta[i][j][k]*tauAlpha[i][j][k][l-1];

				
				sAlpha[i][jnorth][k] = nAlpha[i][j][k];
				
				

							
				//Calculation of coefficients
				if(Fn[i][j][k]>=0)
				{
					Apn = 1.0 - nBeta[i][j][k]; 	  
					Aps = nBeta[i][j][k];
					An  = nBeta[i][j][k];
					As  = 1.0 - nBeta[i][j][k];
				}
					else
				{
					Apn = nBeta[i][j][k];
					Aps = 1.0 - nBeta[i][j][k];
					An  = 1.0 - nBeta[i][j][k];
					As  = nBeta[i][j][k];
				}				
				
				aN[i][j][k][l] 	= An*Fn[i][j][k]/(2*rhon);
				aS[i][jnorth][k][l] = -As*Fn[i][j][k]/(2*rhon);
				aPN[i][j][k][l]	= Apn*Fn[i][j][k]/(2*rhon);
				aPS[i][jnorth][k][l]= -Aps*Fn[i][j][k]/(2*rhon);
	
				//printf("Beta = %e \n",Beta);
				//printf("aN = %e  \n ",aN[i][j][k]);
				//printf("aPN = %e \n",aPN[i][j][k]);
			}	
		}	
	}	
	


	//Constructing alpha fluxes for faces along z direction
	for(k=2;k<=nzm-1;k++)
	{
		ktop = k+1;
		kbottom = k-1;
		ktoptop = k+2;
		dzpt = zc[ktop]-zc[k];
		fzt = (zf[k] - zc[k])/dzpt;
		fzp = 1.0 - fzt;

		for(i=2;i<=nxm;i++)
		{
			for(j=2;j<=nym;j++)
			{


				/* Arithmetic Interpolation of density*/
				rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];


				/* Harmonic Interpolation of density*/
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/

				//Calculation of present and previous Donor cell alpha
				if(Ft[i][j][k]>=0)
				{
					alphaDc = tauAlpha[i][j][k][l];
					alphaAc = tauAlpha[i][j][ktop][l];
					alphaD = tauAlpha[i][j][k][l-1];
					alphaA = tauAlpha[i][j][ktop][l-1];
				}
				else			
				{
					alphaDc = tauAlpha[i][j][ktop][l];
					alphaAc = tauAlpha[i][j][k][l];
					alphaD = tauAlpha[i][j][ktop][l-1];
					alphaA = tauAlpha[i][j][k][l-1];
				}


				//Calculation of Error in each cell 	
				negError = std::max(-alphaDc,0.0);
				posError = std::max(alphaDc-1.0,0.0);


				delAlpha = (alphaA + alphaAc)/2.0 - (alphaD + alphaDc)/2.0; 

				//Calculation of error in Beta
				if(negError!=0.0)
				{
					if(delAlpha>negError)		Betac = std::min(negError*(2.0 + Ct[i][j][k] - 2.0*Ct[i][j][k]*tBeta[i][j][k])/(2.0*Ct[i][j][k]*(delAlpha-negError)),tBeta[i][j][k]);
					else if(delAlpha<=negError)	Betac = 0.0;
				}
			
				else if(posError!=0.0)
				{
					if(delAlpha<-posError)		Betac = std::min(posError*(2.0 + Ct[i][j][k] - 2.0*Ct[i][j][k]*tBeta[i][j][k])/(2.0*Ct[i][j][k]*(-delAlpha-posError)),tBeta[i][j][k]);
					else if(delAlpha>=-posError)	Betac = 0.0;
					else 				Betac = std::min(posError*(2.0 + Ct[i][j][k] - 2.0*Ct[i][j][k]*tBeta[i][j][k])/(2.0*Ct[i][j][k]*(-delAlpha-posError)),tBeta[i][j][k]);
				}

				else	continue;

				//Calculation of new Beta 
				tBeta[i][j][k] = tBeta[i][j][k] - Betac;
				
				//Correction of Alpha at Cell face
				if(Ft[i][j][k]>=0)	tAlpha[i][j][k] = (1.0 - tBeta[i][j][k])*tauAlpha[i][j][k][l] + tBeta[i][j][k]*tauAlpha[i][j][ktop][l];
				else			tAlpha[i][j][k] = (1.0 - tBeta[i][j][k])*tauAlpha[i][j][ktop][l] + tBeta[i][j][k]*tauAlpha[i][j][k][l];

				bAlpha[i][j][ktop] = tAlpha[i][j][k];

				
				//Calculation of coefficients
				if(Ft[i][j][k]>=0)
				{
					Apt = 1.0 - tBeta[i][j][k]; 	  
					Apb = tBeta[i][j][k];
					At  = tBeta[i][j][k];
					Ab  = 1.0 - tBeta[i][j][k];
				}

				else
				{
					Apt = tBeta[i][j][k];
					Apb = 1.0 - tBeta[i][j][k];
					At  = 1.0 - tBeta[i][j][k];
					Ab  = tBeta[i][j][k];
				}				
					
				aT[i][j][k][l] 	= At*Ft[i][j][k]/(2*rhot);
				aB[i][j][ktop][l] = -Ab*Ft[i][j][k]/(2*rhot);
	
				aPT[i][j][k][l]	= Apt*Ft[i][j][k]/(2*rhot);
				aPB[i][j][ktop][l]= -Apb*Ft[i][j][k]/(2*rhot);
	
				//printf("aT = %e  \n ",aT[i][j][k]);
				//printf("aPT = %e \n",aPT[i][j][k]);
					
	
			}	
		}
	}
	
for  (i=2; i<=nxm; i++)
 	{
	 	for (j=2; j<=nym; j++)
   		{
	  		for(k=2; k<=nzm; k++)
			{
				
				ieast = i+1;
				iwest = i-1;
				jnorth = j+1;
				jsouth = j-1;
				ktop = k+1;
				kbottom = k-1;
		
				Ap[i][j][k][l] = aPE[i][j][k][l] + aPW[i][j][k][l] + aPN[i][j][k][l] + aPS[i][j][k][l] + aPT[i][j][k][l] + aPB[i][j][k][l]; 
			
				taP[i][j][k] = ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dt - ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dtau + Ap[i][j][k][l];
				naP[i][j][k] = Ap[i][j][k][l-1] - ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dt;
				ttaP[i][j][k] = ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dtau;

				Sp[i][j][k] = - taP[i][j][k]*tauAlpha[i][j][k][l-1] - naP[i][j][k]*alpha[i][j][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aE[i][j][k][l-1]*alpha[ieast][j][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aW[i][j][k][l-1]*alpha[iwest][j][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aN[i][j][k][l-1]*alpha[i][jnorth][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aS[i][j][k][l-1]*alpha[i][jsouth][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aT[i][j][k][l-1]*alpha[i][j][ktop][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aB[i][j][k][l-1]*alpha[i][j][kbottom][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aE[i][j][k][l]*tauAlpha[ieast][j][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aW[i][j][k][l]*tauAlpha[iwest][j][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aN[i][j][k][l]*tauAlpha[i][jnorth][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aS[i][j][k][l]*tauAlpha[i][jsouth][k][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aT[i][j][k][l]*tauAlpha[i][j][ktop][l-1];
				Sp[i][j][k] = Sp[i][j][k] - aB[i][j][k][l]*tauAlpha[i][j][kbottom][l-1];

/*				printf("Ap = %e , Sp = %e, aP = %e \n",Ap[i][j][k],Sp[i][j][k],aP[i][j][k]);*/

			}
		}
	}


			
}
			
				 



				
