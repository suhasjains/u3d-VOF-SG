void curvatureK8Convo()
{


	double eNorm[105][105][10],wNorm[105][105][10],nNorm[105][105][10],sNorm[105][105][10],tNorm[105][105][10],bNorm[105][105][10];
	
	double N,D;


	for(i=2;i<=nxm-1;i++)	//Calculation of curvature
  	{
		ieast 	= i + 1;
		iwest 	= i-1;
		ieasteast = i+2;
		dx = xf[i]-xf[iwest];

    		for(j=2;j<=nym;j++)
		{
			jnorth	= j + 1;
			jsouth = j-1;
			jnorthnorth = j+2;
			dy = yf[j]-yf[jsouth];
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;
				ktoptop = k+2;
				dz = zf[k]-zf[kbottom];
				
				
				//Calculation of normals along x
				N = -xGradAlphaS[ieast][j][k][l];
				D = pow(xGradAlphaS[ieast][j][k][l],2.0);
				D = D + pow(yGradAlphaS[ieast][j][k][l],2.0);
				D = D + pow(zGradAlphaS[ieast][j][k][l],2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				eNorm[i][j][k] = N/D;
				
				N = -xGradAlphaS[iwest][j][k][l];
				D = pow(xGradAlphaS[iwest][j][k][l],2.0);
				D = D + pow(yGradAlphaS[iwest][j][k][l],2.0);
				D = D + pow(zGradAlphaS[iwest][j][k][l],2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				wNorm[i][j][k] = N/D;
				
				
			}
		}
	}
	
	
	
	for(i=2;i<=nxm;i++)	//Calculation of curvature
  	{
		ieast 	= i + 1;
		iwest 	= i-1;
		ieasteast = i+2;
		dx = xf[i]-xf[iwest];

    		for(j=2;j<=nym-1;j++)
		{
			jnorth	= j + 1;
			jsouth = j-1;
			jnorthnorth = j+2;
			dy = yf[j]-yf[jsouth];
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;
				ktoptop = k+2;
				dz = zf[k]-zf[kbottom];
				
				//Calculation of normals along y
				N = -yGradAlphaS[i][jnorth][k][l];
				D = pow(xGradAlphaS[i][jnorth][k][l],2.0);
				D = D + pow(yGradAlphaS[i][jnorth][k][l],2.0);
				D = D + pow(zGradAlphaS[i][jnorth][k][l],2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				nNorm[i][j][k] = N/D;
				
				N = -yGradAlphaS[i][jsouth][k][l];
				D = pow(xGradAlphaS[i][jsouth][k][l],2.0);
				D = D + pow(yGradAlphaS[i][jsouth][k][l],2.0);
				D = D + pow(zGradAlphaS[i][jsouth][k][l],2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				sNorm[i][j][k] = N/D;
				
			}
		}
	}
	
	
	for(i=2;i<=nxm;i++)	//Calculation of curvature
  	{
		ieast 	= i + 1;
		iwest 	= i-1;
		ieasteast = i+2;
		dx = xf[i]-xf[iwest];

    		for(j=2;j<=nym;j++)
		{
			jnorth	= j + 1;
			jsouth = j-1;
			jnorthnorth = j+2;
			dy = yf[j]-yf[jsouth];
	
	  		for(k=2;k<=nzm-1;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;
				ktoptop = k+2;
				dz = zf[k]-zf[kbottom];
				
				//Calcualtion of normals along z
				N = -zGradAlphaS[i][j][ktop][l];
				D = pow(xGradAlphaS[i][j][ktop][l],2.0);
				D = D + pow(yGradAlphaS[i][j][ktop][l],2.0);
				D = D + pow(zGradAlphaS[i][j][ktop][l],2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				tNorm[i][j][k] = N/D;
				
				N = -zGradAlphaS[i][j][kbottom][l];
				D = pow(xGradAlphaS[i][j][kbottom][l],2.0);
				D = D + pow(yGradAlphaS[i][j][kbottom][l],2.0);
				D = D + pow(zGradAlphaS[i][j][kbottom][l],2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				bNorm[i][j][k] = N/D;
				
			}
		}
	}
	
	
	for(i=2;i<=nxm;i++)	//Calculation of curvature
  	{
		ieast 	= i + 1;
		iwest 	= i-1;
		ieasteast = i+2;
		dx = xf[i]-xf[iwest];

    		for(j=2;j<=nym;j++)
		{
			jnorth	= j + 1;
			jsouth = j-1;
			jnorthnorth = j+2;
			dy = yf[j]-yf[jsouth];
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;
				ktoptop = k+2;
				dz = zf[k]-zf[kbottom];
				
				
				//Curvature calculation
				Kappa[i][j][k] = -(((eNorm[i][j][k]-wNorm[i][j][k])/(2.0*dx)) + ((nNorm[i][j][k]-sNorm[i][j][k])/(2.0*dy)) + ((tNorm[i][j][k]-bNorm[i][j][k])/(2.0*dz)));
				
				//printf("Curvature = %e \n",Kappa[i][j][k]);
				
				
			}
		}
	}
	
	
	for(i=2;i<=nxm;i++)	//Calculation of curvature
  	{
  		ieast = i+1;
  		iwest = i-1;
		dx = xf[i]-xf[iwest];

    		for(j=2;j<=nym;j++)
		{
			jsouth = j-1;
			jnorth = j+1;
			dy = yf[j]-yf[jsouth];
	
	  		for(k=2;k<=nzm;k++)
			{
				kbottom = k-1;
				ktop = k+1;
				dz = zf[k]-zf[kbottom];
				
	
/*				if(KappaHF[i][j][k]!=0.0)*/
/*				{		*/
/*				 	if(KappaHF[ieast][j][k]==0.0)	eKappaHF[i][j][k] = KappaHF[i][j][k];*/
/*				  	else 				eKappaHF[i][j][k] = (KappaHF[i][j][k] + KappaHF[ieast][j][k])/2.0;*/
/*				  	*/
/*				  	if(KappaHF[iwest][j][k]==0.0)	eKappaHF[iwest][j][k] = KappaHF[i][j][k];*/
/*				  	else 				eKappaHF[iwest][j][k] = (KappaHF[iwest][j][k] + KappaHF[i][j][k])/2.0;*/
/*				  	*/
/*				*/
/*				*/
/*				 	if(KappaHF[i][jnorth][k]==0.0)	nKappaHF[i][j][k] = KappaHF[i][j][k];*/
/*					else 				nKappaHF[i][j][k] = (KappaHF[i][j][k] + KappaHF[i][jnorth][k])/2.0;*/
/*					*/
/*					if(KappaHF[i][jsouth][k]==0.0)	nKappaHF[i][jsouth][k] = KappaHF[i][j][k];*/
/*					else 				nKappaHF[i][jsouth][k] = (KappaHF[i][jsouth][k] + KappaHF[i][j][k])/2.0;*/
/*					*/
/*					*/
/*				*/
/*					if(KappaHF[i][j][ktop]==0.0)	tKappaHF[i][j][k] = KappaHF[i][j][k];*/
/*					else				tKappaHF[i][j][k] = (KappaHF[i][j][k] + KappaHF[i][j][ktop])/2.0;*/
/*					*/
/*					if(KappaHF[i][j][kbottom]==0.0)	tKappaHF[i][j][kbottom] = KappaHF[i][j][k];*/
/*					else				tKappaHF[i][j][kbottom] = (KappaHF[i][j][kbottom] + KappaHF[i][j][k])/2.0;*/
/*				}*/
/*				*/
				eKappaHF[i][j][k] = (Kappa[i][j][k] + Kappa[ieast][j][k])/2.0;
				nKappaHF[i][j][k] = (Kappa[i][j][k] + Kappa[i][jnorth][k])/2.0;
				tKappaHF[i][j][k] = (Kappa[i][j][k] + Kappa[i][j][ktop])/2.0;
				
				
			}
		}
	}
	
	
	
	

}
