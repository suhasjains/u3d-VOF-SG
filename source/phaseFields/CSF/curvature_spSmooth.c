void curvatureSmooth()
{
	
	
	double N,D;
	double neAlpha,seAlpha,teAlpha,beAlpha;
	double nwAlpha,swAlpha,twAlpha,bwAlpha;
	double ntAlpha,nbAlpha;
	double stAlpha,sbAlpha;
	
	
	for(i=2;i<=nxm-1;i++)	//Calculation of curvature
  	{
		ieast 	= i + 1;
		iwest 	= i-1;
		ieasteast = i+2;
		dx = xf[i]-xf[iwest];
		dxpe = xc[ieast] - xc[i];
		dxpw = xc[i] - xc[iwest];

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
				
				
				//Interpoalting alpha at unknowm points
				neAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][jnorth][k][l]  + alphaS[ieast][jnorth][k][l]);
				seAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][jsouth][k][l]  + alphaS[ieast][jsouth][k][l]);
				teAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][j][ktop][l]    + alphaS[ieast][j][ktop][l]);
				beAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][j][kbottom][l] + alphaS[ieast][j][kbottom][l]);
				nwAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][jnorth][k][l]  + alphaS[iwest][jnorth][k][l]);
				swAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][jsouth][k][l]  + alphaS[iwest][jsouth][k][l]);
				twAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][j][ktop][l]    + alphaS[iwest][j][ktop][l]);
				bwAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][j][kbottom][l] + alphaS[iwest][j][kbottom][l]);
				ntAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][ktop][l]    + alphaS[i][jnorth][k][l]  + alphaS[i][jnorth][ktop][l]);
				nbAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][kbottom][l] + alphaS[i][jnorth][k][l]  + alphaS[i][jnorth][kbottom][l]);
				stAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][ktop][l]    + alphaS[i][jsouth][k][l]  + alphaS[i][jsouth][ktop][l]);
				sbAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][kbottom][l] + alphaS[i][jsouth][k][l]  + alphaS[i][jsouth][kbottom][l]);
				
				
				//printf("j = %d \n",j);
				
				
				//Calculation of normals along x
				N = -(alpha[ieast][j][k][l] - alpha[i][j][k][l])/dxpe;
				D = pow((alpha[ieast][j][k][l] - alpha[i][j][k][l])/dxpe,2.0);
				D = D + pow((neAlpha-seAlpha)/dy,2.0);
				D = D + pow((teAlpha-beAlpha)/dz,2.0);
				D = pow(D,0.5);
				if(D==0.0)	eNorm[i][j][k] = 0.0;
				else            eNorm[i][j][k] = N/D;
				
				N = -(alpha[i][j][k][l] - alpha[iwest][j][k][l])/dxpw;
				D = pow((alpha[i][j][k][l] - alpha[iwest][j][k][l])/dxpw,2.0);
				D = D + pow((nwAlpha-swAlpha)/dy,2.0);
				D = D + pow((twAlpha-bwAlpha)/dz,2.0);
				D = pow(D,0.5);
				if(D==0.0)	wNorm[i][j][k] = 0.0;
				else            wNorm[i][j][k] = N/D;
				
				
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
			dypn = yc[jnorth] - yc[j];
			dyps = yc[j] - yc[jsouth];
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;
				ktoptop = k+2;
				dz = zf[k]-zf[kbottom];
				
				//Interpoalting alpha at unknowm points
				neAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][jnorth][k][l]  + alphaS[ieast][jnorth][k][l]);
				seAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][jsouth][k][l]  + alphaS[ieast][jsouth][k][l]);
				teAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][j][ktop][l]    + alphaS[ieast][j][ktop][l]);
				beAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][j][kbottom][l] + alphaS[ieast][j][kbottom][l]);
				nwAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][jnorth][k][l]  + alphaS[iwest][jnorth][k][l]);
				swAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][jsouth][k][l]  + alphaS[iwest][jsouth][k][l]);
				twAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][j][ktop][l]    + alphaS[iwest][j][ktop][l]);
				bwAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][j][kbottom][l] + alphaS[iwest][j][kbottom][l]);
				ntAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][ktop][l]    + alphaS[i][jnorth][k][l]  + alphaS[i][jnorth][ktop][l]);
				nbAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][kbottom][l] + alphaS[i][jnorth][k][l]  + alphaS[i][jnorth][kbottom][l]);
				stAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][ktop][l]    + alphaS[i][jsouth][k][l]  + alphaS[i][jsouth][ktop][l]);
				sbAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][kbottom][l] + alphaS[i][jsouth][k][l]  + alphaS[i][jsouth][kbottom][l]);
				
				
				
				
				//Calculation of normals along y
				N = -(alpha[i][jnorth][k][l] - alpha[i][j][k][l])/dypn;
				D = pow((alpha[i][jnorth][k][l] - alpha[i][j][k][l])/dypn,2.0);
				D = D + pow((neAlpha-nwAlpha)/dx,2.0);
				D = D + pow((ntAlpha-nbAlpha)/dz,2.0);
				D = pow(D,0.5);
				if(D==0.0)	nNorm[i][j][k] = 0.0;
				else            nNorm[i][j][k] = N/D;
				
				N = -(alpha[i][j][k][l] - alpha[i][jsouth][k][l])/dyps;				
				D = pow((alpha[i][j][k][l] - alpha[i][jsouth][k][l])/dyps,2.0);
				D = D + pow((seAlpha-swAlpha)/dx,2.0);
				D = D + pow((stAlpha-sbAlpha)/dz,2.0);
				D = pow(D,0.5);
				if(D==0.0)	sNorm[i][j][k] = 0.0;
				else            sNorm[i][j][k] = N/D;
				
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
				dzpt = zc[ktop]-zc[k];
				dzpb = zc[k]-zc[kbottom];
				
				//Interpoalting alpha at unknowm points
				neAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][jnorth][k][l]  + alphaS[ieast][jnorth][k][l]);
				seAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][jsouth][k][l]  + alphaS[ieast][jsouth][k][l]);
				teAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][j][ktop][l]    + alphaS[ieast][j][ktop][l]);
				beAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[ieast][j][k][l]   + alphaS[i][j][kbottom][l] + alphaS[ieast][j][kbottom][l]);
				nwAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][jnorth][k][l]  + alphaS[iwest][jnorth][k][l]);
				swAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][jsouth][k][l]  + alphaS[iwest][jsouth][k][l]);
				twAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][j][ktop][l]    + alphaS[iwest][j][ktop][l]);
				bwAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[iwest][j][k][l]   + alphaS[i][j][kbottom][l] + alphaS[iwest][j][kbottom][l]);
				ntAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][ktop][l]    + alphaS[i][jnorth][k][l]  + alphaS[i][jnorth][ktop][l]);
				nbAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][kbottom][l] + alphaS[i][jnorth][k][l]  + alphaS[i][jnorth][kbottom][l]);
				stAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][ktop][l]    + alphaS[i][jsouth][k][l]  + alphaS[i][jsouth][ktop][l]);
				sbAlpha = 0.25*(alphaS[i][j][k][l] + alphaS[i][j][kbottom][l] + alphaS[i][jsouth][k][l]  + alphaS[i][jsouth][kbottom][l]);
				
				
				//Calcualtion of normals along z
				N = -(alpha[i][j][ktop][l] - alpha[i][j][k][l])/dzpt;
				D = pow((alpha[i][j][ktop][l] - alpha[i][j][k][l])/dzpt,2.0);
				D = D + pow((teAlpha-twAlpha)/dx,2.0);
				D = D + pow((ntAlpha-stAlpha)/dy,2.0);
				D = pow(D,0.5);
				if(D==0.0)	tNorm[i][j][k] = 0.0;
				else            tNorm[i][j][k] = N/D;
				
				N = -(alpha[i][j][k][l] - alpha[i][j][kbottom][l])/dzpb;
				D = pow((alpha[i][j][k][l] - alpha[i][j][kbottom][l])/dzpb,2.0);
				D = D + pow((beAlpha-bwAlpha)/dx,2.0);
				D = D + pow((nbAlpha-sbAlpha)/dy,2.0);
				D = pow(D,0.5);
				if(D==0.0)	bNorm[i][j][k] = 0.0;
				else            bNorm[i][j][k] = N/D;
				
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
				Kappa[i][j][k] = -(((eNorm[i][j][k]-wNorm[i][j][k])/dx) + ((nNorm[i][j][k]-sNorm[i][j][k])/dy) + ((tNorm[i][j][k]-bNorm[i][j][k])/dz));
				
/*				if(Kappa[i][j][k]!=0.0)	printf("Curvature = %e \n",Kappa[i][j][k]);*/
				
				
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
				
				
