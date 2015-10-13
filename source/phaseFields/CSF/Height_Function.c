void heightFunction()
{


	double xNorm[105][105][10],yNorm[105][105][10],zNorm[105][105][10];
	double HF[105][105][10];
	double N,D;
	
	int iHF,jHF,kHF;
	
	double signN;
	double zHF,zzHF,yHF,yyHF,xHF,xxHF,xyHF,xzHF,yzHF;

	for(i=2;i<=nxm-1;i++)	//Calculation of curvature
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
				
				
				//Calculation of normals along x
				N = -(alpha[ieast][j][k][l] - alpha[iwest][j][k][l])/(2.0*dx);
				D = 	pow((alpha[ieast][j][k][l] - alpha[iwest][j][k][l])/(2.0*dx),2.0);
				D = D + pow((alpha[i][jnorth][k][l] - alpha[i][jsouth][k][l])/(2.0*dy),2.0);
				D = D + pow((alpha[i][j][ktop][l] - alpha[i][j][kbottom][l])/(2.0*dz),2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				xNorm[i][j][k] = N/D;
				
				
/*				printf("xN = %e \n",xNorm[i][j][k]);*/
			}
		}
	}
	
	
	
	for(i=2;i<=nxm;i++)	//Calculation of curvature
  	{
  		ieast = i+1;
  		iwest = i-1;
		dx = xf[i]-xf[iwest];

    		for(j=2;j<=nym-1;j++)
		{
			jsouth = j-1;
			jnorth = j+1;
			dy = yf[j]-yf[jsouth];
	
	  		for(k=2;k<=nzm;k++)
			{
				kbottom = k-1;
				ktop = k+1;
				dz = zf[k]-zf[kbottom];
				
				//Calculation of normals along y
				N = -(alpha[i][jnorth][k][l] - alpha[i][jsouth][k][l])/(2.0*dy);
				D = 	pow((alpha[ieast][j][k][l] - alpha[iwest][j][k][l])/(2.0*dx),2.0);
				D = D + pow((alpha[i][jnorth][k][l] - alpha[i][jsouth][k][l])/(2.0*dy),2.0);
				D = D + pow((alpha[i][j][ktop][l] - alpha[i][j][kbottom][l])/(2.0*dz),2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				yNorm[i][j][k] = N/D;
				
				
				//printf("N = %e \n",N);
/*				printf("D = %e \n",D);*/
				
				
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
	
	  		for(k=2;k<=nzm-1;k++)
			{
				kbottom = k-1;
				ktop = k+1;
				dz = zf[k]-zf[kbottom];
				
				//Calcualtion of normals along z
				N = -(alpha[i][j][ktop][l] - alpha[i][j][kbottom][l])/(2.0*dz);
				D = 	pow((alpha[ieast][j][k][l] - alpha[iwest][j][k][l])/(2.0*dx),2.0);
				D = D + pow((alpha[i][jnorth][k][l] - alpha[i][jsouth][k][l])/(2.0*dy),2.0);
				D = D + pow((alpha[i][j][ktop][l] - alpha[i][j][kbottom][l])/(2.0*dz),2.0);
				D = pow(D,0.5);
				if(D==0.0)	D = 1e-100;
				zNorm[i][j][k] = N/D;
				
				
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
				
				if(fabs(xNorm[i][j][k])>=std::max(fabs(yNorm[i][j][k]),fabs(zNorm[i][j][k])))	
				{
				
					
					for(jHF=jsouth;jHF<=jnorth;jHF++)		//Defining Kernels
					{
						for(kHF=kbottom;kHF<=ktop;kHF++)
						{
						
							HF[i][jHF][kHF] = 0.0;		//Initializing HF to zero
						
							for(iHF=i-3;iHF<=i+3;iHF++)
							{
							
								//Condition for going out of boundary
								if(iHF<1||iHF>nx||jHF<1||jHF>ny||kHF<1||kHF>nz)	continue;
							
								//Integration of color functions along the length to find Height function
								HF[i][jHF][kHF] = HF[i][jHF][kHF] + alpha[iHF][jHF][kHF][l]*dx;
								
														
							}
						}
					}
					
					yHF  = (HF[i][jnorth][k] - HF[i][jsouth][k])/(2.0*dy);	
					yyHF = (HF[i][jnorth][k] + HF[i][jsouth][k] - 2.0*HF[i][j][k])/pow(dy,2.0);
					yzHF = (HF[i][jnorth][ktop] - HF[i][jnorth][kbottom] - HF[i][jsouth][ktop] + HF[i][jsouth][kbottom])/(4.0*dy*dz);
					zHF  = (HF[i][j][ktop] - HF[i][j][kbottom])/(2.0*dz);  
					zzHF = (HF[i][j][ktop] + HF[i][j][kbottom] - 2.0*HF[i][j][k])/pow(dz,2.0);
					
/*					N = (alpha[ieast][j][k][l] - alpha[iwest][j][k][l]);*/
/*					D = fabs((alpha[ieast][j][k][l] - alpha[iwest][j][k][l])); */
/*					*/
/*					if(D==0.0)	signN = 1.0;	*/
/*					else 		signN = N/D;*/
/*					*/

/*					if(xNorm[i][j][k]==0)	signN=1.0;*/
/*					else 			signN = xNorm[i][j][k]/fabs(xNorm[i][j][k]);*/

					signN = 1.0;
					
					N = signN * (yyHF + zzHF + yyHF*pow(zHF,2.0) + zzHF*pow(yHF,2.0) - 2.0*yzHF*zHF*yHF);
					D = pow((1.0 + pow(yHF,2.0) + pow(zHF,2.0)),1.5); 

					KappaHF[i][j][k] = N/D;

					if(zNumberOfCells>1)
					{
						if(HF[i][j][k]<=4.0*dx&&HF[i][j][k]>=3.0*dx)	KappaHF[i][j][k] = N/D;
						else 						KappaHF[i][j][k] = 0.0;
					}	
				}


				else if(fabs(yNorm[i][j][k])>=std::max(fabs(xNorm[i][j][k]),fabs(zNorm[i][j][k])))
				{
				
					for(iHF=iwest;iHF<=ieast;iHF++)		//Defining Kernels
					{
						for(kHF=kbottom;kHF<=ktop;kHF++)
						{
						
							HF[iHF][j][kHF] = 0.0;		//Initializing HF to zero
						
							for(jHF=j-3;jHF<=j+3;jHF++)
							{
							
								//Condition for going out of boundary
								if(iHF<1||iHF>nx||jHF<1||jHF>ny||kHF<1||kHF>nz)	continue;
							
								//Integration of color functions along the length to find Height function
								HF[iHF][j][kHF] = HF[iHF][j][kHF] + alpha[iHF][jHF][kHF][l]*dy;						
							}
						}
					}
					
					//if(HF[i][j][k]<=4.0*dy&&HF[i][j][k]>=3.0*dy) printf("HF = %e \n",HF[i][j][k]/dy);
					
					xHF  = (HF[ieast][j][k] - HF[iwest][j][k])/(2.0*dx);	 
					xxHF = (HF[ieast][j][k] + HF[iwest][j][k] - 2.0*HF[i][j][k])/pow(dx,2.0);
					xzHF = (HF[ieast][j][ktop] - HF[ieast][j][kbottom] - HF[iwest][j][ktop] + HF[iwest][j][kbottom])/(4.0*dx*dz);
					zHF  = (HF[i][j][ktop] - HF[i][j][kbottom])/(2.0*dz);  
					zzHF = (HF[i][j][ktop] + HF[i][j][kbottom] - 2.0*HF[i][j][k])/pow(dz,2.0);
					
/*					N = (alpha[i][jnorth][k][l] - alpha[i][jsouth][k][l]);*/
/*					D = fabs((alpha[i][jnorth][k][l] - alpha[i][jsouth][k][l])); */
/*					*/
/*					if(D==0.0)	signN = 1.0;	*/
/*					else 		signN = N/D;*/


/*					if(yNorm[i][j][k]==0)	signN=1.0;*/
/*					else 			signN = yNorm[i][j][k]/fabs(yNorm[i][j][k]);*/


					signN = 1.0;

					
					N = signN * (xxHF + zzHF + xxHF*pow(zHF,2.0) + zzHF*pow(xHF,2.0) - 2.0*xzHF*zHF*xHF);
					D = pow((1.0 + pow(xHF,2.0) + pow(zHF,2.0)),1.5);
					

					KappaHF[i][j][k] = N/D;

					if(zNumberOfCells>1)
					{
						if(HF[i][j][k]<=4.0*dy&&HF[i][j][k]>=3.0*dy)	KappaHF[i][j][k] = N/D;
						else 						KappaHF[i][j][k] = 0.0; 
					}
				
				
				
				}
				
				else if(fabs(zNorm[i][j][k])>=std::max(fabs(yNorm[i][j][k]),fabs(xNorm[i][j][k])))
				{
				
					for(iHF=iwest;iHF<=ieast;iHF++)		//Defining Kernels
					{
						for(jHF=jsouth;jHF<=jnorth;jHF++)
						{
						
							HF[iHF][jHF][k] = 0.0;		//Initializing HF to zero
						
							for(kHF=k-3;kHF<=k+3;kHF++)
							{
							
								//Condition for going out of boundary
								if(iHF<1||iHF>nx||jHF<1||jHF>ny||kHF<1||kHF>nz)	continue;
							
								//Integration of color functions along the length to find Height function
								HF[iHF][jHF][k] = HF[iHF][jHF][k] + alpha[iHF][jHF][kHF][l]*dz;						
							}
						}
					}
					
					
/*					printf("HF = %e \n",HF[i][j][k]);*/
					
					xHF  = (HF[ieast][j][k] - HF[iwest][j][k])/(2.0*dx);	
					xxHF = (HF[ieast][j][k] + HF[iwest][j][k] - 2.0*HF[i][j][k])/pow(dx,2.0);
					xyHF = (HF[ieast][jnorth][k] - HF[ieast][jsouth][k] - HF[iwest][jnorth][k] + HF[iwest][jsouth][k])/(4.0*dx*dy);
					yHF  = (HF[i][jnorth][k] - HF[i][jsouth][k])/(2.0*dy);	
					yyHF = (HF[i][jnorth][k] + HF[i][jsouth][k] - 2.0*HF[i][j][k])/pow(dy,2.0);
					
/*					N = (alpha[i][j][ktop][l] - alpha[i][j][kbottom][l]);*/
/*					D = fabs((alpha[i][j][ktop][l] - alpha[i][j][kbottom][l])); */
/*					*/
/*					if(D==0.0)	signN = 1.0;	*/
/*					else 		signN = N/D;*/

/*					if(zNorm[i][j][k]==0)	signN=1.0;*/
/*					else 			signN = zNorm[i][j][k]/fabs(zNorm[i][j][k]);*/


					signN = 1.0;

					
					N = signN * (xxHF + yyHF + xxHF*pow(yHF,2.0) + yyHF*pow(xHF,2.0) - 2.0*xyHF*yHF*xHF);
					D = pow((1.0 + pow(xHF,2.0) + pow(yHF,2.0)),1.5); 
					

					KappaHF[i][j][k] = N/D;

					if(zNumberOfCells>1)
					{	
						if(HF[i][j][k]<=4.0*dz&&HF[i][j][k]>=3.0*dz)	KappaHF[i][j][k] = N/D;
						else 						KappaHF[i][j][k] = 0.0;
					}
				
				}
				
				
				
/*				if(KappaHF[i][j][k]!=0.0)	printf("K = %e  at %d %d %d\n",KappaHF[i][j][k],i,j,k);*/
/*				printf("N = %e \n",N);*/
/*				printf("D = %e \n",D);*/
				
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
				eKappaHF[i][j][k] = (KappaHF[i][j][k] + KappaHF[ieast][j][k])/2.0;
				nKappaHF[i][j][k] = (KappaHF[i][j][k] + KappaHF[i][jnorth][k])/2.0;
				tKappaHF[i][j][k] = (KappaHF[i][j][k] + KappaHF[i][j][ktop])/2.0;
				
				
/*				printf("Kappa = %e\n",eKappaHF[i][j][k] );*/
				
				
			}
		}
	}
	
	
	
	
}
