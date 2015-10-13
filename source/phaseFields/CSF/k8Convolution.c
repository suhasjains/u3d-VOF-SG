void k8Convo()
{


	int iConv,jConv,kConv;
	double delX,delY,delZ;
	double nF,nCF;
	double Eps;
	double r,A;
	
	double r1;
	double xGradNCF,yGradNCF,zGradNCF;

	for(i=2;i<=nxm;i++)	
  	{
		ieast 	= i + 1;
		iwest 	= i-1;
		dx = xf[i]-xf[iwest];

    		for(j=2;j<=nym;j++)
		{
			jnorth	= j + 1;
			jsouth = j-1;
			dy = yf[j]-yf[jsouth];
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;
				dz = zf[k]-zf[kbottom];
				
				
				
				//Calculation of Normalization factor and smooth alpha Field
				nF = 0.0;
				nCF = 0.0;
				
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				Eps = 2.0*dx;	//Considering equal sides
/*				Eps = 3.0*dx;	//Bigger Kernel*/
				
				for(iConv=i-1;iConv<=i+1;iConv++)	
  				{

			    		for(jConv=j-1;jConv<=j+1;jConv++)
					{
	
				  		for(kConv=k-1;kConv<=k+1;kConv++)
						{
						
							//Condition for going out of boundary
							if(iConv<1||iConv>nx||jConv<1||jConv>ny||kConv<1||kConv>nz)	continue;
						
							
							delX = xc[i] - xc[iConv];
							delY = yc[j] - yc[jConv];
							delZ = zc[k] - zc[kConv];
							
							
							r = vol * pow((1.0 - (pow(delX,2.0) + pow(delY,2.0) + pow(delZ,2.0))/pow(Eps,2.0)),4.0);
							
							//Summation to calculate Normalization factor
							nF = nF + r; 		
							
							//Summation to calculate smooth alpha Field
							nCF = nCF + r * alpha[iConv][jConv][kConv][l]; 	
							
							r1 = alpha[iConv][jConv][kConv][l] * vol * pow((1.0 - (pow(delX,2.0) + pow(delY,2.0) + pow(delZ,2.0))/pow(Eps,2.0)),3.0);
							
							
							//Summation to calculate smoothed alpha Field
							xGradNCF = xGradNCF + r1 * delX/pow(Eps,2.0); 
							yGradNCF = yGradNCF + r1 * delY/pow(Eps,2.0);
							zGradNCF = zGradNCF + r1 * delZ/pow(Eps,2.0);
							
							//printf("r = %e at %d,%d,%d \n",r,iConv,jConv,kConv);
 							 
						}
					}
				}
				
				//Normalization factor			
				A = 1.0/nF;		
				
				
				//Smoothed alpha Field
				alphaS[i][j][k][l] = nCF * A;
				
				//Gradients of smoothed alpha Field at cell centers
				xGradAlphaS[i][j][k][l] = -8.0 * A * xGradNCF;
				yGradAlphaS[i][j][k][l] = -8.0 * A * yGradNCF;
				zGradAlphaS[i][j][k][l] = -8.0 * A * zGradNCF;
				
				//printf("alphaS = %e \n",alphaS[i][j][k][l]);
				
				
			}
		}
	}
	
	
	curvatureK8Convo();
}
				
				
				
				
				
				
				
								
				
				
				
				
				
					
											



