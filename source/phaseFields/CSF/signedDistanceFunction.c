void signDisFun()
{

	double cond;
	int iKer,jKer,kKer;
	double radius;
	double delX,delY,delZ;

	
	for(i=2;i<=nxm;i++)	
  	{

    		for(j=2;j<=nym;j++)
		{
	
	  		for(k=2;k<=nzm;k++)
			{

				//Initializing a large value for sdf outside the kernel
				//sdf[i][j][k] = std::max(std::max(xDomainLength,yDomainLength),zDomainLength); 
				sdf[i][j][k] = 10.0;
				
				for(iKer=i-3;iKer<=i+3;iKer++)	
  				{

			    		for(jKer=j-3;jKer<=j+3;jKer++)
					{
	
				  		for(kKer=k-3;kKer<=k+3;kKer++)
						{

							//Condition for going out of boundary
							if(iKer<1||iKer>nx||jKer<1||jKer>ny||kKer<1||kKer>nz)	continue;
							
							
							//Condition for presence of interface							
							cond = pow(xGradAlphaS[iKer][jKer][kKer][l],2.0);
							cond = cond + pow(yGradAlphaS[iKer][jKer][kKer][l],2.0);
							cond = cond + pow(zGradAlphaS[iKer][jKer][kKer][l],2.0);
							cond = pow(cond,0.5);
							cond = 6.0*alphaS[iKer][jKer][kKer][l]*(1.0 - alphaS[iKer][jKer][kKer][l]);
							
							
							//Calculation of the distance between the cells
							delX = xc[i] - xc[iKer];
							delY = yc[j] - yc[jKer];
							delZ = zc[k] - zc[kKer];
							radius = pow((pow(delX,2.0) + pow(delY,2.0) + pow(delZ,2.0)),0.5); 

							
							if(cond!=0.0)
							{
						
						
						//Assigning the minimum distance as the signed distance function							
								if(fabs(sdf[i][j][k])>radius)
								{
								
									if(alphaS[i][j][k][l]>=0.5)	sdf[i][j][k] = radius;
									else  				sdf[i][j][k] = -radius;
									
								}
							}
							
						}
					}
				}
				
				if(sdf[i][j][k]==10.0)	sdf[i][j][k] = 0;
				
				
				//if(sdf[i][j][k]!=1)	printf("sdf = %e \n",sdf[i][j][k]);
							
			}
		}
	}
}							
							









