//calculates the value of advected alpha
void alphaEqn()
{

	nAlphaLoops = 5;
	
	for(nonLinear=0;nonLinear<=nAlphaLoops;nonLinear++)
  {
	
	
 	for (j=2;j<=nym;j++)
	{
		for (k=2;k<=nzm;k++)
		{
			for(i=2;i<=nxm;i++)
			{
				
                                tauAlpha[i][j][k][l] = Sp[i][j][k]/ttaP[i][j][k];  
								
			}
    	        }
        }
  
  alphaMaxRes = 0;
  alphaTotalRes = 0; 
  for (j=2;j<=nym;j++)
	{
		for (k=2;k<=nzm;k++)
		{
			for(i=2;i<=nxm;i++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
	  			alphaRes[i][j][k] 	= Sp[i][j][k] - ttaP[i][j][k]*tauAlpha[i][j][k][l];
				alphaTotalRes += fabs(ttaP[i][j][k]*tauAlpha[i][j][k][l]);
				
				if(fabs(alphaRes[i][j][k])>alphaMaxRes)	alphaMaxRes = fabs(alphaRes[i][j][k]);
				

				//if(alpha[i][j][k][l]!=1.0&&alpha[i][j][k][l]!=1.0)		printf("alpha[%d][%d][%d] = %e \n",i,j,k,alpha[i][j][k][l]);

/*				if(alphaRes[20][20][2]!=0)	printf(" alphaRes[%d][%d][%d] = %e\n",20,20,2,alphaRes[20][20][2]);*/

	  			
			}
    		}
    	}
    	
/*	printf("alpha Max Res = %e \n",alphaMaxRes);*/
    	
    	alphaMeanRes = alphaTotalRes/((double)((nxm-1)*(nym-1)*(nzm-1)));
    	if(alphaMeanRes==0)		alphaNormRes = 0;
    	else			alphaNormRes = alphaMaxRes/alphaMeanRes;

	
    
    	if(fabs(alphaNormRes)>alphaAccuracy)		nAlphaLoops = nAlphaLoops + 1;
/*    	printf("	Global Alpha Residual = %e\n",alphaMeanRes);*/
  }

/*	printf("	Global Alpha Norm Residual = %e\n",alphaNormRes);*/
	//printf("	No of Alpha Equation Loops = %d\n",nAlphaLoops);

	if(nAlphaLoops>maxnAlphaLoops)	maxnAlphaLoops = nAlphaLoops;

}

