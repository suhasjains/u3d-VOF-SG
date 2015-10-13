void correctAlpha()
{
	
	nAlphaCorrectorLoops = 1;

	printf("Starting VOF Correction Loops\n");

	for(say=0;say<=nAlphaCorrectorLoops;say++)
	{

		maxnAlphaLoops = 0;
		
		maxAlpha = 0;
		minAlpha = 0;

/*		alphaEqnCorrCoeff();*/
		
/*		alphaEqn();*/

		for (i=2;i<=nxm;i++)
	 	{
	  		for (j=2;j<=nym;j++)
		  	{
		  		for(k=2;k<=nzm;k++)
		  	 	{
					if(alpha[i][j][k][l]>maxAlpha)	maxAlpha = alpha[i][j][k][l];
					if(alpha[i][j][k][l]<minAlpha)	minAlpha = alpha[i][j][k][l];
				}
			}
		}

	


/*		if((minAlpha<-1e-5))	nAlphaCorrectorLoops = nAlphaCorrectorLoops + 1; */

		printf("maxAlpha = %e, minAlpha = %e \n",maxAlpha,minAlpha);

		//printf("maxnegError = %e \n",maxnegError);

	}
	
	
		//Correction step of alpha( Hirts way) 
		 for (j=2;j<=nym;j++)
		{
			for (k=2;k<=nzm;k++)
			{
				for(i=2;i<=nxm;i++)
				{
	
/*					if(alpha[i][j][k][l]>1.0)	alpha[i][j][k][l]=1.0;*/
/*					if(alpha[i][j][k][l]<0.0)	alpha[i][j][k][l]=0.0;*/
	
/*					if(alpha[i][j][k][l]>=(1.0-1e-6))	alpha[i][j][k][l]=1.0;*/
/*					if(alpha[i][j][k][l]<=1e-6)		alpha[i][j][k][l]=0.0;*/
				}
			}
		}
	
	
	
	
	

		printf("	Max no of Alpha Equation Loops = %d\n",maxnAlphaLoops);
		printf("	Global Alpha Residual = %e\n",alphaNormRes);
		printf("	No of alpha corrector loops = %d \n",nAlphaCorrectorLoops);
}
