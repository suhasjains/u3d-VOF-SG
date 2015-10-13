void TEqn()
{
	
	for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				ap[i][j][k]  = (- ae[i][j][k] - aw[i][j][k] - an[i][j][k] - as[i][j][k] - at[i][j][k] - ab[i][j][k] + apT[i][j][k])*urecrt;
      				sT[i][j][k]  = sT[i][j][k] + (1 - tUnderRelaxCoeff) * ap[i][j][k] * T[i][j][k][l];
      				
      				apT[i][j][k] = 1 / ap[i][j][k];
      				
			}
    		}
  	}


for(nonLinear=0;nonLinear<=nTemperatureLoops;nonLinear++)
{

	// Creating matrices for ADI-TDMA algorithm  
  	for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
	{
		for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
		{
			for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[i] 	= sT[i][j][k];

				c[i] = c[i] - an[i][j][k]*T[i][jnorth][k][l];
				c[i] = c[i] - as[i][j][k]*T[i][jsouth][k][l];
				c[i] = c[i] - at[i][j][k]*T[i][j][ktop][l];
				c[i] = c[i] - ab[i][j][k]*T[i][j][kbottom][l];
					
				
				c[i] = c[i]*apT[i][j][k];
				
				
				lower[i]= aw[i][j][k]*apT[i][j][k];
				upper[i]= ae[i][j][k]*apT[i][j][k];
				
				mThomas(wcID[e][f][g],ecID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=wcID[e][f][g];thomasI<=ecID[e][f][g];thomasI++)
	  				T[thomasI][j][k][l]=x[thomasI];				
			}
    	}
  }
  
  
  
  	TMaxRes = 0;
  	TTotalRes = 0;
  	for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
	{
		for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
		{
			for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				TRes[i][j][k] 	= sT[i][j][k] - (ap[i][j][k]*T[i][j][k][l]);
				
				TRes[i][j][k] = TRes[i][j][k] - an[i][j][k]*T[i][jnorth][k][l];
				TRes[i][j][k] = TRes[i][j][k] - as[i][j][k]*T[i][jsouth][k][l];
				TRes[i][j][k] = TRes[i][j][k] - at[i][j][k]*T[i][j][ktop][l];
				TRes[i][j][k] = TRes[i][j][k] - ab[i][j][k]*T[i][j][kbottom][l];
				TRes[i][j][k] = TRes[i][j][k] - aw[i][j][k]*T[iwest][j][k][l];
				TRes[i][j][k] = TRes[i][j][k] - ae[i][j][k]*T[ieast][j][k][l];
				
				
				//uTotalRes += fabs(uRes[i][j][k]);
				TTotalRes += fabs(ap[i][j][k]*T[i][j][k][l]);
				if(fabs(TRes[i][j][k])>TMaxRes)	TMaxRes	= fabs(TRes[i][j][k]);
			}
    		}
   	}
	
	TMeanRes = TTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
	if(TMeanRes==0)		TNormRes = 0;
    	else			TNormRes = TMaxRes/TMeanRes;
   
	
	//sync();
	
/*	printf(" T Res = %e\n",TNormRes);*/
	//printf(" u Max Res = %e\n",uTotalRes);
   	if(fabs(TNormRes)>TAccuracy)		nTemperatureLoops = nTemperatureLoops + 1;
   	
   	
}
/*		printf(" T Res = %e\n",TNormRes);*/
  
}
