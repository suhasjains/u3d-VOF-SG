#include"velocity.h"

//calculates the value of initial w velocity
void wEqn()
{
	for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
		 		at[i][j][k]  = (att[i][j][k] + ab[i][j][k] + atE[i][j][k] + atW[i][j][k] + aNt[i][j][k] + aSt[i][j][k] + apw[i][j][k])*urecrw;
      				sw[i][j][k]  = sw[i][j][k] + (1 - wUnderRelaxCoeff) * at[i][j][k] * w[i][j][k][l] ;
      				apw[i][j][k] = 1 / at[i][j][k];	/* For implementation in the pressure correction approach */
		 		
		 	}
 	     }
  	 }
  
   for(nonLinear=0;nonLinear<=wVelocityLoops;nonLinear++)
  {
  
//  Creating matrices for ADI-TDMA algorithm	 
  	 for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[k] 	= sw[i][j][k];
				
				c[k] = c[k] + atW[i][j][k]*w[iwest][j][k][l];
				c[k] = c[k] + atE[i][j][k]*w[ieast][j][k][l];
				c[k] = c[k] + aNt[i][j][k]*w[i][jnorth][k][l];
				c[k] = c[k] + aSt[i][j][k]*w[i][jsouth][k][l];
				
				if(k==2)	c[k] = c[k] + ab[i][j][k]*w[i][j][kbottom][l];
				if(k==nzm)	c[k] = c[k] + att[i][j][k]*w[i][j][ktop][l];
				
				c[k] = c[k]*apw[i][j][k];
				
				
				
				lower[k]= -ab[i][j][k]*apw[i][j][k];
				upper[k]= -att[i][j][k]*apw[i][j][k];
			}	
				mThomas(bcID[e][f][g],tcID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=bcID[e][f][g];thomasI<=tcID[e][f][g];thomasI++)
	  				w[i][j][thomasI][l]=x[thomasI];				
			
    	}
  }
  
  wMaxRes = 0;
  wTotalRes = 0;
  for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
	  			wRes[i][j][k] 	= sw[i][j][k] - at[i][j][k]*w[i][j][k][l];
	  			
				wRes[i][j][k] = wRes[i][j][k] + aNt[i][j][k]*w[i][jnorth][k][l];
				wRes[i][j][k] = wRes[i][j][k] + aSt[i][j][k]*w[i][jsouth][k][l];
				wRes[i][j][k] = wRes[i][j][k] + att[i][j][k]*w[i][j][ktop][l];
				wRes[i][j][k] = wRes[i][j][k] + ab[i][j][k]*w[i][j][kbottom][l];
				wRes[i][j][k] = wRes[i][j][k] + atW[i][j][k]*w[iwest][j][k][l];
				wRes[i][j][k] = wRes[i][j][k] + atE[i][j][k]*w[ieast][j][k][l];
	  			
	  			
				wTotalRes += fabs(at[i][j][k]*w[i][j][k][l]);
				if(fabs(wRes[i][j][k])>wMaxRes)	wMaxRes	= fabs(wRes[i][j][k]);
	  			
			}
    	}
   }
   
   	wMeanRes = wTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
    	if(wMeanRes==0)		wNormRes = 0;
    	else			wNormRes = wMaxRes/wMeanRes;
   
   
   	velocityBoundaryConditions();
   
   	//sync();
   
/*	printf(" w Res = %e\n",wNormRes);*/
   	if(fabs(wNormRes)>vAccuracy)		wVelocityLoops = wVelocityLoops + 1;
 }
  

}

