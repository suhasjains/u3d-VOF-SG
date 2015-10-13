//calculates the value of initial u velocity
void uEqn()
{
	
	for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				ae[i][j][k]  = (aee[i][j][k] + aw[i][j][k] + aNe[i][j][k] + aSe[i][j][k] + aTe[i][j][k] + aBe[i][j][k] + apu[i][j][k])*urecru;
      				su[i][j][k]  = su[i][j][k] + (1 - uUnderRelaxCoeff) * ae[i][j][k] * u[i][j][k][l];
      				apu[i][j][k] = 1 / ae[i][j][k];	/* For implementation in the pressure correction approach */
				
			}
    		}
  	}


for(nonLinear=0;nonLinear<=uVelocityLoops;nonLinear++)
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
	  			
	  			
				c[i] 	= su[i][j][k];

				c[i] = c[i] + aNe[i][j][k]*u[i][jnorth][k][l];
				c[i] = c[i] + aSe[i][j][k]*u[i][jsouth][k][l];
				c[i] = c[i] + aTe[i][j][k]*u[i][j][ktop][l];
				c[i] = c[i] + aBe[i][j][k]*u[i][j][kbottom][l];
				
				if(i==2)	c[i] = c[i] + aw[i][j][k]*u[iwest][j][k][l];
				if(i==nxm)	c[i] = c[i] + aee[i][j][k]*u[ieast][j][k][l];
					
				
				c[i] = c[i]*apu[i][j][k];
				
				
				lower[i]= -aw[i][j][k]*apu[i][j][k];
				upper[i]= -aee[i][j][k]*apu[i][j][k];
			}
				
				mThomas(wcID[e][f][g],ecID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=wcID[e][f][g];thomasI<=ecID[e][f][g];thomasI++)
	  				u[thomasI][j][k][l]=x[thomasI];//*uUnderRelaxCoeff + (1.0-uUnderRelaxCoeff)*u[thomasI][j][k][l];				
			
    	}
  }
  
  
  
  uMaxRes = 0;
  uTotalRes = 0;
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
	  			
				uRes[i][j][k] 	= su[i][j][k] - (ae[i][j][k]*u[i][j][k][l]);
				
				uRes[i][j][k] = uRes[i][j][k] + aNe[i][j][k]*u[i][jnorth][k][l];
				uRes[i][j][k] = uRes[i][j][k] + aSe[i][j][k]*u[i][jsouth][k][l];
				uRes[i][j][k] = uRes[i][j][k] + aTe[i][j][k]*u[i][j][ktop][l];
				uRes[i][j][k] = uRes[i][j][k] + aBe[i][j][k]*u[i][j][kbottom][l];
				uRes[i][j][k] = uRes[i][j][k] + aw[i][j][k]*u[iwest][j][k][l];
				uRes[i][j][k] = uRes[i][j][k] + aee[i][j][k]*u[ieast][j][k][l];
				
				
				//uTotalRes += fabs(uRes[i][j][k]);
				uTotalRes += fabs(ae[i][j][k]*u[i][j][k][l]);
				if(fabs(uRes[i][j][k])>uMaxRes)	uMaxRes	= fabs(uRes[i][j][k]);
			}
    		}
   	}
	
	uMeanRes = uTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
	if(uMeanRes==0)		uNormRes = 0;
    	else			uNormRes = uMaxRes/uMeanRes;
   
	
	
	velocityBoundaryConditions();
	//sync();
	
/*	printf(" u Res = %e\n",uNormRes);*/
	//printf(" u Max Res = %e\n",uTotalRes);
   	if(fabs(uNormRes)>vAccuracy)		uVelocityLoops = uVelocityLoops + 1;
   	
   	
}
  
}

