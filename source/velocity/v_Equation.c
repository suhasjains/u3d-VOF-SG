#include"velocity.h"

//calculates the value of initial v velocity
void vEqn()
{
	 for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
		 		an[i][j][k]  = (anE[i][j][k] + anW[i][j][k] + ann[i][j][k] + as[i][j][k] + anT[i][j][k] + anB[i][j][k] + apv[i][j][k])*urecrv;
      				sv[i][j][k]  = sv[i][j][k] + (1 - vUnderRelaxCoeff) * an[i][j][k] * v[i][j][k][l] ;
      				apv[i][j][k] = 1 / an[i][j][k];	/* For implementation in the pressure correction approach */
		 	
		 	}
 	     }
  	 }
  	 
  for(nonLinear=0;nonLinear<=vVelocityLoops;nonLinear++)
  {
  	
//  Creating matrices for ADI-TDMA algorithm
	for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
	{
		for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[j] 	= sv[i][j][k];
				
				
				c[j] = c[j] + anT[i][j][k]*v[i][j][ktop][l];
				c[j] = c[j] + anB[i][j][k]*v[i][j][kbottom][l];
				c[j] = c[j] + anW[i][j][k]*v[iwest][j][k][l];
				c[j] = c[j] + anE[i][j][k]*v[ieast][j][k][l];

				if(j==2)	c[j] = c[j] + as[i][j][k]*v[i][jsouth][k][l];
				if(j==nym)	c[j] = c[j] + ann[i][j][k]*v[i][jnorth][k][l];				
								
				c[j] = c[j]*apv[i][j][k]; 
				
				
				lower[j]= -as[i][j][k]*apv[i][j][k];
				upper[j]= -ann[i][j][k]*apv[i][j][k];
			}
				mThomas(scID[e][f][g],ncID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=scID[e][f][g];thomasI<=ncID[e][f][g];thomasI++)
	  				v[i][thomasI][k][l]=x[thomasI];				
			
    	}
   }
  
  
	vMaxRes = 0;
	vTotalRes = 0;  
  	for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
	{
		for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
	  			vRes[i][j][k] 	= sv[i][j][k] - an[i][j][k]*v[i][j][k][l];
	  			
				vRes[i][j][k] = vRes[i][j][k] + ann[i][j][k]*v[i][jnorth][k][l];
				vRes[i][j][k] = vRes[i][j][k] + as[i][j][k]*v[i][jsouth][k][l];
				vRes[i][j][k] = vRes[i][j][k] + anT[i][j][k]*v[i][j][ktop][l];
				vRes[i][j][k] = vRes[i][j][k] + anB[i][j][k]*v[i][j][kbottom][l];
				vRes[i][j][k] = vRes[i][j][k] + anW[i][j][k]*v[iwest][j][k][l];
				vRes[i][j][k] = vRes[i][j][k] + anE[i][j][k]*v[ieast][j][k][l];
	  			
	  			
				vTotalRes += fabs(an[i][j][k]*v[i][j][k][l]);
				if(fabs(vRes[i][j][k])>vMaxRes)	vMaxRes	= fabs(vRes[i][j][k]);
	  			
			}
    		}
   	}
	
	vMeanRes = vTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
	if(vMeanRes==0)		vNormRes = 0;
    	else			vNormRes = vMaxRes/vMeanRes;
   
	 velocityBoundaryConditions();  
   	//sync();
/*	printf(" v Res = %e\n",vNormRes);*/
   	if(fabs(vNormRes)>vAccuracy)	vVelocityLoops = vVelocityLoops + 1;
   
 }
   
}
