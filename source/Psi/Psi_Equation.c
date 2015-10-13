void PsiEqn()
{


        double rhob,rhoe,rhon,rhos,rhot,rhow;
        double rr;
        
	
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
			
				ap[i][j][k]  = (- ae[i][j][k] - aw[i][j][k] - an[i][j][k] - as[i][j][k] - at[i][j][k] - ab[i][j][k] + apP[i][j][k])*urecrpsi;
      				sP[i][j][k]  = sP[i][j][k] + (1 - psiUnderRelaxCoeff) * ap[i][j][k] * Psi[i][j][k][l];
      				
      		
							
      				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
                                
                                rr = dt/vol;
                                
/*                                /* Arithmetic Interpolation of density*/
				rhoe = 0.5 * rho[i][j][k] +   0.5 * rho[ieast][j][k];
				rhon = 0.5 * rho[i][j][k] +   0.5 * rho[i][jnorth][k];
				rhot = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][ktop];
				
				rhow = 0.5 * rho[i][j][k] +   0.5 * rho[iwest][j][k];
				rhos = 0.5 * rho[i][j][k] +   0.5 * rho[i][jsouth][k];
				rhob = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][kbottom];
                                

                                sP[i][j][k] += wPsi[i][j][k]*Fe[iwest][j][k]/(2.0*rhow); 
                                sP[i][j][k] -= ePsi[i][j][k]*Fe[i][j][k]/(2.0*rhoe);
                                sP[i][j][k] += sPsi[i][j][k]*Fn[i][jsouth][k]/(2.0*rhos);
                                sP[i][j][k] -= nPsi[i][j][k]*Fn[i][j][k]/(2.0*rhon);
                                sP[i][j][k] += bPsi[i][j][k]*Ft[i][j][kbottom]/(2.0*rhob);
                                sP[i][j][k] -= tPsi[i][j][k]*Ft[i][j][k]/(2.0*rhot);

/*      				*/
/*                                Psi[i][j][k][l] = Psi[i][j][k][l-1];*/
/*                                Psi[i][j][k][l] += wPsi[i][j][k]*Fe[iwest][j][k]*rr/(rhow); */
/*                                Psi[i][j][k][l] -= ePsi[i][j][k]*Fe[i][j][k]*rr/(rhoe);*/
/*                                Psi[i][j][k][l] += sPsi[i][j][k]*Fn[i][jsouth][k]*rr/(rhos);*/
/*                                Psi[i][j][k][l] -= nPsi[i][j][k]*Fn[i][j][k]*rr/(rhon);*/
/*                                Psi[i][j][k][l] += bPsi[i][j][k]*Ft[i][j][kbottom]*rr/(rhob);*/
/*                                Psi[i][j][k][l] -= tPsi[i][j][k]*Ft[i][j][k]*rr/(rhot);*/
                                

      				apP[i][j][k] = 1 / ap[i][j][k];
      				
			}
    		}
  	}


for(nonLinear=0;nonLinear<=nPsiLoops;nonLinear++)
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
	  			


				c[i] 	= sP[i][j][k];

				c[i] = c[i] - an[i][j][k]*Psi[i][jnorth][k][l];
				c[i] = c[i] - as[i][j][k]*Psi[i][jsouth][k][l];
				c[i] = c[i] - at[i][j][k]*Psi[i][j][ktop][l];
				c[i] = c[i] - ab[i][j][k]*Psi[i][j][kbottom][l];
					
				
				c[i] = c[i]*apP[i][j][k];
				
				
				lower[i]= aw[i][j][k]*apP[i][j][k];
				upper[i]= ae[i][j][k]*apP[i][j][k];
/*				*/

/*                            if(Psi[i][j][k][l]>1)   Psi[i][j][k][l] = 1.0;*/
/*                            if(Psi[i][j][k][l]<0)   Psi[i][j][k][l] = 0.0;*/


				
			}
				mThomas(wcID[e][f][g],ecID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=wcID[e][f][g];thomasI<=ecID[e][f][g];thomasI++)
	  				Psi[thomasI][j][k][l]=x[thomasI];				

                    


    	}
  }
  
  
  
  	PsiMaxRes = 0;
  	PsiTotalRes = 0;
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
	  			
				PsiRes[i][j][k] 	= sP[i][j][k] - (ap[i][j][k]*Psi[i][j][k][l]);
				
				PsiRes[i][j][k] = PsiRes[i][j][k] - an[i][j][k]*Psi[i][jnorth][k][l];
				PsiRes[i][j][k] = PsiRes[i][j][k] - as[i][j][k]*Psi[i][jsouth][k][l];
				PsiRes[i][j][k] = PsiRes[i][j][k] - at[i][j][k]*Psi[i][j][ktop][l];
				PsiRes[i][j][k] = PsiRes[i][j][k] - ab[i][j][k]*Psi[i][j][kbottom][l];
				PsiRes[i][j][k] = PsiRes[i][j][k] - aw[i][j][k]*Psi[iwest][j][k][l];
				PsiRes[i][j][k] = PsiRes[i][j][k] - ae[i][j][k]*Psi[ieast][j][k][l];
				
				
				//uTotalRes += fabs(uRes[i][j][k]);
				PsiTotalRes += fabs(ap[i][j][k]*Psi[i][j][k][l]);
				if(fabs(PsiRes[i][j][k])>PsiMaxRes)	PsiMaxRes	= fabs(PsiRes[i][j][k]);
			}
    		}
   	}
	
	PsiMeanRes = PsiTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
	if(PsiMeanRes==0)		PsiNormRes = 0;
    	else			PsiNormRes = PsiMaxRes/PsiMeanRes;
   
	
	//sync();
	
/*	printf(" T Res = %e\n",TMaxRes);*/
	//printf(" u Max Res = %e\n",uTotalRes);
   	if(fabs(PsiNormRes)>TAccuracy)		nPsiLoops = nPsiLoops + 1;
   	
   	
}
  
}
