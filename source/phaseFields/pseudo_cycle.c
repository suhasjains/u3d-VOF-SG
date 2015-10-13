void pseudoCycle()
{


        printf("Starting VOF Loops\n");



/*        Copying values to pseudo variables*/

        for (k=0; k<(zNumberOfCells+5); k++)
	 for (j=0; j<(yNumberOfCells+5); j++)
	  for(i=0; i<(xNumberOfCells+5); i++)
	{
	   
                tauAlpha[i][j][k][l-1]  = alpha[i][j][k][l-1];

	}        
        
        nPseudoLoops = 5;
        for(pseudoLoops=0;pseudoLoops<=nPseudoLoops;pseudoLoops++)
        {

/*                Solving in pseudo time */
        
                alphaEqnCoeff();
	
        	alphaEqn();
        	
        
        
                
        /*        Checking convergence        	*/
        
        	pseudoMaxRes = 0;
        	pseudoTotalRes = 0;
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
        	  			
        	  			pseudoRes[i][j][k] = (tauAlpha[i][j][k][l] - tauAlpha[i][j][k][l-1])/dtau; 
        	  			
/*        	  			printf("tauAlpha = %e  %e\n",tauAlpha[i][j][k][l],tauAlpha[i][j][k][l-1]);*/

                                        pseudoTotalRes += fabs(pseudoRes[i][j][k]);
        				if(fabs(pseudoRes[i][j][k])>pseudoMaxRes)	pseudoMaxRes = fabs(pseudoRes[i][j][k]);
        				
        
				        //if(alpha[i][j][k][l]!=1.0&&alpha[i][j][k][l]!=1.0)		printf("alpha[%d][%d][%d] = %e \n",i,j,k,alpha[i][j][k][l]);
        
        /*				if(alphaRes[20][20][2]!=0)	printf(" alphaRes[%d][%d][%d] = %e\n",20,20,2,alphaRes[20][20][2]);*/
        
        	  			
        			}
            		}
            	}
            	
            	
            	
            	pseudoMeanRes = pseudoTotalRes/((double)((nxm-1)*(nym-1)*(nzm-1)));
	        if(pseudoMeanRes==0)		pseudoNormRes = 0;
	        else			        pseudoNormRes = pseudoMaxRes/pseudoMeanRes;
            	
    	
            	if(fabs(pseudoMaxRes)>pseudoTimeAccuracy)		nPseudoLoops = nPseudoLoops + 1;
            	
            	
            	
/*            	Shift in pseudo time*/
            	
            	shiftPseudo();
            	
/*            	printf("pseudoMaxRes = %e \n",pseudoMaxRes);*/
            	
    	}
    	
    	printf("pseudoMaxRes = %e \n",pseudoMaxRes);
    	printf("No of pseudo loops = %d \n",nPseudoLoops);
    	
    	
/*    	Copying back values from pseudo varibales*/
    	
        for (k=0; k<(zNumberOfCells+5); k++)
	 for (j=0; j<(yNumberOfCells+5); j++)
	  for(i=0; i<(xNumberOfCells+5); i++)
	{
	   
                alpha[i][j][k][l]  = tauAlpha[i][j][k][l];

	} 		
	
	
	correctAlpha();
}
