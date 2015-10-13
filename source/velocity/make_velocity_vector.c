//Interpolation of the velocity vectors back from staggerd cell centers to the pressure cell centers
void makeVelocityVectors()
{

	//Internal Values
	for  (k=2; k<=nzm; k++)
	 for (i=2; i<=nxm; i++)
	  for(j=2; j<=nym; j++)
	  {
	 	finalU[i][j][k][l] = (u[i][j][k][l] + u[i-1][j][k][l])/2.0;
	 	finalV[i][j][k][l] = (v[i][j][k][l] + v[i][j-1][k][l])/2.0;
	 	finalW[i][j][k][l] = (w[i][j][k][l] + w[i][j][k-1][l])/2.0;
	  }
	  
	 //Boundary Values
	 /* Bottom boundary */
  	j = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	
		
			if(BottomBoundType==1)
			{
				//Wall 
    				finalU[i][j-1][k][l] 	= (u[i][j][k][l] + u[i][j-1][k][l])/2.0;
    				finalV[i][j-1][k][l]	= v[i][j-1][k][l];
    				finalW[i][j-1][k][l]	= (w[i][j][k][l] + w[i][j-1][k][l])/2.0;
    				
	    		}
	    		
	    		else if(BottomBoundType==2)
	    		{    	
	    		
	    			//symmetry
	    			finalU[i][j-1][k][l] 	= (u[i][j][k][l] + u[i][j-1][k][l])/2.0;
    				finalV[i][j-1][k][l]	= v[i][j-1][k][l];
    				finalW[i][j-1][k][l]	= (w[i][j][k][l] + w[i][j-1][k][l])/2.0;
			}		
			
			else if(BottomBoundType==3)
			{
				//zero Gradient
				finalU[i][j-1][k][l] 	= (u[i][j][k][l] + u[i][j-1][k][l])/2.0;
    				finalV[i][j-1][k][l]	= (v[i][j][k][l] + v[i][j-1][k][l])/2.0;
    				finalW[i][j-1][k][l]	= (w[i][j][k][l] + w[i][j-1][k][l])/2.0;
			}
			
			else if(BottomBoundType==4)
			{
				//Specified Pressure

			}		
	  	}
  	  }
  	
	    
  	  /* Top boundary */
  	j = nym;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{

		
			if(TopBoundType==1)
			{
			
				//Wall
				finalU[i][j+1][k][l] 	= (u[i][j][k][l] + u[i][j+1][k][l])/2.0;
    				finalV[i][j+1][k][l]	= v[i][j][k][l];
    				finalW[i][j+1][k][l]	= (w[i][j][k][l] + w[i][j+1][k][l])/2.0;
				
	    		}
	    		
	    		else if(TopBoundType==2)
	    		{
	    		    	//symmetry
	    			finalU[i][j+1][k][l] 	= (u[i][j][k][l] + u[i][j+1][k][l])/2.0;
    				finalV[i][j+1][k][l]	= v[i][j][k][l];
    				finalW[i][j+1][k][l]	= (w[i][j][k][l] + w[i][j+1][k][l])/2.0;
			}
			
			else if(TopBoundType==3)
			{
				//zero Gradient
				finalU[i][j+1][k][l] 	= (u[i][j][k][l] + u[i][j+1][k][l])/2.0;
    				finalV[i][j+1][k][l]	= (v[i][j][k][l] + v[i][j-1][k][l])/2.0;
    				finalW[i][j+1][k][l]	= (w[i][j][k][l] + w[i][j+1][k][l])/2.0;
			}
			
			else if(TopBoundType==4)
			{
				//Specified Pressure

 			}
			
	  	}
	  }
  
  	/* West boundary conditions */
	i = 2;
  	for (j=2;j<=nym;j++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{

		
			if(LeftBoundType==1)
			{
				//Wall
    				finalU[i-1][j][k][l] 	= u[i-1][j][k][l];
    				finalV[i-1][j][k][l]	= (v[i][j][k][l] + v[i-1][j][k][l])/2.0;
    				finalW[i-1][j][k][l]	= (w[i][j][k][l] + w[i-1][j][k][l])/2.0;
				
	  		}	
	  		
	  		else if(LeftBoundType==2)
	  		{
			  	//symmetry
	  			finalU[i-1][j][k][l] 	= u[i-1][j][k][l];
    				finalV[i-1][j][k][l]	= (v[i][j][k][l] + v[i-1][j][k][l])/2.0;
    				finalW[i-1][j][k][l]	= (w[i][j][k][l] + w[i-1][j][k][l])/2.0;
			}
			
			else if(LeftBoundType==3)
			{
				//zero Gradient
				finalU[i-1][j][k][l] 	= (u[i][j][k][l] + u[i-1][j][k][l])/2.0;
    				finalV[i-1][j][k][l]	= (v[i][j][k][l] + v[i-1][j][k][l])/2.0;
    				finalW[i-1][j][k][l]	= (w[i][j][k][l] + w[i-1][j][k][l])/2.0;
			}
			
			else if(LeftBoundType==4)
			{
				//specified Pressure
			}
		
		//printf("%e	%e\n",yf[j]-yf[j-1],dy);
		}
    }
          
  	/* East boundary conditions */
 	i = nxm;
  	for (j=2;j<=nym;j++)
  	{
  		for(k=2;k<=nzm;k++)
  		{
		
			if(RightBoundType==1)
			{	
				//Wall
				finalU[i+1][j][k][l] 	= u[i][j][k][l];
    				finalV[i+1][j][k][l]	= (v[i][j][k][l] + v[i+1][j][k][l])/2.0;
    				finalW[i+1][j][k][l]	= (w[i][j][k][l] + w[i+1][j][k][l])/2.0;
			}
	
			else if(RightBoundType==2)
			{
				//symmetry
		  		finalU[i+1][j][k][l] 	= u[i][j][k][l];
    				finalV[i+1][j][k][l]	= (v[i][j][k][l] + v[i+1][j][k][l])/2.0;
    				finalW[i+1][j][k][l]	= (w[i][j][k][l] + w[i+1][j][k][l])/2.0;
			}
	
			else if(RightBoundType==3)
			{		
				//zero Gradient
				finalU[i+1][j][k][l] 	= (u[i][j][k][l] + u[i-1][j][k][l])/2.0;
    				finalV[i+1][j][k][l]	= (v[i][j][k][l] + v[i+1][j][k][l])/2.0;
    				finalW[i+1][j][k][l]	= (w[i][j][k][l] + w[i+1][j][k][l])/2.0;
			}
	
			else if(RightBoundType==4)
			{
				//specified Pressure
			}	
	
  		}
  	}  
  
 	/* Back bounday conditions */
	k = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
		
			if(BackBoundType==1)
			{
				//Wall
    				finalU[i][j][k-1][l] 	= (u[i][j][k][l] + u[i][j][k-1][l])/2.0;
    				finalV[i][j][k-1][l]	= (v[i][j][k][l] + v[i][j][k-1][l])/2.0;
    				finalW[i][j][k-1][l]	= w[i][j][k-1][l];
    				
			}
			
			else if(BackBoundType==2)
			{
				//symmetry
				finalU[i][j][k-1][l] 	= (u[i][j][k][l] + u[i][j][k-1][l])/2.0;
    				finalV[i][j][k-1][l]	= (v[i][j][k][l] + v[i][j][k-1][l])/2.0;
    				finalW[i][j][k-1][l]	= w[i][j][k-1][l];
			}
			
			else if(BackBoundType==3)
			{
				//zero Gradient
				finalU[i][j][k-1][l] 	= (u[i][j][k][l] + u[i][j][k-1][l])/2.0;
    				finalV[i][j][k-1][l]	= (v[i][j][k][l] + v[i][j][k-1][l])/2.0;
    				finalW[i][j][k-1][l]	= (w[i][j][k][l] + w[i][j][k-1][l])/2.0;
			}
			
			else if(BackBoundType==4)
			{
				//Specified Pressure
			}		

		}
    	}
          
  	/* Front boundary conditions */
  	k = nzm;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{

			
			if(FrontBoundType==1)
			{
				//Wall
	    			finalU[i][j][k+1][l] 	= (u[i][j][k][l] + u[i][j][k+1][l])/2.0;
    				finalV[i][j][k+1][l]	= (v[i][j][k][l] + v[i][j][k+1][l])/2.0;
    				finalW[i][j][k+1][l]	= w[i][j][k][l];
	    			
			}	
		
			else if(FrontBoundType==2)
			{
				//symmetry
				finalU[i][j][k+1][l] 	= (u[i][j][k][l] + u[i][j][k+1][l])/2.0;
    				finalV[i][j][k+1][l]	= (v[i][j][k][l] + v[i][j][k+1][l])/2.0;
    				finalW[i][j][k+1][l]	= w[i][j][k][l];
			}
			
			else if(FrontBoundType==3)
			{
				//Zero Gradient
				finalU[i][j][k+1][l] 	= (u[i][j][k][l] + u[i][j][k+1][l])/2.0;
    				finalV[i][j][k+1][l]	= (v[i][j][k][l] + v[i][j][k+1][l])/2.0;
    				finalW[i][j][k+1][l]	= (w[i][j][k][l] + w[i][j][k-1][l])/2.0;
			}
			
			else if(FrontBoundType==4)
			{
				//Specified Pressure
		  		
			}
		}
    	}
	  
	
	  
}
