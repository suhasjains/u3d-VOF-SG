#include"boundaryConditions.h"

//Velocity boundary Conditions
void velocityBoundaryConditions()	
{
	/* Bottom boundary conditions */
  	j = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	
		
			if(BottomBoundType==1)
			{
				//Wall 
    				u[i][j-1][k][l]	= 2.0 * uBottomBoundVel  - u[i][j][k][l];
    				v[i][j-1][k][l]	= vBottomBoundVel;
    				w[i][j-1][k][l]	= 2.0 * wBottomBoundVel  - w[i][j][k][l];
    				
	    		}
	    		
	    		else if(BottomBoundType==2)
	    		{    	
	    		
	    			//symmetry
	    			u[i][j-1][k][l] = u[i][j][k][l];
	    			v[i][j-1][k][l]	= vBottomBoundVel;
	    			w[i][j-1][k][l] = w[i][j][k][l];
			}		
			
			else if(BottomBoundType==3)
			{
				//zero Gradient
				u[i][j-1][k][l] = u[i][j][k][l];
	    			w[i][j-1][k][l] = w[i][j][k][l];
				v[i][j-1][k][l] = v[i][j][k][l];
			}
			
			else if(BottomBoundType==4)
			{
				//Specified Pressure

			}		
	  	}
  	  }
  	
	    
  	  /* Top boundary conditions */
  	j = nym;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{

		
			if(TopBoundType==1)
			{
			
				//Wall
				u[i][j+1][k][l] = 2.0 * uTopBoundVel   - u[i][j][k][l] ;
				v[i][j][k][l] = vTopBoundVel;
				w[i][j+1][k][l] = 2.0 * wTopBoundVel   - w[i][j][k][l];
				
	    		}
	    		
	    		else if(TopBoundType==2)
	    		{
	    		    	//symmetry
	    			u[i][j+1][k][l] = u[i][j][k][l];
	    			v[i][j][k][l] = vTopBoundVel;
	    			w[i][j+1][k][l] = w[i][j][k][l];
			}
			
			else if(TopBoundType==3)
			{
				//zero Gradient
				u[i][j+1][k][l] = u[i][j][k][l];
	    			w[i][j+1][k][l] = w[i][j][k][l];
 				v[i][j][k][l] = v[i][j-1][k][l];
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
    				u[i-1][j][k][l]	= uLeftBoundVel;
    				v[i-1][j][k][l]	= 2.0 * vLeftBoundVel   - v[i][j][k][l];
    				w[i-1][j][k][l]	= 2.0 * wLeftBoundVel   - w[i][j][k][l];
				
	  		}	
	  		
	  		else if(LeftBoundType==2)
	  		{
			  	//symmetry
	  			u[i-1][j][k][l]	= uLeftBoundVel;
	  			v[i-1][j][k][l] = v[i][j][k][l];
				w[i-1][j][k][l] = w[i][j][k][l];
			}
			
			else if(LeftBoundType==3)
			{
				//zero Gradient
				v[i-1][j][k][l] = v[i][j][k][l];
				w[i-1][j][k][l] = w[i][j][k][l];
				u[i-1][j][k][l] = u[i][j][k][l];
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
    				u[i][j][k][l] = uRightBoundVel;
    				v[i+1][j][k][l] = 2.0 * vRightBoundVel  - v[i][j][k][l];
    				w[i+1][j][k][l] = 2.0 * wRightBoundVel  - w[i][j][k][l];
			}
	
			else if(RightBoundType==2)
			{
				//symmetry
		  		u[i][j][k][l] = uRightBoundVel;
		  		v[i+1][j][k][l] = v[i][j][k][l];
				w[i+1][j][k][l] = w[i][j][k][l];
			}
	
			else if(RightBoundType==3)
			{		
				//zero Gradient
				v[i+1][j][k][l] = v[i][j][k][l];
				w[i+1][j][k][l] = w[i][j][k][l];
				u[i][j][k][l] = u[i-1][j][k][l];
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
    				w[i][j][k-1][l]    = wBackBoundVel;
    				u[i][j][k-1][l]    = 2.0 * uBackBoundVel - u[i][j][k][l];
    				v[i][j][k-1][l]    = 2.0 * vBackBoundVel - v[i][j][k][l];
    				
			}
			
			else if(BackBoundType==2)
			{
				//symmetry
				w[i][j][k-1][l]	= wBackBoundVel;
				u[i][j][k-1][l]	= u[i][j][k][l];
				v[i][j][k-1][l] = v[i][j][k][l];
			}
			
			else if(BackBoundType==3)
			{
				//zero Gradient
				u[i][j][k-1][l] = u[i][j][k][l];
				v[i][j][k-1][l] = v[i][j][k][l];
				w[i][j][k-1][l] = w[i][j][k][l];
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
	    			w[i][j][k][l]	= wFrontBoundVel;
	    			u[i][j][k+1][l] = 2.0 * uFrontBoundVel    - u[i][j][k][l];
	    			v[i][j][k+1][l] = 2.0 * vFrontBoundVel    - v[i][j][k][l];
	    			
			}	
		
			else if(FrontBoundType==2)
			{
				//symmetry
				w[i][j][k][l]	= wFrontBoundVel;
				u[i][j][k+1][l] = u[i][j][k][l];
				v[i][j][k+1][l] = v[i][j][k][l];
			}
			
			else if(FrontBoundType==3)
			{
				//Zero Gradient
				u[i][j][k+1][l] = u[i][j][k][l];
				v[i][j][k+1][l] = v[i][j][k][l];
				w[i][j][k][l] = w[i][j][k-1][l]; 
			}
			
			else if(FrontBoundType==4)
			{
				//Specified Pressure
		  		
			}
		}
    	}
}
