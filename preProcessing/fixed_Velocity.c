void fixedVelocity()
{

	double uTheta,vTheta; 	//angle for u and v centres
	double uVel,vVel;	//tangential velocity at u and v cell centers.
	double angVel;		//angular velocity_BC
	double r; 		//radius
	double xcd,ycd;		//centres of domain
	double pi;
	double uN,uD,vN,vD;
	
	double rhoe,rhon,rhot;
	
	
	

	pi = 3.14159265;
	angVel = 2.0*pi/5.0;
	xcd = xDomainLength/2.0;
	ycd = yDomainLength/2.0;


	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			for(k=1;k<=nz;k++)
			{
				r = pow((pow((xcd-xf[i]),2.0) + pow((ycd-yc[j]),2.0)),0.5);
				uVel =  r*angVel;

				
				r = pow((pow((xcd-xc[i]),2.0) + pow((ycd-yf[j]),2.0)),0.5);
				vVel = r*angVel;


				uN = yc[j]-ycd;
				uD = xf[i]-xcd;
				uTheta = atan2(uN,uD);
				
				vN = yf[j]-ycd;
				vD = xc[i]-xcd;
				vTheta = atan2(vN,vD);
				
/*				printf("uTheta = %f \n",uTheta*180/pi);*/
				
			
/*				u[i][j][k][l] = uVel*sin(uTheta);			//Irrotational flow*/
/*				v[i][j][k][l] = -vVel*cos(vTheta);*/
/*				w[i][j][k][l] = 0.0;*/

				
				if(simulationTime<endTime/2.0)
				{
					v[i][j][k][l] = sin(pi*uD)*cos(pi*uN);
					u[i][j][k][l] = -cos(pi*vD)*sin(pi*vN);			//Shear flow
					w[i][j][k][l] = 0.0;
				
				}
				else
				{
					v[i][j][k][l] = -sin(pi*uD)*cos(pi*uN);
					u[i][j][k][l] = cos(pi*vD)*sin(pi*vN);			//Shear flow
					w[i][j][k][l] = 0.0;
				
				}
				
				
			
			}
		}
	}
	
	
	for (i=2;i<=nxm;i++)
	  {
		ieast = i+1;
		iwest = i-1;
		
	  	for (j=2;j<=nym;j++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	{
	      
		      		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
		      		
				rhoe = 0.5 * rho[i][j][k] +   0.5 * rho[ieast][j][k];
				
				/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
	      			
	      			Fe[i][j][k] = rhoe*u[i][j][k][l]*s; 
              		}
	      	  }
	     }
	     
	for (j=2;j<=nym;j++)
  	{
		jnorth = j+1;
		jsouth = j-1;
	
  		for (i=2;i<=nxm;i++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	{
		  		
		  		s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
		  		
		  		rhon = 0.5 * rho[i][j][k] +   0.5 * rho[i][jnorth][k];
			
				/* Harmonic Interpolation of density*/
/*				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/
	
				Fn[i][j][k] = rhon*v[i][j][k][l]*s;
      				
     			}
    		}
  	}
  	
  	for (k=2;k<=nzm;k++)
  	{
		ktop = k+1;
		kbottom = k-1;
		
  		for (i=2;i<=nxm;i++)
		  {
		  	for(j=2;j<=nym;j++)
		  	{
		  		
		  		s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
				
				rhot = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][ktop];
				
				/* Harmonic Interpolation of density*/
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/
	
				Ft[i][j][k] = rhot*w[i][j][k][l]*s;
		      		
     			}
    		}
  	}		
		
	
	


}
