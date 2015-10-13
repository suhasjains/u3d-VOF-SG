void definePhaseField()
{

	double r1,r2;
	
	double Rad,x1,y1,x2,y2;
	
	double alphaR[105][105][10];
	
	double dRRx,dRRy,dRx,dRy;
	
	double xRc[20],yRc[20];
	
	int iRef,jRef;
	

	//Initial Condition of the Phase field
	for   (l=0; l<2; l++)
	 for  (k=0; k<(zNumberOfCells+5); k++)
	  for (j=0; j<(yNumberOfCells+5); j++)
	   for(i=0; i<(xNumberOfCells+5); i++)
	{
		if(nFluids>1)
		{
		
		
			alpha[i][j][k][l] = 0;
		
			//Circle relative to domain
/*			r1 = 0.25*0.25*nx*ny;*/
/*			r2 = pow(r1 - pow((j - ny/2.0),2.0),0.5);*/
/*			if((i>=nx/2.0-r2)&&(i<=nx/2.0+r2))	alpha[i][j][k][l] = 1;	*/
			


/*			//Circle based on dimensions*/
			Rad = 0.025;	//Radius of circle
			x1 = 0.05;
			y1 = 0.05;	//Coordinates of centre of circle 

			r1 = Rad*Rad;
			r2 = pow(r1 - pow((yc[j] - y1),2.0),0.5);
			if((xc[i]>=x1-r2)&&(xc[i]<=x1+r2))			alpha[i][j][k][l] = 1;
			
			
/*			Rad = 0.15;	//Radius of circle*/
/*			x2 = 0.66;*/
/*			y2 = 0.5;	//Coordinates of centre of circle */

/*			r1 = Rad*Rad;*/
/*			r2 = pow(r1 - pow((yc[j] - y2),2.0),0.5);*/
/*			if((xc[i]>=x2-r2)&&(xc[i]<=x2+r2))			alpha[i][j][k][l] = 1;*/


			//zaleski disk
/*			Rad = 0.25;	//Radius of circle*/
/*			x1 = 0.5;*/
/*			y1 = 0.5;	//Coordinates of centre of circle */

/*			r1 = Rad*Rad;*/
/*			r2 = pow(r1 - pow((yc[j] - y1),2.0),0.5);*/
/*			if((xc[i]>=x1-r2)&&(xc[i]<=x1+r2))						alpha[i][j][k][l] = 1;*/
			
				
			//Square
/*			if((j<=(2*ny/3))&&((j>=(ny/3)))&&(i<=(2*nx/3))&&(i>=(nx/3))&&(k<=(2*nz/3))&&(k>=(nz/3)))	alpha[i][j][k][l] = 1;*/
			
			//Dam
/*			if((j<=(ny/2))&&(i<=(nx/3)))	alpha[i][j][k][l] = 1;*/

		
/*			if((j<=(ny/2)))		alpha[i][j][k][l] = 1;*/

			
			//Dam 2
/*			if((i<=(200.0/5.366))&&(j<=(50.0/3.0)))		alpha[i][j][k][l] = 1;*/
			
		
			//Nanofluids
/*			alpha[i][j][k][l] = 0.01;*/


			//Rayleigh Taylor instability
/*			if((j<=(ny/2)))		alpha[i][j][k][l] = 1;*/
		
		
		
					
		}
		
		else	alpha[i][j][k][l] = 1;
		
	}
	
	
	//Refinement of alpha 
	
	l=1;
	 for (i=1;i<=nxm;i++)
	  {
		ieast = i+1;
		iwest = i-1;
	
	  	for (j=2;j<=nym;j++)
		  {
		  	jnorth = j+1;
			jsouth = j-1;
		  	
		  	for(k=2;k<=nzm;k++)
		  	{
		  		ktop = k+1;
				kbottom = k-1;
				
				alphaR[i][j][k] = alpha[i][j][k][l];
	
				if((alpha[i][j][k][l]==alpha[ieast][j][k][l])&&(alpha[i][j][k][l]==alpha[iwest][j][k][l])&&(alpha[i][j][k][l]==alpha[i][jnorth][k][l])&&(alpha[i][j][k][l]==alpha[i][jsouth][k][l]))
				{
					continue;
				}
		
				else
				{
					//printf("Yes\n");
				
					alphaR[i][j][k] = 0;
				
					//Stepsize of refined grid
					dRRx = (xf[i] - xf[iwest])/32.0;
					dRx = (xf[i] - xf[iwest])/16.0;
					dRRy = (yf[j] - yf[jsouth])/32.0;
					dRy = (yf[j] - yf[jsouth])/16.0;
					
					//Giving coordinate values for refined grid
					xRc[8] = xc[i] - dRRx;
					yRc[8] = yc[j] - dRRy;
					xRc[9] = xc[i] + dRRx;
					yRc[9] = yc[j] + dRRy;


					for(iRef=7;iRef>=1;iRef--)
						xRc[iRef] = xRc[iRef+1] - dRx;
						
					for(iRef=9;iRef<=16;iRef++)
						xRc[iRef] = xRc[iRef-1] + dRx;
				
					for(jRef=7;jRef>=1;jRef--)
						yRc[jRef] = yRc[jRef+1] - dRy;
						
					for(jRef=9;jRef<=16;jRef++)
						yRc[jRef] = yRc[jRef-1] + dRy;
						
						
					for(iRef=1;iRef<=16;iRef++)
					{
						for(jRef=1;jRef<=16;jRef++)
						{
						
							r2 = pow(r1 - pow((yRc[jRef] - y1),2.0),0.5);
							if((xRc[iRef]>=x1-r2)&&(xRc[iRef]<=x1+r2))	alphaR[i][j][k] = alphaR[i][j][k] + 1.0/256.0;
							
							r2 = pow(r1 - pow((yRc[jRef] - y2),2.0),0.5);
							if((xRc[iRef]>=x2-r2)&&(xRc[iRef]<=x2+r2))	alphaR[i][j][k] = alphaR[i][j][k] + 1.0/256.0;
							
						
							
						}
					
					}
					
				}	
					alpha[i][j][k][0] = alphaR[i][j][k];
					alpha[i][j][k][1] = alphaR[i][j][k];
			
				
			}
		}
	}
/*	*/
/*	for   (l=0; l<2; l++)*/
/*	 for  (k=0; k<(zNumberOfCells+5); k++)*/
/*	  for (j=0; j<(yNumberOfCells+5); j++)*/
/*	   for(i=0; i<(xNumberOfCells+5); i++)*/
/*	{*/
/*		if(nFluids>1)*/
/*		{*/
/*		*/
/*			//slot for zaleski disk*/
/*			if((xc[i]>=0.45)&&(xc[i]<=0.55)&&(yc[j]>=0.1)&&(yc[j]<=0.45))		alpha[i][j][k][l] = 0;*/
/*			*/
/*					*/
/*		}*/
/*		*/
/*		else	alpha[i][j][k][l] = 1;*/
/*		*/
/*	}*/
	
	
}
