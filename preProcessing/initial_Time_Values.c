#include"preProcessing.h"

//Gives initial values to all the variables and fields at time 0
void zeroTimeValues()
{

	//INPUT OUTPUT VARIABLES
	nx = ny = nz = n = 0;
	deltaT = simulationTime = 0;
	adjustTimeStep = 0;

	//COEFFICIENTS OF VELOCITY EQUATION
	dx = dy = dz = dt = 0;

	//COEFFICIENTS OF PRESSURE CORRECTION EQUATION
	maxMassResidual=0;
	uMaxRes = vMaxRes = wMaxRes = 0;

	//VELOCITIES and PRESSURE
	for   (l=0; l<3; l++)
	 for  (k=0; k<(zNumberOfCells+5); k++)
	  for (j=0; j<(yNumberOfCells+5); j++)
	   for(i=0; i<(xNumberOfCells+5); i++)
	{
		u[i][j][k][l]	  = 0;
		v[i][j][k][l]	  = 0;
		w[i][j][k][l] 	  = 0;
		
		T[i][j][k][l]	 = 0;


		p[i][j][k][l]     = 0;
		pp[i][j][k][l] = 0;
		

		xCoord[i][j][k] = 0;
		yCoord[i][j][k] = 0;
		zCoord[i][j][k] = 0;


       	 	Fe[i][j][k] = Fn[i][j][k] = Ft[i][j][k] = 0;

		ae[i][j][k] = an[i][j][k] = at[i][j][k] = aw[i][j][k] = as[i][j][k] = ab[i][j][k] = ap[i][j][k] = 0;

	    	dpx[i][j][k] = dpy[i][j][k] = dpz[i][j][k] = 0;
	    
	    	rho[i][j][k] = 0;
		mu[i][j][k]  = 0;

	}


	for(i=0; i<(xNumberOfCells+5); i++)
	{
		lower[i] =  0;
		upper[i] = 0;
		c[i] = 0;
		x[i] = 0;
		m[i] = 0;
		cc[i] = 0;
	}

	for(i=0; i<2; i++)
	{
		uMaxResidual[i]=0;
		vMaxResidual[i]=0;
		wMaxResidual[i]=0;
	}

	i = j = k = l = t = thomasI = simpleVar = nonLinear = 0;   //indices

	nx = ny = nz = n = 0;                                  //Number of Cells

	uMax = vMax = wMax = maxCourantNumber = 0;

	check = f1 = f2 = f3 = f4 = f5 = f6 = 0;


	

}
