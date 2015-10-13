#include"solver_Controls.h"

//calculates the max courant number
void courantNumber()
{
	uMax=0;
	vMax=0;
	wMax=0;

	for  (k=1; k<=nz; k++)
	 for (j=1; j<=ny; j++)
	  for(i=1; i<=nx; i++)
	{
		if(u[i][j][k][l]>uMax)	uMax=u[i][j][k][l];
		if(v[i][j][k][l]>vMax)	vMax=v[i][j][k][l];			//finding out max velocity
		if(w[i][j][k][l]>wMax)	wMax=w[i][j][k][l];
	}

	maxCourantNumber = (fabs(uMax/dx) + fabs(vMax/dy) + fabs(wMax/dz))*dt;

	printf("Max courant number: %lf \n",maxCourantNumber);


}

//Increments the time variables
void incrementTime()
{

	if(adjustableCourantNumber==1)
	{
	
	
		if(maxCourantNumber>allowableCourantNumber)
		{
			dt = allowableCourantNumber/(1.2*(uMax/dx + vMax/dy + wMax/dz));
		}
		
		if(maxCourantNumber<leastCourantNumber)
		{
			dt = 1.2*dt;
		}
		
		
		if(dt>=capillaryTimeStep)	dt = capillaryTimeStep;		//Capillary time constraint
	}
	
	
        dtau = dt*0.02/allowableCourantNumber;

	simulationTime = simulationTime + dt;	//time increment

	t=t+1;			//output file name parameter

}



//shifts the values to l=0
void shift()
{

	for  (k=0; k<(zNumberOfCells+5); k++)
	 for (j=0; j<(yNumberOfCells+5); j++)
	  for(i=0; i<(xNumberOfCells+5); i++)
	{
	  u[i][j][k][0]=u[i][j][k][1];
	  v[i][j][k][0]=v[i][j][k][1];
	  w[i][j][k][0]=w[i][j][k][1];
	  p[i][j][k][0]=p[i][j][k][1];			//copying values to l=0
	  pp[i][j][k][0]=pp[i][j][k][1];
	  T[i][j][k][0]=T[i][j][k][1];
	  Psi[i][j][k][0]=Psi[i][j][k][1];
	  alpha[i][j][k][0]=alpha[i][j][k][1];
	  ue[i][j][k][0]=ue[i][j][k][1];
	  vn[i][j][k][0]=vn[i][j][k][1];
	  wt[i][j][k][0]=wt[i][j][k][1];
	  xBFe[i][j][k][0]=xBFe[i][j][k][1];
	  yBFn[i][j][k][0]=yBFn[i][j][k][1];
	  zBFt[i][j][k][0]=zBFt[i][j][k][1];
	  eAlpha[i][j][k][0]=eAlpha[i][j][k][1];
	  wAlpha[i][j][k][0]=wAlpha[i][j][k][1];
	  nAlpha[i][j][k][0]=nAlpha[i][j][k][1];
	  sAlpha[i][j][k][0]=sAlpha[i][j][k][1];
	  tAlpha[i][j][k][0]=tAlpha[i][j][k][1];
	  bAlpha[i][j][k][0]=bAlpha[i][j][k][1];


	  u[i][j][k][1]=u[i][j][k][2];
	  v[i][j][k][1]=v[i][j][k][2];
	  w[i][j][k][1]=w[i][j][k][2];
	  p[i][j][k][1]=p[i][j][k][2];			//copying values to l=0
	  pp[i][j][k][1]=pp[i][j][k][2];
	  T[i][j][k][1]=T[i][j][k][2];
	  Psi[i][j][k][1]=Psi[i][j][k][2];
	  alpha[i][j][k][1]=alpha[i][j][k][2];
	  ue[i][j][k][1]=ue[i][j][k][2];
	  vn[i][j][k][1]=vn[i][j][k][2];
	  wt[i][j][k][1]=wt[i][j][k][2];
	  xBFe[i][j][k][1]=xBFe[i][j][k][2];
	  yBFn[i][j][k][1]=yBFn[i][j][k][2];
	  zBFt[i][j][k][1]=zBFt[i][j][k][2];
	  eAlpha[i][j][k][1]=eAlpha[i][j][k][2];
	  wAlpha[i][j][k][1]=wAlpha[i][j][k][2];
	  nAlpha[i][j][k][1]=nAlpha[i][j][k][2];
	  sAlpha[i][j][k][1]=sAlpha[i][j][k][2];
	  tAlpha[i][j][k][1]=tAlpha[i][j][k][2];
	  bAlpha[i][j][k][1]=bAlpha[i][j][k][2];

	  pp[i][j][k][2]=0;	                          //equating values in l=1 to zero
	  
/*	   Shifting Coefficients*/

        aE[i][j][k][0]=aE[i][j][k][1];
        aW[i][j][k][0]=aW[i][j][k][1];
        aN[i][j][k][0]=aN[i][j][k][1];
        aS[i][j][k][0]=aS[i][j][k][1];
        aT[i][j][k][0]=aT[i][j][k][1];
        aB[i][j][k][0]=aB[i][j][k][1];
        
        aPB[i][j][k][0]=aPB[i][j][k][1];
        aPE[i][j][k][0]=aPE[i][j][k][1];
        aPN[i][j][k][0]=aPN[i][j][k][1];
        aPS[i][j][k][0]=aPS[i][j][k][1];
        aPT[i][j][k][0]=aPT[i][j][k][1];
        aPW[i][j][k][0]=aPW[i][j][k][1];

        Ap[i][j][k][0]=Ap[i][j][k][1];
        
        aE[i][j][k][1]=aE[i][j][k][2];
        aW[i][j][k][1]=aW[i][j][k][2];
        aN[i][j][k][1]=aN[i][j][k][2];
        aS[i][j][k][1]=aS[i][j][k][2];
        aT[i][j][k][1]=aT[i][j][k][2];
        aB[i][j][k][1]=aB[i][j][k][2];
        
        aPB[i][j][k][1]=aPB[i][j][k][2];
        aPE[i][j][k][1]=aPE[i][j][k][2];
        aPN[i][j][k][1]=aPN[i][j][k][2];
        aPS[i][j][k][1]=aPS[i][j][k][2];
        aPT[i][j][k][1]=aPT[i][j][k][2];
        aPW[i][j][k][1]=aPW[i][j][k][2];
        
        Ap[i][j][k][1]=Ap[i][j][k][2];

	}

}

//Starts CPU clock
void startCPUClock()
{

	t=1;			//output file starts from 1st time step

	l = 2;   //making the simulation run at l=1 so that previous time step values are stored in l=0

	simulationTime=dt;		//time counting starts from dt

	assert((start = clock())!=-1);			//starting real time

}

//Calculation of the spacial and temporal step sizes
void stepSize()
{

	dx = xDomainLength/xNumberOfCells;
	dy = yDomainLength/yNumberOfCells;
	dz = zDomainLength/zNumberOfCells;
	deltaT=(endTime-startTime)/nTimeSteps;  //initial time step

	dt=deltaT;  				//equating initial time step to time step

	nx=xNumberOfCells + 2;
	ny=yNumberOfCells + 2;
	nz=zNumberOfCells + 2;
	
	nxm = nx - 1;
	nym = ny - 1;
	nzm = nz - 1;


}


