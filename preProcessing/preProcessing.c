#include"preProcessing.h"

//Finalizes things inside simple loop
void finalize()
{

	  stop = clock();
	  realTime = (double) (stop-start)/CLOCKS_PER_SEC;
          printf("Run time: %2.2f seconds \n", realTime);				//printing real time

}

//Initializes things inside simple loop
void initialize()
{

	l = 2;   //making the simulation run at l=1 so that previous time step values are stored in l=0

	nSimpleLoops = 2;  //making initial no of simple loops = 1
	nSimplerLoops = 1; //making initial no of simpler loops = 1

	printf("\n\n\nsimulationTime :%2.10lf seconds\n",simulationTime);
	
	printf("dt = %lf\n",dt);
	
	urecru	= 1 / uUnderRelaxCoeff;
	urecrv	= 1 / vUnderRelaxCoeff;
	urecrw	= 1 / wUnderRelaxCoeff;
	
	urecrp	= 1 / pUnderRelaxCoeff;
	urecrt  = 1 / tUnderRelaxCoeff;
	urecrpsi = 1 / psiUnderRelaxCoeff; 

	
}

//Initializing number of loops
void nLoops()
{

	uVelocityLoops = 1;
	vVelocityLoops = 1;
	wVelocityLoops = 1;
	nPressureLoops = 1;
	nTemperatureLoops = 1;
	nPsiLoops = 1;
	 
	
	
	
  /* Initializing temporary variables to 0*/
  for (i=1;i<=nx;i++)
  {
    for (j=1;j<=ny;j++)
	{
		for(k=1;k<=nz;k++)
		{
			apu[i][j][k] = 0;
      			apv[i][j][k] = 0;
       			apw[i][j][k] = 0;
       			apT[i][j][k] = 0;
       			apP[i][j][k] = 0;
     
      			su[i][j][k] = 0;
      			sv[i][j][k] = 0;		
			sw[i][j][k] = 0;
			sT[i][j][k] = 0;
			sP[i][j][k] = 0;
			
			//ap[i][j][k] = 0;		
			
		}
      
    }
  }
	

}







