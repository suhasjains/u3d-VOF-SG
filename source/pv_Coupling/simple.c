#include"pv_Coupling.h"
#include"../velocity/velocity.h"
#include"../Temperature/Temperature.h"
#include"../Psi/Psi.h"
#include"../update/update.h"


//SIMPLE algorithm loop
void simple()
{

	printf("Starting SIMPLE loop\n");
	
	int stst=0;
	for (simplerVar = 1;simplerVar <= nSimpleLoops;simplerVar++)
 {

		nLoops();	//Initializes no of velocity  and pressure inner loops
		
		boundary_pressure(p);
		
		blockVelocity();

		pressure_correction_east();
		
		pressure_correction_north();
		
		pressure_correction_top();
		
		pressure_correction_source();
		
		ppBoundary_Conditions();
		
		pCorrEqn();		//finding pressure correction based on initial velocities

		boundary_pressure(pp);
				
		update_Pressure_Velocity();       //updating pressure
		
/*		blockTemperature();*/
		
/*		blockPsi();*/
		
		
		if(fabs(maxMassResidual)>massAccuracy)		nSimpleLoops = nSimpleLoops + 1;
		
		printf("Global mass residual = %e \n",maxMassResidual);

		//printf("u = %e\n",u[10][10][4][l]);
		
 }

	printf("	No of SIMPLE loops: %d\n",nSimpleLoops);
	printf("	No of pressure loops: %d\n",nPressureLoops);
	printf("	No of U Velocity loops: %d\n",uVelocityLoops);
	printf("	No of V Velocity loops: %d\n",vVelocityLoops);
	printf("	No of W Velocity loops: %d\n",wVelocityLoops);
/*	printf("	No of T loops: %d\n",nTemperatureLoops);*/
/*	printf("	No of Psi loops: %d\n",nPsiLoops);*/
	printf("	Max of u,v,w,p Normal residuals = %e\n",std::max(std::max(std::max(std::max(uNormRes,vNormRes),wNormRes),pNormRes),TNormRes));
/*	printf("	T Normal residual %e\n",TNormRes);*/
/*	printf("	Psi Normal residual %e\n",PsiNormRes);*/
	printf("	Global mass residual = %e\n",maxMassResidual);

	//check = check + nSimpleLoops;
	//printf("Total loops= %d\n",check);	
}

