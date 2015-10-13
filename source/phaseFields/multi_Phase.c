void multiPhase()
{

	if(nFluids>1)
	{
	
	
	        pseudoCycle();
	
		createFields();
		
		if(curvature==1||curvature==3||curvature==2)		spSmoothing();		//Seven point smoothing and curvature calculation
		
		if(curvature==2)					k8Convo();		//Kernel Convolution with K8 and curvature calculation 
		
		if(curvature==3) 					heightFunction();	//second order height function to calculate curvature
		
		createFields();
	
		surfaceTension();
		
/*		signDisFun();*/
	}
}
