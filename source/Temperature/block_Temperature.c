void blockTemperature()
{

	
	b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			TEqnCoeff();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  }
	  


	  TemperatureBoundaryConditions();
	  

	  b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			TEqn();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  }
	  		
	  	 
}
	  			
	  				
		
	

