void blockVelocity()
{

		
		
	velocityBoundaryConditions();
	
		
	b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			uEqnCoeff();
	  			
	  			vEqnCoeff();
	  			
	  			wEqnCoeff();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  }
	  
	  b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			uEqn();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  }
	  		
	    b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			vEqn();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  }			
	 
	 
	 b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			wEqn();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  } 
	  			
	  			
	  	 
}
	  			
	  				
		
	

