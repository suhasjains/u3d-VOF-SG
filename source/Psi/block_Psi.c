void blockPsi()
{

	
	b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			PsiEqnCoeff();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  }
	  


	  PsiBoundaryConditions();
	  

	  b = 1;
	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
	  		
	  			PsiEqn();
	  			
	  			b += 1;
	  			
	  		}
	  	}
	  }
	  		
	  	 
}
	  			
	  				
		
	

