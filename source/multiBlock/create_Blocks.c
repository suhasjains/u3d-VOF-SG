void createBlocks()
{

	//Assigning starting and ending cell indexes
	for(e=1;e<=xNB;e++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(g=1;g<=zNB;g++)
	  		{
	  			if(e==1)
	  			{
	  				wcID[e][f][g] = 2;
					ecID[e][f][g] = (nxm-1)/xNB + 1;
				}
				else
				{
					wcID[e][f][g] = ecID[e-1][f][g] + 1;
					ecID[e][f][g] = ecID[e-1][f][g] + (nxm-1)/xNB;
				}
				
				if(f==1)
				{
					scID[e][f][g] = 2;
					ncID[e][f][g] = (nym-1)/yNB + 1;
				}
				else
				{
					scID[e][f][g] = ncID[e][f-1][g] + 1;
					ncID[e][f][g] = ncID[e][f-1][g] + (nym-1)/yNB;
				}
				
				if(g==1)
	  			{
	  				bcID[e][f][g] = 2;
					tcID[e][f][g] = (nzm-1)/zNB + 1;
				}
				else
				{
					bcID[e][f][g] = tcID[e][f][g-1] + 1;
					tcID[e][f][g] = tcID[e][f][g-1] + (nzm-1)/zNB;
				}
				
				//printf("wcID[%d][%d][%d] = %d \n",e,f,g,wcID[e][f][g]);
				//printf("ecID[%d][%d][%d] = %d \n",e,f,g,ecID[e][f][g]);
				//printf("ncID[%d][%d][%d] = %d \n",e,f,g,ncID[e][f][g]);
				//printf("scID[%d][%d][%d] = %d \n",e,f,g,scID[e][f][g]);
				//printf("tcID[%d][%d][%d] = %d \n",e,f,g,tcID[e][f][g]);
				//printf("bcID[%d][%d][%d] = %d \n",e,f,g,bcID[e][f][g]);
			}
		}
	}


	nBlocks = xNB*yNB*zNB;	//No of Blocks
}



