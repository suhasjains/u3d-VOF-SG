void syncBufferArrays()
{

	b = 1;

	for(g=1;g<=zNB;g++)
	{
		for (f=1;f<=yNB;f++)
	  	{
	  		for(e=1;e<=xNB;e++)
	  		{
			

				sync();

				b += 1; 	//increment in block number

				
				//printf("u = %e \n",u[28][10][4][l]);
			}
		}
	}
	
	
	
	
	
/*	for  (i=2; i<=nxm; i++)*/
/* 	{*/
/*	 	for (j=2; j<=nym; j++)*/
/*   		{*/
/*	  		for(k=2; k<=nzm; k++)*/
/*			{*/
/*				*/
/*				uB[i][j][k] = u[i][j][k][l];*/
/*				vB[i][j][k] = v[i][j][k][l];*/
/*				wB[i][j][k] = w[i][j][k][l];*/
/*			}*/
/*		}*/
/*	}*/

}



