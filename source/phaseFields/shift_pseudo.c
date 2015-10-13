void shiftPseudo()
{

        for  (k=0; k<(zNumberOfCells+5); k++)
	 for (j=0; j<(yNumberOfCells+5); j++)
	  for(i=0; i<(xNumberOfCells+5); i++)
	{
	   
	        tauAlpha[i][j][k][0]=tauAlpha[i][j][k][1];

	        tauAlpha[i][j][k][1]=tauAlpha[i][j][k][2];

	}

}
