void spSmoothing()
{

	for(i=2;i<=nxm;i++)	//Calculation of components of Surface Tension
	{
		ieast 	= i + 1;
		iwest 	= i-1;

    		for(j=2;j<=nym;j++)
		{
			jnorth	= j + 1;
			jsouth = j-1;
	
	  		for(k=2;k<=nzm;k++)
			{
				ktop  	= k + 1;
				kbottom = k-1;

				alphaS[i][j][k][l] = 0.5*alpha[i][j][k][l];
				alphaS[i][j][k][l] = alphaS[i][j][k][l] + (alpha[ieast][j][k][l]/12.0);
				alphaS[i][j][k][l] = alphaS[i][j][k][l] + (alpha[iwest][j][k][l]/12.0);
				alphaS[i][j][k][l] = alphaS[i][j][k][l] + (alpha[i][jnorth][k][l]/12.0);
				alphaS[i][j][k][l] = alphaS[i][j][k][l] + (alpha[i][jsouth][k][l]/12.0);
				alphaS[i][j][k][l] = alphaS[i][j][k][l] + (alpha[i][j][ktop][l]/12.0);
				alphaS[i][j][k][l] = alphaS[i][j][k][l] + (alpha[i][j][kbottom][l]/12.0);
				
			}
		}
	}
	
	
	curvatureSmooth();

}
