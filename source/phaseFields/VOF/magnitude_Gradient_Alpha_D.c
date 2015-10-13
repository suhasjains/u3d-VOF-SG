double magGradAlphaD()
{

	double A;
	
	dx = xf[i]-xf[i-1];
	dy = yf[j]-yf[j-1];
	dz = zf[k]-zf[k-1];

	ieast = i+1;
	jnorth = j+1;
	ktop = k+1;

				A = 0;
	
	if(curvature==1||curvature==3)
	{
		if(Fe[i][j][k]>=0)	A = A + pow((eAlpha[i][j][k][l]	 - wAlpha[i][j][k][l])/dx,2);
		else			A = A + pow((eAlpha[ieast][j][k][l] - wAlpha[ieast][j][k][l])/(xf[ieast]-xf[i]),2);
		if(Fn[i][j][k]>=0)	A = A + pow((nAlpha[i][j][k][l]	 - sAlpha[i][j][k][l])/dy,2);
		else			A = A + pow((nAlpha[i][jnorth][k][l]- sAlpha[i][jnorth][k][l])/(yf[jnorth]-yf[j]),2);
		if(Ft[i][j][k]>=0)	A = A + pow((tAlpha[i][j][k][l]	 - bAlpha[i][j][k][l])/dz,2);
		else 			A = A + pow((tAlpha[i][j][ktop][l]	 - bAlpha[i][j][ktop][l])/(zf[ktop]-zf[k]),2);
	}


	if(curvature==2)
	{
	
		if(Fe[i][j][k]>=0)	A = A + pow(xGradAlphaS[i][j][k][l],2);
		else			A = A + pow(xGradAlphaS[ieast][j][k][l],2);
		if(Fn[i][j][k]>=0)	A = A + pow(yGradAlphaS[i][j][k][l],2);
		else			A = A + pow(yGradAlphaS[i][jnorth][k][l],2);
		if(Ft[i][j][k]>=0)	A = A + pow(zGradAlphaS[i][j][k][l],2);
		else 			A = A + pow(zGradAlphaS[i][j][ktop][l],2);
	}

				A = pow(A,0.5);

	return A;

}
