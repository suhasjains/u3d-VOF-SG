void createFields()
{


	/*Calculate internal density and viscosity fields*/
	for(i=2;i<=nxm;i++)
	{
		for(j=2;j<=nym;j++)
		{
			for(k=2;k<=nzm;k++)
			{
			
/*				rho[i][j][k] = alphaS[i][j][k][l]*rho1 + (1.0 - alphaS[i][j][k][l])*rho2;*/
/*	*/
/*				*/
/*				mu[i][j][k] = alphaS[i][j][k][l]*mu1 + (1.0 - alphaS[i][j][k][l])*mu2; */
/*				*/
/*				Kcond[i][j][k] = alphaS[i][j][k][l]*Kcond1 + (1.0 - alphaS[i][j][k][l])*Kcond2;*/
/*				*/
/*				Beta[i][j][k] = alphaS[i][j][k][l]*Beta1 + (1.0 - alphaS[i][j][k][l])*Beta2;*/
/*				*/
/*				Cp[i][j][k] = alphaS[i][j][k][l]*Cp1 + (1.0 - alphaS[i][j][k][l])*Cp2;*/


                                rho[i][j][k] = Psi[i][j][k][l]*rho1 + (1.0 - Psi[i][j][k][l])*rho2;
	
				
				mu[i][j][k] = Psi[i][j][k][l]*mu1 + (1.0 - Psi[i][j][k][l])*mu2; 
				
				Kcond[i][j][k] = Psi[i][j][k][l]*Kcond1 + (1.0 - Psi[i][j][k][l])*Kcond2;
				
				Beta[i][j][k] = Psi[i][j][k][l]*Beta1 + (1.0 - Psi[i][j][k][l])*Beta2;
				
				Cp[i][j][k] = Psi[i][j][k][l]*Cp1 + (1.0 - Psi[i][j][k][l])*Cp2;




			}
		}
	}

	/*Calculate boundary alpha, density and viscosity fields*/
	for(i=2;i<=nxm;i++)
	{
		for(j=2;j<=nym;j++)
		{
		
			rho[i][j][1] = rho[i][j][2];
			rho[i][j][nz] = rho[i][j][nzm];
			
			mu[i][j][1] = mu[i][j][2];
			mu[i][j][nz] = mu[i][j][nzm];
			
			Kcond[i][j][1] = Kcond[i][j][2];
			Kcond[i][j][nz] = Kcond[i][j][nzm];
			
			Beta[i][j][1] = Beta[i][j][2];
			Beta[i][j][nz] = Beta[i][j][nzm];
			
			Cp[i][j][1] = Cp[i][j][2];
			Cp[i][j][nz] = Cp[i][j][nzm];

			alpha[i][j][1][l] = alpha[i][j][2][l];
			alpha[i][j][nz][l] = alpha[i][j][nzm][l];
			
			alphaS[i][j][1][l] = alphaS[i][j][2][l];
			alphaS[i][j][nz][l] = alphaS[i][j][nzm][l];
			
		}
	}

	for(i=2;i<=nxm;i++)
	{
		for(k=2;k<=nzm;k++)
		{
		
			rho[i][1][k] = rho[i][2][k];
			rho[i][ny][k] = rho[i][nym][k];

			mu[i][1][k] = mu[i][2][k];
			mu[i][ny][k] = mu[i][nym][k];
			
			Kcond[i][1][k] = Kcond[i][2][k];
			Kcond[i][ny][k] = Kcond[i][nym][k];
			
			Beta[i][1][k] = Beta[i][2][k];
			Beta[i][ny][k] = Beta[i][nym][k];
			
			Cp[i][1][k] = Cp[i][2][k];
			Cp[i][ny][k] = Cp[i][nym][k];

			alpha[i][1][k][l] = alpha[i][2][k][l];
			alpha[i][ny][k][l] = alpha[i][nym][k][l];
			
			alphaS[i][1][k][l] = alphaS[i][2][k][l];
			alphaS[i][ny][k][l] = alphaS[i][nym][k][l];	
		}
	}

	for(j=2;j<=nym;j++)
	{
		for(k=2;k<=nzm;k++)
		{
		
			rho[1][j][k] = rho[2][j][k];
			rho[nx][j][k] = rho[nxm][j][k];

			mu[1][j][k] = mu[2][j][k];
			mu[nx][j][k] = mu[nxm][j][k];
			
			Kcond[1][j][k] = Kcond[2][j][k];
			Kcond[nx][j][k] = Kcond[nxm][j][k];
			
			Cp[1][j][k] = Cp[2][j][k];
			Cp[nx][j][k] = Cp[nxm][j][k];
			
			Beta[1][j][k] = Beta[2][j][k];
			Beta[nx][j][k] = Beta[nxm][j][k];

			alpha[1][j][k][l] = alpha[2][j][k][l];
			alpha[nx][j][k][l] = alpha[nxm][j][k][l];
			
			alphaS[1][j][k][l] = alphaS[2][j][k][l];
			alphaS[nx][j][k][l] = alphaS[nxm][j][k][l];
		}
	}


}
