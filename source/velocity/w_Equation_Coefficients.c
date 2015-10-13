#include"velocity.h"
#include"../pv_Coupling/pv_Coupling.h"

//Updates u Equation Coefficients using hybrid scheme
void wEqnCoeff()
{
	


	double rhot,rhob,rhoc;
	
	double Tt,Tb;
	
	double Betat,Betab;
	
	double SB;
	
	double Fted,Ftw,Fntd,Fst,FT,FPdd;
	double DT,DPdd,Dted,Dtw,Dntd,Dst;
	double rhont,rhost,rhoT,rhoP,rhote,rhotw;
	double munt,must,muT,muP,mute,mutw;


	//Constructing convective fluxes along top face
	for(k=bcID[e][f][g]-1;k<=tcID[e][f][g];k++)
	{
		ktop = k+1;
		kbottom = k-1;
		dzpt = zc[ktop]-zc[k];
		fzt = (zf[k] - zc[k])/dzpt;
		fzp = 1.0 - fzt;
		
		if(k>=nzm)	continue; 
		if(k<2)		continue;
		
		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			ieast = i+1;
			iwest = i-1;
			
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
			
				jnorth = j+1;
   				jsouth = j-1;
			
				/*Interpolating viscosity field*/				
				munt = 0.25*(mu[i][j][k]+mu[i][jnorth][k]+mu[i][j][ktop]+mu[i][jnorth][ktop]);
				must = 0.25*(mu[i][j][k]+mu[i][jsouth][k]+mu[i][j][ktop]+mu[i][jsouth][ktop]);			
      				muT  = mu[i][j][ktop];
      				muP  = mu[i][j][k];
      				mute = 0.25*(mu[i][j][k]+mu[ieast][j][k]+mu[ieast][j][ktop]+mu[i][j][ktop]);
				mutw = 0.25*(mu[i][j][k]+mu[iwest][j][k]+mu[iwest][j][ktop]+mu[i][j][ktop]);
      
      				/* Arithmetic Interpolation of density*/
				rhont = 0.25*(rho[i][j][k]+rho[i][jnorth][k]+rho[i][j][ktop]+rho[i][jnorth][ktop]);
				rhost = 0.25*(rho[i][j][k]+rho[i][jsouth][k]+rho[i][j][ktop]+rho[i][jsouth][ktop]);			
      				rhoT  = rho[i][j][ktop];
      				rhoP  = rho[i][j][k];
      				rhote = 0.25*(rho[i][j][k]+rho[ieast][j][k]+rho[ieast][j][ktop]+rho[i][j][ktop]);
				rhotw = 0.25*(rho[i][j][k]+rho[iwest][j][k]+rho[iwest][j][ktop]+rho[i][j][ktop]);
      
      				//Convective Fluxes
				Fted = rhote * (0.5*u[i][j][k+1][l] + 0.5*u[i][j][k][l])*dy*dz; //weighted mean of ue and uTe
				Ftw  = rhotw * (0.5*u[i-1][j][k+1][l] + 0.5*u[i-1][j][k][l])*dy*dz; //weighted mean of uw and uTw
				Fntd = rhont * (0.5*v[i][j][k+1][l] + 0.5*v[i][j][k][l])*dx*dz; //weighted mean of vn and vnT
				Fst  = rhost * (0.5*v[i][j-1][k+1][l] + 0.5*v[i][j-1][k][l])*dx*dz; //weighted mean of vs and vsT
				FT   = rhoT * (0.5*w[i][j][k+1][l] + 0.5*w[i][j][k][l])*dx*dy; //weighted mean of wt and wtt
				FPdd = rhoP * (0.5*w[i][j][k-1][l] + 0.5*w[i][j][k][l])*dx*dy; //weighted mean of wt and wb	

				//Diffusion Fluxes
	 			DT   =     muT * dx*dy/dz; //mu at T
	 			DPdd =     muP * dx*dy/dz; //mu at P
	 			Dted =     mute * dy*dz/dx; //mu at te
	 			Dtw  =     mutw * dy*dz/dx; //mu at tw
	 			Dntd =     munt * dx*dz/dy; //mu at nt
	 			Dst  =     must * dx*dz/dy; //mu at st
	 			
	 			if((i==2)||(i==(nx))||(j==2)||(j==(ny))) 			DT = DPdd = 0;
	 			
	 			//Coefficients of Equation
	 			att[i][j][k] = DT   + std::max(-FT,0.0);	//upwind
	 			ab[i][j][k]  = DPdd + std::max(FPdd,0.0);	//upwind
	 			atE[i][j][k] = Dted   + std::max(-Fted,0.0);				//central differencing
	 			atW[i][j][k] = Dtw   + std::max(Ftw,0.0);				//central differencing
	 			aNt[i][j][k] = Dntd   + std::max(-Fntd,0.0);				//central differencing
	 			aSt[i][j][k] = Dst   + std::max(Fst,0.0);				//central differencing				
	
				fwuds = w[i][j][k][l] * std::max(FT,0.0) - w[i][j][k+1][l] * std::max(-FT,0.0) - w[i][j][k-1][l] * std::max(FPdd,0.0) + w[i][j][k][l] * std::max(-FPdd,0.0);

	 			fwcds = FT*(w[i][j][k][l] + w[i][j][k+1][l])/2 - FPdd*(w[i][j][k][l] + w[i][j][k-1][l])/2;

	 			sw[i][j][k] = DCvalue*(fwuds-fwcds);
						
				
			}
		}
	}
	
	
	//Constructing Source Terms
	for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		dx = xf[i] - xf[i-1];
		ieast 	= i + 1;
      		iwest 	= i - 1;


		for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			dy = yf[j] - yf[j-1];
			jnorth	= j + 1;
      			jsouth	= j - 1;

			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				dz = zf[k] - zf[k-1];
				ktop  	= k + 1;
		  		kbottom	= k - 1;
				
				
      				/* Updating source terms with pressure gradient */
				sw[i][j][k] = sw[i][j][k] + (p[i][j][k][l] - p[i][j][ktop][l])*dx*dy;
								
				
			}
		}
	}
	
	
//	Adding unsteady contribution to source terms
	for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
			
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				rhot = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][ktop];
				
				/* Harmonic Interpolation of density*/
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/
				
				apt	= rhot*vol/dt;
				
				if(simulationTime==dt)
				{
					sw[i][j][k]	= sw[i][j][k] + apt*w[i][j][k][l-1];
					apw[i][j][k]	= apw[i][j][k] + apt;
				}
				else
				{
					sw[i][j][k]	= sw[i][j][k] + (1.0 + Gamma)*apt*w[i][j][k][l-1] - 0.5*Gamma*apt*w[i][j][k][l-2];
					apw[i][j][k]	= apw[i][j][k] + (1.0 + 0.5*Gamma)*apt;	
				}			
			}
		}
	}
	

	//Adding gravity contribution to source terms
	for(i=2;i<=nxm;i++)
	{
		for(j=2;j<=nym;j++)
		{
			for(k=2;k<=nzm;k++)
			{
				dz = zf[k] - zf[k-1];
				ktop  	= k + 1;
		  		kbottom	= k - 1;
				dzpt = zc[ktop]-zc[k];
				fzt = (zf[k] - zc[k])/dzpt;
				fzp = 1.0 - fzt;
				dzpb = zc[k]-zc[kbottom];
				fzb = (zf[kbottom] - zc[kbottom])/dzpb;
			
			
				rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];
			

				/* Harmonic Interpolation of density*/
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/
				
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				/* Arithmetic Interpolation of Temperature*/
				Tt = fzp * T[i][j][k][l] +   fzt * T[i][j][ktop][l];
				Tb = (1.0-fzb) * T[i][j][k][l] +   fzb * T[i][j][kbottom][l];
				
				/* Arithmetic Interpolation of Temperature*/
				Betat = fzp * Beta[i][j][k] +   fzt * Beta[i][j][ktop];
				Betab = (1.0-fzb) * Beta[i][j][k] +   fzb * Beta[i][j][kbottom];
				
				zBFt[i][j][k][l] = rhot*vol*zGrav - rhot*vol*zTempGrav*Betat*(Tt-Tref);				//Volume effect not yet included
				
				//Adding Surface Tension Body force
				zBFt[i][j][k][l] = zBFt[i][j][k][l] + zSTt[i][j][k]*vol;   
				
				sw[i][j][k]	= sw[i][j][k] + zBFt[i][j][k][l];
			}
		}
	}


}

