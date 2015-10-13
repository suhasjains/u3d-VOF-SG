#include"velocity.h"
#include"../pv_Coupling/pv_Coupling.h"

//Updates u Equation Coefficients using hybrid scheme
void uEqnCoeff()
{
	
	
	double rhoe,rhow,rhoc;
	
	double Te,Tw;
	
	double Betae,Betaw;
	
	double SB;
	
	double Fte,Fbe,FE,FP,Fne,Fse;
	double DE,DP,Dne,Dse,Dte,Dbe;
	double rhoE,rhoP,rhobe,rhone,rhose,rhote;
	double muE,muP,mube,mune,muse,mute;
	
	//Constructing convective fluxes along east face 
	for  (i=wcID[e][f][g]-1; i<=ecID[e][f][g]; i++)
 	{
 		dx = xf[i] - xf[i-1];
		ieast = i+1;
		iwest = i-1;
		
		if(i>=nxm)	continue;
		if(i<2)		continue;

	 	for (j=scID[e][f][g]; j<=ncID[e][f][g]; j++)
   		{
   			dy = yf[j] - yf[j-1];
   			jnorth = j+1;
   			jsouth = j-1;
   			
	  		for(k=bcID[e][f][g]; k<=tcID[e][f][g]; k++)
			{
				dz = zf[k] - zf[k-1];
				ktop = k+1;
				kbottom = k-1;
		
		
				/* Evaluation of diffusive term */
				/*Interpolating viscosity field*/	
				mute = 0.25*(mu[i][j][k]+mu[ieast][j][k]+mu[ieast][j][ktop]+mu[i][j][ktop]);
				mube = 0.25*(mu[i][j][k]+mu[ieast][j][k]+mu[ieast][j][kbottom]+mu[i][j][kbottom]);			
      				muE  = mu[ieast][j][k];
      				muP  = mu[i][j][k];
      				mune = 0.25*(mu[i][j][k]+mu[ieast][j][k]+mu[ieast][jnorth][k]+mu[i][jnorth][k]);
				muse = 0.25*(mu[i][j][k]+mu[ieast][j][k]+mu[ieast][jsouth][k]+mu[i][jsouth][k]);	
				
				/* Arithmetic Interpolation of density*/
				rhote = 0.25*(rho[i][j][k]+rho[ieast][j][k]+rho[ieast][j][ktop]+rho[i][j][ktop]);
				rhobe = 0.25*(rho[i][j][k]+rho[ieast][j][k]+rho[ieast][j][kbottom]+rho[i][j][kbottom]);			
      				rhoE  = rho[ieast][j][k];
      				rhoP  = rho[i][j][k];
      				rhone = 0.25*(rho[i][j][k]+rho[ieast][j][k]+rho[ieast][jnorth][k]+rho[i][jnorth][k]);
				rhose = 0.25*(rho[i][j][k]+rho[ieast][j][k]+rho[ieast][jsouth][k]+rho[i][jsouth][k]);
				
      				
      				//Convective Fluxes
	 			Fte = rhote * (0.5*w[i][j][k][l]   + 0.5*w[i+1][j][k][l]) *dx*dy ;  //mean of wt and wtE
				Fbe = rhobe * (0.5*w[i][j][k-1][l] + 0.5*w[i+1][j][k-1][l]) *dx*dy ;  //weighted mean of wb and wbE
				FE  = rhoE * (0.5*u[i][j][k][l]   + 0.5*u[i+1][j][k][l]) *dy*dz ;  //weighted mean of ue and uee
				FP  = rhoP * (0.5*u[i][j][k][l]   + 0.5*u[i-1][j][k][l]) *dy*dz ;  //weighted mean of ue and uw
				Fne = rhone * (0.5*v[i][j][k][l]   + 0.5*v[i+1][j][k][l]) *dx*dz ;  //weighted mean of vn and vnE
				Fse = rhose * (0.5*v[i][j-1][k][l] + 0.5*v[i+1][j-1][k][l]) *dx*dz ;  //weighted mean of vs and vsE
				
				//Diffusion Fluxes
				DE  =	muE * dy * dz / dx; //mu at E
				DP  =	muP * dy * dz / dx; //mu at P
				Dne =   mune * dx * dz / dy; //mu at ne
				Dse =   muse * dx * dz / dy; //mu at se
				Dte =   mute * dx * dy / dz; //mu at te
				Dbe =   mube * dx * dy / dz; //mu at be

				if((i==2)||(i==(nx))||(j==2)||(j==(ny))) 			DE = DP = 0;


				//Coeff of the Equation
				aee[i][j][k] = std::max(-FE,0.0) + DE;  //upwind Scheme
				aw[i][j][k]  = std::max(FP,0.0)  + DP;  //upwind Scheme
				aNe[i][j][k] = std::max(-Fne,0.0)+ Dne;	//central difference
				aSe[i][j][k] = std::max(Fse,0.0) + Dse;	//central difference
				aTe[i][j][k] = std::max(-Fte,0.0)+ Dte;	//central difference
				aBe[i][j][k] = std::max(Fbe,0.0) + Dbe; //central difference
				

				fuuds = u[i][j][k][l] * std::max(FE,0.0) - u[i+1][j][k][l] * std::max(-FE,0.0) - u[i-1][j][k][l] * std::max(FP,0.0) + u[i][j][k][l] * std::max(-FP,0.0);
				fucds = FE*(u[i][j][k][l] + u[i+1][j][k][l])/2 - FP*(u[i][j][k][l] + u[i-1][j][k][l])/2;

				
				su[i][j][k] = DCvalue*(fuuds - fucds);

			}
  		}
	}
	
	

	
	
	//Constructing Source Terms
	for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		dx = xf[i] - xf[i-1];
		ieast 	= i + 1;


		for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			dy = yf[j] - yf[j-1];

			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				dz = zf[k] - zf[k-1];
				
      			
      				/* Updating source terms with pressure gradient */
      				su[i][j][k] = su[i][j][k] + (p[i][j][k][l] - p[ieast][j][k][l])*dy*dz;
      				
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
				
				rhoe = 0.5 * rho[i][j][k] +   0.5 * rho[ieast][j][k];
				
				/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
				
				apt	= rhoe*vol/dt; //rho at e
				
				if(simulationTime==dt)
				{
					su[i][j][k]	= su[i][j][k] + apt*u[i][j][k][l-1];
					apu[i][j][k]	= apt;  
				}
				else
				{
					su[i][j][k]	= su[i][j][k] + (1.0 + Gamma)*apt*u[i][j][k][l-1] - 0.5*Gamma*apt*u[i][j][k][l-2];
					apu[i][j][k]	= (1.0 + 0.5*Gamma)*apt;  
				}			
			}
		}
	}

	

	//Adding gravity contribution to source terms
	for(i=2;i<=nxm;i++)
	{
		dx = xf[i] - xf[i-1];
		ieast 	= i + 1;
      		iwest 	= i - 1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;
		dxpw = xc[i] - xc[iwest];
		fxw = (xf[iwest]-xc[iwest])/dxpw;
	
		for(j=2;j<=nym;j++)
		{
			for(k=2;k<=nzm;k++)
			{
			
			
				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];
				
				/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
				
				
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				/* Arithmetic Interpolation of Temperature*/
				Te = fxp * T[i][j][k][l] +   fxe * T[ieast][j][k][l];
				Tw = (1.0-fxw) * T[i][j][k][l] +   fxw * T[iwest][j][k][l];
				
				/* Arithmetic Interpolation of Temperature*/
				Betae = fxp * Beta[i][j][k] +   fxe * Beta[ieast][j][k];
				Betaw = (1.0-fxw) * Beta[i][j][k] +   fxw * Beta[iwest][j][k];
				
				
				xBFe[i][j][k][l] = rhoe*vol*xGrav  - rhoe*vol*xTempGrav*Betae*(Te-Tref); 
				
				//Adding Surface Tension Body force
				xBFe[i][j][k][l] = xBFe[i][j][k][l] + xSTe[i][j][k]* vol;  		
				
				su[i][j][k]	= su[i][j][k] + xBFe[i][j][k][l];
			}
		}
	}


}

