#include"velocity.h"
#include"../pv_Coupling/pv_Coupling.h"

//Updates u Equation Coefficients using hybrid scheme
void vEqnCoeff()
{
	
	
	double rhon,rhos,rhoc;
	
	double Tn,Ts;
	
	double Betan,Betas;

	double SB;
	
	double Fned,Fnw,FN,FPd,Fnt,Fnb;
	double Dned,Dnw,DN,DPd,Dnt,Dnb;
	double rhont,rhonb,rhoN,rhoP,rhone,rhonw;
	double munt,munb,muN,muP,mune,munw;

	//Constructing convective fluxes along north face
	for(j=scID[e][f][g]-1;j<=ncID[e][f][g];j++)
	{
		jnorth = j+1;
   		jsouth = j-1;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;
	    	fyp = 1.0 - fyn;
	    	
	    	if(j>=nym)	continue;
	    	if(j<2)		continue;

		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			ieast = i+1;
			iwest = i-1;
		
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
			
				ktop = k+1;
				kbottom = k-1;	
      				
				/*Interpolating viscosity field*/
				munt = 0.25*(mu[i][j][k]+mu[i][jnorth][k]+mu[i][j][ktop]+mu[i][jnorth][ktop]);
				munb = 0.25*(mu[i][j][k]+mu[i][jnorth][k]+mu[i][j][kbottom]+mu[i][jnorth][kbottom]);			
      				muN  = mu[i][jnorth][k];
      				muP  = mu[i][j][k];
      				mune = 0.25*(mu[i][j][k]+mu[ieast][j][k]+mu[ieast][jnorth][k]+mu[i][jnorth][k]);
				munw = 0.25*(mu[i][j][k]+mu[iwest][j][k]+mu[iwest][jnorth][k]+mu[i][jnorth][k]);
				
				
				/* Arithmetic Interpolation of density*/
				rhont = 0.25*(rho[i][j][k]+rho[i][jnorth][k]+rho[i][j][ktop]+rho[i][jnorth][ktop]);
				rhonb = 0.25*(rho[i][j][k]+rho[i][jnorth][k]+rho[i][j][kbottom]+rho[i][jnorth][kbottom]);			
      				rhoN  = rho[i][jnorth][k];
      				rhoP  = rho[i][j][k];
      				rhone = 0.25*(rho[i][j][k]+rho[ieast][j][k]+rho[ieast][jnorth][k]+rho[i][jnorth][k]);
				rhonw = 0.25*(rho[i][j][k]+rho[iwest][j][k]+rho[iwest][jnorth][k]+rho[i][jnorth][k]);


				//Convective Fluxes
        			Fned = rhone * (0.5*u[i][j+1][k][l] + 0.5*u[i][j][k][l])*dy*dz; //mean of ue and uNe
				Fnw  = rhonw * (0.5*u[i-1][j+1][k][l] + 0.5*u[i-1][j][k][l])*dy*dz; //weighted mean of uw and uNw
				FN   = rhoN * (0.5*v[i][j+1][k][l] + 0.5*v[i][j][k][l])*dx*dz; //weighted mean of vn and vnn
				FPd  = rhoP * (0.5*v[i][j-1][k][l] + 0.5*v[i][j][k][l])*dx*dz; //weighted mean of vn and vs
				Fnt  = rhont * (0.5*w[i][j+1][k][l] + 0.5*w[i][j][k][l])*dx*dy; //weighted mean of wt and wNt
				Fnb  = rhonb * (0.5*w[i][j+1][k-1][l] + 0.5*w[i][j][k-1][l])*dx*dy; //weighted mean of wb and wNb 
				
				//Diffusion Fluxes
				Dned=     mune * dy * dz/dx; //mu at ne
				Dnw =     munw * dy * dz/dx; //mu at nw
				DN  =     muN * dx * dz/dy; //mu at N
				DPd =     muP * dx * dz/dy; //mu at P
				Dnt =     munt * dx * dy/dz; //mu at nt
				Dnb =     munb * dx * dy/dz; //mu at nb				
				     	
      				if((i==2)||(i==(nx))||(j==2)||(j==(ny))) 			DN = DPd = 0;
      			
      				//Coeff of the Equation
				anE[i][j][k] = Dned  + std::max(-Fned,0.0);				//central differencing
				anW[i][j][k] = Dnw  + std::max(Fnw,0.0);				//central differencing
				ann[i][j][k] = DN  + std::max(-FN,0.0);  	//upwind Scheme
				as[i][j][k]  = DPd + std::max(FPd,0.0);  	//upwind Scheme
				anT[i][j][k] = Dnt  + std::max(-Fnt,0.0);				//central differencing
				anB[i][j][k] = Dnb  + std::max(Fnb,0.0);				//central differencing
      			
      			
      				fvuds = v[i][j][k][l] * std::max(FN,0.0) - v[i][j+1][k][l] * std::max(-FN,0.0) - v[i][j-1][k][l] * std::max(FPd,0.0) + v[i][j][k][l] * std::max(-FPd,0.0);
      				fvcds = FN*(v[i][j][k][l] + v[i][j+1][k][l])/2 - FPd*(v[i][j][k][l] + v[i][j-1][k][l])/2;
      				
      				sv[i][j][k] = DCvalue*(fvuds-fvcds);
      			
      
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
				
				vol	= dx*dy*dz;
				
      				/* Updating source terms with pressure gradient */
      				sv[i][j][k] = sv[i][j][k] + (p[i][j][k][l] - p[i][jnorth][k][l])*dx*dz;
				
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
				
				rhon = ( rho[i][j][k] +  rho[i][jnorth][k])/2.0;
				
				/* Harmonic Interpolation of density*/
/*				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/
				
				apt	= rhon*vol/dt;
				
				if(simulationTime==dt)
				{
					sv[i][j][k]	= sv[i][j][k] + apt*v[i][j][k][l-1]; 
					apv[i][j][k]	= apv[i][j][k] + apt;
				}
				else
				{
					sv[i][j][k]	= sv[i][j][k] + (1.0 + Gamma)*apt*v[i][j][k][l-1] - 0.5*Gamma*apt*v[i][j][k][l-2]; 
					apv[i][j][k]	= apv[i][j][k] + (1.0 + 0.5*Gamma)*apt;
				}			
			}
		}
	}
	

	//Adding gravity contribution to source terms
	for(i=2;i<=nxm;i++)
	{
		for(j=2;j<=nym;j++)
		{
			dy = yf[j] - yf[j-1];
			jnorth	= j + 1;
      			jsouth	= j - 1;
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
			dyps = yc[j] - yc[jsouth];
			fys = (yf[jsouth] - yc[jsouth])/dyps;
		
		
			for(k=2;k<=nzm;k++)
			{
			

				rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];

				/* Harmonic Interpolation of density*/
				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);

				
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				/* Arithmetic Interpolation of Temperature*/
				Tn = fyp * T[i][j][k][l] +   fyn * T[i][jnorth][k][l];
				Ts = (1.0-fys) * T[i][j][k][l] +   fys * T[i][jsouth][k][l];
				
				/* Arithmetic Interpolation of Temperature*/
				Betan = fyp * Beta[i][j][k] +   fyn * Beta[i][jnorth][k];
				Betas = (1.0-fys) * Beta[i][j][k] +   fys * Beta[i][jsouth][k];
				
				
				yBFn[i][j][k][l] = rhon*vol*yGrav - rhon*vol*yTempGrav*Betan*(Tn-Tref);			//Volume effect not yet included
				
/*				printf("yBFn = %e \n",yBFn[i][j][k][l]);*/
				
				//Adding Surface Tension Body force
				yBFn[i][j][k][l] = yBFn[i][j][k][l] + ySTn[i][j][k]*vol;   
				
				sv[i][j][k]	= sv[i][j][k] + yBFn[i][j][k][l]; 
			}
		}
	}


}

