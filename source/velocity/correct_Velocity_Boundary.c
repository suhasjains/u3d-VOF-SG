 //Velocity correction at boundary due to specified pressure boundary conditions
 
void correct_WestVelBoundary()
{

	double apue,rhoe,vole,s,d;

	//Along west	
	i=1;
	for (j=2;j<=nym;j++)
	{
	  	for(k=2;k<=nzm;k++)
	  	{
	  	
		  	
	  		ieast = i+1;
			iwest = i-1;
			ieasteast = i+2;
			dxpe = xc[ieast] - xc[i];
			fxe = (xf[i]-xc[i])/dxpe;		
			fxp = 1.0 - fxe;
	  	
	  	
	  		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
			vole = dxpe * s;
	     			
				/* Arithmetic Interpolation of density*/
/*				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];*/
				
				/* Harmonic Interpolation of density*/
				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);
			d = rhoe*s;
	     			
			apue  = 1.5*apu[ieast][j][k]  - 0.5*apu[ieasteast][j][k];  
	    			
			ue[i][j][k][l] = ue[i][j][k][l] - vole*apue*(pp[ieast][j][k][l]-pp[i][j][k][l])/dxpe;
	     			
			Fe[i][j][k] = d*ue[i][j][k][l];
	      			
	      	}
	      	
	}


}


void correct_EastVelBoundary()
{

	double apue,rhoe,vole,s,d;

	//Along east
	i=nxm;
	for (j=2;j<=nym;j++)
	{
	  	for(k=2;k<=nzm;k++)
	  	{
		  	
	  		ieast = i+1;
			iwest = i-1;
			ieasteast = i+2;
			dxpe = xc[ieast] - xc[i];
			fxe = (xf[i]-xc[i])/dxpe;		
			fxp = 1.0 - fxe;
		  	
		  	
	  		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
			vole = dxpe * s;
	      			
			/* Arithmetic Interpolation of density*/
/*				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];*/
				
				/* Harmonic Interpolation of density*/
				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);
			d = rhoe*s;
	     			
			apue  = 1.5*apu[i][j][k]  - 0.5*apu[iwest][j][k];  
	     			
			ue[i][j][k][l] = ue[i][j][k][l] - vole*apue*(pp[ieast][j][k][l]-pp[i][j][k][l])/dxpe;
	     			
			//printf("%e \n",pp[ieast][j][k][l]-pp[i][j][k][l]);
	     			
			Fe[i][j][k] = d*ue[i][j][k][l];
	      			
		}
	      	
	}


}

void correct_NorthVelBoundary()
{

	double apvn,rhon,voln,s,d;

	//Along North
	j=nym;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{

			jnorth = j+1;
			jnorthnorth = j+2;
			jsouth = j-1;
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
		
			s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
			voln = s * dypn;
			
			/* Arithmetic Interpolation of density*/
/*			rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];*/

			/* Harmonic Interpolation of density*/
			rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);
			d = rhon * s;
			
			apvn  = 1.5 * apv[i][j][k]  - 0.5 * apv[i][jsouth][k];
			
			vn[i][j][k][l] = vn[i][j][k][l] - voln*apvn*(pp[i][jnorth][k][l]-pp[i][j][k][l])/dypn;
			
			Fn[i][j][k] = d*vn[i][j][k][l];
			
		}
	}

}

void correct_SouthVelBoundary()
{

	double apvn,rhon,voln,s,d;

	//Along South
	j=1;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{

			jnorth = j+1;
			jnorthnorth = j+2;
			jsouth = j-1;
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
		
			s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
			voln = s * dypn;
			
			/* Arithmetic Interpolation of density*/
/*			rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];*/

			/* Harmonic Interpolation of density*/
			rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);
			d = rhon * s;
			
			apvn  = 1.5 * apv[i][jnorth][k]  - 0.5 * apv[i][jnorthnorth][k];
			
			vn[i][j][k][l] = vn[i][j][k][l] - voln*apvn*(pp[i][jnorth][k][l]-pp[i][j][k][l])/dypn;
			
			Fn[i][j][k] = d*vn[i][j][k][l];
			
		}
	}

}

void correct_TopVelBoundary()
{

	double apwt,rhot,volt,s,d;

	//Along Top
	k=nzm;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
	  	
	  		ktop = k+1;
			ktoptop = k+2;
			kbottom = k-1;
			dzpt = zc[ktop]-zc[k];
			fzt = (zf[k] - zc[k])/dzpt;
			fzp = 1.0 - fzt;
			
			s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
			volt = s * dzpt;
			
			/* Arithmetic Interpolation of density*/
/*			rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];*/

			/* Harmonic Interpolation of density*/
			rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);
			d = rhot * s; 
	
			apwt  = 1.5 * apw[i][j][k]  - 0.5 * apw[i][j][kbottom];
			
			wt[i][j][k][l] = wt[i][j][k][l] - volt*apwt*(pp[i][j][ktop][l]-pp[i][j][k][l])/dzpt;
			
			Ft[i][j][k] = d*wt[i][j][k][l];
			
		}
	}
}

void correct_BottomVelBoundary()
{

	double apwt,rhot,volt,s,d;


	//Along Bottom
	k=1;
  	for(i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
	  	
	  		ktop = k+1;
			ktoptop = k+2;
			kbottom = k-1;
			dzpt = zc[ktop]-zc[k];
			fzt = (zf[k] - zc[k])/dzpt;
			fzp = 1.0 - fzt;
			
			s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
			volt = s * dzpt;
			
			/* Arithmetic Interpolation of density*/
/*			rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];*/

			/* Harmonic Interpolation of density*/
			rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);
			d = rhot * s; 
	
			apwt  = 1.5 * apw[i][j][ktop]  - 0.5 * apw[i][j][ktoptop];
			
			wt[i][j][k][l] = wt[i][j][k][l] - volt*apwt*(pp[i][j][ktop][l]-pp[i][j][k][l])/dzpt;
			
			Ft[i][j][k] = d*wt[i][j][k][l];
			
		}
	}


}


