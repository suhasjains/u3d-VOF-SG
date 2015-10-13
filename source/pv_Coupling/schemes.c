#include"pv_Coupling.h"

void upwind_implicit(double ***flux)
{
  /* Contribution of convective terms to the coefficients */
  convf = std::min(flux[i][j][k],0.0);	/* f - indicates that convective fluxes are evaluated at the face which can be east face or the north face */
  convp = std::max(flux[i][j][k],0.0);	/* p - indicates for the adjacent cell */
  
}

void v_east_source_dc(double f1,double f2,double f3,double f4,double f5,double f6)
{
	/* Constructing source terms for the u v w momentum equation */
      	/* Here, DCvalue denotes the blending factor */ 
      	su[i][j][k]     = su[i][j][k]     + DCvalue*(f1 - f4);
      	su[ieast][j][k] = su[ieast][j][k] - DCvalue*(f1 - f4);
      	sv[i][j][k]     = sv[i][j][k]     + DCvalue*(f2 - f5);
      	sv[ieast][j][k] = sv[ieast][j][k] - DCvalue*(f2 - f5);
      	sw[i][j][k]	= sw[i][j][k]	  + DCvalue*(f3 - f6);
      	sw[ieast][j][k]	= sw[ieast][j][k] - DCvalue*(f3 - f6);
}

void v_north_source_dc(double f1,double f2,double f3,double f4,double f5,double f6)
{
	/* Constructing source terms for the u v w momentum equation */
      	/* through deferred correction approach */
      	/* Here, DCvalue denotes the blending factor */ 
      	su[i][j][k]      = su[i][j][k]      + DCvalue*(f1 - f4);
      	su[i][jnorth][k] = su[i][jnorth][k] - DCvalue*(f1 - f4);
      	sv[i][j][k]      = sv[i][j][k]      + DCvalue*(f2 - f5);
      	sv[i][jnorth][k] = sv[i][jnorth][k] - DCvalue*(f2 - f5);
      	sw[i][j][k]	= sw[i][j][k] 	    + DCvalue*(f3 - f6);
      	sw[i][jnorth][k] = sw[i][jnorth][k] - DCvalue*(f3 - f6); 
}

void v_top_source_dc(double f1,double f2,double f3,double f4,double f5,double f6)
{
	/* Constructing source terms for the u v w momentum equation */
      	/* through deferred correction approach */
      	/* Here, DCvalue denotes the blending factor */ 
      	su[i][j][k]      = su[i][j][k]      + DCvalue*(f1 - f4);
      	su[i][j][ktop] 	 = su[i][j][ktop]   - DCvalue*(f1 - f4);
      	sv[i][j][k]      = sv[i][j][k]      + DCvalue*(f2 - f5);
      	sv[i][j][ktop]   = sv[i][j][ktop]   - DCvalue*(f2 - f5);
      	sw[i][j][k]	 = sw[i][j][k] 	    + DCvalue*(f3 - f6);
      	sw[i][j][ktop]   = sw[i][j][ktop]   - DCvalue*(f3 - f6); 
}

double DC_QUICK_flux_east(double f0,double ****vel)
{
	double del1 = xf[i] - xc[i];
	double del2 = xc[ieast] - xf[i];
	
	double del3;
	if(Fe[i][j][k]==0)	del3 = xc[iwest] - xf[i];
	else			del3 = (xc[iwest] - xf[i])*convp/Fe[i][j][k] + (xc[ieasteast] - xf[i])*convf/Fe[i][j][k];	
	
	double B1 = del2*del3*(del2 - del3);
	double B2 = del1*del3*(del1 + del3);
	double B3 = del1*del2*(del1 + del2);
	double B = B1 - B2 + B3;
	
	if(i==ecID[e][f][g])
	{
		f0 = 0;
      		if(i==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*ueBA[1][j][k][l] + (2*B3/B)*vel[iwest][j][k][l])*convp;	
      		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*ueBA[1][j][k][l] + (B3/B)*vel[iwest][j][k][l])*convp;
		if(i==(nxm-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*ueBA[1][j][k][l] + (2*B3/B)*ueBA[2][j][k][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*ueBA[1][j][k][l] + (B3/B)*ueBA[2][j][k][l])*convf;	
		return f0;
	}
	else if(i==wcID[e][f][g])
	{
		f0 = 0;
      		if(i==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*vel[ieast][j][k][l] + (2*B3/B)*uwBA[1][j][k][l])*convp;	
      		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[ieast][j][k][l] + (B3/B)*uwBA[1][j][k][l])*convp;
		if(i==(nxm-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*vel[ieast][j][k][l] + (2*B3/B)*vel[ieasteast][j][k][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[ieast][j][k][l] + (B3/B)*vel[ieasteast][j][k][l])*convf;	
		return f0;
	}
	else
	{
		f0 = 0;
      		if(i==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*vel[ieast][j][k][l] + (2*B3/B)*vel[iwest][j][k][l])*convp;	
      		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[ieast][j][k][l] + (B3/B)*vel[iwest][j][k][l])*convp;
		if(i==(nxm-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*vel[ieast][j][k][l] + (2*B3/B)*vel[ieasteast][j][k][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[ieast][j][k][l] + (B3/B)*vel[ieasteast][j][k][l])*convf;	
		return f0;
	}
}

double DC_QUICK_flux_north(double f0,double ****vel)
{
	double del1 = yf[j] - yc[j];
	double del2 = yc[jnorth] - yf[j];
	
	double del3;
	if(Fn[i][j][k]==0)	del3 = yc[jsouth] - yf[j];
	else			del3 = (yc[jsouth] - yf[j])*convp/Fn[i][j][k] + (yc[jnorthnorth] - yf[j])*convf/Fn[i][j][k];	
	
	double B1 = del2*del3*(del2 - del3);
	double B2 = del1*del3*(del1 + del3);
	double B3 = del1*del2*(del1 + del2);
	double B = B1 - B2 + B3;
	
	
	if(j==ncID[e][f][g])
	{
		f0 = 0;
		if(j==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*vnBA[i][1][k][l] + (2*B3/B)*vel[i][jsouth][k][l])*convp;
		else 			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vnBA[i][1][k][l] + (B3/B)*vel[i][jsouth][k][l])*convp;
		if(j==(nym-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*vnBA[i][1][k][l] + (2*B3/B)*vnBA[i][2][k][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vnBA[i][1][k][l] + (B3/B)*vnBA[i][2][k][l])*convf;
		return f0;
	}
	else if(j==scID[e][f][g])
	{
		f0 = 0;
		if(j==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*vel[i][jnorth][k][l] + (2*B3/B)*vsBA[i][1][k][l])*convp;
		else 			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][jnorth][k][l] + (B3/B)*vsBA[i][1][k][l])*convp;
		if(j==(nym-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*vel[i][jnorth][k][l] + (2*B3/B)*vel[i][jnorthnorth][k][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][jnorth][k][l] + (B3/B)*vel[i][jnorthnorth][k][l])*convf;
		return f0;
	}
	else
	{
		f0 = 0;
		if(j==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*vel[i][jnorth][k][l] + (2*B3/B)*vel[i][jsouth][k][l])*convp;
		else 			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][jnorth][k][l] + (B3/B)*vel[i][jsouth][k][l])*convp;
		if(j==(nym-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*vel[i][jnorth][k][l] + (2*B3/B)*vel[i][jnorthnorth][k][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][jnorth][k][l] + (B3/B)*vel[i][jnorthnorth][k][l])*convf;
		return f0;
	}
		

}

double DC_QUICK_flux_top(double f0,double ****vel)
{
	double del1 = zf[k] - zc[k];
	double del2 = zc[ktop] - zf[k];
	
	double del3;
	if(Ft[i][j][k]==0)	del3 = zc[kbottom] - zf[k];
	else			del3 = (zc[kbottom] - zf[k])*convp/Ft[i][j][k] + (zc[ktoptop] - zf[k])*convf/Ft[i][j][k];	
	
	double B1 = del2*del3*(del2 - del3);
	double B2 = del1*del3*(del1 + del3);
	double B3 = del1*del2*(del1 + del2);
	double B = B1 - B2 + B3;

	if(k==tcID[e][f][g])
	{
		f0 = 0;
		if(k==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*wtBA[i][j][1][l] + (2*B3/B)*vel[i][j][kbottom][l])*convp;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*wtBA[i][j][1][l] + (B3/B)*vel[i][j][kbottom][l])*convp;
		if(k==(nzm-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*wtBA[i][j][1][l] + (2*B3/B)*wtBA[i][j][2][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*wtBA[i][j][1][l] + (B3/B)*wtBA[i][j][2][l])*convf;
		return f0;
	}
	else if(k==bcID[e][f][g])
	{
		f0 = 0;
		if(k==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*vel[i][j][ktop][l] + (2*B3/B)*wbBA[i][j][1][l])*convp;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][j][ktop][l] + (B3/B)*wbBA[i][j][1][l])*convp;
		if(k==(nzm-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*vel[i][j][ktop][l] + (2*B3/B)*vel[i][j][ktoptop][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][j][ktop][l] + (B3/B)*vel[i][j][ktoptop][l])*convf;
		return f0;
	}
	else
	{
		f0 = 0;
		if(k==2)		f0 = f0 + (((B1-B3)/B)*vel[i][j][k][l] - (B2/B)*vel[i][j][ktop][l] + (2*B3/B)*vel[i][j][kbottom][l])*convp;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][j][ktop][l] + (B3/B)*vel[i][j][kbottom][l])*convp;
		if(k==(nzm-1))		f0 = f0 + ((B1/B)*vel[i][j][k][l] - ((B2+B3)/B)*vel[i][j][ktop][l] + (2*B3/B)*vel[i][j][ktoptop][l])*convf;
		else			f0 = f0 + ((B1/B)*vel[i][j][k][l] - (B2/B)*vel[i][j][ktop][l] + (B3/B)*vel[i][j][ktoptop][l])*convf;
		return f0;
	}
}
