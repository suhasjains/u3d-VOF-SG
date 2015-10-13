#include"create_Mesh.h"

void uniformGrid()
{
	double const_stretch_factor = 0.5;


	for(i=1;i<=nx;i++)
		fx[i] = const_stretch_factor;

	for(j=1;j<=ny;j++) 
		fy[j] = const_stretch_factor;

	for(k=1;k<=nz;k++) 
		fz[k] = const_stretch_factor;


	/* Generating xc, yc and zc */

	for(i=1;i<=nx;i++) 
	{
		if(i==1)		xc[i] = 0;
		else if(i==2||i==nx)	xc[i] = xc[i-1] + 0.5*dx;
		else			xc[i] = xc[i-1] + dx;
	}

	for(j=1;j<=ny;j++) 
	{
		if(j==1)		yc[j] = 0;
		else if(j==2||j==ny)	yc[j] = yc[j-1] + 0.5*dy;
		else			yc[j] = yc[j-1] + dy;
	}

	for(k=1;k<=nz;k++) 
	{
		if(k==1)		zc[k] = 0;
		else if(k==2||k==nz)	zc[k] = zc[k-1] + 0.5*dz;
		else			zc[k] = zc[k-1] + dz;
	}

	/* Creating xf, yf and zf */

	for(i=2;i<=nx-2;i++)	xf[i] = xc[i] + 0.5*dx;
	for(j=2;j<=ny-2;j++) 	yf[j] = yc[j] + 0.5*dy;
	for(k=2;k<=nz-2;k++) 	zf[k] = zc[k] + 0.5*dz;

	

}
