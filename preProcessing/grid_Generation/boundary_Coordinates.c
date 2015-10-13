#include"create_Mesh.h"

void boundaryCoordinates()
{
	
	xc[1] = 0;
	xc[nx] = xDomainLength;
	yc[1] = 0;
	yc[ny] = yDomainLength;
	zc[1] = 0;
	zc[nz] = zDomainLength;

	xf[1] = 0;
	yf[1] = 0;
	zf[1] = 0;

	xf[nxm] = xDomainLength;
	yf[nym] = yDomainLength;
	zf[nzm] = zDomainLength;

}
