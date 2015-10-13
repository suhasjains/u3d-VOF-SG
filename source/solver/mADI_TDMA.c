#include"solver.h"

//Solving AX=B using Thomas algorithm
void mThomas(int lb,int ub)
{

	m[lb]=1;  //changed principal diagonal of A matrix
	cc[lb]=c[lb];  //changed B matrix

	for(thomasI=(lb+1);thomasI<=ub;thomasI++)
   {

  		m[thomasI]  = 1 - (lower[thomasI] * upper[thomasI-1] / m[thomasI-1]);
  		cc[thomasI] = c[thomasI] - (cc[thomasI-1] * lower[thomasI] / m[thomasI-1]);
   }


	x[ub] = cc[ub] / m[ub];


	for(thomasI=(ub-1);thomasI>=lb;thomasI--)
	{
  	 x[thomasI] = (cc[thomasI] - upper[thomasI] * x[thomasI+1]) / m[thomasI];
	}

}
