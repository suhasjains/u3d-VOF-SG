#include"solver.h"

//Solving AX=B using Thomas algorithm
void Thomas(int n)
{

	m[2]=1;  //changed principal diagonal of A matrix
	cc[2]=c[2];  //changed B matrix

	for(thomasI=3;thomasI<=n;thomasI++)
   {

  		m[thomasI]  = 1 - (lower[thomasI] * upper[thomasI-1] / m[thomasI-1]);
  		cc[thomasI] = c[thomasI] - (cc[thomasI-1] * lower[thomasI] / m[thomasI-1]);
   }


	x[n] = cc[n] / m[n];


	for(thomasI=(n-1);thomasI>=2;thomasI--)
	{
  	 x[thomasI] = (cc[thomasI] - upper[thomasI] * x[thomasI+1]) / m[thomasI];
	}

}
