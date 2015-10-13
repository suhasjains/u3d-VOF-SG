#include"create_Mesh.h"

void writeMesh()
{
	fp=NULL;
	fp=fopen("./preProcessing/grid_Generation/mesh/mesh.csv","w");  //opening the mesh file

	if(fp!=NULL)
		fprintf(fp,"x Coordinate,y Coordinate,z Coordinate\n");
	else
		printf("fopen ERROR!!\n");


	for  (k=1; k<=nz; k++)
	 for (j=1; j<=ny; j++)
	  for(i=1; i<=nx; i++)
		fprintf(fp,"%e,%e,%e\n",xc[i],yc[j],zc[k]);

	fclose(fp);



}
