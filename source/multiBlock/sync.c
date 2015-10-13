void sync()
{

				//Along x
				for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
				{
					for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
					{


						if(ecID[e][f][g]!= nxm)	
						{
					
							ueBA[1][j][k][b] = u[ecID[e][f][g] + 1][j][k][l];
							//printf("ueBA[1][%d][%d] = %e \n",j,k,ueBA[1][j][k][b]);
							veBA[1][j][k][b] = v[ecID[e][f][g] + 1][j][k][l];
							weBA[1][j][k][b] = w[ecID[e][f][g] + 1][j][k][l];
							ppeBA[1][j][k][b] = pp[ecID[e][f][g] + 1][j][k][l];
							alphaeBA[1][j][k][b] = alpha[ecID[e][f][g] + 1][j][k][l];
		
							ueBA[2][j][k][b] 	  = u[ecID[e][f][g] + 2][j][k][l];
							//printf("ueBA[2][%d][%d][%d] = %e \n",j,k,b,ueBA[2][j][k][b]);
							veBA[2][j][k][b] 	  = v[ecID[e][f][g] + 2][j][k][l];
							weBA[2][j][k][b] 	  = w[ecID[e][f][g] + 2][j][k][l];
							ppeBA[2][j][k][b]	  = pp[ecID[e][f][g] + 2][j][k][l];
							alphaeBA[2][j][k][b] = alpha[ecID[e][f][g] + 2][j][k][l];
						}

						if(wcID[e][f][g]!= 2)	
						{
					
							uwBA[1][j][k][b]     = u[wcID[e][f][g] - 1][j][k][l];
							vwBA[1][j][k][b]     = v[wcID[e][f][g] - 1][j][k][l];
							wwBA[1][j][k][b]     = w[wcID[e][f][g] - 1][j][k][l];
							ppwBA[1][j][k][b]    = pp[wcID[e][f][g] - 1][j][k][l];
							alphawBA[1][j][k][b] = alpha[wcID[e][f][g] - 1][j][k][l];
		
							uwBA[2][j][k][b] 	  = u[wcID[e][f][g] - 2][j][k][l];
							vwBA[2][j][k][b] 	  = v[wcID[e][f][g] - 2][j][k][l];
							wwBA[2][j][k][b] 	  = w[wcID[e][f][g] - 2][j][k][l];
							ppwBA[2][j][k][b]	  = pp[wcID[e][f][g] - 2][j][k][l];
							alphawBA[2][j][k][b] 	  = alpha[wcID[e][f][g] - 2][j][k][l];
						}
					}
				}



				//Along y
				for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
				{
					for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
					{


						if(ncID[e][f][g]!= nym)	
						{
					
							unBA[i][1][k][b] 	= u[i][ncID[e][f][g] + 1][k][l];
							vnBA[i][1][k][b] 	= v[i][ncID[e][f][g] + 1][k][l];
							wnBA[i][1][k][b] 	= w[i][ncID[e][f][g] + 1][k][l];
							ppnBA[i][1][k][b] 	= pp[i][ncID[e][f][g] + 1][k][l];
							alphanBA[i][1][k][b] 	= alpha[i][ncID[e][f][g] + 1][k][l];
		
							unBA[i][2][k][b] 	  = u[i][ncID[e][f][g] + 2][k][l];
							vnBA[i][2][k][b] 	  = v[i][ncID[e][f][g] + 2][k][l];
							wnBA[i][2][k][b] 	  = w[i][ncID[e][f][g] + 2][k][l];
							ppnBA[i][2][k][b]	  = pp[i][ncID[e][f][g] + 2][k][l];
							alphanBA[i][2][k][b] 	  = alpha[i][ncID[e][f][g] + 2][k][l];
						}

						if(scID[e][f][g]!= 2)	
						{
					
							usBA[i][1][k][b]     = u[i][scID[e][f][g] - 1][k][l];
							vsBA[i][1][k][b]     = v[i][scID[e][f][g] - 1][k][l];
							wsBA[i][1][k][b]     = w[i][scID[e][f][g] - 1][k][l];
							ppsBA[i][1][k][b]    = pp[i][scID[e][f][g] - 1][k][l];
							alphasBA[i][1][k][b] = alpha[i][scID[e][f][g] - 1][k][l];
		
							usBA[i][2][k][b] 	  = u[i][scID[e][f][g] - 2][k][l];
							vsBA[i][2][k][b] 	  = v[i][scID[e][f][g] - 2][k][l];
							wsBA[i][2][k][b] 	  = w[i][scID[e][f][g] - 2][k][l];
							ppsBA[i][2][k][b]	  = pp[i][scID[e][f][g] - 2][k][l];
							alphasBA[i][2][k][b] = alpha[i][scID[e][f][g] - 2][k][l];
						}
					}
				}

				//Along z
				for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
				{
					for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
					{


						if(tcID[e][f][g]!= nzm)	
						{
					
							utBA[i][j][1][b] = u[i][j][tcID[e][f][g] + 1][l];
							vtBA[i][j][1][b] = v[i][j][tcID[e][f][g] + 1][l];
							wtBA[i][j][1][b] = w[i][j][tcID[e][f][g] + 1][l];
							pptBA[i][j][1][b] = pp[i][j][tcID[e][f][g] + 1][l];
							alphatBA[i][j][1][b] = alpha[i][j][tcID[e][f][g] + 1][l];
		
							utBA[i][j][2][b] 	  = u[i][j][tcID[e][f][g] + 2][l];
							vtBA[i][j][2][b] 	  = v[i][j][tcID[e][f][g] + 2][l];
							wtBA[i][j][2][b] 	  = w[i][j][tcID[e][f][g] + 2][l];
							pptBA[i][j][2][b]	  = pp[i][j][tcID[e][f][g] + 2][l];
							alphatBA[i][j][2][b] = alpha[i][j][tcID[e][f][g] + 2][l];
						}

						if(bcID[e][f][g]!= 2)	
						{
					
							ubBA[i][j][1][b]     = u[i][j][bcID[e][f][g] - 1][l];
							vbBA[i][j][1][b]     = v[i][j][bcID[e][f][g] - 1][l];
							wbBA[i][j][1][b]     = w[i][j][bcID[e][f][g] - 1][l];
							ppbBA[i][j][1][b]    = pp[i][j][bcID[e][f][g] - 1][l];
							alphabBA[i][j][1][b] = alpha[i][j][bcID[e][f][g] - 1][l];
		
							ubBA[i][j][2][b] 	  = u[i][j][bcID[e][f][g] - 2][l];
							vbBA[i][j][2][b] 	  = v[i][j][bcID[e][f][g] - 2][l];
							wbBA[i][j][2][b] 	  = w[i][j][bcID[e][f][g] - 2][l];
							ppbBA[i][j][2][b]	  = pp[i][j][bcID[e][f][g] - 2][l];
							alphabBA[i][j][2][b] = alpha[i][j][bcID[e][f][g] - 2][l];
						}
					}
				}

				//printf("ueBA = %e \n",ueBA[1][10][4][b]);







}
