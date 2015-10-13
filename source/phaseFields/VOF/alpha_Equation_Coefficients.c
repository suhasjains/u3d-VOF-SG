void alphaEqnCoeff()
{

	double N,D;	//Numerator, Denominator
	double Theta,Psi;
	
	int q;
	double r;

	double rhoe,rhon,rhot;
	double rhow,rhos,rhob;
	
	double ue,uw,un,us,ut,ub;
	double ve,vw,vn,vs,vt,vb;
	double we,ww,wn,ws,wt,wb;
	


	//Calcualtion of local Courant Number
	localCourantNumber();


	

	//Constructing alpha fluxes for faces along x direction
	for  (i=2; i<=nxm-1; i++)
 	{
		ieast 	= i+1;
		iwest 	= i-1;
		ieasteast = i+2;
		dx = xf[i]-xf[iwest];
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;

	 	for (j=2; j<=nym; j++)
   		{
	  		for(k=2; k<=nzm; k++)
			{
				

				/* Arithmetic Interpolation of density*/
				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];


				/* Harmonic Interpolation of density*/
/*				rhoe = 2.0 * rho[i][j][k] * rho[ieast][j][k]/( rho[i][j][k] + rho[ieast][j][k]);*/
	
				if(Fe[i][j][k]>=0)
				{

					alphaD = tauAlpha[i][j][k][l-1];
					alphaA = tauAlpha[ieast][j][k][l-1];
					alphaU = tauAlpha[iwest][j][k][l-1];
				}
				else			
				{
					alphaD = tauAlpha[ieast][j][k][l-1];
					alphaA = tauAlpha[i][j][k][l-1];
					alphaU = tauAlpha[ieasteast][j][k][l-1];
				}

				N = alphaA - alphaU;
				D = alphaD - alphaU;
				
				if(D==0)	rS = 1.0;
				else 		rS = N/D;
				
				kS = N*(1.0-Ce[i][j][k])*3.0/8.0 + D*(3.0+Ce[i][j][k])/4.0 + alphaU;


				//VOF scheme
				HiRAC();
									
				
				//Calculation of angle between interface normal and cell face
				//Numerator
				if(curvature==1||curvature==3)	
				{
					if(Fe[i][j][k]>=0)	N = (eAlpha[i][j][k][l]	- wAlpha[i][j][k][l])/dx;
					else			N = (eAlpha[ieast][j][k][l]- wAlpha[ieast][j][k][l])/(xf[ieast]-xf[i]);
				}
	
				if(curvature==2)
				{
											
					if(Fe[i][j][k]>=0)	N = xGradAlphaS[i][j][k][l];
					else			N = xGradAlphaS[ieast][j][k][l];
				}

				//Denominator
				D = magGradAlphaD();

				if(D==0)	D=1e-100;
	
				r = N/D;
					
				if(r>1.0)	r=1.0;
				if(r<-1.0)	r=-1.0;
	
	
				//Blending Factor
				Psi = BF_HiRAC(r);
					
					
				//Calculation of Alpha at Cell face
				eAlpha[i][j][k][l] = Psi*alphaBD + (1.0-Psi)*alphaHR;
                                
                                if(eAlpha[i][j][k][l]>1)        eAlpha[i][j][k][l] = 1.0;
                                if(eAlpha[i][j][k][l]<0)        eAlpha[i][j][k][l] = 0;
                                 
                                
				wAlpha[ieast][j][k][l] = eAlpha[i][j][k][l];
				
/*				printf("eAlpha = %e\n",eAlpha[i][j][k][l]);*/
                                	
	                        
	                        ue = u[i][j][k][l];
	                        ve = 0.25*(v[i][j][k][l] + v[ieast][j][k][l] + v[i][jsouth][k][l] + v[ieast][jsouth][k][l]);
	                        we = 0.25*(w[i][j][k][l] + w[ieast][j][k][l] + w[i][j][kbottom][l] + w[ieast][j][kbottom][l]);
	                        
	                        uw = u[iwest][j][k][l];
	                        vw = 0.25*(v[i][j][k][l] + v[iwest][j][k][l] + v[i][jsouth][k][l] + v[iwest][jsouth][k][l]);
	                        ww = 0.25*(w[i][j][k][l] + w[iwest][j][k][l] + w[i][j][kbottom][l] + w[iwest][j][kbottom][l]);
	                        
	                        magVe[i][j][k] = pow((pow(ue,2)+pow(ve,2)+pow(we,2)),0.5); 
	                        magVw[i][j][k] = pow((pow(uw,2)+pow(vw,2)+pow(ww,2)),0.5); 
	
				//printf("nAlphaD=%e,eAlpha = %e,%d,%d,%d \n",nAlphaD,eAlpha[i][j][k],i,j,k);  


				
					
			}	
		}	
	}


	//Constructing alpha fluxes for faces along y direction
	for(j=2;j<=nym-1;j++)
	{
		jnorth = j+1;
		jsouth = j-1;
		jnorthnorth = j+2;
		dy = yf[j]-yf[jsouth];
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;	
    		fyp = 1.0 - fyn;

		for(i=2;i<=nxm;i++)
		{
			for(k=2;k<=nzm;k++)
			{
				
				/* Arithmetic Interpolation of density*/
				rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];


				/* Harmonic Interpolation of density*/
/*				rhon = 2.0 * rho[i][j][k] * rho[i][jnorth][k]/( rho[i][j][k] + rho[i][jnorth][k]);*/
				
				if(Fn[i][j][k]>=0)
				{

					alphaD = tauAlpha[i][j][k][l-1];
					alphaA = tauAlpha[i][jnorth][k][l-1];
					alphaU = tauAlpha[i][jsouth][k][l-1];
				}
				else			
				{
					alphaD = tauAlpha[i][jnorth][k][l-1];
					alphaA = tauAlpha[i][j][k][l-1];
					alphaU = tauAlpha[i][jnorthnorth][k][l-1];
				}

				N = alphaA - alphaU;
				D = alphaD - alphaU;
				
				if(D==0)	rS = 1.0;
				else 		rS = N/D;
				
				kS = N*(1.0-Cn[i][j][k])*3.0/8.0 + D*(3.0+Cn[i][j][k])/4.0 + alphaU;


				//VOF scheme
				HiRAC();	
						
				//Calculation of angle between interface normal and cell face
				//Numerator
				if(curvature==1||curvature==3)
				{
					if(Fn[i][j][k]>=0)	N = (nAlpha[i][j][k][l]	- sAlpha[i][j][k][l])/dy;
					else			N = (nAlpha[i][jnorth][k][l]- sAlpha[i][jnorth][k][l])/(yf[jnorth]-yf[j]);
				}
						
				if(curvature==2)
				{
					if(Fn[i][j][k]>=0)	N = yGradAlphaS[i][j][k][l];
					else			N = yGradAlphaS[i][jnorth][k][l];
				}
						
				//printf("N=%e \n",N);
			
				//Denominator
				D = magGradAlphaD();

				//printf("D=%e \n",D);
				if(D==0)	D=1e-100;
	
				r = N/D;
/*				printf("r = %lf \n",r);*/
					
				if(r>1.0)	r=1.0;
				if(r<-1.0)	r=-1.0;
					
	
				//Blending Factor
				Psi = BF_HiRAC(r);
					
				//Calculation of Alpha at Cell face
				nAlpha[i][j][k][l] = Psi*alphaBD + (1.0-Psi)*alphaHR;
                                
                                if(nAlpha[i][j][k][l]>1)        nAlpha[i][j][k][l] = 1.0;
                                if(nAlpha[i][j][k][l]<0)        nAlpha[i][j][k][l] = 0;
                                
				sAlpha[i][jnorth][k][l] = nAlpha[i][j][k][l];
	

                                un = 0.25*(u[i][j][k][l] + u[iwest][j][k][l] + u[i][jnorth][k][l] + u[iwest][jnorth][k][l]);
	                        vn = v[i][j][k][l];
	                        wn = 0.25*(w[i][j][k][l] + w[i][j][kbottom][l] + w[i][jnorth][kbottom][l] + w[i][jnorth][k][l]);
	                        
	                        us = 0.25*(u[i][j][k][l] + u[iwest][j][k][l] + u[i][jsouth][k][l] + u[iwest][jsouth][k][l]);
	                        vs = v[i][jsouth][k][l];
	                        ws = 0.25*(w[i][j][k][l] + w[i][j][kbottom][l] + w[i][jsouth][kbottom][l] + w[i][jsouth][k][l]);
	                        
	                        magVn[i][j][k] = pow((pow(un,2)+pow(vn,2)+pow(wn,2)),0.5); 
	                        magVs[i][j][k] = pow((pow(us,2)+pow(vs,2)+pow(ws,2)),0.5); 


			}	
		}	
	}	
	


	//Constructing alpha fluxes for faces along z direction
	for(k=2;k<=nzm-1;k++)
	{
		ktop = k+1;
		kbottom = k-1;
		ktoptop = k+2;
		dz = zf[k]-zf[kbottom];
		dzpt = zc[ktop]-zc[k];
		fzt = (zf[k] - zc[k])/dzpt;
		fzp = 1.0 - fzt;

		for(i=2;i<=nxm;i++)
		{
			for(j=2;j<=nym;j++)
			{

				//Calculating sign function using convective fluxes 				
				upwind_implicit(Ft);


				/* Arithmetic Interpolation of density*/
				rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];

				/* Harmonic Interpolation of density*/
/*				rhot = 2.0 * rho[i][j][k] * rho[i][j][ktop]/( rho[i][j][k] + rho[i][j][ktop]);*/

				if(Ft[i][j][k]>=0)
				{

					alphaD = tauAlpha[i][j][k][l-1];
					alphaA = tauAlpha[i][j][ktop][l-1];
					alphaU = tauAlpha[i][j][kbottom][l-1];
				}
				else			
				{
					alphaD = tauAlpha[i][j][ktop][l-1];
					alphaA = tauAlpha[i][j][k][l-1];
					alphaU = tauAlpha[i][j][ktoptop][l-1];
				}

				N = alphaA - alphaU;
				D = alphaD - alphaU;
				
				if(D==0)	rS = 1.0;
				else 		rS = N/D;
				
				kS = N*(1.0-Ct[i][j][k])*3.0/8.0 + D*(3.0+Ct[i][j][k])/4.0 + alphaU;


				//VOF scheme
				HiRAC();
					
				//Calculation of angle between interface normal and cell face
				//Numerator
				if(curvature==1||curvature==3)
				{
					if(Ft[i][j][k]>=0)	N = (tAlpha[i][j][k][l]	- bAlpha[i][j][k][l])/dz;
					else			N = (tAlpha[i][j][ktop][l] - bAlpha[i][j][ktop][l])/(zf[ktop]-zf[k]);
				}
						
				if(curvature==2)
				{
					if(Ft[i][j][k]>=0)	N = zGradAlphaS[i][j][k][l];
					else			N = zGradAlphaS[i][j][ktop][l];
				}
						
						//printf("N=%e \n",N);
			
				//Denominator
				D = magGradAlphaD();

				//printf("D=%e \n",D);
				if(D==0)	D=1e-100;
	
				r = N/D;
/*				printf("r = %lf \n",r);*/
					
				if(r>1.0)	r=1.0;
				if(r<-1.0)	r=-1.0;
					
	
				//Blending Factor
				Psi = BF_HiRAC(r);
					
					
				//Calculation of Alpha at Cell face
				tAlpha[i][j][k][l] = Psi*alphaBD + (1.0-Psi)*alphaHR;	
				
				if(tAlpha[i][j][k][l]>1)        tAlpha[i][j][k][l] = 1.0;
                                if(tAlpha[i][j][k][l]<0)        tAlpha[i][j][k][l] = 0;
                                				
				bAlpha[i][j][ktop][l] = tAlpha[i][j][k][l];

                                ut = 0.25*(u[i][j][k][l] + u[iwest][j][k][l] + u[i][j][ktop][l] + u[iwest][j][ktop][l]);
	                        vt = 0.25*(v[i][j][k][l] + v[i][jsouth][k][l] + v[i][j][ktop][l] + v[i][jsouth][ktop][l]);
	                        wt = 0.25*w[i][j][k][l];
	                        
	                        ub = 0.25*(u[i][j][k][l] + u[iwest][j][k][l] + u[i][j][kbottom][l] + u[iwest][j][kbottom][l]);
	                        vb = 0.25*(v[i][j][k][l] + v[i][jsouth][k][l] + v[i][j][kbottom][l] + v[i][jsouth][kbottom][l]);
	                        wb = 0.25*w[i][j][kbottom][l];
	                        
	                        magVt[i][j][k] = pow((pow(ut,2)+pow(vt,2)+pow(wt,2)),0.5); 
	                        magVb[i][j][k] = pow((pow(ub,2)+pow(vb,2)+pow(wb,2)),0.5); 
	
			}	
		}
	}
	

	for  (i=2; i<=nxm; i++)
 	{
	 	for (j=2; j<=nym; j++)
   		{
	  		for(k=2; k<=nzm; k++)
			{
				
				ieast = i+1;
				iwest = i-1;
				jnorth = j+1;
				jsouth = j-1;
				ktop = k+1;
				kbottom = k-1;
				
				/* Arithmetic Interpolation of density*/
				rhoe = 0.5 * rho[i][j][k] +   0.5 * rho[ieast][j][k];
				rhon = 0.5 * rho[i][j][k] +   0.5 * rho[i][jnorth][k];
				rhot = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][ktop];
				
				rhow = 0.5 * rho[i][j][k] +   0.5 * rho[iwest][j][k];
				rhos = 0.5 * rho[i][j][k] +   0.5 * rho[i][jsouth][k];
				rhob = 0.5 * rho[i][j][k] +   0.5 * rho[i][j][kbottom];
				
		
			
				taP[i][j][k] = ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dt - ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dtau;
				naP[i][j][k] =  - ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dt;
				ttaP[i][j][k] = ((xf[i]-xf[i-1])*(yf[j]-yf[j-1])*(zf[k]-zf[k-1]))/dtau;

                                Sp[i][j][k] = 0;
				Sp[i][j][k] = - taP[i][j][k]*tauAlpha[i][j][k][l-1] - naP[i][j][k]*alpha[i][j][k][l-1];
				Sp[i][j][k] -= Fe[i][j][k]         *(eAlpha[i][j][k][l] + eAlpha[i][j][k][l-1])/(2.0*rhoe);
				Sp[i][j][k] += Fe[iwest][j][k]     *(wAlpha[i][j][k][l] + wAlpha[i][j][k][l-1])/(2.0*rhow);
				Sp[i][j][k] -= Fn[i][j][k]         *(nAlpha[i][j][k][l] + nAlpha[i][j][k][l-1])/(2.0*rhon);
				Sp[i][j][k] += Fn[i][jsouth][k]    *(sAlpha[i][j][k][l] + sAlpha[i][j][k][l-1])/(2.0*rhos);
				Sp[i][j][k] -= Ft[i][j][k]         *(tAlpha[i][j][k][l] + tAlpha[i][j][k][l-1])/(2.0*rhot);
				Sp[i][j][k] += Ft[i][j][kbottom]   *(bAlpha[i][j][k][l] + bAlpha[i][j][k][l-1])/(2.0*rhob);
/*				Sp[i][j][k] += 0.5*eNorm[i][j][k]*magVe[i][j][k]*eAlpha[i][j][k][l]*(1.0-eAlpha[i][j][k][l])*dy*dz;*/
/*				Sp[i][j][k] -= 0.5*wNorm[i][j][k]*magVw[i][j][k]*wAlpha[i][j][k][l]*(1.0-wAlpha[i][j][k][l])*dy*dz;*/
/*				Sp[i][j][k] += 0.5*nNorm[i][j][k]*magVn[i][j][k]*nAlpha[i][j][k][l]*(1.0-nAlpha[i][j][k][l])*dx*dz;*/
/*				Sp[i][j][k] -= 0.5*sNorm[i][j][k]*magVs[i][j][k]*sAlpha[i][j][k][l]*(1.0-sAlpha[i][j][k][l])*dx*dz;*/
/*				Sp[i][j][k] += 0.5*tNorm[i][j][k]*magVt[i][j][k]*tAlpha[i][j][k][l]*(1.0-tAlpha[i][j][k][l])*dx*dy;*/
/*				Sp[i][j][k] -= 0.5*bNorm[i][j][k]*magVb[i][j][k]*bAlpha[i][j][k][l]*(1.0-bAlpha[i][j][k][l])*dx*dy;*/
				
				
/*				printf("eNorm=%e,wNorm=%e,nNorm=%e,sNorm=%e,tNorm=%e,bNorm=%e  \n",eNorm[i][j][k],wNorm[i][j][k],nNorm[i][j][k],sNorm[i][j][k],tNorm[i][j][k],bNorm[i][j][k]);*/
				
				
				

/*				printf("Sp=%e taP=%e naP=%e ttaP=%e nAlpha = %e nAlpha = %e eAlpha = %e eAlpha = %e wAlpha = %e wAlpha = %e sAlpha = %e sAlpha = %e res = %e  tauAlpha = %e alpha = %e \n",Sp[i][j][k],taP[i][j][k],naP[i][j][k],ttaP[i][j][k],nAlpha[i][j][k][l],nAlpha[i][j][k][l-1],eAlpha[i][j][k][l],eAlpha[i][j][k][l-1],wAlpha[i][j][k][l],wAlpha[i][j][k][l-1],sAlpha[i][j][k][l],sAlpha[i][j][k][l-1],pseudoMaxRes,tauAlpha[i][j][k][l-1],alpha[i][j][k][l-1]);*/

			}
		}
	}


			
}
			
				 



				
