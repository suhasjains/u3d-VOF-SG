//Declaration of all variables


#ifndef VARIABLES_H
#define VARIABLES_H


	//INPUT OUTPUT VARIABLES
	double xDomainLength,yDomainLength,zDomainLength;
	int   xNumberOfCells,yNumberOfCells,zNumberOfCells,nx,ny,nz,n,nxm,nym,nzm;   //number of cells
	int   nTimeSteps;
	double startTime,endTime,deltaT,simulationTime,realTime,start,stop;
	char  adjustTimeStep;
	double uTopBoundVel,uBottomBoundVel,uLeftBoundVel,uRightBoundVel,uFrontBoundVel,uBackBoundVel;
	double vTopBoundVel,vBottomBoundVel,vLeftBoundVel,vRightBoundVel,vFrontBoundVel,vBackBoundVel;
	double wTopBoundVel,wBottomBoundVel,wLeftBoundVel,wRightBoundVel,wFrontBoundVel,wBackBoundVel;
	double TopBoundTemp,BottomBoundTemp,LeftBoundTemp,RightBoundTemp,FrontBoundTemp,BackBoundTemp;
	int TopBoundType,BottomBoundType,LeftBoundType,RightBoundType,FrontBoundType,BackBoundType;


	//COEFFICIENTS OF VELOCITY EQUATION
	double dx,dy,dz,dt;
	double xGrav,yGrav,zGrav;
	double ***Fe,***Fn,***Ft;
	
	double ***aNe,***aBe,***aTe,***aSe,***aee;
	double ***anB,***anE,***anT,***anW,***ann;
	double ***aNt,***aSt,***atE,***atW,***att;
	
	double ****finalU,****finalV,****finalW;


	//COEFFICIENTS OF PRESSURE CORRECTION EQUATION
	double ***ae,***aw,***an,***as,***at,***ab,***ap;
	double ***su,***sv,***sw,maxMassResidual;
	double ***apu,***apv,***apw,apt;
	double uMaxRes,vMaxRes,wMaxRes,pMaxRes,TMaxRes,***uRes,***vRes,***wRes,***pRes,uMaxResidual[2],vMaxResidual[2],wMaxResidual[2],***TRes;
	double uMeanRes,vMeanRes,wMeanRes,pMeanRes,TMeanRes;
	double uTotalRes,vTotalRes,wTotalRes,pTotalRes,TTotalRes;
	double uNormRes,vNormRes,wNormRes,pNormRes,TNormRes;
	int pMassResidual,nMassResidual,nsMassResidual,psMassResidual,nuMaxRes,nvMaxRes,nwMaxRes,puMaxRes,pvMaxRes,pwMaxRes;

	//VELOCITIES and PRESSURE
	double ****u,****v,****w;
	double pe,pw,pn,ps,pt,pb,ppo;
	double ***dpx,***dpy,***dpz;													//pressure gradient
	double ****pp,****p;
	double ***xCoord,***yCoord,***zCoord;
	double convf,convp,***flux;														//convective mass fluxes
	double diff;
	
	double vol,vole,voln,volt;
	
	//Rhie Chow vriables
	double ****ue,****vn,****wt;
	
	
	double uMax,vMax,wMax;															//max velocities
	double maxCourantNumber,allowableCourantNumber,leastCourantNumber;									//max courant number
	int adjustableCourantNumber;
	int nSimpleLoops,uVelocityLoops,vVelocityLoops,wVelocityLoops,nPressureLoops,nSimplerLoops,nAlphaLoops,nTemperatureLoops;

	int i,j,k,l,t,courantNumberControl,thomasI,simpleVar,simplerVar,nonLinear;      //indices
	int ieast,iwest,jnorth,jsouth,ktop,kbottom,ieasteast,jnorthnorth,ktoptop;								
	
	double fxe,fxp,fyn,fyp,fzt,fzp,fxw,fys,fzb;

	double dxpe,dypn,dzpt,dxpw,dyps,dzpb;													//transoformation
	
	double lower[250],upper[250],c[250],x[250],m[250],cc[250],middle[250];       	//MATRIX A,X AND B

	double pUnderRelaxCoeff,vUnderRelaxCoeff,uUnderRelaxCoeff,wUnderRelaxCoeff,tUnderRelaxCoeff,urecru,urecrv,urecrw,urecrp,urecrt; // under relaxation coefficient for pressure

	double ***rho,***mu,rho1,rho2,mu1,mu2;										//density and viscosity
	
	//area
	double s;

	//Output variables	
	FILE *fp; 													//file pointer
	char str[100];													//file name
	int writeOutputSteps;


	//Deferred correction
	double DCvalue,fuuds,fvuds,fwuds,fucds,fvcds,fwcds,fuQ,fvQ,fwQ,fuH,fvH,fwH,fuQHP,fvQHP,fwQHP,fuP,fvP,fwP,f1,f2,f3,f4,f5,f6;

	//grid
	double xc[250],yc[250],zc[250],xf[250],yf[250],zf[250],fx[250],fy[250],fz[250];

	//accuracy of outer iterations and p and v loops
	double vAccuracy,pAccuracy,alphaAccuracy,massAccuracy,TAccuracy;
	
	//scheme
	double scheme;
	
	//counts no of outer iterations
	int check;

	//Time stepping
	double Gamma;

	//Phase field
	double ****alpha;
	double ****alphaS;
	double ****xGradAlphaS,****yGradAlphaS,****zGradAlphaS;

	//Viscosity at cell faces	
	double mue,mun,mut;
	
	//Conductivity 
	double ***Kcond,Kconde,Kcondn,Kcondt,Kcond1,Kcond2;
	
	//Specific heat 
	double Cp1,Cp2,***Cp,Cpe,Cpn,Cpt;
	
	//Temperature
	double ****T;
	
	//reference Temperature
	double Tref;
	
	
	
	//volumetric expansion coefficient
	double Beta1,Beta2,***Beta;
	
	//acceleration due to temperature difference
	double xTempGrav,yTempGrav,zTempGrav;
	
	//Inverse of Prandtl number
	double Prr;
	
	//Coefficients of Temperature Equation
	double ***sT, fTcds, fTuds,***apT;


	//Coefficients of Phase Fields
	double ****eAlpha,****wAlpha,****nAlpha,****sAlpha,****tAlpha,****bAlpha;
	double ****aE,****aW,****aN,****aS,****aT,****aB;
	double ****aPE,****aPW,****aPN,****aPS,****aPT,****aPB,****Ap;
	double ***Sp;

//	Local Courant Number
	double ***COutD,***Ce,***Cn,***Ct;

	double ***eBeta,***nBeta,***tBeta;

	//Residual of Alpha
	double alphaMaxRes,alphaTotalRes,alphaMeanRes,alphaNormRes;	
	double ***alphaRes;
	int maxnAlphaLoops,nAlphaCorrectorLoops,say;
	double maxAlpha,minAlpha;

	double maxnegError;
	
	//Curvature
	double ***Kappa;
	double ***KappaHF,***eKappaHF,***nKappaHF,***tKappaHF;
	int curvature;
	
	//Surface Tension
	double ***xST,***yST,***zST,***xSTe,***ySTn,***zSTt;
	double sigma;
	
	//No of FLuids
	int nFluids;
	
	//No of Blocks
	int xNB,yNB,zNB;
	
	
	//multi Block variables
	int e,f,g,b;					//Block indices
	int BlockID[10][10][10];
	int ecID[10][10][10],wcID[10][10][10],tcID[10][10][10],bcID[10][10][10],ncID[10][10][10],scID[10][10][10];			//Cell indices
	double ueBA[5][55][10][10],uwBA[5][55][10][10],unBA[55][5][10][10],usBA[55][5][10][10],utBA[55][55][5][10],ubBA[55][55][5][10];		//Buffer Array or Extension cells
	double veBA[5][55][10][10],vwBA[5][55][10][10],vnBA[55][5][10][10],vsBA[55][5][10][10],vtBA[55][55][5][10],vbBA[55][55][5][10];		//Buffer Array or Extension cells
	double weBA[5][55][10][10],wwBA[5][55][10][10],wnBA[55][5][10][10],wsBA[55][5][10][10],wtBA[55][55][5][10],wbBA[55][55][5][10];		//Buffer Array or Extension cells
	double ppeBA[5][55][10][10],ppwBA[5][55][10][10],ppnBA[55][5][10][10],ppsBA[55][5][10][10],pptBA[55][55][5][10],ppbBA[55][55][5][10];	//Buffer Array or Extension cells
	double alphaeBA[5][55][10][10],alphawBA[5][55][10][10],alphanBA[55][5][10][10],alphasBA[55][5][10][10],alphatBA[55][55][5][10],alphabBA[55][55][5][10];	//Buffer Array or Extension cells
	int nBlocks;
	
	double ***uB,***vB,***wB;
	
	
	//Improved Rhie Chow for large Body Forces
	double ****xBFe,****yBFn,****zBFt,***xBFw,***yBFs,***zBFb,***xBF,***yBF,***zBF;
	
	
	//Signed Distance Function
	double ***sdf;
	
	//VOFschemes
	double nAlphaBD,nAlphaHR;		//Bounded Downwind and High Resolution 
	int VOFScheme;
	
	
	double capillaryTimeStep;
	
	int fixedAdvectionVelocity;
	
//	Pseudo time variables
        double ***ttaP,***taP,***naP,****tauAlpha,pseudoMaxRes,pseudoMeanRes,pseudoTotalRes,pseudoNormRes,***pseudoRes;
        double pseudoTimeAccuracy;
        double dtau;
        int pseudoLoops,nPseudoLoops;
        
//        Succesive Over Relaxation
        double SOR;
        
//        HiRAC variables
	double rS,kS,alphaBD,alphaHR,alphaA,alphaD,alphaU;    
	
//	Artificial Compressinility
	double ***magVe,***magVb,***magVn,***magVs,***magVt,***magVw;    

//        interface normals                
        double ***eNorm,***wNorm,***nNorm,***sNorm,***tNorm,***bNorm;
        
//        new scalar
        double ****Psi,***sP,***apP,***PsiRes;
        double urecrpsi,PsiMaxRes,PsiMeanRes,PsiNormRes,PsiTotalRes,psiUnderRelaxCoeff;
        int nPsiLoops;
        double ***ePsi,***wPsi,***nPsi,***sPsi,***tPsi,***bPsi;
        
#endif

