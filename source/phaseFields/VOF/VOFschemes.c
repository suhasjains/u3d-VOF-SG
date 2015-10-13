
//calculation of CICSAM scheme
void CICSAM(double nAlphaD)
{
	//Calculation of CBC in CICSAM scheme
	if(nAlphaD>0.0&&nAlphaD<=1.0)
	{
		if(COutD[i][j][k]==0.0)	nAlphaBD = 1.0;			
		else			nAlphaBD = std::min(nAlphaD/COutD[i][j][k],1.0);
	}
		else			nAlphaBD = nAlphaD;

	//Calculation of UQ in CICSAM scheme
	if(nAlphaD>=0||nAlphaD<=1)			nAlphaHR = std::min(((8.0*COutD[i][j][k]*nAlphaD + (1.0 - COutD[i][j][k])*(6.0*nAlphaD + 3.0))/8.0),nAlphaBD);
	else 						nAlphaHR = nAlphaD;
	
/*	printf("nAlphaBD=%e\n",nAlphaBD);*/
/*	printf("nAlphaD=%lf\n",nAlphaD);*/
/*	printf("nAlphaHR=%lf\n",nAlphaHR);*/
/*	printf("COutD=%e\n",COutD[i][j][k]);*/
		
}


//Blending Factor for CICSAM scheme
double BF_CICSAM(double Theta)
{

	double BF;
	
	BF = std::min(((cos(2.0*Theta)+1.0)/2.0),1.0);
	
	return BF;

}


void STACS(double nAlphaD)
{
	//Calculation of SUPERBEE in STACS scheme 
	if(nAlphaD<=0.0)	nAlphaBD = nAlphaD;
	else if(nAlphaD>=1.0)	nAlphaBD = nAlphaD;
	else 			nAlphaBD = 1.0;	

	//Calculation of STOIC in STACS scheme
	if(nAlphaD<=0.0)				nAlphaHR = nAlphaD;
	else if((nAlphaD>0.0)&&(nAlphaD<=0.5))		nAlphaHR = 0.5 + 0.5*nAlphaD;
	else if((nAlphaD>0.5)&&(nAlphaD<=(5.0/6.0)))	nAlphaHR = (3.0/8.0) + (3.0/4.0)*nAlphaD;
	else if((nAlphaD>(5.0/6.0))&&(nAlphaD<=1.0))	nAlphaHR = 1.0;
	else 						nAlphaHR = nAlphaD;

/*	printf("nAlphaBD=%e\n",nAlphaBD);*/
/*	printf("nAlphaD=%lf\n",nAlphaD);*/
/*	printf("nAlphaHR=%lf\n",nAlphaHR);*/

}

double BF_STACS(double Theta)
{

	double BF;
	
	BF = pow(cos(Theta),4.0);
	
	return BF;	

}

void HiRAC()
{
	//Calculation of CBC in HiRAC scheme
	if(rS>1.0&&alphaD>alphaU)
	{
		if(COutD[i][j][k]==0.0)	alphaBD = 1.0;			
		else			alphaBD = std::min(((alphaD-alphaU)/COutD[i][j][k] + alphaU),alphaA);
	}
	else if(rS>1.0&&alphaD<alphaU)
	{
		if(COutD[i][j][k]==0.0)	alphaBD = 1.0;			
		else			alphaBD = std::max(((alphaD-alphaU)/COutD[i][j][k] + alphaU),alphaA);
	}
	else
	{
		alphaBD = alphaD;
	}
	
	
	
	//Calculation of UQ in HiRAC scheme
	if(rS>1.0&&alphaD>alphaU)
	{
		alphaHR = std::min(kS,alphaBD);
	}
	else if(rS>1.0&&alphaD<alphaU)
	{
		alphaHR = std::max(kS,alphaBD);
	}
	else
	{
		alphaHR = alphaD;
	}
	

}

double BF_HiRAC(double r)
{

	double BF;
	
	BF = std::min(pow(r,2.0),1.0);
	
	return BF;	

}

