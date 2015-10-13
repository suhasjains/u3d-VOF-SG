#include"preProcessing.h"


//dynamically allocates the size of array
void dynamicAllocate()
{


/*	4D arrays*/
	
	u = create4DdoubleArray();
	v = create4DdoubleArray();
	w = create4DdoubleArray();
	p = create4DdoubleArray();
	pp = create4DdoubleArray();
	alpha = create4DdoubleArray();
	finalU = create4DdoubleArray();
	finalV = create4DdoubleArray();
	finalW = create4DdoubleArray();
	ue = create4DdoubleArray();
	vn = create4DdoubleArray();
	wt = create4DdoubleArray();
	alphaS = create4DdoubleArray();
	xGradAlphaS = create4DdoubleArray();
	yGradAlphaS = create4DdoubleArray();
	zGradAlphaS = create4DdoubleArray();
	T = create4DdoubleArray();
	Psi = create4DdoubleArray();
	aE = create4DdoubleArray();
	aW = create4DdoubleArray();
	aN = create4DdoubleArray();
	aS = create4DdoubleArray();
	aT = create4DdoubleArray();
	aB = create4DdoubleArray();
	aPE = create4DdoubleArray();
	aPW = create4DdoubleArray();
	aPN = create4DdoubleArray();
	aPS = create4DdoubleArray();
	aPT = create4DdoubleArray();
	aPB = create4DdoubleArray();
	Ap = create4DdoubleArray();
	xBFe = create4DdoubleArray();
	yBFn = create4DdoubleArray();
	zBFt = create4DdoubleArray();
	tauAlpha = create4DdoubleArray();
	eAlpha = create4DdoubleArray();
        wAlpha = create4DdoubleArray();
        nAlpha = create4DdoubleArray();
        sAlpha = create4DdoubleArray();
        tAlpha = create4DdoubleArray();
        bAlpha = create4DdoubleArray();
	
	
	
	
	

/*        3D arrays	*/

        xCoord = create3DdoubleArray();
        yCoord = create3DdoubleArray();
        zCoord = create3DdoubleArray();
        Ft = create3DdoubleArray();
        Fn = create3DdoubleArray();
        Fe = create3DdoubleArray();
        su = create3DdoubleArray();
        sv = create3DdoubleArray();
        sw = create3DdoubleArray();
        ae = create3DdoubleArray();
        aw = create3DdoubleArray();
        an = create3DdoubleArray();
        as = create3DdoubleArray();
	at = create3DdoubleArray();
        ab = create3DdoubleArray();
        ap = create3DdoubleArray();
        apu = create3DdoubleArray();
        apv = create3DdoubleArray();
        apw = create3DdoubleArray();
        dpx = create3DdoubleArray();
        dpy = create3DdoubleArray();
        dpz = create3DdoubleArray();
        uRes = create3DdoubleArray();
        vRes = create3DdoubleArray();
        wRes = create3DdoubleArray();
        pRes = create3DdoubleArray();
        mu = create3DdoubleArray();
        rho = create3DdoubleArray();
        aNe = create3DdoubleArray();
        aBe = create3DdoubleArray();
        aTe = create3DdoubleArray();
        aSe = create3DdoubleArray();
        aee = create3DdoubleArray();
        anB = create3DdoubleArray();
        anE = create3DdoubleArray();
        anT = create3DdoubleArray();
        anW = create3DdoubleArray();
        ann = create3DdoubleArray();
        aNt = create3DdoubleArray();
        aSt = create3DdoubleArray();
        atE = create3DdoubleArray();
        atW = create3DdoubleArray();
        att = create3DdoubleArray();
        Sp = create3DdoubleArray();
        COutD = create3DdoubleArray();
        Ce = create3DdoubleArray();
        Cn = create3DdoubleArray();
        Ct = create3DdoubleArray();
        eBeta = create3DdoubleArray();
        nBeta = create3DdoubleArray();
        tBeta = create3DdoubleArray();
        alphaRes = create3DdoubleArray();
        Kappa = create3DdoubleArray();
        KappaHF = create3DdoubleArray();
        eKappaHF = create3DdoubleArray();
        nKappaHF = create3DdoubleArray();
        tKappaHF = create3DdoubleArray();
        xST = create3DdoubleArray();
        yST = create3DdoubleArray();
        zST = create3DdoubleArray();
        xSTe = create3DdoubleArray();
        ySTn = create3DdoubleArray();
        zSTt = create3DdoubleArray();
        uB = create3DdoubleArray();
        vB = create3DdoubleArray();
        wB = create3DdoubleArray();
        sdf = create3DdoubleArray();
        ttaP = create3DdoubleArray();
        taP = create3DdoubleArray();
        naP = create3DdoubleArray();
        pseudoRes = create3DdoubleArray();
        TRes = create3DdoubleArray();
        Kcond = create3DdoubleArray();
        Cp = create3DdoubleArray();
        Beta = create3DdoubleArray();
        sT = create3DdoubleArray();
        apT = create3DdoubleArray();
        xBFw = create3DdoubleArray();
        yBFs = create3DdoubleArray();
        zBFb = create3DdoubleArray();
        xBF = create3DdoubleArray();
        yBF = create3DdoubleArray();
        zBF = create3DdoubleArray();
        magVb = create3DdoubleArray();  
        magVe = create3DdoubleArray();
        magVn = create3DdoubleArray();
        magVs = create3DdoubleArray();
        magVt = create3DdoubleArray();
        magVw = create3DdoubleArray();
        eNorm = create3DdoubleArray();
        wNorm = create3DdoubleArray();
        nNorm = create3DdoubleArray();
        sNorm = create3DdoubleArray();
        tNorm = create3DdoubleArray();
        bNorm = create3DdoubleArray();
        apP = create3DdoubleArray();
        sP = create3DdoubleArray();
        PsiRes = create3DdoubleArray();
        ePsi = create3DdoubleArray();
        wPsi = create3DdoubleArray();
        nPsi = create3DdoubleArray();
        sPsi = create3DdoubleArray();
        tPsi = create3DdoubleArray();
        bPsi = create3DdoubleArray();
        
        
        
        
        
}


double ****create4DdoubleArray()
{

        double ****array4D;

        array4D	 = new double***[xNumberOfCells+5];
     
        for (int i=0; i<(xNumberOfCells+5); ++i)
	{
	        array4D[i]	 = new double**[yNumberOfCells+5];
	        
	        for (int j=0; j<(yNumberOfCells+5); ++j)
	        {
                        array4D[i][j]	 = new double*[zNumberOfCells+5];
                        
                        for (int k=0; k<(zNumberOfCells+5); ++k)
		        {
                                array4D[i][j][k] 	 = new double[3];
                        }
                }
        }
        
        return array4D;

}


double ***create3DdoubleArray()
{

        double ***array3D;

        array3D	 = new double**[xNumberOfCells+5];
     
        for (int i=0; i<(xNumberOfCells+5); ++i)
	{
	        array3D[i]	 = new double*[yNumberOfCells+5];
	        
	        for (int j=0; j<(yNumberOfCells+5); ++j)
	        {
                        array3D[i][j]	 = new double[zNumberOfCells+5];
                }
        }
        
        return array3D;
        
        
}
