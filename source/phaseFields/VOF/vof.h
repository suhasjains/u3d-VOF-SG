#ifndef VOF
#define VOF


void alphaEqnCoeff();
void alphaEqnCorrectionCoeff();
double magGradAlphaD();
void alphaEqn();
void CICSAM(double nAlphaD);
double BF_CICSAM(double Theta);
void STACS(double nAlphaD);
double BF_STACS(double Theta);

#endif
