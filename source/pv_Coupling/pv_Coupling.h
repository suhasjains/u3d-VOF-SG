#ifndef PV_COUPLING
#define PV_COUPLING

void simple();
void upwind_implicit(double ***flux);
void v_east_source_dc(double f1,double f2,double f3,double f4,double f5,double f6);
void v_north_source_dc(double f1,double f2,double f3,double f4,double f5,double f6);
void v_top_source_dc(double f1,double f2,double f3,double f4,double f5,double f6);

#endif
