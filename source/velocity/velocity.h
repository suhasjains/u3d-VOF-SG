#ifndef VELOCITY
#define VELOCITY

void uEqn();
void vEqn();
void wEqn();
void uEqnCoeff();
void vEqnCoeff();
void wEqnCoeff();
void blockVelocity();
void correct_BottomVelBoundary();
void correct_EastVelBoundary();
void correct_NorthVelBoundary();
void correct_SouthVelBoundary();
void correct_TopVelBoundary();
void correct_WestVelBoundary();
void makeVelocityVectors();

#endif
