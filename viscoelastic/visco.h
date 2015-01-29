#ifndef VISCO_H
#define VISCO_H

#include "drying.h"

double strain(double, double, double, double, double, double, double, double);
double maxstrain(double, double, drydat);

double strainpc(double, double, drydat);
double MaxwellStrainPc(double, double, drydat);
double ZhuMaxwellStrain(double, double, drydat);
double EqStrainPc(double, double, drydat);
double displacement(int, vector*, double);
vector* displacementV(vector*, double);

#endif

