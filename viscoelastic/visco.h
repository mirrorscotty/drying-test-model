#ifndef VISCO_H
#define VISCO_H

#include "drying.h"

double strain(double, double, double, double, double, double, double, double);
double maxstrain(double, double, drydat);
double etastrain(double, double, drydat, double);

double strainpc(double, double, drydat,
        double (*)(double, double, drydat),
        double (*)(double, double, double, double, int));

double MaxwellStrainPc(double, double, drydat);
double ZhuMaxwellStrain(double, double, drydat);
double EqStrainPc(double, double, drydat);
double displacement(int, vector*, double);
vector* displacementV(vector*, double);

double CreepGina(double, double, double, double, int);
double CreepLaura(double, double, double, double, int);
double CreepZhu(double, double, double, double, int);

#endif

