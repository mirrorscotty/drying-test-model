#ifndef VISCO_H
#define VISCO_H

#include "drying.h"
#include "matrix.h"


double strain(double, double, double, double, double, double, double, double);
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
double CreepLaura2(double, double, double, double, int);
double CreepZhu(double, double, double, double, int);

double RelaxCummings(double, double, double, int);
double RelaxLaura(double, double, double, int);
double RelaxGina(double, double, double, int);
double RelaxZhu(double, double, double, int);
double VEStress(double, double, drydat,
        double (*E)(double, double, drydat, double),
        double (G)(double, double, double, int));

#endif

