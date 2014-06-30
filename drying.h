#ifndef DRYING_H
#define DRYING_H

#define CONSTX0 0
#define CONSTXe 18.261700
#define CONSTnterms 50 
#define BETA0 1e-4

#define SLABWIDTH 6e-3
#define SLABLENGTH 8e-3

double CrankEquationFx(double, double, double, double, double, double, int);

double Esurf(double, double);

#endif

