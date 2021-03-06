#ifndef DRYING_H
#define DRYING_H

#define CONSTX0 0
#define CONSTXe 18.261700
#define CONSTnterms 50 
#define BETA0 1e-4
#define ETA .28 // 40 C
//#define ETA .15 // 80 C

#define POISSON .45

#define SLABWIDTH 6e-3
#define SLABLENGTH 8e-3

#define NTERMS 100
#define NTSTEPS 100

typedef struct {
    double L;
    double D;
    double X0;
    double Xe;
    double T;
    int nterms;
} drydat;


double CrankEquationFx(double, double, drydat);
double CrankEquation(double, drydat);

double Esurf(double, double);
double GradEsurf(double, double, drydat);
double PoreP(double, double, drydat);
double EffPoreP(double, double, drydat);

double maxstrain(double, double, drydat);
double etastrain(double, double, drydat, double);
double solidfrac(double, double, drydat, double);


#endif

