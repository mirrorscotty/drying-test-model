#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "pasta.h"
#include "visco.h"
#include "drying.h"

#define NX 20
#define NTERMS 40

double SurfDisplace(double t, double RH, double D, double X0, double Xe, double L, double T)
{
    int nx = NX;
    double dx = L/nx,
           d;
    int i;

    vector *x, *Xdb, *str, *u;
    x = CreateVector(nx);
    Xdb = CreateVector(nx);
    str = CreateVector(nx);
    u = CreateVector(nx);

    for(i=0; i<nx; i++)
        setvalV(x, i, dx*i);
    for(i=0; i<nx; i++)
        setvalV(str, i, strainpc(t, valV(x,i), RH, D, X0, Xe, L, T));
    for(i=0; i<nx; i++)
        setvalV(u, i, displacement(i, str, L));
    
    d = valV(u, nx-1);
    DestroyVector(x);
    DestroyVector(Xdb);
    DestroyVector(str);
    DestroyVector(u);

    return d;
}

double SurfMoistureFlux(double t, double RH, double D, double X0, double Xe, double L, double T)
{
    double J, h = 1e-7, Xdbs, Xdbsmh;

    Xdbs = CrankEquationFx(L, t, L, D, Xe, X0, NTERMS);
    Xdbsmh = CrankEquationFx(L-h, t, L, D, Xe, X0, NTERMS);

    J = D*(Xdbs-Xdbsmh) / h;

    return J;
}

int main(int argc, char *argv[])
{
    double L = 1e-3, /* Length [m] */
           D, /* Diffusivity [m^2/s] */
           X0, /* Initial moisture content */
           Xe, /* Equilibrium moisture content */
           T, /* Drying temperature */
           h = 2,
           t,
           nt,
           RH,
           V,
           dt,
           kf;

    oswin *d; /* Isotherm data */

    int i; /* Loop index */
    
    vector *Xdb, /* Slab moisture content */
           *Vs, /* Surface Velocity */
           *Js, /* Moisture flux at the surface */
           *tv, /* Time vector */
           *Disp;
    matrix *out;

    /* Print out a usage statement if needed */
    if(argc != 6) {
        puts("Usage:");
        puts("dry <T> <X0> <aw> <t> <nt>");
        puts("T: Drying temperature [C]");
        puts("X0: Initial moisture content of the pasta. [kg/kg db]");
        puts("aw: Final water activity to dry to. [-]");
        puts("t: Length of time to simulate. [s]");
        puts("nt: Number of time steps to use.");
        exit(0);
    }

    /* Store command line arguments */
    T = atof(argv[1]) + 273.15;
    X0 = atof(argv[2]);
    RH = atof(argv[3]);
    t = atof(argv[4]);
    nt = atof(argv[5]);

    dt = t/nt;

    /* Calculate equilibrium moisture content and average diffusivity */
    d = CreateOswinData();
    Xe = OswinIsotherm(d, RH, T);
    D = DiffCh10((Xe+X0)/2, T);

    tv = CreateVector(nt);
    Xdb = CreateVector(nt);
    Vs = CreateVector(nt);
    Js = CreateVector(nt);
    Disp = CreateVector(nt);

    for(i=0; i<nt; i++) {
        setvalV(tv, i, dt*i);

        kf = M_PI*M_PI*D/(L*L);
        setvalV(Xdb, i, CrankEquation(kf, i*dt, X0, Xe, NTERMS));

        V = (SurfDisplace(nt*i+h, RH, D, X0, Xe, L, T) - SurfDisplace(nt*i-h, RH, D, X0, Xe, L, T))/(2*h);
        setvalV(Vs, i, V);

        setvalV(Js, i, SurfMoistureFlux(i*nt, RH, D, X0, Xe, L, T));

        setvalV(Disp, i, SurfDisplace(nt*i, RH, D, X0, Xe, L, T));
    }

    out = CatColVector(5, tv, Xdb, Disp, Vs, Js);

    mtxprntfilehdr(out, "out.csv", "Time [s],Surface Displacement [m],Surface Velocity [m/s],Surface Mass Flux [kg/m^2]\n");

    return 0;
}

