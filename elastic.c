#include "pasta.h"
#include "matrix.h"
#include "drying.h"

#include <stdlib.h>
#include <stdio.h>

#define NTERMS 40

double strain(double t, double x, double RH, double D, double X0, double Xe, double L, double T)
{
    int i,
        nt = 100;
    double Xdb1, Xdb2,
           Eb1, Eb2,
           GradEb,
           h = 1e-5,
           e = 0,
           E = 5e9; /* From Guinea et al. 2004 */

    oswin *d;

    /* Point 1: x-h, t
     * Point 2: x+h, t
     */

    d = CreateOswinData();
    
    for(i=0; i<nt; i++) {
        Xdb1 = CrankEquationFx(x-h, t, L, D, Xe, X0, NTERMS);
        Xdb2 = CrankEquationFx(x+h, t, L, D, Xe, X0, NTERMS);

        Eb1 = BindingEnergyOswin(d, Xdb1, T);
        Eb2 = BindingEnergyOswin(d, Xdb2, T);

        GradEb = (Eb2 - Eb1)/(2*h);
        e = -GradEb / E;
    }

    return e;
}

double strainpc(double t, double x, double RH, double D, double X0, double Xe, double L, double T)
{
    int i,
        nt = 100;
    double Xdb, Pc,
           e = 0,
           E = 5e9; /* From Guinea et al. 2004 */

    for(i=0; i<nt; i++) {
        Xdb = CrankEquationFx(x, t, L, D, Xe, X0, NTERMS);

        Pc = pore_press(Xdb, T);
        e = Pc / E;
    }

    return e;
}

double displacement(int xi, vector* strain, double L)
{
    int i;
    double dx = L/len(strain),
           u = 0;
    for(i=0; i<xi; i++) {
        u += valV(strain, i) * dx;
    }

    return u;
}

int main(int argc, char *argv[])
{
    double L = 1e-3, /* Length [m] */
           D, /* Diffusivity [m^2/s] */
           X0 = .3, /* Initial moisture content */
           Xe, /* Equilibrium moisture content */
           T = 60+273.15, /* Drying temperature */
           h = 1e-7,
           Xdbph,
           t,
           RH;

    oswin *d; /* Isotherm data */

    int npts = 50,
        i;
    
    vector *z, *Xdb, *Eb, *Es, *GradEb, *GradEs, *str, *strpc, *u;
    matrix *output;

    if(argc != 3) {
        puts("Usage:");
        puts("dry <aw> <t>");
        puts("aw: Final water activity to dry to.");
        puts("t: Time to output a profile for.");
        exit(0);
    }

    RH = atof(argv[1]);
    t = atof(argv[2]);
    d = CreateOswinData();
    Xe = OswinIsotherm(d, RH, T);
    D = DiffCh10((Xe+X0)/2, T);

    z = CreateVector(npts);
    Xdb = CreateVector(npts);
    Eb = CreateVector(npts);
    Es = CreateVector(npts);
    GradEb = CreateVector(npts);
    GradEs = CreateVector(npts);
    str = CreateVector(npts);
    strpc = CreateVector(npts);
    u = CreateVector(npts);

//    CrankEquationFx(L/2, t, L, D, Xe, X0, NTERMS);
    printf("Xe = %g\n", Xe);
    
    for(i=0; i<npts; i++)
        setvalV(z, i, L/npts * i);

    for(i=0; i<npts; i++) {
        setvalV(Xdb, i, CrankEquationFx(valV(z, i), t, L, D, Xe, X0, NTERMS));
        setvalV(Eb, i, BindingEnergyOswin(d, valV(Xdb, i), T));
        //setvalV(Eb, i, BindingEnergyOswin(d, h, T));
        setvalV(Es, i, Esurf(valV(Xdb, i), T));
        Xdbph = CrankEquationFx(valV(z,i)+h, t, L, D, Xe, X0, NTERMS);
        setvalV(GradEb, i, (BindingEnergyOswin(d, Xdbph, T) - BindingEnergyOswin(d, valV(Xdb, i), T))/h);
        setvalV(GradEs, i, (Esurf(Xdbph, T) - Esurf(valV(Xdb, i), T))/h);
        setvalV(str, i, strain(t, valV(z, i), RH, D, X0, Xe, L, T));
        setvalV(strpc, i, strainpc(t, valV(z, i), RH, D, X0, Xe, L, T));
    }
    for(i=0;i<npts;i++)
        setvalV(u, i, displacement(i, str, L));

    output = CatColVector(6, z, Xdb, GradEb, str, strpc, u);
    mtxprntfile(output, "out.csv");

    return 0;
}

