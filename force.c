#include "material-data.h"
#include "drying.h"

/**
 * Calculate "surface energy" by integrating the binding energy with respect to
 * moisture content at constant temperature. Since binding energy has units of
 * J/mol, the moisture content is converted to units of moles per unit volume.
 * @param Xdb Dry basis moisture content [-]
 * @param T Temperature [K]
 * @returns [J/m^3]
 */
double Esurf(double Xdb, double T)
{
    double Mw = 0.0180153, /* kg/mol */
           rhos,
           Xmin = 1e-9, /* Since binding energy is undefined at Xdb = 0, use a
                           larger minimum moisture content instead. */
           dX, /* Size of each integration slice */
           X, /* Current moisture content for calculations */
           E = 0; /* Set the energy to zero initially. */
    int nslice = 100, /* Number of slices to integrate with */
        i; /* Loop index */
    oswin *d;
    choi_okos *co;

    /* Get the property data we need */
    d = CreateOswinData();
    co = CreateChoiOkos(PASTACOMP);
    
    /* Calculate slice width */
    dX = (Xdb-Xmin)/nslice;
    /* Integrate from Xmin to Xdb */
    for(i=0; i<nslice; i++) {
        X = Xmin+i*dX;
        E += BindingEnergyOswin(d, X, T) * dX;
    }

    /* Calculate solid density for unit conversion */
    rhos = rho(co, T);
    /* Ensure the units are consitent and return the result. */
    return E*Xdb*rhos/Mw;
}

double GradEsurf(double x, double t, drydat cond)
{
    double h = 1e-6, El, Eh, gE, Xl, Xh;

    if(x-h < 0) {
        Xl = CrankEquationFx(x, t, cond);
        Xh = CrankEquationFx(x+h, t, cond);
        gE = (Esurf(Xh, cond.T) - Esurf(Xl, cond.T))/h;
    } else if(x+h > cond.L) {
        Xl = CrankEquationFx(x-h, t, cond);
        Xh = CrankEquationFx(x, t, cond);
        gE = (Esurf(Xh, cond.T) - Esurf(Xl, cond.T))/h;
    } else {
        Xl = CrankEquationFx(x-h, t, cond);
        Xh = CrankEquationFx(x+h, t, cond);
        gE = (Esurf(Xh, cond.T) - Esurf(Xl, cond.T))/(2*h);
    } 

    return gE;
}

double PoreP(double x, double t, drydat cond)
{
    double Xdb;
    Xdb = CrankEquationFx(x, t, cond);
    return pore_press(Xdb, cond.T);
}

