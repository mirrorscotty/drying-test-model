#include <math.h>

#include "matrix.h"
#include "drying.h"
#include "material-data.h"

#include <stdlib.h>
#include <stdio.h>

/* Number of terms to use for the Crank equation when solving for moisture
 * profile. */

/**
 * Calculate the maximum possible strain assuming that none of the water leaving
 * the pores is replaced by air. The strain returned here is the volumetric
 * strain, and the densities are based on the Choi-Okos equations.
 * @param t Time [s]
 * @param x Position [m]
 * @param drydat Set of drying parameters
 * @returns Volumetric strain [-]
 */
double maxstrain(double t, double x, drydat d)
{
    double Xdb = CrankEquationFx(x, t, d),
           rhow, rhos,
           vs, vw;
    choi_okos *co;

        co = CreateChoiOkos(WATERCOMP);
        rhow = rho(co, d.T);
        DestroyChoiOkos(co);

        co = CreateChoiOkos(PASTACOMP);
        rhos = rho(co, d.T);
        DestroyChoiOkos(co);

        vs = 1/(1+rhos*d.X0/rhow);
        vw = Xdb*rhos/rhow * vs;

        return (vs+vw)-1;
}

/**
 * Calculate volumetric shrinkage due to water loss, assuming that a portion of
 * the water is replaced by air. The value of eta determines exactly how much
 * shrinkage occurs. This equation is from Mercier et al. 2011.
 * eta > 1: swelling
 *     = 1: no volume change
 *     < 1: shrinkage
 *     = 0: total shrinkage
 *     < 0: collapse
 * @param t Time [s]
 * @param x Position [m]
 * @param eta Shrinkage parameter [-]
 * @returns Volumetric strain [-]
 */
double etastrain(double t, double x, drydat d, double eta)
{
    double Xdb = CrankEquationFx(x, t, d),
           rhow, rho0,
           X0 = d.X0;
    choi_okos *co, *co0;

    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, d.T);
    DestroyChoiOkos(co);

    co = CreateChoiOkos(PASTACOMP);
    co0 = AddDryBasis(co, X0);
    rho0 = rho(co0, d.T);
    DestroyChoiOkos(co);
    DestroyChoiOkos(co0);
    
    return ((1-eta)*(X0-Xdb))/((rhow/rho0)*(1+X0));
}

double solidfrac(double t, double x, drydat d, double eta)
{
    double rhos, rhoapp, rhoapp0,
           X = CrankEquationFx(x, t, d),
           X0 = d.X0,
           chit = etastrain(t, x, d, eta);
    choi_okos *co, *co0;

    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, d.T);
    co0 = AddDryBasis(co, X0);
    rhoapp0 = rho(co0, d.T);
    DestroyChoiOkos(co);
    DestroyChoiOkos(co0);

    rhoapp = rhoapp0*(1+X)/((1+X0)*(1-chit));

    return rhoapp/(1+X) * 1/rhos;
}

