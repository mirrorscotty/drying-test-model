#include "drying.h"
#include "matrix.h"
#include <math.h>

/**
 * Crank equation for a diffusion in a sheet.
 * \f[
 * \frac{X_{db}-X_e}{X_0-X_e}
 *     = \frac{8}{\pi^2}\sum_{n=0}^\infty\frac{1}{(2n+1)^2}
 *             \exp\left\{-k_F t (2n+1)^2\right\}
 * \f]
 * \f[
 * k_F = \frac{\pi^2 D}{l^2}
 * \f]
 * (Crank 1956)
 * @param kf Diffusivity constant (D*pi^2/l^2) where D is diffusivity and l is
 *      the slab thickness [1/s]
 * @param t Time [sec]
 * @param d Drying conditions
 * @returns Moisture content [kg/kg db]
 */
double CrankEquation(double t, drydat d)
{
    double kF; /* Dimensionless diffusion coefficient */
    double value = 0; /* Variable for summing up all the terms */
    int n; /* Current term */

    kF = M_PI*M_PI*d.D/(d.L*d.L);

    /* Calculate the value for each of the n terms and add them up */
    for(n=0; n<d.nterms; n++) {
        value += 8/(pow(2*n+1, 2) * M_PI*M_PI) *
            exp(-kF * t * pow(2*n+1, 2));
    }

    /* Solve for Xdb */
    value = value * (d.X0-d.Xe) + d.Xe;

    return value;
}

/**
 * Equation for sorption/desorption by a membrane
 * @param x X-coordinate in the membrane [m]
 * @param t Time [s]
 * @param d Drying parameters
 *
 * @returns Moisture content at the specified point in the slab [kg/kg db]
 */
double CrankEquationFx(double x, double t, drydat d)
{
    double value = 0; /* Variable for summing up all the terms */
    int n; /* Current term */
    double kF = d.D*M_PI*M_PI/(d.L*d.L),
           X1 = d.Xe,
           X0 = d.X0,
           term;

    /* Calculate the value for each of the n terms and add them up */
    for(n=0; n<d.nterms; n++) {
        term = pow(-1, n)/(2*n+1) * exp(-kF*pow(2*n+1, 2)*t/4)
            * cos( ((2*n+1)*M_PI*x)/(2*d.L) );
        value += term;
    }

    return (1 - 4/M_PI*value) * (X1-X0) + X0;
}

/**
 * Solve the Crank equation for kF using Newton's method. The kF value has the
 * following form:
 * \f[
 * k_F = \frac{\pi^2 D}{l^2}
 * \f]
 * @param t Time [s]
 * @param X Moisture content [kg/kg db]
 * @param X0 Initial moisture content [kg/kg db]
 * @param Xe Equilibrium moisture content [kg/kg db]
 * @returns kF value [1/s]
 */
//double CrankkF(double t, double X, double X0, double Xe, double beta0)
//{
//    double kf = beta0, /* Initial guess for kF */
//           kfp = 0, /* kF from the previous loop iteration. */
//           f, /* Value of the residual */
//           df, /* Derivative of f */
//           h = 1e-10, /* dx value used for differentiation */
//           tol = 1e-10; /* Tolerance for Newton's method */
//    int nterms = CONSTnterms; /* Number of terms to use */
//
//    /* Newton's method */
//    do {
//        /* Calculate f and df */
//        f = CrankEquation(kf, t, X0, Xe, nterms) - X;
//        df = (CrankEquation(kf+h, t, X0, Xe, nterms) 
//                - CrankEquation(kf-h, t, X0, Xe, nterms))/(2*h);
//        /* Set old kF value to kfp for checking convergence */
//        kfp = kf;
//        /* Calculate the new value of kF */
//        kf = kf - f/df;
//    } while(fabs(kfp - kf) > tol); /* Check convergence */
//
//    return kf;
//}

/**
 * Function to allow the Crank equation to be used in the fitnlm function
 * @param t Time [s]
 * @param beta 1x1 matrix containing the value for kF
 * @returns Moisture content [kg/kg db]
 */
//double CrankModel(double t, matrix *beta)
//{
//    double X0 = CONSTX0, /* Initial moisture content */
//           Xe = CONSTXe,  /* Equilibrium moisture content */
//           kf = val(beta, 0, 0); /* Get kF from the beta matrix */
//    int nterms = CONSTnterms; /* Number of terms to use */
//
//    return CrankEquation(kf, t, X0, Xe, nterms);
//}
//
