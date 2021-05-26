#ifndef CHEVALIER_H
#define CHEVALIER_H

#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#define parsec_in_cm 3.085678e+18 //parsec in cm
#define parsec_in_km 3.085678e+13 //parsec in km
#define msun_in_g    1.988920e+33	//msun in g
#define year_in_sec  3.155760e+07	//year in s


typedef struct Chevalier
{
    double M_dot;               //mass input rate
    double E_dot;               //energy input rate
    double R;                   //limit of mass and energy input
    double gamma;               //adiabatic index
    double V;                   //Volume of input area
    double q;                   //M_dot / V
    double Q;                   //E_dot / V
    gsl_root_fsolver *grsolve;  //GSL root solver
} Chevalier;

// Function declarations for Chevalier struct
void SetChevalier(Chevalier *c, double M_dot, double E_dot, double gamma, double R);
double MachNumber(Chevalier *c, double r);      //mach number vs. radius
double MomentumDensity(Chevalier *c, double r); //momentum density vs. radius in Msun/pc^3 * km/s
double Pressure(Chevalier *c, double r);        //Pressure in dyn cm^-2
double Density(Chevalier *c, double r);         //density in g cm^-3
double Density_Msunpc3(Chevalier *c, double r); //density in Msun/pc3
double WindVelocity(Chevalier *c, double r);    //wind velocity in km/s
double EnergyIntegral(Chevalier *c, double r);  //msun/yr (km/s)^2 / pc^3
double u_star(Chevalier *c, double r);          //normalized wind velocity
double rho_star(Chevalier *c, double r);        //normalized density
double P_star(Chevalier *c, double r);          //normalized pressure

double mach_crossing_A(double M, void *fp);     //functions for Mach number root finding
double mach_crossing_B(double M, void *fp);

// Function definitions
void SetChevalier(Chevalier *c, double M_dot, double E_dot, double gamma, double R) {
    // M_dot = Mass input in Msun/yr
    // E_dot = Energy input in erg / s
    // gamma = Adiabatic index
    // R     = Radius of the wind region in parsec
    
    //Set model parameters
    c->M_dot = M_dot;
    c->E_dot = E_dot;
    c->gamma = gamma;
    c->R     = R;
    c->V     = (4.*M_PI/3. * R*R*R); //volume in pc^3
    //compute input density rate (msun / yr / pc^3)
    c->q     = c->M_dot / c->V; 
    //compute input energy density rate (msun/yr (km/s)^2 / pc^3)
    c->Q     = (c->E_dot/msun_in_g*year_in_sec / 1.0e10) / c->V;
    //define the GSL root solver to use
    const gsl_root_fsolver_type *grsolve_T = gsl_root_fsolver_brent;
    c->grsolve = gsl_root_fsolver_alloc(grsolve_T);
}


double u_star(Chevalier *c, double r) {
    // Returns dimensionless wind velocity
    double u = WindVelocity(c, r);
    double u_prime = sqrt(c->E_dot)/sqrt(c->M_dot*msun_in_g/year_in_sec)/1.0e5;
    return u/u_prime;
}


double rho_star(Chevalier *c, double r) {
    // Returns dimensionless density
    double rho = Density(c,r); //g cm^-3
    double M_dot_prime = c->M_dot * msun_in_g/year_in_sec; // g/s
    double R_prime = c->R * parsec_in_cm; //cm
    double rho_prime = pow(M_dot_prime,1.5)/(sqrt(c->E_dot)*R_prime*R_prime);
    return rho/rho_prime;
}


double P_star(Chevalier *c, double r) {
    // Returns dimensionless pressure
    double P = Pressure(c,r); //dyn cm^-2
    double M_dot_prime = c->M_dot * msun_in_g/year_in_sec; // g/s
    double R_prime = c->R * parsec_in_cm; //cm
    double P_prime = sqrt(M_dot_prime)*sqrt(c->E_dot)/(R_prime*R_prime);
    return P/P_prime;
}


double Pressure(Chevalier *c, double r) {
    // Returns pressure in dyn cm^-2
    double Mach = MachNumber(c,r);          //mach number
    double u    = WindVelocity(c,r);        //km/s
    double cs   = u/Mach * 1.0e5;           //sound speed in cm/s
    u *= year_in_sec / parsec_in_km;        //to pc/yr
    double rho  = MomentumDensity(c,r)/u;   //Msun / pc^3
    rho *= msun_in_g;                       //g / pc^3
    rho /= pow(parsec_in_cm, 3);            //g / cm^3
    double P = cs*cs*rho/c->gamma;          //pressure in cgs
    return P;                               //in dyn cm^-2
}


double Density(Chevalier *c, double r) {
    // Returns density in g cm^-3
    double rho = Density_Msunpc3(c, r);     //Msun/pc^3
    rho *= msun_in_g;                       //g / pc^3
    rho /= pow(parsec_in_cm, 3);            //g / cm^3
    return rho;                             // in g/cm^3
}


double Density_Msunpc3(Chevalier *c, double r) {
    // Returns density in g cm^-3
    double u = WindVelocity(c, r);          //km/s
    u *= year_in_sec / parsec_in_km;        //to pc/yr
    return MomentumDensity(c,r)/u;          //Msun/pc^3;
}


double WindVelocity(Chevalier *c, double r) {
    // Returns wind velocity in km/s

    //First, find integrated energy input density
    double Qint = EnergyIntegral(c, r);     //msun/yr (km/s)^2 / pc^3
    //Second, find momentum density
    double rhou = MomentumDensity(c, r);     //msun/yr/pc^2
    //find sq of velocity (modulo gamma + Mach correction)
    double usq = Qint/rhou;                 //1/2 u^2 + (gamma/(gamma-1))*P/rho
                                            //1/2 u^2 + c^2 / (gamma-1)
                                            //1/2 u^2 + u^2 / (M^2 (gamma-1))
                                            //u^2 * ( 1/2 + 1/(M^2 (gamma-1)) )
    //get the mach number
    double Mach = MachNumber(c,r);
    //find the adjustment factor
    double fac = (0.5 + 1./(Mach*Mach*(c->gamma-1.)));
    return sqrt(usq/fac);	//km/s
}


double EnergyIntegral(Chevalier *c, double r) {
    double Qint = r<c->R ? 1./3.*c->Q*r : 1./3.*c->Q*c->R*c->R*c->R/(r*r);
    return Qint;    //msun/yr (km/s)^2 / pc^3
}


double MomentumDensity(Chevalier *c, double r) {
    double rhou = r<c->R ? 1./3.*c->q*r : 1./3.*c->q*c->R*c->R*c->R/(r*r) ;
    return rhou;    // in Msun/yr/pc^2
}


double MachNumber(Chevalier *c, double r) {
    
    int    status, iter=0, max_iter=100;
    double M_lo = 1.0e-5, M_hi = 5.;
    double Mx, answer = 0;
    double x = r/c->R;
    gsl_function func;
    double fp[2] = {c->gamma, x};

    //choose which solution to use
    if(x<=1.0) {
        M_lo = 1.0e-5;
        M_hi = 1.0;
        func.function = &mach_crossing_A;
    }
    else {
        M_lo = 1.0;
        M_hi = 10000.0;
        func.function = &mach_crossing_B;
    }
    func.params   = &fp[0];

    gsl_root_fsolver_set(c->grsolve, &func, M_lo, M_hi);

    while(++iter<max_iter){
        status = gsl_root_fsolver_iterate(c->grsolve);
        Mx = gsl_root_fsolver_root(c->grsolve);
        M_lo = gsl_root_fsolver_x_lower(c->grsolve);
        M_hi = gsl_root_fsolver_x_upper(c->grsolve);
        status = gsl_root_test_interval(M_lo,M_hi,0,1.0e-5);
        if(status==GSL_SUCCESS) {answer = Mx; break;}
    }

    return answer;
}


double mach_crossing_A(double M, void *fp) {
    double *g     = (double *) fp;
    double gamma  = g[0];
    double x      = g[1];
    double alpha  = -1.*(3.*gamma+1.)/(5.*gamma+1.);
    double beta   = (gamma+1.)/(2.*(5.*gamma+1.));
    double A      = pow( (3.*gamma + 1./(M*M))/(1.+3.*gamma), alpha);
    double B      = pow( (gamma-1.+2./(M*M))/(1.+gamma), beta);
    return A*B - x;
}


double mach_crossing_B(double M, void *fp) {
    double *g     = (double *) fp;
    double gamma = g[0];
    double x     = g[1];
    double alpha = 2./(gamma-1.);
    double beta  = (gamma+1.)/(2.*(gamma-1.));
    double A     = pow( M, alpha);
    double B     = pow( (gamma-1.+2./(M*M))/(1.+gamma), beta);
    return A*B - x*x;
}

#endif //CHEVALIER_H