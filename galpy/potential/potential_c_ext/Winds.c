#include <galpy_potentials.h>
#include <math.h>
#include <stdio.h>
#include <Chevalier.h>

int isCached(double *args, double R,double z,double phi,double t,
             double vR,double vT, double vz) {

    return ( R==args[0] &&  z==args[1] && phi==args[2] && t==args[3] \
             && vR==args[4] && vT==args[5] && vz==args[6] );
}

//////////////////////////////////////////////////////////////////////////////
//  CONSTANT VERTICAL FORCE
//////////////////////////////////////////////////////////////////////////////
double ConstantVerticalForceAmplitude(double R,double z,double phi,double t,
                                      struct potentialArg * potentialArgs,
                                      double vR,double vT, double vz){
                                          
  double *args= potentialArgs->args;
  //Get args
  double forceAmplitude = args[8];
    
  // Caching
  args[0] = R;
  args[1] = z;
  args[2] = phi;
  args[3] = t;
  args[4] = vR;
  args[5] = vT;
  args[6] = vz;
  args[7] = forceAmplitude;
  return forceAmplitude;
}

double ConstantVerticalForceRforce(double R,double z,double phi,double t,
                                   struct potentialArg * potentialArgs,
                                   double vR,double vT,double vz){
                                       
  //double forceAmplitude = potentialArgs->args[7];
  //if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
  //  forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return 0;//forceAmplitude * vR;
}

double ConstantVerticalForcezforce(double R,double z,double phi,double t,
                                   struct potentialArg * potentialArgs,
                                   double vR,double vT,double vz){
                                       
  double forceAmplitude = potentialArgs->args[7];
  
  if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
    forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return forceAmplitude;
}

double ConstantVerticalForcephiforce(double R,double z,double phi,double t,
                                     struct potentialArg * potentialArgs,
                                     double vR,double vT,double vz){
                                       
  //double forceAmplitude = potentialArgs->args[7];
  
  //if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
  //forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return 0;
}



//////////////////////////////////////////////////////////////////////////////
//  CONSTANT WIND
//////////////////////////////////////////////////////////////////////////////
double* calcWindVelocityComponents(double R, double z, double *args){
    
    double vwind = args[11];
    int isRadial = (int)(args[12]);
    static double vw[2];   // vR and vz components

    if (isRadial) {  // vwind is radial, e.g. V(r)
        double theta = atan2(R,z);
        vw[0] = vwind*sin(theta);
        vw[1] = vwind*cos(theta);
    }
    else {           // vwind is vertical, e.g. V(z)
        vw[0] = 0.;
        vw[1] = z >= 0 ? vwind : -vwind;
    }
    return vw;
}

double ConstantWindAmplitude(double R,double z,double phi,double t,
                             double *args, double vR, double vT,
                             double vz, double vw_R, double vw_z){

  //Get args
  double Mcl   = args[8];
  double Rcl   = args[9];
  double rhoh  = args[10];
  
  double vs = sqrt((vR-vw_R)*(vR-vw_R)+vT*vT+(vw_z-vz)*(vw_z-vz));
  double C = M_PI*Rcl*Rcl*rhoh/Mcl;
  double forceAmplitude = C*vs;
  
  // Caching
  args[0] = R;
  args[1] = z;
  args[2] = phi;
  args[3] = t;
  args[4] = vR;
  args[5] = vT;
  args[6] = vz;
  args[7] = forceAmplitude;
  return forceAmplitude;

}

double ConstantWindRforce(double R, double z, double phi, double t,
                          struct potentialArg *pArgs,
                          double vR, double vT, double vz){
                    
    double *vwp = calcWindVelocityComponents(R,z,pArgs->args);
    double forceAmplitude = pArgs->args[7];
    if (!isCached(pArgs->args,R,z,phi,t,vR,vT,vz)) 
        forceAmplitude = ConstantWindAmplitude(R,z,phi,t,pArgs->args,vR,vT,vz,vwp[0],vwp[1]);
    double forceR = forceAmplitude*(vwp[0]-vR);
    
    return forceR;
}

double ConstantWindzforce(double R, double z, double phi, double t,
                          struct potentialArg *pArgs,
                          double vR, double vT, double vz){

  double *vwp = calcWindVelocityComponents(R,z,pArgs->args);
  double forceAmplitude = pArgs->args[7];
  if (!isCached(pArgs->args,R,z,phi,t,vR,vT,vz))
    forceAmplitude = ConstantWindAmplitude(R,z,phi,t,pArgs->args,vR,vT,vz,vwp[0],vwp[1]);
  return forceAmplitude*(vwp[1]-vz);
}


double ConstantWindphiforce(double R, double z, double phi, double t,
                            struct potentialArg *pArgs,
                            double vR, double vT, double vz){

  double *vwp = calcWindVelocityComponents(R,z,pArgs->args);
  double forceAmplitude = pArgs->args[7];
  if (!isCached(pArgs->args,R,z,phi,t,vR,vT,vz))
    forceAmplitude = ConstantWindAmplitude(R,z,phi,t,pArgs->args,vR,vT,vz,vwp[0],vwp[1]);
  return forceAmplitude*(-vT*R);
}


//////////////////////////////////////////////////////////////////////////////
// CC85 WIND
//////////////////////////////////////////////////////////////////////////////

Chevalier CC85; 

void defineCC85model(double *args) {
    
    double gamma = args[10];
    double scale_radius = args[11]*1000;
    double Edot = args[12];
    double Mdot = args[13];
    
    int alreadyset = CC85.gamma==gamma && CC85.R==scale_radius && 
                     CC85.E_dot==Edot && CC85.M_dot==Mdot;
    
    if (!alreadyset) SetChevalier(&CC85,Mdot,Edot,gamma,scale_radius);
}


double WindDensity(double R, double z, double *args){
    
    double R0 = args[15];
    double V0 = args[16];
    double r = sqrt(R*R+z*z)*R0*1000;
    double density = Density_Msunpc3(&CC85, r);  // Msun/pc3
    // Converting to galpy units
    double G = 4.302*1E-03;                      // pc / Msun (km/s)^2
    double dens_in_msolpc3 = V0*V0/R0/R0/G*1E-06;
    return density / dens_in_msolpc3;
}


double* WindVelocityComponents(double R, double z, double *args){
    
    double R0 = args[15]*1000;
    double V0 = args[16];
    double r = sqrt(R*R + z*z)*R0;
    double vwind = WindVelocity(&CC85,r)/V0;
    int isRadial = (int)(args[14]);
    static double vw[2];   // vR and vz components

    if (isRadial) {  // vwind is radial, e.g. V(r)
        double theta = atan2(R,z);
        vw[0] = vwind*sin(theta);
        vw[1] = vwind*cos(theta);
    }
    else {           // vwind is vertical, e.g. V(z)
        vw[0] = 0.;
        vw[1] = z >= 0 ? vwind : -vwind;
    }
    return vw;
}


double CC85WindAmplitude(double R,double z,double phi,double t,
                         double *args, double vR, double vT,
                         double vz, double vw_R, double vw_z){

  //Get args
  double cloud_radius = args[8];
  double cloud_mass = args[9];
  double density = WindDensity(R,z,args);
  double vs = sqrt((vR-vw_R)*(vR-vw_R)+vT*vT+(vw_z-vz)*(vw_z-vz));  
  double C = M_PI*cloud_radius*cloud_radius*density/cloud_mass;
  double forceAmplitude = C*vs;
  
  // Caching
  args[0] = R;
  args[1] = z;
  args[2] = phi;
  args[3] = t;
  args[4] = vR;
  args[5] = vT;
  args[6] = vz;
  args[7] = forceAmplitude;

  return forceAmplitude;

}


double CC85WindRforce(double R, double z, double phi, double t, struct potentialArg *pArgs,
                      double vR, double vT, double vz){
        
    defineCC85model(pArgs->args);
    double *vwp = WindVelocityComponents(R,z,pArgs->args);
    double forceAmplitude = pArgs->args[7];
    if (!isCached(pArgs->args,R,z,phi,t,vR,vT,vz)) 
        forceAmplitude = CC85WindAmplitude(R,z,phi,t,pArgs->args,vR,vT,vz,vwp[0],vwp[1]);
    double forceR = forceAmplitude*(vwp[0]-vR);
    return forceR;
}


double CC85Windzforce(double R, double z, double phi, double t, struct potentialArg *pArgs,
                      double vR, double vT, double vz){
    defineCC85model(pArgs->args);
    double *vwp = WindVelocityComponents(R,z,pArgs->args);
    double forceAmplitude = pArgs->args[7];
    if (!isCached(pArgs->args,R,z,phi,t,vR,vT,vz))
        forceAmplitude = CC85WindAmplitude(R,z,phi,t,pArgs->args,vR,vT,vz,vwp[0],vwp[1]);
    double forceZ = forceAmplitude*(vwp[1]-vz);
    return forceZ;
}


double CC85Windphiforce(double R, double z, double phi, double t, struct potentialArg *pArgs,
                        double vR, double vT, double vz){
    defineCC85model(pArgs->args);
    double *vwp = WindVelocityComponents(R,z,pArgs->args);
    double forceAmplitude = pArgs->args[7];
    if (!isCached(pArgs->args,R,z,phi,t,vR,vT,vz))
        forceAmplitude = CC85WindAmplitude(R,z,phi,t,pArgs->args,vR,vT,vz,vwp[0],vwp[1]);
    double forcephi = forceAmplitude*(-vT*R);

  return forcephi;
}





