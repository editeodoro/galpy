#include <galpy_potentials.h>

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
                                       
  double forceAmplitude = potentialArgs->args[7];
  if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
      forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return 0;//forceAmplitude * vR;
}

double ConstantVerticalForcezforce(double R,double z,double phi,double t,
                                   struct potentialArg * potentialArgs,
                                   double vR,double vT,double vz){
                                       
  double forceAmplitude = potentialArgs->args[7];
  
  if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
    forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return forceAmplitude;//forceAmplitude * vz;
}

double ConstantVerticalForcephiforce(double R,double z,double phi,double t,
                                   struct potentialArg * potentialArgs,
                                   double vR,double vT,double vz){
                                       
  double forceAmplitude = potentialArgs->args[7];
  
  if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
    forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return 0;//forceAmplitude*vT*R;
}


/*
//////////////////////////////////////////////////////////////////////////////
//  CONSTANT WIND
//////////////////////////////////////////////////////////////////////////////
double ConstantWindAmplitude(double R,double z,double phi,double t,
                             struct potentialArg * potentialArgs,
                             double vR,double vT, double vz){
                                          
  double *args= potentialArgs->args;
  //Get args
  double Mcl   = args[8];
  double Rcl   = args[9];
  double rhoh  = args[10];
  double vwind = args[11];
  
  double vs = sqrt(vR**2.+vT**2.+(vwind-vz)*(vwind-vz));
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

double ConstantWindRforce(double R,double z,double phi,double t,
                          struct potentialArg * potentialArgs,
                          double vR,double vT,double vz){
                                       
  double forceAmplitude = potentialArgs->args[7];
  if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
      forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return 0;//forceAmplitude * vR;
}

double ConstantWindzforce(double R,double z,double phi,double t,
                                   struct potentialArg * potentialArgs,
                                   double vR,double vT,double vz){
                                       
  double forceAmplitude = potentialArgs->args[7];
  double vwind = potentialArgs->args[11];
  
  if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
    forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return forceAmplitude*(vwind-vz);
}


double ConstantWindphiforce(double R,double z,double phi,double t,
                            struct potentialArg * potentialArgs,
                            double vR,double vT,double vz){
                                       
  double forceAmplitude = potentialArgs->args[7];
  
  if (!isCached(potentialArgs->args,R,z,phi,t,vR,vT,vz))
    forceAmplitude = ConstantVerticalForceAmplitude(R,z,phi,t,potentialArgs,vR,vT,vz);
  
  return 0;//forceAmplitude*vT*R;
}
*/




