###############################################################################
#  DragForce: Class that implements the coronal drag force (by AA)
#  Link this file (with full path!) to galpy/potential:
#  
#  > ln -s $PWD/DragForce.py /usr/local/lib/python3.9/site-packages/galpy/potential
#
###############################################################################
import copy
import hashlib
import numpy as np
from scipy import special, interpolate
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
from .DissipativeForce import DissipativeForce
from .Potential import _APY_LOADED, evaluateDensities
from .Potential import flatten as flatten_pot
from ..util import conversion

if _APY_LOADED:
    import astropy.units as u


class ConstantVerticalForce(DissipativeForce):
    """Class that implements a constant force """
    def __init__(self,acc,ro=None,vo=None, amp_units=None):
        
        DissipativeForce.__init__(self,amp=1,ro=ro,vo=vo,amp_units=amp_units)
        
        if _APY_LOADED and isinstance(acc,u.Quantity):
            conv = self._vo**2/self._ro
            acc = acc.to(u.km**2/u.s**2/u.kpc).value/conv
        
        self.acc = acc
        self._cached_force = None
        self.hasC = True
        
    def _calc_force(self,R,phi,z,v,t): 
        self._cached_force = self.acc

    def _Rforce(self,R,z,phi=0.,t=0.,v=None):
        return 0

    def _phiforce(self,R,z,phi=0.,t=0.,v=None): 
        return 0
        
    def _zforce(self,R,z,phi=0.,t=0.,v=None):
        self._calc_force(R,phi,z,v,t)
        force_z = self._cached_force
        return force_z


class ConstantRadialForce(DissipativeForce):
    """Class that implements a constant force """
    def __init__(self,acc,ro=None,vo=None, amp_units=None):
        
        DissipativeForce.__init__(self,amp=1,ro=ro,vo=vo,amp_units=amp_units)
        
        if _APY_LOADED and isinstance(acc,u.Quantity):
            conv = self._vo**2/self._ro
            acc = acc.to(u.km**2/u.s**2/u.kpc).value/conv
        
        self.acc = acc
        self._cached_force = None
        self._force_hash= None

    
    def _calc_cylindrical(self,R,phi,z,Fr):
        theta = np.arctan2(R,z)
        Fz    = Fr*np.cos(theta)
        FR    = Fr*np.sin(theta)
        return (FR,0,Fz)
    
    
    def _calc_force(self,R,phi,z,v,t): 
        self._cached_force = self.acc
        self._force_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()


    def _Rforce(self,R,z,phi=0.,t=0.,v=None):
        new_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash: 
            self._calc_force(R,phi,z,v,t)
        return self._calc_cylindrical(R,phi,z,self._cached_force)[0]


    def _phiforce(self,R,z,phi=0.,t=0.,v=None):
        new_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash: 
            self._calc_force(R,phi,z,v,t)
        return self._calc_cylindrical(R,phi,z,self._cached_force)[1]

    def _zforce(self,R,z,phi=0.,t=0.,v=None):
        new_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash: 
            self._calc_force(R,phi,z,v,t)
        return self._calc_cylindrical(R,phi,z,self._cached_force)[2]



class ConstantWind(DissipativeForce):
    """Class that implements a constant vertical wind driven by drag"""
    def __init__(self,Rcl,Mcl,rhoh,vwind,isRadial=True,ro=None,vo=None):
        
        DissipativeForce.__init__(self,amp=1,ro=ro,vo=vo,amp_units=None)
        
        if _APY_LOADED and isinstance(Mcl,u.Quantity):
            Mcl = Mcl.to(u.Msun).value/conversion.mass_in_msol(self._vo,self._ro)
        if _APY_LOADED and isinstance(Rcl,u.Quantity):
            Rcl = Rcl.to(u.kpc).value/self._ro
        if _APY_LOADED and isinstance(rhoh,u.Quantity):
            rhoh = rhoh.to(u.Msun/u.pc**3).value/conversion.dens_in_msolpc3(self._vo,self._ro)
        if _APY_LOADED and isinstance(vwind,u.Quantity):
            vwind = vwind.to(u.km/u.s).value/self._vo

        self.Mcl = Mcl
        self.Rcl = Rcl
        self.rhoh = rhoh
        self.vwind = vwind
        self.isRadial = isRadial
        self._cached_force = None
        self._force_hash = None
        self.hasC = True


    def _calc_cylindrical(self,R,phi,z,Fr):
        
        if self.isRadial:
            # Vwind is radial, e.g. V(r)
            theta = np.arctan2(R,z)
            Fz    = Fr*np.cos(theta)
            FR    = Fr*np.sin(theta)
        else:
            # Vwind is vertical, e.g. V(z)
            Fz    = Fr if z>0 else -Fr
            FR    = 0
        return (FR,0,Fz)


    def _calc_force(self,R,phi,z,v,t):
        
        vw_R, _, vw_z = self._calc_cylindrical(R,phi,z,self.vwind)
        vs = np.sqrt((vw_R-v[0])**2.+v[1]**2.+(vw_z-v[2])**2.)
        c = np.pi*self.Rcl**2*self.rhoh/self.Mcl
        self._cached_force = c*vs
        self._force_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        

    def _Rforce(self,R,z,phi=0.,t=0.,v=None):

        vw_R, _, vw_z = self._calc_cylindrical(R,phi,z,self.vwind)
        new_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        force_R = self._cached_force*(vw_R-v[0])
        return force_R


    def _phiforce(self,R,z,phi=0.,t=0.,v=None): 
        
        new_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        force_phi = self._cached_force*(-v[1]*R)
        return force_phi
        
        
    def _zforce(self,R,z,phi=0.,t=0.,v=None):
        
        vw_R, _, vw_z = self._calc_cylindrical(R,phi,z,self.vwind)
        new_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        force_z = self._cached_force*(vw_z-v[2])
        return force_z



class CC85Wind(DissipativeForce):
    """ Class that implements a drag driven by a Chevalier&Clegg85 wind """
    
    
    # NB PYTHON VERSION NEEDS TO BE FIXED!

    def __init__(self,Rcl=1.,Mcl=1.,gamma=5/3,Rs=200*u.pc,Edot=1E43*u.erg/u.s, \
                 Mdot=1*u.Msun/u.yr,isRadial=True,ro=None,vo=None):

        DissipativeForce.__init__(self,amp=1,ro=ro,vo=vo)
        
        if _APY_LOADED and isinstance(Rcl,u.Quantity):
            Rcl = Rcl.to(u.kpc).value/(self._ro)
        if _APY_LOADED and isinstance(Mcl,u.Quantity):
            Mcl = Mcl.to(u.Msun).value/conversion.mass_in_msol(self._vo,self._ro)
        if _APY_LOADED and isinstance(Rs,u.Quantity):
            Rs = Rs.to(u.kpc).value/self._ro
        
        Edot_dim = Edot
        Mdot_dim = Mdot
        if _APY_LOADED and isinstance(Edot,u.Quantity):
            Edot_dim = Edot.to('erg/s').value
            Edot = Edot.to(u.Msun * (u.pc/u.Myr**2) * (u.kpc/u.Gyr)).value \
                /(conversion.mass_in_msol(self._vo,self._ro)*conversion.force_in_pcMyr2(self._vo,self._ro) \
                * conversion.velocity_in_kpcGyr(self._vo,self._ro))
        if _APY_LOADED and isinstance(Mdot,u.Quantity):
            Mdot_dim = Mdot.to('Msun/yr').value
            Mdot = Mdot.to(u.Msun/u.Gyr).value\
                    /(conversion.mass_in_msol(self._vo,self._ro) * conversion.time_in_Gyr(self._vo,self._ro))

        self._Rcl = Rcl
        self._Mcl = Mcl
        
        # Dimensional parameters of CC85 wind
        self._Edot_dim = Edot_dim           # erg/s
        self._Mdot_dim = Mdot_dim           # Msun/yr
        self._Rs_dim   = Rs*self._ro        # kpc
        # Dimensionless parameters of CC85 wind
        self._Edot = Edot
        self._Mdot = Mdot
        self._Rs = Rs
        self._gamma = gamma

        self.isRadial = isRadial

        self._force_hash = None
        self._cached_force = None
        self.hasC = True

        return None

    def Mach_estimator(self,M,r,R,gamma):

        if r<R:
            A = (3*gamma + (1/M**2))/(1 + 3*gamma)
            B = (gamma - 1 + (2/M**2))/(1 + gamma)
            alpha = (-3*gamma - 1)/(5*gamma + 1)
            beta = (gamma + 1)/(10*gamma + 2)
            value = A**(alpha) * B**(beta) - r/R
        else:
            alpha = 2/(gamma-1)
            B = (gamma-1 + (2/M**2))/(1+gamma)
            beta = (gamma+1)/(2*gamma - 2)
            value = M**(alpha) * B**(beta) - (r**2/R**2)
        return value 


    def Mach_number(self,r,R,gamma):
    
        if r<R:
            mach_value = optimize.fsolve(self.Mach_estimator,0.001,args=(r,R,gamma))
        else:
            mach_value = optimize.fsolve(self.Mach_estimator,2,args=(r,R,gamma))
    
        return mach_value


    def _wind_velocity(self,r,R,gamma):
    
        mach_value = self.Mach_number(r,R,gamma)[0]
        #print("Mach number is ",mach_value)
        velocity = ((2*mach_value**2)/(mach_value**2 + (2/(gamma-1))))**(1/2) * self._Mdot**(-1/2) * self._Edot**(1/2)
        return velocity


    def _calc_cylindrical(self,R,phi,z,Fr):
        
        if self.isRadial:
            # Vwind is radial, e.g. V(r)
            theta = np.arctan2(R,z)
            Fz    = Fr*np.cos(theta)
            FR    = Fr*np.sin(theta)
        else:
            # Vwind is vertical, e.g. V(z)
            Fz    = Fr if z>0 else -Fr
            FR    = 0
        return (FR,0,Fz)


    def _wind_density(self,r,R,gamma):
    
        mach_value = self.Mach_number(r,R,gamma)[0]
        if r < R:
            density = (r/R) * 1/(4*np.pi) * np.sqrt((mach_value**2 + 2/(gamma-1))/(2*mach_value**2)) * self._Mdot**(3/2) * self._Edot**(-1/2) * (self._scale_radius)**(-2)
        else:
            density = (4*np.pi*(r/R)**2)**(-1) * np.sqrt((mach_value**2 + 2/(gamma-1))/(2*mach_value**2)) * self._Mdot**(3/2) * self._Edot**(-1/2) * self._scale_radius**(-2)
        return density


    def _calc_force(self,R,phi,z,v,t):

        r_position = np.sqrt(R**2 + z**2)
        vw_r = self._wind_velocity(r_position,self._scale_radius,self._gamma)
        vw_R, _, vw_z = self._calc_cylindrical(R,phi,z,vw_r)
        
        print ("THIS IS PYTHON")
        print("r fraction is",r_position/self._scale_radius)
        print("time is ",t)
        print("vR, vT and vz are",v[0],v[1],v[2])

        vs = np.sqrt((vw_R-v[0])**2.+v[1]**2.+(vw_z-v[2])**2.)

        drag_force = np.pi *self._Rcl**2*self._wind_density(r_position,self._Rs,self._gamma)/self._Mcl

        self._cached_force = drag_force * vs
        self._force_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        


    def _Rforce(self,R,z,phi=0.,t=0.,v=None):

        r = np.sqrt(R**2 + z**2)
        vw_r = self._wind_velocity(r,self._scale_radius,self._gamma)
        vw_R, _, vw_z = self._calc_cylindrical(R,phi,z,vw_r)
        new_hash = hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        force_R = self._cached_force*(vw_R-v[0])
        return force_R

    def _phiforce(self,R,z,phi=0.,t=0.,v=None):

        new_hash= hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        force_phi = self._cached_force*(-v[1]*R)
        return force_phi

    def _zforce(self,R,z,phi=0.,t=0.,v=None):

        r = np.sqrt(R**2 + z**2)
        vw_r = self._wind_velocity(r,self._scale_radius,self._gamma)
        vw_R, _, vw_z = self._calc_cylindrical(R,phi,z,vw_r)
        new_hash= hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t])).hexdigest()
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        force_z = self._cached_force*(vw_z-v[2])
        return force_z
