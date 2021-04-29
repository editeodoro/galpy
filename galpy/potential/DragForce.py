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
    from astropy import units


#####GENERAL CONSTANTS
Msolar=1.989e33
km=1.e5
pc=3.086e18
Kpc=3.086e21
yr=3.154e7
Gyr=1.e9*yr
kmpers=1.e5

mp=1.673e-24/Msolar
kb=1.381e-16/Msolar/kmpers**2
G=6.673e-8/pc/kmpers**2*Msolar

####CONSTANTS
mu=0.6
T_cloud=2.e4



class DragForce(DissipativeForce):
    """Class that implements the Drag force
    """
    def __init__(self,amp=1.,Mcl=None, dens=None, mp=mp, G=G, rzero=None,r_s=None,kb=kb,n_zero=None,rozero=None, alfa=0., temp=None, T0=None, rvir=None, gamma=None,
                 minr=0.0001*units.kpc, ro=None,vo=None):
        """
        NAME:

           __init__

        PURPOSE:

           initialize drag force


        OUTPUT:

           (none)

        HISTORY:

           2019-10-21 - Started - Afruni (UofGr)

        """

        DissipativeForce.__init__(self,amp=amp,ro=ro,vo=vo,
                                  amp_units=None)

        ########CONVERSION TO GALPY UNITS

        if _APY_LOADED and isinstance(minr,units.Quantity):
            minr= minr.to(units.kpc).value/self._ro
        if _APY_LOADED and isinstance(Mcl,units.Quantity):
            Mcl= Mcl.to(units.Msun).value/conversion.mass_in_msol(self._vo,self._ro)

        self._ms= Mcl        
        self._minr= minr
        self._mp=mp/conversion.mass_in_msol(self._vo,self._ro)
        self._G=G/conversion._G
        self._rzero=rzero/self._ro
        self._r_s=r_s/self._ro
        self._rvir=rvir/self._ro
        self._kb=kb/conversion.mass_in_msol(self._vo,self._ro)/self._vo**2
        self._n_zero=n_zero/conversion.dens_in_msolpc3(self._vo,self._ro)
        self._rozero=rozero/conversion.dens_in_msolpc3(self._vo,self._ro)
        self._alfa=alfa/conversion.freq_in_Gyr(self._vo,self._ro)
        self._T0=T0
        self._gamma=gamma
        
        ####################################
        
        ##########CORONA TEMPERATURE AND DENSITY##############
        self._temp= lambda R,z,phi=0.,t=0.: self._T0*(1+(self._gamma-1)/self._gamma*(mu*self._mp/(self._kb*self._T0)*(4*np.pi*self._G*self._rozero*self._r_s**2*(np.log(1+np.sqrt((R)**2+(z)**2)/self._r_s)\
                                            /(np.sqrt((R)**2+(z)**2)/self._r_s)\
                                            -np.log(1+self._rzero/self._r_s)/(self._rzero/self._r_s)))))
        self._dens= lambda R,z,phi=0.,t=0.: (self._temp(R,z,phi=phi,t=t)/self._T0)**(1/(self._gamma-1))*2.1*self._temp(R,z,phi=phi,t=t)**(-2)


        self._amp*= 1
        self._force_hash= None
        return None
 

    def _calc_force(self,R,phi,z,v,t): 

		        
        r= np.sqrt(R**2.+z**2.)
        
                     
        if r < self._minr:
            self._cached_force= 0.
        elif r >= 1.5*self._rvir:  
            self._cached_force=0.
        else:
            vs= np.sqrt(v[0]**2.+(v[1])**2.+v[2]**2.)   
            mcloud= self._ms
            C_drag=np.pi*(3/4.)**2.  ###drag coefficient
            #############DRAG FORCE##########################
            self._cached_force=\
                -C_drag**(1/3.)*(T_cloud)**(2/3.)*(mu*self._n_zero*np.exp(self._alfa*t))**(1/3.)*mcloud**(-1./3.)*(self._dens(R,z,phi=phi,t=t)**(1./3.)*vs)
            ###################################################################


    def _Rforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _Rforce
        PURPOSE:
           evaluate the radial force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v= current velocity in cylindrical coordinates
        OUTPUT:
           the radial force
        HISTORY:
           2018-03-18 - Started - Bovy (UofT)
        """
        new_hash= hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t]))\
            .hexdigest()
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        
        force_R= self._cached_force*v[0]
        return force_R


    def _phiforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _phiforce
        PURPOSE:
           evaluate the azimuthal force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v= current velocity in cylindrical coordinates
        OUTPUT:
           the azimuthal force
        HISTORY:
           2018-03-18 - Started - Bovy (UofT)
        """
        new_hash= hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t]))\
            .hexdigest()
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        

        force_phi= self._cached_force*(v[1]*R)

        return force_phi
        
        
    def _zforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _zforce
        PURPOSE:
           evaluate the vertical force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v= current velocity in cylindrical coordinates
        OUTPUT:
           the vertical force
        HISTORY:
           2018-03-18 - Started - Bovy (UofT)
        """
        new_hash= hashlib.md5(np.array([R,phi,z,v[0],v[1],v[2],t]))\
            .hexdigest()
        
        
        if new_hash != self._force_hash:
            self._calc_force(R,phi,z,v,t)
        
        ####ACCELERATE INTEGRATION IF CLOUDS ESCAPE 1.5R_VIR OR COME BACK TO THE DISC 
        
        if z>=np.sqrt((1.5*self._rvir)**2-R**2):
            force_z=1e5 
        elif (z>=0.01/self._ro and z<np.sqrt((1.5*self._rvir)**2-R**2)):  
            force_z= self._cached_force*v[2]
        else:
            force_z=-1e5 
        
        ############################################################################
        
        return force_z


class ConstantVerticalForce(DissipativeForce):
    """Class that implements a constant force """
    def __init__(self,acc,ro=None,vo=None, amp_units=None):
        
        DissipativeForce.__init__(self,amp=1,ro=ro,vo=vo,amp_units=amp_units)
        
        if _APY_LOADED and isinstance(acc,units.Quantity):
            conv = self._vo**2/self._ro
            acc = acc.to(units.km**2/units.s**2/units.kpc).value/conv
        
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
        
        if _APY_LOADED and isinstance(acc,units.Quantity):
            conv = self._vo**2/self._ro
            acc = acc.to(units.km**2/units.s**2/units.kpc).value/conv
        
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
        
        if _APY_LOADED and isinstance(Mcl,units.Quantity):
            Mcl = Mcl.to(units.Msun).value/conversion.mass_in_msol(self._vo,self._ro)
        if _APY_LOADED and isinstance(Rcl,units.Quantity):
            Rcl = Rcl.to(units.kpc).value/self._ro
        if _APY_LOADED and isinstance(rhoh,units.Quantity):
            rhoh = rhoh.to(units.Msun/units.pc**3).value/conversion.dens_in_msolpc3(self._vo,self._ro)
        if _APY_LOADED and isinstance(vwind,units.Quantity):
            vwind = vwind.to(units.km/units.s).value/self._vo

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