!! Binary
!! Copyright (C)  2013  Marco Tazzari, Giuseppe Lodato
!! Version: 13.0
!! Last Modified: 22 May 2013
!! 
!! Aim: This program computes the evolution of an accretion disc interacting with a binary system. 
!!      It also computes the evolution of an accretion disc orbiting a single central object.
!! 
!! This program is free software: you can redistribute it and/or modify it under the terms of the 
!! GNU General Public License as published by the Free Software Foundation, either version 3 
!! of the License, or (at your option) any later version.
!! 
!! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!! See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with this program.  
!! If not, see <http://www.gnu.org/licenses/>.
!! 
!! If you make use of this code, please quote:
!!     Tazzari, Lodato, Estimating the fossil disc mass during supermassive black hole mergers: 
!!                      the importance of torque implementation, MNRAS, Volume 449, Issue 1, p.1118-1128
!!     DOI: 10.1093/mnras/stv352


module constants

	double precision, parameter :: pi = acos(-1.)	! greek pi !already defined pi_d double precision pi from Numerical recipes
	double precision, parameter :: twopi	= 2.0d0*pi
	double precision, parameter :: threepi	= 3.0d0*pi
	double precision, parameter :: fourpi	= 4.0d0*pi
	double precision, parameter :: sixpi	= 6.0d0*pi

	double precision, parameter :: tiny 				= 1.d-16	
	double precision, parameter :: onethird				= 1./3.0d0
	double precision, parameter :: twothirds			= 2./3.0d0
	double precision, parameter :: two_at_twothirds		= 2.**(2./3.)
	double precision, parameter :: onehalf_at_twothirds	= 2.**(-2./3.)
	double precision, parameter :: oneseventh 			= 1./7.0d0
	double precision, parameter :: twothirteenths 		= 2./13.0d0
	
	
	double precision :: k_b,mu,m_p,Grav,M_sun,AU,sig_sb,h_planck,c_light,c_light2,eta,sigma_T
	double precision :: R_sun,year,L_sun, M_Jup
	double precision :: k_T,k_ff_0,a_rad,M_sun_year,km_s
	double precision :: parsec
	
	! first some general constants
	parameter(k_b      = 1.380658d-16)			! Boltzmann constant in erg/K
	parameter(m_p      = 1.6726231d-24)			! proton mass in g
	parameter(Grav     = 6.67259d-8)			! gravitational constant in cm^3 g^-1 s^-2
	parameter(M_sun    = 1.989d33)				! mass of the sun in g
	parameter(M_sun_year    = 6.3312d25)		! mdot: M_sun/year in g/s
	parameter(M_Jup    = 1.8986d30)             ! mass of jupiter in g
	parameter(AU       = 1.496d13)				! astronomical unit in cm
	parameter(sig_sb   = 5.670512d-5)			! Stefan-Boltzmann constant in g s^-3 K^-4
	parameter(h_planck = 6.626068d-27)			! Planck's constant in erg s
	parameter(c_light  = 2.99792458d10)			! speed of light in cm/s
	parameter(c_light2 = 8.98755179d20)			! square of the speed of light in cm^2/s^2
	parameter(year     = 31558149.54d0)			! year in s
	parameter(R_sun    = 0.0047d0*AU)         	! radius of the sun in cm
	parameter(L_sun    = 3.839d33)              ! solar luminosity in erg/s
	parameter(parsec   = 3.0856d18)				! parsec in cm
	parameter(k_T	   = 0.3977d0)				! Thomson opacity cm^2 g^-1
	parameter(sigma_T  = 6.6525d-25)			! Thomson cross-section cm^2
	parameter(k_ff_0   = 6.6e22)				! cm^2 g^-1	3.68d22 or 6.6e22 ? check!
	parameter(R_gas	   = 8.3145d7)				! erg K^-1 mol^-1
	parameter(mu 	   = 0.67d0)				! mean molecular mass in proton masses
	parameter(a_rad    = 7.565917d-15)			! a constant sig_sb=ac/4. 
	parameter(km_s	   = 1.0d3)					! km/s in m/s
	parameter(eta	   = 0.1d0)					! BH accretion efficiency
end module constants
