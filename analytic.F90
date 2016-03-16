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


! This file contains analytic solutions of the accretion discs.

!------------------------------------------------
!	SPREADING RING
!------------------------------------------------
double precision function spreading_ring(xin,x,tau)
	! Analytical solution taken from:
	!  ---> Lodato, G. Classical disc physics. New Astronomy Reviews 52, 21–41 (2008).
	! The total mass is zero because it is adimensionalized through Sigma(R,t) = Sigma0 * sigma(x,tau) and R=R0*x, t=R0^2/nu*tau
	
	! NOTE: the modified Bessel function I1/4(x) is a function that increases very rapidly.
	! 		for a good computation the arument x should not exceed 300 (I1/4(300)=1.e143 !!)
	!		for these reasons, a working regime is: rin=0.01, rout=2. xin=1. t0=0.001
	
	! Numerical Recipes modules; needed to compute the modified fractional Bessel Function of first kind.
	use nrtype
	use nrutil
	use nr, only : bessik
	!----------------------------------
	
	double precision,intent(in) :: tau		! adimensional time
	double precision,intent(in) :: x		! adimensional radial coordinate
	double precision,intent(in) :: xin		! initial (adimensional) position of the ring
	double precision :: argbessel			! argument of modified I_nu(x) function
	double precision :: xnu					! index nu: it can be either real or integer
	double precision :: bessel_factor		! in the sigma expression, the factor due to Bessel function
	double precision :: trash				! trash variable: bessik() returns other unwanted results.

	argbessel=xin*x/(6.*tau)
	xnu=1/4.
	call bessik(argbessel,xnu,bessel_factor,trash,trash,trash)
!	print *, xin,x,tau,argbessel,bessel_factor ! check !OK
	spreading_ring=(1./(12.*pi_d*tau))*((xin/x)**(1/4.))*exp(-(xin**2.+x**2.)/(12.*tau))*bessel_factor

end function spreading_ring