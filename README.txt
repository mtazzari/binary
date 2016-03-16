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


Compile with: 
    make
Default option: __PRINGLE__ : uses J. Pringle prescription (as in Lodato+2009) to solve the energy balance equation

Other options:
1) __ONLYPGAS__ : considers only gas pressure (neglects radiation pressure)
2) __ENERGYCHANG__: considers only gas pressure and assuemes fixed opacity (as in Chang+2010)
To compile:
Option 1:
	make gas
Option 2:
	make chang


Below, an annotated sample input file:

&INPUTVARIABLES
 SIM_MODE='chang09sm1hbis            ',			# simulation name (choose freely) 
 INR=        200,								# radial grid, number of cells
 GRIDTYPE=          1,							# grid type: 0:linear  1: log 
 T0=  0.0000000000000000     ,					# initial computation time
 TOT_TIME=  3000000.0000     ,					# final computation time (the code stops when this time is reached)
 DTAUSNAP=  100000.00000000000     ,			# time between two snapshots [yr] (files in snap/ directory containing the disk structure)
 DTAUSNAP_PARAM=  10000.0000000000000     ,		# time between two snapshots of the disk parameters [yr] (in the mass.dat file; usually this is smaller than DTAUSNAP)
 MP=  10000000.000000000     ,					# mass of the primay BH [Msun] (or of the central object) 
 Q= 0.10000000000000001     ,					# mass ratio q=M2/Mp
 ROUT=  100000.00000000000     ,				# radial grid, outer radius [ROUT_UNIT]
 ROUT_UNIT=          2,							# unit of the radial grid: 1:parsec   2:Rg=GM/c^2   3:see routines.F90, dimensional_variables()  
 A0=  10000.000000000000     ,					# initial binary separation A0 = a(t=0)
 F=  1.0000000000000000E-002,					# torque parameter f
 INITIAL_SIGMA_TYPE=          8,				# initial surface density profile: see routines.F90, set_initialsigma()
 INITIAL_SIGMA_FILE='steady.dat     ',			# input file for the initial surface density
 SIGMADOT_TYPE=          2,						# type of inflow rate at the grid outer radius: 0:no accretion 2:constant inflow at a rate MDOTIN (see routines.F90, set_sigmadot())
 MDOTIN= 1.607904d-07     ,						# inflow rate at the grid outer radius [Msun/year]
 M_DISC0=  1000.0000000000000     ,				# initial disk mass [Msun] (note that this is relevant only for some choices of INITIAL_SIGMA_TYPE)
 ISPLANET=          1,							# 0:no secondary BH  1:secondary BH is present
 TINJECTPLANET=  0.0000000000000000     ,		# time at which the secondary BH is injected in the disk [yr]
 TMIGRATION=  0.0000000000000000     ,			# time after which the migration BH is switched ON [yr] (see routines.F90, dimensional_variables())
 ALPHA= 0.10000000000000001     ,				# alpha parameter from Shakura-Sunyaev prescripion
 /

Note:
For INITIAL_SIGMA_TYPE = 5, 8, 9, M_DISC0 is used to re-normalize the disk surface density to match M_DISC0 initial mass.
This is useful, e.g. to reproduce Lodato+2009, Chang+2009, Armitage+2002 simulations.
It is also useful to start the simulation from a steady-state surface density profile (previously obtained with a very long simulation at a fixed MDOTIN, letting it evlove until it reached steady state).