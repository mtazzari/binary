Documentation
=============

How to compile
--------------
Compile with:
``` 
    make
```

Default option: __PRINGLE__ : uses J. Pringle prescription (as in Lodato+2009) to solve the energy balance equation

Other options:

1) __ONLYPGAS__ : considers only gas pressure (neglects radiation pressure). Compile with:
```
	make gas
```

2) __ENERGYCHANG__: considers only gas pressure and assuemes fixed opacity (as in Chang+2010). Compile with:
```
	make chang
```

Input file
----------
Below, an annotated sample input file:

```fortran
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
```

Note:
For INITIAL_SIGMA_TYPE = 5, 8, 9, the value of M_DISC0 is used to re-normalize the disk surface density to match M_DISC0 initial mass.

This is useful, e.g. to reproduce Lodato+2009, Chang+2009, Armitage+2002 simulations.

It is also useful to start the simulation from a steady-state surface density profile (previously obtained with a very long simulation at a fixed MDOTIN, letting it evlove until it reached steady state).

Output files
-----
Each execution creates a directory named as SIM_MODE with the following content:

* code/ 		  		# directory containing the code and the input file used for the run (for reproducibility)
* snap/					# directory containing the snapshots of the disk, one snapshot every DTAUSNAP
    * snapshot_01.dat # snapshot of the disk (quantities as a function of disc radius)
* snap_latest/			# directory containing the snapshots of the disk, one snapshot every time the binary separation decreases by 1 Rg.
    * snapshot_01.dat # same as above
* mass.dat				# snapshot of the disk structure throughout the run (one line per each snapshot)
* sigma_neg.dat: 		# file containing a trace of points where sigma < 0 (just a sanity check for consistency; should be empty)
* splash.columns        # splash file with column names
* splash.defaults		# splash file with default parameters
* splash.filenames		# splash file with filenames of the snapshots to be shown together
* splash.limits			# splash file with plot limits 

To have a look of the output snapshots, just change directory into SIM_MODE/ and do:
```bash
	splash
```

Since the splash.* files are already provided in the directory, it should be already set up with proper labels and limits.

Here more details on the output:
File mass.dat:
```fortran
	open(unit=u_massfile, file=massfilename,access='append')
	write(u_massfile,f_massfile) 		isncount		,&				! 1. snap number
									&	t/year			,&				! 2. time
									&	mass/M_sun		,&				! 3. disc mass
									&	mass_inner/M_sun,&				! 4. inner disc mass
									&	tot_acc_mass/M_sun		,&		! 5. total mass accreted
									&	tot_sim_mass/M_sun		,&		! 6. total mass of the simulation (must be constant)
									& 	mdot(1)/(M_sun_year)	,&		! 7. mdot at R_in
									&	mdot(inr)/(M_sun_year)	,&		! 8. mdot at R_out
									& 	a				,&				! 9. secondary's position in cm
									&	a/parsec		,&				!10. secondary's position in parsec
									&	a/rg			,&				!11. secondary's position in rg
									&	adot1			,&				!12. acceleration due to viscous torque
									&	adot2			,&				!13. acceleration due to gw torque
									&	Acc_L/L_Edd						!14. Luminosity / Eddington Luminosity
	close(u_massfile)	
```
	
File snapshot_###.dat:
```fortran
	write(u_snapshots_latest,f_snapshots) 	xrg(i)				,&		!  1. radial coordinate in rg
									&	max(prec_out,sigma(i))		,&		!  2. sigma
									&	max(prec_out,nu(i))			,&		!  3. viscosity
									&	max(prec_out,T_eff(i))		,&		!  4. effective temperature
									&	max(prec_out,T_c(i))		,&		!  5. central temperature
									&	max(prec_out,H(i)/x(2*i))	,&		!  6. thickness/R
									& 	lambda(i)					,&		!  7. torque
									&	xpc(i)						,&		!  8. radial coordinate in parsec
									&	prad						,&		!  9. radiation pressure
									&	pgas						,&		! 10. gas pressure
									&	prad_pgas					,&		! 11. radiation to gas pressure ratio
									&	tnu									! 12. viscous timescale									
	enddo
	write(u_snapshots_latest,f_snapshots) 	xrg(inr),prec_out,prec_out,prec_out,& 
									& prec_out,prec_out,prec_out,xpc(inr),prec_out,prec_out,prec_out,prec_out		
	close(u_snapshots_latest)
```
	