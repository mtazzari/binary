%.o : %.F90	
	gfortran -c -D__PRINGLE__ $< $
	
binary: variables.o constants.o routines.o binary.o
	gfortran $^ $ -D__PRINGLE__ -o $@ 

gas:	variables.F90 constants.F90 routines.F90 binary.F90
	gfortran -c -D__ONLYPGAS__ variables.F90
	gfortran -c -D__ONLYPGAS__ constants.F90
	gfortran -c -D__ONLYPGAS__ routines.F90
	gfortran -c -D__ONLYPGAS__ binary.F90
	gfortran variables.o constants.o routines.o binary.o -D__ONLYPGAS__ -o gasbinary
	
	
chang:	variables.F90 constants.F90 routines.F90 binary.F90
	gfortran -c -D__ENERGYCHANG__ variables.F90
	gfortran -c -D__ENERGYCHANG__ constants.F90
	gfortran -c -D__ENERGYCHANG__ routines.F90
	gfortran -c -D__ENERGYCHANG__ binary.F90
	gfortran variables.o constants.o routines.o binary.o -D__ENERGYCHANG__ -o changbinary

clean:
	rm variables.o constants.o routines.o binary.o variables.mod units.mod formats.mod constants.mod

cleannr:
	rm variables.o constants.o routines.o analytic.o binary.o nrtype.o nr.o nrutil.o variables.mod units.mod formats.mod constants.mod nrtype.mod nr.mod nrutil.mod chebev.o bessik.o beschb.o
	
cleansnap:
	rm Snap*
	
cleandat:
	rm *.dat	
	
launch:	launch.o variables.o
	gfortran $^ $ -o $@ 