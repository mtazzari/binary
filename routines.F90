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


!------------------------------------------------
!	this subroutine sets the defaults values
! 	of variables reading from a file if exists,
! 	otherwise sets default values.								
!------------------------------------------------
subroutine read_variables()

	use units,		only :	u_mainoutput
	
	use variables,	only :	inputfile, nu0,INITIAL_SIGMA_FILE,								&
						 &	sim_mode, inr, gridtype, t0, tot_time, dtausnap, dtausnap_param,&
						 &	mp, q, rout, rout_unit,a0, f, initial_sigma_type, sigmadot_type,&
						 & 	mdotin, m_disc0, isplanet, tinjectplanet, tmigration,alpha
							
	implicit none
	integer :: ierror
	integer :: findsigmacell
	
	! defines the namelist of the input variables
namelist / inputvariables / &
		&	 SIM_MODE			, &
		&	 INR				, &
		&	 GRIDTYPE			, &
		&	 T0					, &
		&	 TOT_TIME			, &
		&	 DTAUSNAP			, &
		&	 DTAUSNAP_PARAM		, &
		&	 MP					, &
		&	 Q					, &
		&	 ROUT				, &
		&	 ROUT_UNIT			, &	
		&	 A0					, &
		&	 F					, &
		&	 INITIAL_SIGMA_TYPE	, &
		&	 INITIAL_SIGMA_FILE , &
		&	 SIGMADOT_TYPE		, &
		&	 MDOTIN				, &
		&	 M_DISC0			, &
		&	 ISPLANET			, &
		&	 TINJECTPLANET		, &
		&	 TMIGRATION			, &
		&	 ALPHA
		
	! set the default values
	 SIM_MODE	= 'defaults'	
	 INR		= 300				
	 GRIDTYPE	= 1					
	 T0			= 0.0d0     		
	 TOT_TIME	= 3.0d2     	
	 DTAUSNAP	= 1.0d3	
	 DTAUSNAP_PARAM	 	= 1.0d3
	 MP			= 1.0d-3 			
	 Q			= 0.1d0				
	 ROUT		= 20.0d0  
	 ROUT_UNIT	= 1	
	 A0			= 1.0d0     		
	 F			= 1.0d0				
	 INITIAL_SIGMA_TYPE	= 5	
	 INITIAL_SIGMA_FILE = 'restart.dat'
	 SIGMADOT_TYPE		= 0 		
	 MDOTIN		= 1.0d0    			
	 M_DISC0	= 1.0d1		
	 ISPLANET 	= 1
	 TINJECTPLANET		= 0.0d0     
	 TMIGRATION	= 1.0d-2		
	 ALPHA		= 0.1d0
	! set switches...to do
	
	
	! note for input.nml:
	!   1. put the name of the namelist at the beginning, e.g.:					&inputvariables
	!	2. write the variables, one per line, with a comma at the end, e.g.: 	rin=0.1,
	! 	3. end the values with / on the same line of the last variable, e.g.:	mdotin=150./		
	!	4. press enter to end the file
	! see for an example: sample_input.nml
	
	! reads from input.nml, if it exists
	open(unit=u_mainoutput,file=trim(inputfile),delim='apostrophe',status='old',action='read',iostat=ierror)
	if (ierror==0) then
		write(*,'(a,$)') ' - reading from '//trim(inputfile)//' ...'
		read(u_mainoutput,nml=inputvariables)
		write(*,'(a)') 'done'
	else
		write(*,'(a)') ' - input.nml doesn''t exist, using defaults'
	endif
	close(u_mainoutput)
	
	if (sim_mode=='') sim_mode='manual'
	
	write(*,nml=inputvariables)
	! writes the used input variables
	!open(unit=u_mainoutput,file=trim(homedir)//'/used_input.nml',delim='apostrophe')
	!write(u_mainoutput,nml=inputvariables)
	!close(u_mainoutput)
	
	if (INITIAL_SIGMA_TYPE	== 6) then
		! if I want to restart the simulation
		write (*,*) 'Restart from : ', 	 INITIAL_SIGMA_FILE
	endif
	
	! checks dtausnap and dtausnap_param
	if ((DTAUSNAP<DTAUSNAP_PARAM).or.((DTAUSNAP/DTAUSNAP_PARAM)/=int(DTAUSNAP/DTAUSNAP_PARAM))) then
		write (*,*) 'ERROR: DTAUSNAP MUST BE A MULTIPLE OF DTAUSNAP_PARAM'
		stop
	endif

end subroutine read_variables
!************************************************


!------------------------------------------------
!	this subroutine sets the defaults values
! 	of variables reading from a file if exists,
! 	otherwise sets default values.								
!------------------------------------------------
subroutine dimensional_variables()

	use variables,	only :	mp,GMp,Ms,q,rg,rISCO,rin,rout,rout_unit,mdotin,a0,m_disc0,gwprefactor,tot_time,dtausnap,tmigration,&
						 &	t0,T_c_prefactor,alpha,Rovermu,T_eff_prefactor,acrit,isncount_latest,dtausnap_param,rg10,&
						 &	initial_dtausnap_param,L_Edd
	
	use constants, 	only :	c_light,c_light2,M_sun,parsec,year,Grav,R_gas,mu,sig_sb,twopi,fourpi,eta,sigma_T,m_p
							
	implicit none
	
	mp		=	mp		*	M_sun
	GMp		=	Grav	*	mp
	Ms		=	q		*	mp
	
	! note : rg with 2* is the Schwarschild radius. I made this choice for a more direct comparison with Chang & al. 2010.
	rg		=	2.0d0	* GMp / c_light**2.
	rg10 	= 	10.0d0	* rg
	rISCO	=	3.0d0	* rg
	acrit 	=	500		* rg
	isncount_latest = 500			! (decreasing) index for the snapshots prior to merger
	
	rin		=	rISCO		
	if (rout_unit==1) then
		! a0 and rout are given in parsec
		a0		=	a0		* 	parsec
		rout	=	rout	*	parsec
	elseif (rout_unit==2) then
		! a0 and rout are given in rg
		! for Chang et al. 2009 simulations
		a0		=	a0		* 	rg
		rout	=	rout	*	rg
	elseif (rout_unit==3) then
		! a0 is given in parsec
		! rout is given in a0
		! for Lodato et al. 2009 simulations
		a0		=	a0		* 	parsec
		rout	=	rout	*	a0
	endif
	
	mdotin	=	mdotin	*	M_sun / year
	m_disc0 =	m_disc0	*	M_sun

	gwprefactor	= 64./5. * ((rg/2)**3.) * c_light * (1+q) * q
	
	t0						=	t0			*	year
	tot_time				=	tot_time	*	year
	dtausnap				=	dtausnap	*	year
	dtausnap_param 			= dtausnap_param * year
	initial_dtausnap_param 	= dtausnap_param
	tmigration 				= tmigration	*	year
	if (tmigration == 0.0d0 ) tmigration = 300.0d0*twopi/sqrt(GMp/(a0**3.))		! migration starts after 300 orbits

	Rovermu= R_gas / mu
	T_c_prefactor	=	9./16. * alpha * Rovermu / sig_sb
	T_eff_prefactor =	9./8.  / sig_sb

	!Eddington Luminosity
	L_Edd = fourpi * Grav * c_light * m_p * Mp / sigma_T

	print *, 'check dimensional variables'
	!write(*,'(g13.6,g13.6,g13.6,g13.6,g13.6,g13.6)') mp,GMp,Ms,rg,rISCO,rin
	!write(*,'(g13.6,g13.6,g13.6,g13.6)') a0,rout,mdotin,m_disc0,gwprefactor
	!write(*,'(g13.6,g13.6)') gwprefactor, tmigration/year
	write (*,'(a,ES14.6)') 'Orbit period [year]: ',tmigration/300.0d0/year
end subroutine dimensional_variables
!************************************************


!------------------------------------------------
!	prepares files
!------------------------------------------------
subroutine initialize_files()

	use variables,	only:	massfilename,planetfilename,splashfilenames,homedir,today,sim_mode,sigma_negfile,inputfile
	
	use units, 		only: 	u_snapshots,u_massfile,&
							& u_splashfilenames,u_sigma_neg

	
	implicit none
	
	! sets the current day/month/year in the vector today(3)
	call idate(today)
	
	write(homedir,'(a,i4,i2.2,i2.2,a,a)') 'sim_',today(3),today(2),today(1),'_',trim(sim_mode)
	call system('mkdir '//trim(homedir))
	call system('mkdir '//trim(homedir)//'/snap')
	call system('mkdir '//trim(homedir)//'/snap_latest')
	call system('mkdir '//trim(homedir)//'/code')
	
	
	!if already exists, deletes mass.dat
	write(massfilename,'(a,a)') trim(homedir),'/mass.dat'
	open(unit=u_massfile, file=massfilename,status='unknown')
	close(u_massfile,status='delete')


	
	! if already exists, deletes splash.filenames
	write(splashfilenames,'(a,a)') trim(homedir),'/splash.filenames'
	open(unit=u_splashfilenames, file=splashfilenames,status='unknown')
	close(u_splashfilenames,status='delete')
	

	! save the code and the parameters used for the current simulation
	call system('cp '//trim(inputfile)//' '//trim(homedir)//'/code/usedinput.nml')
	call system('cp binary.F90 '//trim(homedir)//'/code/binary.F90')
	call system('cp routines.F90 '//trim(homedir)//'/code/routines.F90')
	call system('cp variables.F90 '//trim(homedir)//'/code/variables.F90')
	call system('cp constants.F90 '//trim(homedir)//'/code/constants.F90')
	!call system('cp analytic.F90 '//trim(homedir)//'/code/analytic.F90')
	call system('cp makefile '//trim(homedir)//'/code/makefile')
	call system('cp splash.limits '//trim(homedir)//'/splash.limits')
	call system('cp splash.defaults '//trim(homedir)//'/splash.defaults')
	call system('cp splash.columns '//trim(homedir)//'/splash.columns')
	
	! if already exists, deletes file sigma_neg.dat
	write(sigma_negfile,'(a,a)') trim(homedir),'/sigma_neg.dat'
	open(unit=u_sigma_neg, file=sigma_negfile,status='unknown')
	close(u_sigma_neg,status='delete')	
	! reopens sigma_neg.dat to prepare it for writing
	open(unit=u_sigma_neg, file=sigma_negfile,status='unknown')
	
end subroutine initialize_files
!************************************************



!------------------------------------------------
!	allocates all the arrays
!------------------------------------------------
subroutine allocate_arrays()
	
	use variables, 	only:	inr,inr2,x,xpc,xrg,x12,x2,dx,dx2,xdx,D,sigma,oldsigma,mdot,&
						&	lambda,advvel,nu,H,exactsigma,gg,area,omega,omega2,inside_rhill,T_c,T_eff
	
	implicit none
	
	!write(*,'(a,$)') ' - allocating arrays...'
	inr2=2*inr	
	
	allocate(x(inr2),xpc(inr),xrg(inr),x12(inr2),x2(inr2))
	allocate(dx(inr2),dx2(inr2),xdx(inr),area(inr),D(inr))
	allocate(sigma(inr),oldsigma(inr),exactsigma(inr),mdot(inr))
	allocate(lambda(inr),advvel(inr),nu(inr),H(inr),gg(inr),omega(inr),omega2(inr))
	allocate(T_c(inr),T_eff(inr))
	allocate(inside_rhill(inr))
	!write(*,*) 'done'	
	
end subroutine allocate_arrays
!************************************************

!------------------------------------------------
!	allocates all the arrays
!------------------------------------------------
subroutine deallocate_arrays()
	
	use variables, 	only:	x,xpc,xrg,x12,x2,dx,dx2,xdx,D,sigma,oldsigma,mdot,nu,omega,omega2,inside_rhill,&
						&	lambda,advvel,nu,H,exactsigma,gg,area,T_c,T_eff
	
	implicit none
	
	
	deallocate(x,xpc,xrg,x12,x2,dx,dx2,xdx,area,D,sigma,oldsigma,exactsigma,mdot,lambda,advvel,H,gg)
	deallocate(nu,omega,omega2,T_c,T_eff)
	deallocate(inside_rhill)
end subroutine deallocate_arrays
!************************************************



!------------------------------------------------
!	sets the spatial grid: x, dx
!------------------------------------------------
subroutine set_grid()
 	
 	use variables, 	only: 	i,inr,inr2,x,xpc,xrg,x12,x2,dx,dx2,xdx,rin,rout,dx_lin,gridtype,mass,area,omega,omega2,GMp,rg

 	use constants,	only:	twopi,parsec
 	
	implicit none
	
	integer :: findxcell
	
	write(*,'(a,$)') ' - setting radial grid...'
	
	!Sets radial limits. Warning: rin=0. gives error because D(i) is proportional to x(i)**-1
	
	select case (gridtype)
		case (0)
			! LINEAR GRID
			!Warning: x(1) is not defined
			x(2)=rin
			dx_lin=(rout-rin)/real(inr2-2)
			do i=3,inr2
				x(i)=x(2)+(i-2)*dx_lin	
			enddo
			!print *, x(i),dx_lin 	! check OK
			! x(2) and x(inr2) are the position of ghost cells
		case (1)
			! LOG GRID
			! Warning: x(1) is not defined
			! in a log grid: - delta_x/x is constant
			! 				 - x(i+1)/x(i) is constant
			!				 - in the log space, the points are equally spaced
			! 10-based logarithm
			x(2)=log10(rin)
			dx_lin=(log10(rout)-log10(rin))/(inr2-2)
			do i=3,inr2
				x(i)=10.**(x(2)+dx_lin*(i-2))
			enddo
			x(2)=10.**x(2)			
	end select

	!Sets dx of the spacing grid.
	!Warning: dx(1) is not defined
	do i=3,inr2-1
		dx(i)=x(i+1)-x(i-1)
		dx2(i)=dx(i)**2.
		x12(i)=sqrt(x(i))
		x2(i)=x(i)**2.
	enddo
	x(1)=0.0d0
	
	x12(1)=0.0d0
	x12(2)=sqrt(x(2))
	x12(inr2)=sqrt(x(inr2))
	
	x2(1)=0.0d0
	x2(2)=0.0d0
	x2(inr2)=x(inr2)**2.
	
	dx(2)=dx(3)
	dx(inr2)=dx(inr2-1)
	
	dx2(2)=dx(3)**2.
	dx2(inr2)=dx(inr2-1)**2.	
	

	! area(i) contains annuli area
	do i=1,inr
		xdx(i)=x(2*i)*dx(2*i)
		area(i)=twopi*xdx(i)
		omega2(i)=GMp/(x(2*i)**3.)
		omega(i)=sqrt(omega2(i))
		xpc(i)=x(2*i)/parsec
		xrg(i)=x(2*i)/rg
		!print *, x(2*i),dx(2*i),area(i)
		!print *, i,xdx(i),area(i),omega2(i),omega(i)
	enddo

	write(*,*) 'done'
	
	
	! major check of the grid initialization
	!do i=1,inr2
		!print *, i,x(i),(x(i+1)-x(i))/x(i),dx(i)
		!mass=mass+area(i/2)
		!print *, i,area(i),mass
	!enddo
	
end subroutine set_grid
!************************************************



!------------------------------------------------
!	sets the time related variables
!------------------------------------------------
subroutine set_time()

	use variables, 	only: t,t0,dtau,tot_time,dtausnap,isnapnum,timenextsnap,initial_sigma_type,timenextsnap_param,dtausnap_param
	
	implicit none
	
	!write(*,'(a,$)') ' - setting time-related variables and snapshots...'
	t=t0
	
	! sets the first dtau in order to correctly initialize the sigma
	if (initial_sigma_type/=6) dtau=1.d-3
	
	! sets the time interval between two snapshots and the total number of snapshots
	!!dtausnap=0.1
	isnapnum=nint((tot_time-t0)/dtausnap)+1
	
	! sets the time of the next (the first) snapshot after the initial one
	timenextsnap=t0+dtausnap
	timenextsnap_param = t0+dtausnap_param
	
	!write(*,*) 'done'
	!write(*,'(a,i7,a)') ' - this simulation will compute:',isnapnum-1,' snapshots'
end subroutine set_time
!************************************************


!------------------------------------------------
!	sets the initial distribution of sigma
!------------------------------------------------
subroutine set_initialsigma()
	use variables, 	only: 	i,inr,inr2,oldsigma,sigma,sigmadot_rin,sigmadot_rout,&
						& 	rin,rout,x,nu,tau0,t0,xin,initial_sigma_type,mdotin,nu,nu0,a,ia,&
						&	a0,sigma0,m_disc0,INITIAL_SIGMA_FILE,rout,nu,T_eff,T_c,H,mass,sigmadot_type,prec_out,area
							
	use formats, 	only:	f_snapshots
	! use constants, only: pi_d !it's commented because pi_d is already defined in nrtype; pi_d is the double precision version.
	use constants, 	only:	threepi,pi,parsec,M_Sun
	
	! Numerical Recipes modules
	!use nrtype
	!use nrutil
	!use nr, only : bessik
	!----------------------------------
	
	implicit none

	! functions declaration
	!double precision :: spreading_ring
	double precision :: trash
	integer :: i0
	integer :: findxcell
	double precision :: mass_normalization
	double precision :: r0
	
	nu0=1.5d20
	
	!Sets the sigma equal to zero everywhere.
	write(*,'(a,$)') ' - setting the initial sigma distribution...'
	do i=1,inr
		oldsigma(i)=0.0d0
		sigma(i)=0.0d0
		!nu(i)=nu0*(x(2*i)/(0.1*parsec))**(1.5d0)
	enddo		
	!print *, 'nu(inr)',inr, nu(inr)	! check OK
	sigmadot_rout=0.0d0
	
	select case(initial_sigma_type)
		case(0)
			! no initial sigma
	
			if (sigmadot_type == 2) then
				! correction to avoid errors in self-consistent calculation of nu, T_c, H
				sigmadot_rout=1.0d-4

				sigma(inr-1)=1.0d-8
				sigma(inr)=1.0d-8
				oldsigma(inr-1)=sigma(inr-1)
				oldsigma(inr)=sigma(inr)
				nu(inr-1)=1.0d0
				nu(inr)=1.0d0
			endif
		case(1)
			! Dirac-delta function centered in r0
			r0=10.*a0
			i0=findxcell(x,inr2,r0)/2
			sigma(i0)=100.
			
		case(2)
			! linear+exponential decay
			do i=1,inr
				sigma(i)=(rout-x(2*i))*exp(-x(2*i)/(rout/10.))
			enddo
			
		case(3)
			! exponential
			mass_normalization=pi*rout**2*(exp(-rin/rout)*(rin/rout+1)-2.*exp(-1.))
			do i=1,inr
				sigma(i)=10.*exp(-(x(2*i)/(rout/10.)))/mass_normalization
			enddo
			
		case(4)
			! spreading ring
			! needs analytic.o to be included
			do i=1,inr
				!sigma(i)=spreading_ring(xin,x(2*i),t0)	
!				print *,sigma(i),xin,x(2*i),t0
			enddo
			
		case(5)
			! stationary disc
			!ia=findxcell(x,inr2,a)/2
			!mass_normalization=mdotin/(threepi*nu(ia))*(1.-sqrt(rin/x(2*ia)))*(a**2.)
			mass = 0
			!do i=1,findxcell(x,inr2,rout*0.99)/2
			do i = 1, inr-1
				!sigma(i)=a**(-2.)*(nu(ia)/nu(i))*((1.-sqrt(rin/x(2*i)))/(1.-sqrt(rin/x(2*ia))))
				!sigma(i)=mdotin/(threepi*nu(i))*(1.-sqrt((rin/x(2*i))))
				sigma(i)=x(2*i)**(-3/5.)*(1.-sqrt((rin/x(2*i))))
				mass = mass + area(i)*sigma(i)
				oldsigma(i)=sigma(i)
				!print *,x(2*i),sigma(i),nu(i)
			enddo	
			
			! re-normalize the sigma distribution in order to obtain a disc of mass M_DISC0
			do i = 1, inr
				sigma(i)=sigma(i)/mass*m_disc0
			enddo
			
			if (sigmadot_type == 2) then
				! correction to avoid errors in self-consistent calculation of nu, T_c, H
				sigmadot_rout=1.0d-4
				sigma(inr)=1.0d-8
				oldsigma(inr)=sigma(inr)
				nu(inr-1)=1.0d0
				nu(inr)=1.0d0
			endif
		case(6)
			! initialize sigma from a file
			! i.e. : it is used to restart simulation
			open(unit=30,file=INITIAL_SIGMA_FILE,status='old')
			do i=1,inr
				read(30,f_snapshots) trash,sigma(i),nu(i),T_eff(i),T_c(i),trash,trash
				!write(*,'(5(ES14.6,1x))'), x(2*i),sigma(i),nu(i),T_eff(i),T_c(i)
			enddo
			close(30)

			trash=0.0d0
			
		case(7)
			! initial condition Lodato et al. 2009
			sigma0=m_disc0 / (pi*((2.0*a0)**2.-(0.8*a0)**2.))
			!print*, m_disc0,2.0*a0,*a0
			do i = findxcell(x,inr2,0.8*a0)/2, findxcell(x,inr2,2.0*a0)/2
				sigma(i)=sigma0
				!print *, i , sigma(i)
			enddo
			
			if (sigmadot_type == 2) then
					! correction to avoid errors in self-consistent calculation of nu, T_c, H
					sigmadot_rout=1.0d-4

					sigma(inr-1)=1.0d-8
					sigma(inr)=1.0d-8
					oldsigma(inr-1)=sigma(inr-1)
					oldsigma(inr)=sigma(inr)
					nu(inr-1)=1.0d0
					nu(inr)=1.0d0
			endif
			
		case(8)
			! initial steady state of a b) region (gas pressure dominated)
			! Chang initial condition
			! initialize sigma from a file
			mass = 0.0d0
			open(unit=30,file=INITIAL_SIGMA_FILE,status='old')
			do i=1,inr
				read(30,f_snapshots) trash,sigma(i),nu(i),T_eff(i),T_c(i),trash,trash
				!write(*,'(5(ES14.6,1x))'), x(2*i),sigma(i),nu(i),T_eff(i),T_c(i)
				mass = mass + area(i)*sigma(i)
			enddo
			close(30)
			trash=0.0d0
			
			
			! re-normalize the sigma distribution in order to obtain a disc of mass M_DISC0
			do i = 1, inr
				sigma(i)=sigma(i)/mass*m_disc0
			enddo
				
			if (sigmadot_type == 2) then
					! correction to avoid errors in self-consistent calculation of nu, T_c, H
					sigmadot_rout=1.0d-4

					sigma(inr)=1.0d-8
					oldsigma(inr)=sigma(inr)
					nu(inr-1)=1.0d0
					nu(inr)=1.0d0
			endif
		
		case(9)
			! same initial steady state of case(8), i.e. of Chang et al.
			! here I clear a region around the secondary's initial radius 0.5 < r/a0 < 2.0
			
						mass = 0.0d0
			open(unit=30,file=INITIAL_SIGMA_FILE,status='old')
			do i=1,inr
				read(30,f_snapshots) trash,sigma(i),nu(i),T_eff(i),T_c(i),trash,trash
				!write(*,'(5(ES14.6,1x))'), x(2*i),sigma(i),nu(i),T_eff(i),T_c(i)
				mass = mass + area(i)*sigma(i)
			enddo
			close(30)
			trash=0.0d0
			
			
			! re-normalize the sigma distribution in order to obtain a disc of mass M_DISC0
			do i = 1, inr
				sigma(i)=sigma(i)/mass*m_disc0
			enddo
				
			if (sigmadot_type == 2) then
					! correction to avoid errors in self-consistent calculation of nu, T_c, H
					sigmadot_rout=1.0d-4

					sigma(inr)=1.0d-8
					oldsigma(inr)=sigma(inr)
					nu(inr-1)=1.0d0
					nu(inr)=1.0d0
			endif
			
			! now I clear the gap
			do i = findxcell(x,inr2,0.5*a0)/2 , findxcell(x,inr2,2.0*a0)/2
				sigma(i) = 0.0d0
			enddo
		end select
	
	
	write(*,*) 'done'
	call calculate_mass()
	write (*,'(a,ES14.6)') 'Initial disc mass [M_Sun]: ',mass/M_sun
	
	!Sets the initial boundary condition
	!Zero-torque boundary conditions at r=rin and r=rout
	sigmadot_rin=0.
	sigma(1)=sigmadot_rin
	sigma(inr)=sigmadot_rout

	! check
!	do i=1,inr
!		write(*,'(5(ES14.6,1x))'), x(2*i),sigma(i),nu(i),T_eff(i),T_c(i)
!	enddo
!	stop
	
end subroutine set_initialsigma
!************************************************


!------------------------------------------------
!	sets dtau according to 
!	courant stability (necessary condition)	 
!------------------------------------------------
subroutine set_dtau()
	
	use variables,	only:	i,inr2,dtau,nu,dx,dx2,xdx,area,advvel,oldsigma,dtaumin,sigmacell,CFL,a,t,x,adot,iax,acc_factor
	use constants,	only:	year
	implicit none

	
	dtaumin=1.d10

	do i=2,inr2-1
		sigmacell=i/2
		
		dtaumin=min(dtaumin,dx2(i)/nu(sigmacell),dx2(i+1)/nu(sigmacell),&
						& xdx(sigmacell)/abs(advvel(sigmacell)),xdx(sigmacell)/abs(advvel(sigmacell)))

		!print *,dtaumin,dx2(i)/nu(sigmacell),dx2(i+1)/nu(sigmacell),&
		!				& area(i)/abs(advvel(sigmacell)),area(i)/abs(advvel(sigmacell))

	enddo
	
	! Lodato uses accuracy acc=0.1 (line 632 main.f). I call it acc_factor
	dtaumin=min(dtaumin,acc_factor*dx(iax)/abs(adot))
!	if (dtaumin >= 0.1*dx(iax)/abs(adot)) then
!		!print *,'**** ECCO ****', dtaumin/year,a,iax,dx(iax)
!		dtaumin=min(dtaumin,0.1*dx(iax)/abs(adot))
!	endif
	dtau=dtaumin*CFL

end subroutine set_dtau
!************************************************



!------------------------------------------------
!	sigma dot   								
!------------------------------------------------
subroutine set_sigmadot()

	use variables, only: sigma,oldsigma,idrip,inr,inr2,mdotin,dtau,dx,nu,x,sigmadotin,sigmadot_type
	use constants, only: sixpi
	
	implicit none
	
	select case(sigmadot_type)
		case(0)
			oldsigma(inr)=0.
			return
!		case(1)
!			! matter drips in at a certain radius at every timestep
!			oldsigma(idrip)=oldsigma(idrip)+sigmadotin*dtau
		case(2)
			! performs constant inflow of matter at inr
			oldsigma(inr)=(sqrt(x(inr2-2))*nu(inr-1)*oldsigma(inr-1)+(mdotin*(dx(inr2-2)+dx(inr2-1)))/&
					  &(sixpi*sqrt(x(inr2-1))))/(sqrt(x(inr2))*nu(inr))
			sigma(inr)=oldsigma(inr)		 
			!print *, 'added',nu(inr),oldsigma(inr)		! check 
			return
	end select
	
	! Major check of sigma initialization
	!do i=1,inr
	!		print *,i,x(2*i),sigma(i)
	!enddo

end subroutine set_sigmadot
!************************************************


!------------------------------------------------
!	computes the torque and the advective velocity
!	uses torque clipping
!------------------------------------------------
subroutine set_torque_advvel_clip

	use variables, 	only: 	i,inr,inr2,lambda,advvel,H,x,x12,x2,dx,f,q,&
						&	r_hill,a,ia,iax,ILR,OLR,iILR,iOLR,delta,deltamin,iIlRx,iOLRx,iax,&
						&	Rey,factor,omega,omega2,sigma,dtau,inside_rhill,r_hill,clipped_torque,rg
						
	
	use constants, only : twothirds,two_at_twothirds,onehalf_at_twothirds,parsec

	implicit none
	
	integer :: findxcell

	! SMOOTHING PRESCRIPTION
	! From: Lodato et al. 2009, Armitage & Natarajan 2002, Goldreich & Tremaine 1980, Armitage "Astrophysics of planet formation", page 236
	! Statement: - smoothing where |r-a| < max(H,r_hill), where H is the scale height
	! Besides, we impose an exponential decay of lambda at radius outside the first outer Lindblad resonance.
	

	! finds the Lindblad resonances
	OLR=two_at_twothirds*a
	ILR=onehalf_at_twothirds*a
	!iILRx=findxcell(x,inr2,ILR)
	!iOLRx=findxcell(x,inr2,OLR)
	!iILR=iILRx/2
	!iOLR=iOLRx/2
	
	! gap size
	deltamin=max(H(ia),r_hill)
	! maybe H(ia) is the semi-thickness of the disc : chiedere a Lodato
	!write(*,'(3(ES14.6))') ilr/rg,a/rg,olr/rg
	do i=2,inr-1

		
			factor=1.0d0
			delta=max(deltamin,abs(x(2*i)-a),2.*dx(2*i))

			if (x(2*i) <= a) then 
				if (x(2*i) <= ILR) factor = exp(-((x(2*i)-ILR)/(20.*H(ia)))**2.)
				factor = -factor * (x(2*i)/delta)**4.
			else
				if (x(2*i) >= OLR) factor = exp(-((x(2*i)-OLR)/(20.*H(ia)))**2.)
				factor =  factor * (a/delta)**4.
			endif
!			if (abs(x(2*i)-a)<= r_hill) then
!				inside_rhill(i)=1.
!			else	
!				inside_rhill(i)=0.
!			endif
			!write(*,'(2(i3.3),7(ES14.6))') i,ia,H(ia-30),r_hill,deltamin,x(2*i)-a,2.*dx(2*i),delta,factor
		lambda(i)=factor*0.5*f*(q**2.)*omega2(i)*x2(2*i)
			
			! torque clipping
			! if the local value of the torque is greater than a certain amount, then clip the torque
			! I use the prescription in Alexander & Armitage 2009, page 3
			clipped_torque=0.1*omega2(i)*x(2*i)*H(ia)
			if (abs(lambda(i))>=clipped_torque) then
				!write(*,'(2(ES14.6))') lambda(i), clipped_torque
				lambda(i)=sign(clipped_torque,lambda(i))
			endif
			advvel(i)=2.*lambda(i)/omega(i)
	enddo
	
	
	! mass conservation condition: the torque has to be null in the first near cells
	lambda(ia)=0.
	lambda(ia+1)=0.
	advvel(ia)=0.
	advvel(ia+1)=0.
	
	! checks
	!print *, lambda(iILR+1),lambda(iOLR-2),OLR,x(2*iOLR-1),x(2*iOLR),x(2*iOLR+1)
	!do i = 1,inr	
	!	print *,x(2*i),lambda(i),advvel(i) !check OK	
	!enddo
	
end subroutine set_torque_advvel_clip
!************************************************

!------------------------------------------------
!	computes the torque 
!	and the advective velocity	 
!------------------------------------------------
subroutine set_torque_advvel_hill()

	use variables, 	only: 	i,inr,inr2,lambda,advvel,H,x,x12,x2,dx,f,q,&
						&	r_hill,a,ia,iax,ILR,OLR,iILR,iOLR,delta,deltamin,iIlRx,iOLRx,iax,&
						&	Rey,factor,omega,omega2,sigma,dtau,inside_rhill,r_hill,rg,d_smoothing_in,d_smoothing_out
						
	
	use constants, only : twothirds,two_at_twothirds,onehalf_at_twothirds,parsec

	implicit none
	
	integer :: findxcell

	! SMOOTHING PRESCRIPTION
	! From: Lodato et al. 2009, Armitage & Natarajan 2002, Goldreich & Tremaine 1980, Armitage "Astrophysics of planet formation", page 236
	! Statement: - smoothing where |r-a| < max(H,r_hill), where H is the scale height
	! Besides, we impose an exponential decay of lambda at radius outside the first outer Lindblad resonance.
	

	! finds the Lindblad resonances
	OLR=two_at_twothirds*a
	ILR=onehalf_at_twothirds*a
	!iILRx=findxcell(x,inr2,ILR)
	!iOLRx=findxcell(x,inr2,OLR)
	!iILR=iILRx/2
	!iOLR=iOLRx/2
	
	! gap size
	deltamin=max(2.*H(ia),r_hill)
	! maybe H(ia) is the semi-thickness of the disc : chiedere a Lodato
	
	do i=2,inr-1
			! added on 21 May 2013
			d_smoothing_in=370.0d0*H(i)
			d_smoothing_out=75.0d0*H(i)
			if (H(i)==0.) then 
				d_smoothing_in = dx(2*i)
				d_smoothing_out = dx(2*i)
			endif	
			factor=1.0d0
			delta=max(deltamin,abs(x(2*i)-a),2.*dx(2*i))
			
			if (x(2*i) <= a) then 
				!if (x(2*i) <= ILR) factor = exp(-((x(2*i)-ILR)/(H(i)))**2.)
				if (x(2*i) <= ILR) factor = exp(-((x(2*i)-ILR)/d_smoothing_in)**2.)
				factor = -factor * (x(2*i)/delta)**4.
			else
				!if (x(2*i) >= OLR) factor = exp(-((x(2*i)-OLR)/(H(i)))**2.)
				if (x(2*i) >= OLR) factor = exp(-((x(2*i)-OLR)/d_smoothing_out)**2.)				
				factor =  factor * (a/delta)**4.
			endif
!			if (abs(x(2*i)-a)<= r_hill) then
!				inside_rhill(i)=1.
!			else	
!				inside_rhill(i)=0.
!			endif
			!write(*,'(2(i3.3),7(ES14.6))') i,ia,H(ia-30),r_hill,deltamin,x(2*i)-a,2.*dx(2*i),delta,factor
		lambda(i)=factor*0.5*f*(q**2.)*omega2(i)*x2(2*i)
		advvel(i)=2.*lambda(i)/omega(i)
		
	enddo
	
	
	! mass conservation condition: the torque has to be null in the first near cells
	lambda(ia)=0.
	lambda(ia+1)=0.
	advvel(ia)=0.
	advvel(ia+1)=0.
	
	! checks
	!print *, lambda(iILR+1),lambda(iOLR-2),OLR,x(2*iOLR-1),x(2*iOLR),x(2*iOLR+1)
	!do i = 1,inr	
	!	print *,x(2*i),lambda(i),advvel(i) !check OK	
	!enddo
	
end subroutine set_torque_advvel_hill
!************************************************

!------------------------------------------------
!	computes the torque 
!	and the advective velocity	 
!------------------------------------------------
subroutine set_torque_advvel()

	use variables, 	only: 	i,inr,inr2,lambda,advvel,H,x,x12,x2,dx,f,q,&
						&	r_hill,a,ia,iax,ILR,OLR,iILR,iOLR,delta,deltamin,iIlRx,iOLRx,iax,&
						&	Rey,factor,omega,omega2,sigma,dtau,inside_rhill
						
	
	use constants, only : twothirds,two_at_twothirds,onehalf_at_twothirds

	implicit none
	
	integer :: findxcell

	! SMOOTHING PRESCRIPTION
	! From: Lodato et al. 2009, Armitage & Natarajan 2002, Goldreich & Tremaine 1980, Armitage "Astrophysics of planet formation", page 236
	! Statement: - smoothing where |r-a| < max(H,r_hill), where H is the scale height
	! Besides, we impose an exponential decay of lambda at radius outside the first outer Lindblad resonance.
	

	! finds the Lindblad resonances
	!OLR=two_at_twothirds*a
	!ILR=onehalf_at_twothirds*a
	!iILRx=findxcell(x,inr2,ILR)
	!iOLRx=findxcell(x,inr2,OLR)
	!iILR=iILRx/2
	!iOLR=iOLRx/2
	
	! set initial oldsigma
	do i=1,inr
		! smoothing of the torque outside the Lindblad resonances.
		! the resonances are not dirac delta functions; they're instead gaussian-like of width H(i)
		! result from : Artymowicz, P., ApJ, v.419 419, 155 (1993).
		
		! Marco Tazzari Torque
			factor=1.0d0
			delta=max(abs(x(2*i)-a),H(i))
			
			if (abs(x(2*i)-a) < H(i)) then
				factor = (x(2*i)-a)/H(i)		
			else
				factor = sign(1.0d0,x(2*i)-a)	
			endif

			if (x(2*i) <= a) then 
				!if (x(2*i) <= ILR) factor = exp(-((x(2*i)-ILR)/H(i))**2.)
				lambda(i)= (x(2*i)/delta)**4. * factor
			else
				!if (x(2*i) > OLR) factor = exp(-((x(2*i)-OLR)/H(i))**2.)
				lambda(i)=  (a/delta)**4. * factor
			endif
!			if (abs(x(2*i)-a)<= r_hill) then
!				inside_rhill(i)=1.
!			else	
!				inside_rhill(i)=0.
!			endif



		lambda(i)=lambda(i)*0.5*f*(q**2.)*omega2(i)*x2(2*i)
		advvel(i)=2.*lambda(i)/omega(i)
	enddo
	
	
	! mass conservation condition: the torque has to be null in the first near cells
	lambda(ia)=0.
	lambda(ia+1)=0.
	advvel(ia)=0.
	advvel(ia+1)=0.
	
	! checks
	!print *, lambda(iILR+1),lambda(iOLR-2),OLR,x(2*iOLR-1),x(2*iOLR),x(2*iOLR+1)
	!do i = 1,inr	
	!	print *,x(2*i),lambda(i),advvel(i) !check OK	
	!enddo
	
end subroutine set_torque_advvel
!************************************************


!------------------------------------------------
!	Disc thickness with Chang et al. 2010 prescription
!	computes the disc thickness taking into account the gas pressure only and fixed opacity.
!------------------------------------------------
subroutine set_thickness_chang()

	use variables, only : H,i,inr,x,T_c,alpha,T_c_prefactor,omega,omega2,oldsigma,Rovermu,nu,aconst,bconst,cconst,thin_corr,tau,&
						& a,i_t,accT,t,T_guess,k_ff,errT,T_eff,T_eff_prefactor,T_new
	
	use constants, only : onethird,twothirds,a_rad,k_T,k_ff_0,oneseventh,twothirteenths,parsec
	
	implicit none
	
	
	do i=2,inr-1
		! central temperature
		T_c(i)	= (T_c_prefactor * oldsigma(i)**2. * omega(i) * k_T)**onethird
		
		! local viscosity
		nu(i) 	= alpha * Rovermu * T_c(i) / omega(i)
		
		! local thickness
		H(i)	= sqrt(T_c(i) * Rovermu / omega2(i))
		
		! for a fixed thickness
		!H(i) = 0.05 * x(2*i)
		
		T_eff(i)= (T_eff_prefactor * nu(i) * oldsigma(i) * omega2(i) )**0.25

		! fixed H/R
		!H(i)=x(2*i)*0.05d0	
		!nu(i)=1.5d20*(x(2*i)/(0.1*parsec))**(1.5d0)
		!write(*,'(i3,g13.6,g13.6)') i, x(2*i),H(i)	! check OK
		!print *, H(i),T_c(i),nu(i)
	enddo
	!print *, oldsigma(inr-1),oldsigma(inr), nu(inr-1), nu(inr)
	return
end subroutine set_thickness_chang
!************************************************




!------------------------------------------------
!	computes the disc thickness taking into account the gas pressure only.
!------------------------------------------------
subroutine set_thickness_pgas()

	use variables, only : H,i,inr,x,T_c,alpha,T_c_prefactor,omega,omega2,oldsigma,Rovermu,nu,aconst,bconst,cconst,thin_corr,tau,&
						& a,i_t,accT,t,T_guess,k_ff,errT,T_eff,T_eff_prefactor,T_new
	
	use constants, only : onethird,twothirds,a_rad,k_T,k_ff_0,oneseventh,twothirteenths,parsec
	
	implicit none
	
	
	do i=2,inr-1
		! central temperature
		T_c(i)	= (T_c_prefactor * oldsigma(i)**2. * omega(i) * k_T)**onethird
		
		! local viscosity
		nu(i) 	= alpha * Rovermu * T_c(i) / omega(i)
		
	if (oldsigma(i) > 1.0d-20) then
		! local thickness
!		bconst 	= twothirds * a_rad * T_c(i)**4. / omega2(i) / oldsigma(i) * thin_corr
		cconst 	= T_c(i) * Rovermu / omega2(i)
		H(i)	= sqrt(cconst)
		
		! optically thin regime correction
		tau		=0.5 * k_T * oldsigma(i)
		thin_corr= tau/(1.+tau)
		!print * , tau,thin_corr ! check 
		! I see that optically thin regime occurs inside the gap and near it; outside this region opt thick regime dominates
		
		k_ff=k_ff_0 * oldsigma(i) / H(i) / T_c(i)**3.5

		if (k_ff >= k_T) then
			! accuracy of the temperature iteration scheme
			accT= 1.0d-3
			i_t = 0
			aconst = T_c_prefactor * omega(i) * k_ff_0
			
			if (aconst==0.) then
				print *, 'Forced stop: aconst = 0.',' time:',t,' a: ',a
				stop
			endif
			T_guess = (aconst * oldsigma(i)**3. / sqrt(Rovermu/omega2(i)))**oneseventh
			T_c(i)=T_guess
			errT=accT+1.
			
			
			T_iteration : do while (errT > accT) 
!				bconst 	= twothirds * a_rad * T_c(i)**4. / omega2(i) / oldsigma(i) * thin_corr
				cconst 	= T_c(i) * Rovermu / omega2(i)
				H(i)	= sqrt(cconst)
				
				T_new 	= ( aconst * oldsigma(i)**3. / H(i) )**twothirteenths
				errT 	= abs(T_new-T_c(i))/T_c(i)
				
				T_c(i)= T_new
				
				
				i_t = i_t+1
				if (i_t > 100) then
					print *, 'Forced stop: i_t = ', i_t,' time:',t,' a: ',a
					stop
				endif
				
			enddo T_iteration
			
!			bconst 	= twothirds * a_rad * T_c(i)**4. / omega2(i) / oldsigma(i) * thin_corr
			cconst 	= T_c(i) * Rovermu / omega2(i)
			H(i)	= sqrt(cconst)
			
		endif
		T_eff(i)= (T_eff_prefactor * nu(i) * oldsigma(i) * omega2(i) )**0.25
	endif	
		! fixed H/R
		!H(i)=x(2*i)*0.05d0	
		!nu(i)=1.5d20*(x(2*i)/(0.1*parsec))**(1.5d0)
		!write(*,'(i3,g13.6,g13.6)') i, x(2*i),H(i)	! check OK
		!print *, H(i),T_c(i),nu(i)
	enddo
	!print *, oldsigma(inr-1),oldsigma(inr), nu(inr-1), nu(inr)
	return
end subroutine set_thickness_pgas
!************************************************



!------------------------------------------------
!	computes the disc thickness with Pringle numerical scheme
!	takes into account both gas pressure and radiation pressure
!------------------------------------------------
subroutine set_thickness()

	use variables, only : H,i,inr,x,T_c,alpha,T_c_prefactor,omega,omega2,oldsigma,Rovermu,nu,aconst,bconst,cconst,thin_corr,tau,&
						& a,i_t,accT,t,T_guess,k_ff,errT,T_eff,T_eff_prefactor,T_new
	
	use constants, only : onethird,twothirds,a_rad,k_T,k_ff_0,oneseventh,twothirteenths,parsec
	
	implicit none
	
	
	do i=2,inr-1
		! central temperature
		T_c(i)	= (T_c_prefactor * oldsigma(i)**2. * omega(i) * k_T)**onethird
		
		! local viscosity
		nu(i) 	= alpha * Rovermu * T_c(i) / omega(i)
		
		!H(i)=x(2*i)*0.05d0	
		!print *, i,T_c(i),nu(i),H(i)
	
		
	if (oldsigma(i) > 1.0d-20) then
		! local thickness
		bconst 	= twothirds * a_rad * T_c(i)**4. / omega2(i) / oldsigma(i) * thin_corr
		cconst 	= T_c(i) * Rovermu / omega2(i)
		H(i)	= 0.5 * (bconst + sqrt(bconst**2.+4.*cconst))
		
		! optically thin regime correction
		tau		=0.5 * k_T * oldsigma(i)
		thin_corr= tau/(1.+tau)
		!print * , tau,thin_corr ! check 
		! I see that optically thin regime occurs inside the gap and near it; outside this region opt thick regime dominates
		
		k_ff=k_ff_0 * oldsigma(i) / H(i) / T_c(i)**3.5

		if (k_ff >= k_T) then
			! accuracy of the temperature iteration scheme
			accT= 1.0d-3
			i_t = 0
			aconst = T_c_prefactor * omega(i) * k_ff_0
			
			if (aconst==0.) then
				print *, 'Forced stop: aconst = 0.',' time:',t,' a: ',a
				stop
			endif
			T_guess = (aconst * oldsigma(i)**3. / sqrt(Rovermu/omega2(i)))**oneseventh
			T_c(i)=T_guess
			errT=accT+1.
			
			
			T_iteration : do while (errT > accT) 
				bconst 	= twothirds * a_rad * T_c(i)**4. / omega2(i) / oldsigma(i) * thin_corr
				cconst 	= T_c(i) * Rovermu / omega2(i)
				H(i)	= 0.5 * (bconst + sqrt(bconst**2.+4.*cconst))
				
				T_new 	= ( aconst * oldsigma(i)**3. / H(i) )**twothirteenths
				errT 	= abs(T_new-T_c(i))/T_c(i)
				
				T_c(i)= T_new
				
				
				i_t = i_t+1
				if (i_t > 100) then
					print *, 'Forced stop: i_t = ', i_t,' time:',t,' a: ',a
					stop
				endif
				
			enddo T_iteration
			
			bconst 	= twothirds * a_rad * T_c(i)**4. / omega2(i) / oldsigma(i) * thin_corr
			cconst 	= T_c(i) * Rovermu / omega2(i)
			H(i)	= 0.5 * (bconst + sqrt(bconst**2.+4.*cconst))
			
		endif
		T_eff(i)= (T_eff_prefactor * nu(i) * oldsigma(i) * omega2(i) )**0.25
	endif	
		! fixed H/R
		!H(i)=x(2*i)*0.05d0	
		!nu(i)=1.5d20*(x(2*i)/(0.1*parsec))**(1.5d0)
		!write(*,'(i3,g13.6,g13.6)') i, x(2*i),H(i)	! check OK
		!print *, H(i),T_c(i),nu(i)
	enddo
	!print *, oldsigma(inr-1),oldsigma(inr), nu(inr-1), nu(inr)
	return
end subroutine set_thickness
!************************************************



!-----------------------------------------!
!	computes the position of the planet	  !
!-----------------------------------------!
subroutine set_planet()

	use variables,	only:	a,ia,iax,x,inr,inr2,i,r_hill,t,tot_time,old_a,dtau,q,&
						&	angular_momentum,adot,BB,Rey,AA,oldsigma,rin,tmigration,isplanet,adot1,adot2,adot,&
						&	GMp,Ms,gwprefactor,lambda,advvel,t,tot_time,timenextsnap,timenextsnap_param,&
						&	tot_sim_mass,mass,tot_acc_mass,tot_out_mass,CFL,initial_dtausnap_param,&
						&	dtausnap_param,dtausnap
						
	use constants,	only: 	fourpi,onethird,year
	
	implicit none
	
	integer :: findxcell
	
	call calculate_angular_momentum()

	
	! computes the derivative of a(t)
	adot1=-(2./Ms)*sqrt(old_a/GMp)*angular_momentum
	
	adot2=-gwprefactor/(old_a**3.)
	
	adot=adot1+adot2

	! computes the new position of the planet
	a=old_a+adot*dtau*max(0.,t-tmigration)/(t-tmigration)
	
	! sets the new Hill's radius
	r_hill = (q**onethird)*a  

	if (a >= rin) then
		isplanet = 1
		
		iax = findxcell(x,inr2,a)
		ia = iax/2
		!call set_torque_advvel_clip()	! routine WITH		Torque Clipping (Alexander & Armitage 2009)
		call set_torque_advvel_hill() 	! routine WITH 		Hill Radius check
		!call set_torque_advvel()		! routine WITHOUT	Hill Radius check
	else
		
		CFL=0.1d0
		!dtausnap_param=initial_dtausnap_param
		dtausnap_param=1.0d3*year
		dtausnap=1.0d3*year
		timenextsnap_param = t+dtausnap_param
		timenextsnap = timenextsnap_param
		! makes the last snapshot
			! calculates mass and disc-planet exchanged angular momentum 
			call calculate_mass()

			tot_sim_mass=mass+tot_acc_mass-tot_out_mass
			! writes the snapshots and the other dump files
			call takesnapshot_parameters()
			
		a 	 = 0.
		adot = 0.
		adot1= 0.
		adot2= 0.
		do i = 1,inr
			lambda(i)=0.
			advvel(i)=0.
		enddo
		isplanet = 3
		
		write(*,'(a,ES14.6,a)'), '*************************  M E R G E R   t =',t/year,'  *************************'
		
		! sets the next snapshot one year after the merger.
		timenextsnap_param = t+year
		timenextsnap=t+year
	endif
		
		

	
	!NBBB: increasing the hill's radius prevents the angular_momentum to become negative.
	! this is because if there is too much mass where lambda is negative(r<a) then 
	! the integral becomes negative.
	! The consequence is that the planet doesn't move inward.
	! it would be interesting to understand why.
	
	!print *, a,ia,x(2*ia),x(2*ia+1),r_hill ! check OK
	!print *,dtau,a,angular_momentum !check
end subroutine set_planet
!************************************************


!------------------------------------------------
!	this is a handle function to interpolate functions 
!------------------------------------------------
double precision function eval(arr,n,xx,ixx)

	use variables, 	only:	x,dx
	
	implicit none
	

	integer,			intent(in) 	:: n		! array size
	double precision, 	intent(in)	:: arr(n)		! array name
	double precision, 	intent(in) 	:: xx
	integer,			intent(in) 	:: ixx
	
	eval=arr(ixx/2)+(arr(ixx/2+1)-arr(ixx/2))/(dx((ixx/2)*2)+dx((ixx/2)*2+1))*(xx-x((ixx/2)*2))

	return
end function eval
!************************************************


!------------------------------------------------
!	computes the integral of angular momentum 
!	exchanged with the planet   
!------------------------------------------------
subroutine calculate_angular_momentum()

	use variables, only: i,inr,angular_momentum,lambda,sigma,dx,area

	implicit none
	 
	! initialize angular momentum
	angular_momentum=0.

	! rectangles
	do i = 1, inr-1
		angular_momentum = angular_momentum + lambda(i)*sigma(i)*area(i)
	enddo

end subroutine calculate_angular_momentum
!************************************************


!------------------------------------------------
!	computes the mass of the disc 
!	and the integral of angular momentum 
!	exchanged with the planet   
!------------------------------------------------
subroutine calculate_mass()

	use variables, 	only: 	i,x,dx,area,gridtype,mass,sigma,inr,mass_inner,ia
	implicit none

		
	! initialize mass
	mass=0.
	mass_inner=0.
	
	! rectangles
	do i= 1, ia
		mass=mass+area(i)*sigma(i)
	enddo
	mass_inner=mass
	do i= ia+1, inr-1
		mass=mass+area(i)*sigma(i)
	enddo
	
end subroutine calculate_mass
!************************************************

!------------------------------------------------
!	writes snapshot files,
!	calculates mass and writes mass file
!------------------------------------------------
subroutine takesnapshot()

	use variables,	only:	i,inr,x,x2,xpc,xrg,mass,mdot,&
							& snapfile,mdotfile, strcount,rin,rout,sigma,&
							& D,H,lambda,nu,advvel,isncount,homedir,&
							& splashfilenames,oldsigma,T_eff,T_c,&
							& prec_out,t,a,mass_inner,prad,pgas,prad_pgas,rovermu,tnu,rg,adot,adot1,adot2,inside_rhill
	
	use units, 		only: 	u_snapshots, u_splashfilenames,u_sigma_neg
	
	use formats, 	only: 	f_snapshots
	
	use constants,	only:	year,M_sun,onethird,a_rad
		
	implicit none
	
	! convert to string the counter of the snapshots
	write(strcount,'(i10.10)') isncount	
	snapfile=trim(homedir)//'/snap/snap'//trim(strcount)//'.dat'

	open(unit=u_snapshots, file=snapfile,status='unknown')
	write(u_snapshots,'(a)') 'time/yr   a/cm   M_in/M_sun   M_disc/M_sun   rg/cm  abs(adotGAS) abs(adotGW) abs(adot) a/abs(adot)'
	write(u_snapshots,'(ES17.11,3x,8(ES14.6,3x))') t/year,a,mass_inner/M_sun,mass/M_sun,rg,&
									&	abs(adot1),abs(adot2),abs(adot),a/abs(adot)
	write(u_snapshots,'(a)') 'r/rg   sigma   nu   Teff   Tc   H/r   lambda   r/pc   Prad   Pgas   Prad/Pgas   tnu'
	do i=1,inr-1

		! computes radiation pressure and gas pressure (if H and T are nonzero)
		if ((T_c(i)>prec_out).and.(H(i)>prec_out)) then
			prad = a_rad * T_c(i)**4. /3.
			pgas = rovermu * T_c(i) * sigma(i) / (2.*H(i))
			prad_pgas = prad/pgas
		else
			prad = 0.0d0
			pgas = 0.0d0
			prad_pgas = 0.0d0
		endif
		
		! computes the viscous timescale (if nu is nonzero)
		if (nu(i)>1.0d-5) then
			tnu = x2(2*i)/nu(i)/year
			!write(*,'(5ES14.6)') x(2*i),x(2*i)/rg,x2(2*i),nu(i),tnu
		else
			tnu = 0.0d0
		endif
		
		! writes the output
		write(u_snapshots,f_snapshots) 	xrg(i)				,&				!  1. radial coordinate in rg
									&	max(prec_out,sigma(i))		,&		!  2. sigma
									&	max(prec_out,nu(i))			,&		!  3. viscosity
									&	max(prec_out,T_eff(i))		,&		!  4. effective temperature
									&	max(prec_out,T_c(i))		,&		!  5. central temperature
									&	max(prec_out,H(i)/x(2*i))	,&		!  6. thickness/R
									& 	lambda(i)					,&		!  7. torque
									&	xpc(i)						,&		!  8. radial coordinate in parsec
!									&	inside_rhill(i)						,&		!  8. inside 1 outside 0									
									&	prad						,&		!  9. radiation pressure
									&	pgas						,&		! 10. gas pressure
									&	prad_pgas					,&		! 11. radiation to gas pressure ratio
									&	tnu									! 12. viscous timescale											
	enddo
	write(u_snapshots,f_snapshots) 	xrg(inr),prec_out,prec_out,prec_out,& 
									& prec_out,prec_out,prec_out,xpc(inr),prec_out,prec_out,prec_out,prec_out					
	close(u_snapshots)							
	
	
	! writes the snapshots' filenames list in splash.filenames
	open(unit=u_splashfilenames, file=splashfilenames,access='append')
	write(u_splashfilenames,'(a)') 'snap/snap'//trim(strcount)//'.dat'
	close(u_splashfilenames)
	
	!close(u_sigma_neg)
	!open(unit=u_sigma_neg, file=sigma_negfile,access='append')
end subroutine takesnapshot
!************************************************

!------------------------------------------------
!	writes snapshot files,
!	calculates mass and writes mass file
!------------------------------------------------
subroutine takesnapshot_latest()

	use variables,	only:	i,inr,x,x2,xpc,xrg,&
							& snapfile_latest,strcount_latest,rin,rout,&
							& D,H,lambda,nu,advvel,homedir,tot_sim_mass,&
							& sigma,sigma_negfile,T_eff,T_c,&
							& prec_out,isncount_latest,mass,&
							& prec_out,t,a,mass_inner,prad,pgas,prad_pgas,rovermu,tnu,rg,adot1,adot2,adot
	
	use units, 		only: 	u_snapshots_latest
	
	use formats, 	only: 	f_snapshots
	
	use constants, 	only:	year,M_sun,onethird,a_rad
		
	implicit none
	
	! convert to string the counter of the snapshots
	write(strcount_latest,'(i10.10)') isncount_latest	
	snapfile_latest=trim(homedir)//'/snap_latest/snap_latest'//trim(strcount_latest)//'.dat'


	open(unit=u_snapshots_latest, file=snapfile_latest,status='unknown')	
	write(u_snapshots_latest,'(a)') 'time/yr a/cm M_in/M_sun M_disc/M_sun  rg/cm  abs(adotGAS) abs(adotGW) abs(adot) a/abs(adot)'
 	write(u_snapshots_latest,'(ES17.11,3x,8(ES14.6,3x))') 	t/year,a,mass_inner/M_sun,mass/M_sun,rg,&
 											&	abs(adot1),abs(adot2),abs(adot),a/abs(adot)
	write(u_snapshots_latest,'(a)') 'r/rg   sigma   nu   Teff   Tc   H/r   lambda   r/pc   Prad   Pgas   Prad/Pgas   tnu'
	
	do i=1,inr-1

		! computes radiation pressure and gas pressure (if H and T are nonzero)
		if ((T_c(i)>prec_out).and.(H(i)>prec_out)) then
			prad = a_rad * T_c(i)**4. /3.
			pgas = rovermu * T_c(i) * sigma(i) / (2.*H(i))
			prad_pgas = prad/pgas
		else
			prad = 0.0d0
			pgas = 0.0d0
			prad_pgas = 0.0d0
		endif
		
		! computes the viscous timescale (if nu is nonzero)
		if (nu(i)>1.0d-5) then
			tnu = x2(2*i)/nu(i)/year
			!write(*,'(5ES14.6)') x(2*i),x(2*i)/rg,x2(2*i),nu(i),tnu
		else
			tnu = 0.0d0
		endif
		
		! computes the viscous timescale (if nu is nonzero)
		if (nu(i)>prec_out) then
			tnu = x2(i)/nu(i)/year	
		endif
	
		! writes the output
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

end subroutine takesnapshot_latest
!************************************************

!------------------------------------------------
!	writes snapshot files,
!	calculates mass and writes mass file
!------------------------------------------------
subroutine takesnapshot_parameters()

	use variables,	only:	inr,t,dtau,mass,mdot,tot_acc_mass,ia,a,adot1,adot2,isncount,tot_sim_mass, &
						&	massfilename,mass_inner,rg,Acc_L,L_Edd
	
	use units, 		only: 	u_massfile

	use formats, 	only: 	f_massfile
	
	use constants,	only:	parsec,year,M_sun,M_sun_year,km_s,eta,c_light2
	
	implicit none

	! calculates the Luminosity
	Acc_L=mdot(1) * eta * c_light2
	
	! writes the file mass.dat
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
	
end subroutine takesnapshot_parameters
!************************************************



! ###### USEFUL GENERAL PURPOSE FUNCTIONS ###### !

!------------------------------------------------
!	finds the cell where is the planet.
!	global a, here r0: planet position in x coordinate (it's a real number).
!	the planet is between x(ia) and x(ia+1).
!	here I use the simplest algorithm, but there are more efficient ones.
!------------------------------------------------
integer function findxcell(r,n,r0)
	
	implicit none
	
	integer, intent(in) :: n
	double precision, intent(in) :: r(n)
	double precision, intent(in) :: r0
	
	integer i

		
	do i=1,n
		if (r(i).gt.r0) then
			findxcell=i-1
			return
		endif
	enddo
	
	return
end function findxcell
!************************************************

!------------------------------------------------
!	finds the cell of sigma where is the planet.
!	global a, here r0: planet position in x coordinate (it's a real number).
!------------------------------------------------
integer function findsigmacell(r,n,rr0)
	implicit none
	integer, intent(in) :: n
	double precision, intent(in) :: r(n)
	double precision, intent(in) :: rr0
	
	integer :: findxcell
	
	findsigmacell=int(findxcell(r,n,rr0)/2)
	!print *,findxcell(r,n,r0),findsigmacell ! check Ok
	
	return
end function findsigmacell
!************************************************

!------------------------------------------------
!	self-made signum function
!------------------------------------------------
double precision function signum(a)
	implicit none
	double precision, intent(in) :: a
	
	if (a.gt.(0.)) then
		signum = 1.
	elseif (a.lt.(0.)) then
		signum = -1.
	else
		signum = 0.
	endif
	
	return
end function signum
!************************************************

    

