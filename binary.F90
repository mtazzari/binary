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


program binary
	
	! defining variables 
    use variables, only: i,D,dtau,dx,x,x12,inr,a,ia,inr2,lambda,advvel, &
& oldsigma,sigma,isncount,t,diff1,diff2,tot_acc_mass,q,adv,&
& sigmadot_rin,sigmadot_rout,j,mdot,nu,r_hill,tot_time,adot,iax,&
& timenextsnap,timenextsnap_param,mass,old_a,massfilename,planetfilename,splashfilenames,angular_momentum,inputfile,&
& sigmadotin,dtausnap,dtausnap_param,tot_sim_mass,tot_out_mass,isplanet,adot1,adot2,thin_corr,isncount_latest,mass_inner,&
& sigmadot_type,initial_sigma_type,sim_mode,tinjectplanet,Rey,sigma_diff,CFL,gg,initial_mass,H,acrit,a0,delta,prec_out,rg,rg10
    
    use formats, 	only: 	f_sigma_neg
    
    use units, 		only: 	u_massfile,u_splashfilenames,u_snapshots,u_sigma_neg
    
    use constants, 	only: 	pi,fourpi,threepi,twothirds,sixpi,year,M_sun,parsec
    

    implicit none
    
	! defining functions
	integer :: findxcell
	
	! use count variable to make a number of main cycles equal to : count
	! look for other occurrences of count variable to understand this check system
	!integer :: count
	!count=0
	
	! read inputfile from the first argument
	call getarg(1,inputfile)
	if (inputfile=='') then
		print *, 'ERROR: invalid input file'
		stop
	endif
	inputfile = trim(adjustl(inputfile))
	!print*,'begin', inputfile,'end'! check

	! initialise this variable for a correct computation of thickness since the first cycle.
	thin_corr=1.
	
	! sets the lower limit for output values. Trick for a correct visualisation with Splash
	! I set prec_out to 3.3d-25 so it is recognisable where the output correction occurred.
	prec_out = 3.3333d-20

	call read_variables()
	call initialize_files()
	
	if (sim_mode=='manual') then
		write(*,*)
		write(*,*) '============================================='
		write(*,*) '|                  BINARY                   |'
		write(*,*) '|             by  Marco Tazzari             |'
		write(*,*) '|          marco.tazzari@gmail.com          |'
		write(*,*) '|                22.05.2013                 |'
		write(*,*) '============================================='
		write(*,*)
	endif
	
	call dimensional_variables()
	
	call allocate_arrays()
	call set_grid()
	call set_time()
	
	! sets the initial sigma distribution with the switch initial_sigma_type
	call set_initialsigma(initial_sigma_type)	
	
	! sets the sigmadot
	!sigmadotin=mdotin/(2.*pi*x(2*idrip)*2.*dx(2*idrip))

	! initialization of the velocity in order to calculate correctly the first dtau
	do i=1,inr-1
		advvel(i)=1.
	enddo

	open(unit=u_splashfilenames, file=splashfilenames)
	
	! writes the initial snapshot and initialise the heading
	open(unit=u_massfile,file=massfilename)
	write(u_massfile,'(a,a)') 'count   time/yr   mass/M_sun   mass_inner/M_sun   tot_acc_mass/M_sun   tot_sim_mass/M_sun',&
		& '   mdot(1)/M_sun/yr   mdot(inr)/M_sun/yr   a/cm   a/pc   a/rg   adotGAS   adotGW   Acc_L/L_Edd'
	isncount=0


	! sets the secondary position in the radial grid. Needed to compute the initial inner disc mass.
	ia = findxcell(x,inr2,a0)/2
		
	! calculates initial mass
	call calculate_mass()
	
	
	! takes the first snapshot of the initial configuration
	a=a0
	tot_sim_mass = mass
	call takesnapshot()
	call takesnapshot_parameters()
	isncount=0

	
	write(*,*)
	write(*,'(a,a,$)') '  SIMULATION MODE: ',trim(sim_mode)
	if (INITIAL_SIGMA_TYPE==6) then
		! if it is a restart
		write(*,'(a)') '   RESTART'
	endif
	write(*,*)
	write(*,'(a,T6,a,T17,a,T31,a,T45,a,T57,a,T70,a)') '| sn ','|  time','|  disc mass','| inner mass  ','|    a/rg',&
		& '|  adotDISC ' ,'|  adotGW  |'
	!write(*,'(a)') '|snap|  time  |   disc mass    |   inner mass  |   a   |   adotDISC   |   adotGW    '
	write(*,'(a)') '---------------------------------------------------------------------------------------------'
	Write(*,'(i4,T6,ES10.3,T17,ES13.6,T31,ES13.6,T45,ES13.6,T59,ES13.6,T73,ES13.6)') isncount,&
		& t/year,mass/M_sun,mass_inner/M_sun,a/rg,0.0d0	! initial snapshot
	isncount=1

	!___________________________________!	
	!									!
	!	BEGINNING OF THE MAIN CYCLE		!
	!___________________________________!
	r_hill=q**(1/3.)*a
	CFL=0.1d0
	! saving values
	do i=2,inr
		oldsigma(i) = sigma(i)
		!print *,i,oldsigma(i) !check
	enddo
	a=a0
	old_a=a

	! sets the first disc thickness
#if defined(__ONLYPGAS__) 
		call set_thickness_pgas()
#endif
#if defined(__ENERGYCHANG__)
		call set_thickness_chang()
#endif
#if defined(__PRINGLE__)
		call set_thickness()
#endif

	
	! sets the first torques
	if (isplanet==1) call set_planet()
	
	! sets the first timestep	
	call set_dtau()
	
	maincycle: do while (t <= tot_time)

		! sets the new disc thickness
#if defined(__ONLYPGAS__) 
		call set_thickness_pgas()
#endif
#if defined(__ENERGYCHANG__)
		call set_thickness_chang()
#endif
#if defined(__PRINGLE__)
		call set_thickness()
#endif
			
		do i = 1, inr
			gg(i)=x12(2*i)*oldsigma(i)*nu(i)
		enddo
		
! Since I am executing all simulations with a planet
!		if (isplanet==1) then
			! in the case there is a planet
			
			do i=2,inr-1				
				! computes FIRST diffusive terms
				diff1=x12(2*i+1)*(gg(i+1)-gg(i))/dx(2*i+1)
				diff2=x12(2*i-1)*(gg(i)-gg(i-1))/dx(2*i-1)
    			
    			! computes THEN advection (tidal) term
    			! performs the convective derivative according to the direction of the fluid local velocity
				if (i <= ia) then
					adv=(oldsigma(i+1)*advvel(i+1)-oldsigma(i)*advvel(i))/dx(2*i+1)
				else
					adv=(oldsigma(i)*advvel(i)-oldsigma(i-1)*advvel(i-1))/dx(2*i-1)
				endif
				
				
				! computes the new sigma
				sigma_diff=(3.*(diff1-diff2)/dx(2*i)-adv)*dtau/x(2*i)
				sigma(i)=oldsigma(i)+sigma_diff
				
			enddo

		! Mass conservation: at any time mass+tot_acc_mass must be constant
		mdot(1)=sixpi*(x12(3)*gg(2))/dx(3)
		mdot(inr)=sixpi*x12(inr2-1)*(gg(inr)-gg(inr-1))/dx(inr2-1)
		!tot_acc_mass=tot_acc_mass+(mdot(1)-mdot(inr))*dtau
		
		! total mass accreted onto the central object
		tot_acc_mass=tot_acc_mass+mdot(1)*dtau
		tot_out_mass=tot_out_mass+mdot(inr)*dtau
		!print*, mdot(1)*dtau,mdot(inr)*dtau,tot_acc_mass,tot_out_mass
		!print *,'mdot(1)',mdot(1),'mdot(inr)',mdot(inr),'tot_acc_mass',tot_acc_mass ! check 

		! when GW-driven dominates, it takes parameters snapshot at every timestep until merger
		! if binary separation is smaller than acrit, then take a snapshot of parameters at every rg
		if ((a <= acrit).and.(isplanet==1)) then 
		
			! calculates mass and disc-planet exchanged angular momentum 
			call calculate_mass()
			tot_sim_mass=mass+tot_acc_mass-tot_out_mass
			Write(*,'(i4,T6,ES10.3,T17,ES13.6,T31,ES13.6,T45,ES13.6,T59,ES13.6,T73,ES13.6)') &
			& isncount,t/year,mass/M_sun,mass_inner/M_sun,a/rg,adot1,adot2
			
			call takesnapshot_parameters()
			call takesnapshot_latest()
			
			isncount_latest = isncount_latest-1
			acrit = acrit-rg
			
			! if the separation is smaller than 2rg, then takes snapshots of parameters every 1 day (86400 s)
			if (acrit <= rg10) then 
				dtausnap_param = 86400.0d0
				! lowers the cfl by 10 times
				CFL=0.01d0
				timenextsnap_param = t
			endif
			!print *, acrit/rg

		endif
		
		! if time is due, write the snapshot
		if (t == timenextsnap_param) then
			
			! calculates mass and disc-planet exchanged angular momentum 
			call calculate_mass()

			tot_sim_mass=mass+tot_acc_mass-tot_out_mass
			! prints in Terminal the runtime control parameters
			Write(*,'(i4,T6,ES10.3,T17,ES13.6,T31,ES13.6,T45,ES13.6,T59,ES13.6,T73,ES13.6)') &
			& isncount,t/year,mass/M_sun,mass_inner/M_sun,a/rg,adot1,adot2

			! writes the snapshots and the other dump files
			call takesnapshot_parameters()
			
			if (t == timenextsnap) then
				call takesnapshot()
				isncount=isncount+1
				timenextsnap=timenextsnap+dtausnap	
			endif
			
			timenextsnap_param = timenextsnap_param + dtausnap_param
			if (isplanet==3) then
				! this stops the simulation smoothly
				t = 2*tot_time
			endif
			!if (abs(oldmass-mass)<=2.0d32) then
				! if the difference in mass between two snapshots is smaller than a certain value (i.e. 0.1 MSun)
				! then the steadiness is considered reached.
			!	t = 2* tot_time
			!endif
		endif
		
		
		! updates boundary conditions
		! it does here the update so in the snapshots are written the boundary values of sigma(1) and sigma(inr)
		sigma(1)=sigmadot_rin
		!sigma(inr)=sigmadot_rout
		
		! saving values at previous time
		do i=2,inr-1
			oldsigma(i) = sigma(i)
			!print *,i,oldsigma(i) !check
		enddo
		old_a=a
	
	
		if (isplanet == 1) call set_planet()
		
		! advancing global time
		call set_dtau() 
		!print *,dtau ! check OK
		
		if (((timenextsnap_param-t) > 0.).and.((timenextsnap_param-t) < dtau)) then
			 dtau=timenextsnap_param-t
			 !print *,t,dtau,timenextsnap_param
		endif
		t=t+dtau
		
		! add sigmadot
		call set_sigmadot(sigmadot_type)
	end do maincycle
	
	close(u_massfile)
	close(u_splashfilenames) 
	close(u_sigma_neg)
	
	call deallocate_arrays()

	write(*,*)
	write(*,'(a)') '=======================  E N D    O F    S I M U L A T I O N  ======================='
	write(*,*)
	
	return
	
end program binary

