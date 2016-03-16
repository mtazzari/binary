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


module variables

implicit none
save


double precision :: AA
double precision :: BB
double precision :: Rey

! arrays
double precision, allocatable :: x(:)			! radial coordinate
double precision, allocatable :: xpc(:)			! radial coordinate in parsec
double precision, allocatable :: xrg(:)			! radial coordinate in rg
double precision, allocatable :: x12(:)			! radial coordinate (square root)
double precision, allocatable :: x2(:)			! radial coordinate (square)
double precision, allocatable :: dx(:)			! radial spacing
double precision, allocatable :: dx2(:)			! radial spacing
double precision, allocatable :: xdx(:)

double precision, allocatable :: D(:)			! diffusion coefficient
double precision, allocatable :: sigma(:)		! surface density at time N+1
double precision, allocatable :: oldsigma(:)	! surface density at time N
double precision, allocatable :: exactsigma(:)	! surface density of an analytical solution
double precision, allocatable :: mdot(:)		! accretion rate
double precision, allocatable :: lambda(:)		! adimensional torque
double precision, allocatable :: advvel(:)		! advective velocity
double precision, allocatable :: nu(:)			! viscosity
double precision, allocatable :: H(:)			! scale height
double precision, allocatable :: area(:)		! annulus area
double precision, allocatable :: gg(:)
double precision, allocatable :: omega(:)
double precision, allocatable :: omega2(:)
double precision, allocatable :: inside_rhill(:)
double precision, allocatable :: T_c(:)			! central temperature
double precision, allocatable :: T_eff(:)		! effective temperature


double precision :: Mp,GMp,Ms
double precision :: rg, rg10
double precision :: rISCO
integer :: rout_unit

! viscosity-related variables
double precision :: T_c_prefactor
double precision :: T_eff_prefactor
double precision :: alpha
double precision :: Rovermu
double precision :: aconst
double precision :: bconst
double precision :: cconst
double precision :: thin_corr
double precision :: tau
double precision :: accT
double precision :: errT
integer :: i_t
double precision :: T_guess
double precision :: k_ff
double precision ::	T_new
double precision :: prad_pgas,prad,pgas

! arrays-related variables
integer :: inr						! number of meshpoints where sigma is defined
integer :: inr2						! number of meshpoints where x is defined (inr*2)
	
! cycle-related variables
integer :: i						! counter
integer :: j						! counter
integer :: in_iter					! number of iteration needed to reach tot_time
	
! output-related variables
integer :: inextsnap				! counter of the next snapshot to be taken	
integer :: istepsnapshot			! number of iteration between two snapshots
integer :: isnapnum					! total number of snapshots
integer :: isncount					! counter of the taken snapshots
integer :: isncount_latest			! counter of the taken snapshots just before merger
	
! position-related variables
double precision :: dx_lin			! radial spacing between two meshpoints
double precision :: rin				! inner radius
double precision :: rout			! outer radius
double precision :: r_hill			! hill radius of the planet
double precision :: rlso			! radius of the innermost last stable orbit
	
! mass-related variables
double precision :: initial_mass	! initial mass of the disc (computed at time t=t0)
double precision :: mass			! mass of the disc
double precision :: mass_exact		! mass of the disc of the analytical solution
double precision :: tot_acc_mass	! total accreted mass at time t from t0
double precision :: sigmadot_rin	! boundary condition at the inner edge r-->rin
double precision :: sigmadot_rout	! boundary condition at the outer edge r-->rout
double precision :: rdrip			! radius at which mass is drop onto the disc
double precision :: mdotin			! mass injected in the disc in mass/s
double precision :: sigmadotin
double precision :: oldmass			! mass at the previous snapshot
integer :: initial_sigma_type		! kind of initial sigma distribution
integer :: sigmadot_type			! kind of sigmadot to perform: none, constant added mass, constat flux
double precision :: nu0				! factor of viscosity (for power-law implementation)
double precision :: tot_sim_mass	! total mass present in the simulation grid
double precision :: tot_out_mass 	! net mass 
double precision :: sigma0			! initial surface density of the disc
double precision :: m_disc0			! initial mass of the disc
double precision :: mass_inner		! inner disc mass
integer :: idrip					! index of rdrip

! planet-related variables
integer :: ia						! index of the cell in which the planet is (it's an integer)
integer :: iax

integer :: iILR						! index of the position of the inner lindblad resonance (m=2) (sigma grid)
integer :: iOLR						! index of the position of the outer lindblad resonance (m=1) (sigma grid)
integer :: iILRx					! index of the position of the inner lindblad resonance (m=2) (x grid)
integer :: iOLRx					! index of the position of the outer lindblad resonance (m=1) (x grid)
double precision :: old_a			! previous radial position of the planet
double precision :: a				! current  radial position of the planet
double precision :: a0				! initial  radial position of the planet
double precision :: tmigration		! time at which the planet starts to migrate
double precision :: f				! adimensional parameter of the torque
double precision :: q 				! Msecondary/Mprimary ratio
double precision :: angular_momentum

double precision :: ILR
double precision :: OLR
double precision :: adot1			! acceleration due to viscous torques
double precision :: adot2			! acceleration due to GW radiation emission
double precision :: adot			! total acceleration 

double precision :: acrit
double precision :: tinjectplanet
double precision :: delta
double precision :: factor

integer :: isplanet					! swithc

double precision :: gwprefactor

!double precision :: to_add
! analytical solution variables
double precision,parameter :: xin=1.! initial position of spreading ring, k=r0(S.R.)/r0(adimensionalization). 
double precision :: tau0			! initial time of spreading ring solution, k=r0(S.R.)/r0(adimensionalization). 	

! calculation variables
double precision :: adv
double precision :: diff1					! diffusive term 1
double precision :: diff2					! diffusive term 2
double precision,parameter :: acc_factor=0.1	! accuracy factor for the movement of the planet in the radial grid
character(len=100) :: inputfile
integer 		:: sigmacell
double precision :: sigma_diff
integer :: gridtype
double precision :: deltamin
double precision :: prec_out
double precision :: tnu
double precision :: clipped_torque

! time-related variables
double precision :: t				! global time of the simulation
double precision :: t0				! global initial time of the simulation
double precision :: tot_time		! time that the simulation has to reach	
double precision :: dtau			! timestep
double precision :: dtausnap		! time interval between two snapshots
double precision :: dtausnap_param	! time interval between two snapshots of parameters
double precision :: initial_dtausnap_param	! time interval between two snapshots of parameters
double precision :: dtaumin
double precision :: timenextsnap
double precision :: timenextsnap_param
double precision :: CFL				! Courant stability parameter
double precision :: d_smoothing_in 	! smoothing width of torque at Lindblad resonances
double precision :: d_smoothing_out 	! smoothing width of torque at Lindblad resonances
double precision :: L_Edd			! Eddington Luminosity
double precision :: Acc_L	! Accretion Luminosity

! filenames
character(len=200) :: splashfilenames
character(len=200) :: planetfilename
character(len=200) :: massfilename
character(len=200) :: homedir							! home directory of the current simulation
character(len=200) :: sim_mode
character(len=200) :: snapfile		! filename of the snapshots
character(len=200) :: snapfile_latest	! filename of the snapshots
character(len=200) :: mdotfile		! filename of the mdot snapshots
character(len=200) :: sigma_negfile	! filename of the sigma_neg.dat
character(len=200) :: initial_sigma_file	! filename of the initial sigma used to restart the simulation
character(len=10) strcount			! string casting of isncount; allowed snapnum<10**10
character(len=10) strcount_latest	! string casting of isncount_latest; allowed snapnum<10**10

integer :: today(3) ! will contain the current system date
integer :: now(3)	! will contain the current system time
		
end module variables

module units

integer,parameter :: u_snapshots		=11		! unit for the snapshots	
integer,parameter :: u_snapshots_latest	=12		! unit for the snapshots latest (prior to merger)
integer,parameter :: u_massfile			=13		! unit for the mass file
integer,parameter :: u_splashfilenames	=15 	! unit for the splash.filenames file
integer,parameter :: u_mainoutput		=20		! unit for the main output: read/write defaults
integer,parameter :: u_sigma_neg		=21		! unit for the main output: read/write defaults

end module units


module formats

character(len=30),parameter :: f_snapshots	='(12(ES14.7,1x))'					! snapshot###.dat format
character(len=30),parameter :: f_snapfiles	=''									! snapfile variable format
character(len=60),parameter :: f_massfile	= '(i6,1x,ES17.11,1x,13(ES14.7,1x))'	! mass.dat file format
character(len=60),parameter :: f_sigma_neg	='(ES14.7,1x,i5,1x,3(ES14.7,2x))'	! sigma_neg.dat format
	
end module formats