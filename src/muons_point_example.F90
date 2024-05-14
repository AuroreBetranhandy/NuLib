!-*-f90-*-
program muons_point_example

  use nulib
  use inputparser
  
  use eosmodule, ONLY  : energy_shift
#if NUCLEI_HEMPEL
  use nuclei_hempel, only : set_up_Hempel
#endif
  implicit none

  !many people use different number of species, this is to denote how they are devided up.
  ! mytable_neutrino_scheme = 1 (three output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino+antineutrino
  ! neutrino_scheme = 2 (four output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino             #4 mu and tau antineutrino
  ! neutrino_scheme = 3 (six output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon neutrino                 #4 mu antineutrino
  !         #5: tau neutrino                  #6 tau antineutrino
  integer :: mypoint_neutrino_scheme = 3

  !number of species for nulib to calculate interactions for, must
  !be six currently, average is done via above parameter
  integer :: mypoint_number_species = 6

  !number of species for nulib to output interactions for, must
  !be commensurate with neutrino_scheme above
  integer :: mypoint_number_output_species = 6

  !number of energy groups
  integer :: mypoint_number_groups = 18 !18

  character*200 :: parameters = "./parameters"

  !local variables
  real*8, allocatable,dimension(:,:) :: local_emissivity
  real*8, allocatable,dimension(:,:) :: local_absopacity
  real*8, allocatable,dimension(:,:) :: local_scatopacity
  real*8, allocatable,dimension(:,:) :: local_delta
  real*8, allocatable,dimension(:,:) :: local_Phi0, local_Phi1
  real*8, allocatable,dimension(:,:) :: local_Phi0_muons, local_Phi1_muons
  real*8, allocatable,dimension(:,:) :: blackbody_spectra
  real*8, allocatable,dimension(:) :: eos_variables
  real*8, allocatable, dimension(:) :: Phi0_brem,Phi0_brem2
  real*8,allocatable, dimension(:) :: Phi0_brem_ann,Phi0_brem2_ann
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  real*8 :: NES_Phi0_ThompsonBruenn
  real*8 :: proton_number_density !proton number density, # neutron/cm^3
  real*8 :: neutron_number_density !neutron number density, # neutron/cm^3

	real*8 :: matter_muhat0 !
  real*8 :: nue_absorption_on_n !function declaration
  real*8 :: anue_absorption_on_p !function declaration
  real*8 :: numu_absorption_on_n !function declaration
  real*8 :: anumu_absorption_on_p !function declaration
  real*8 :: nue_absorption_on_A !function declaration
  real*8 :: nux_absorption_on_n_and_p !function declaration

  real*8 :: Q,Q2
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
  real*8 :: xrho, xtemp, xye
  real*8 :: n_N,eta_star
  integer :: i,j,inde,in_N,zone
  real*8 :: dxfac,mindx
  real*8 :: fermi
!~   real*8 :: find_s
  real*8,dimension(100) :: integral1,integral2,integral3,temp_array,dens_array,ener_array
  real*8:: integral
  real*8 :: J_1,J_1_bar,phi_a,phi_p,energy,energy_2
  real*8,dimension(100) :: M1_mom,M1_mom_inv
  real*8 :: ye_array(13) = (/(i,i = 43,559,43)/)/1000.0d0  
  real*8 :: rho_array(143) = (/(i,i = 6,148,1)/)/10.0d0  
  real*8,dimension(4) :: rho_figure 
!~   real*8,dimension(100) :: rho_array
  real*8,dimension(13) :: eps_array
  real*8,dimension(796) :: rho_test,R_test,T_test,Ye_test
  character(len=13),dimension(4) :: name_list
  
  
  !fermi function
  real*8 :: fermidirac,fermidirac_dimensionless
  real*8 :: x
  real*8 :: find_s,find_s2
  real*8 :: find_g

  

  real*8,parameter :: nulib_emissivity_gf = 5.59424238d-55/(6.77140812d-06**3*2.03001708d+05)
  real*8,parameter :: nulib_opacity_gf = 1.0d0/6.77140812d-06
  real*8,parameter :: nulib_energy_gf = 1.60217733d-6*5.59424238d-55
  real*8,parameter :: nulib_kernel_gf = 6.77140812d-06**3/2.03001708d+05
  

  call input_parser(parameters)

  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)
#if NUCLEI_HEMPEL
  call set_up_Hempel ! set's up EOS for nuclear abundances
#endif
  !read in EOS table & set reference mass
  call read_eos_table(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)

  !example point
  xrho = 2d14
  xtemp = 46
  xye = 0.2
  
  
  
  !set up energies bins
  do_integrated_BB_and_emissivity = .true.
  mindx = 2.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 2.0d0 !MeV
  bin_bottom(3) = bin_bottom(2)+mindx
  bin_bottom(number_groups) = 200.0d0 ! MeV
  
  call nulib_series2(number_groups-1,bin_bottom(2),bin_bottom(number_groups),mindx,dxfac)
  do i=4,number_groups
     bin_bottom(i) = bin_bottom(i-1)+(bin_bottom(i-1)-bin_bottom(i-2))*dxfac
  enddo
  
  !calculate bin widths & energies from the bottom of the bin & energy at top on bin
  do i=1,number_groups-1
     energies(i) = (bin_bottom(i)+bin_bottom(i+1))/2.0d0
     bin_widths(i) = bin_bottom(i+1)-bin_bottom(i)
     bin_top(i) = bin_bottom(i+1)
  enddo
  
  energies(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups-1)*dxfac/2.0d0
  bin_widths(number_groups) = 2.0*(energies(number_groups)-bin_bottom(number_groups))
  bin_top(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups)
  
  
  
  

  allocate(eos_variables(total_eos_variables))
  eos_variables(:) = 0.0d0
  eos_variables(rhoindex) = xrho
  eos_variables(tempindex) = xtemp
  eos_variables(yeindex) =xye
  
  
  !! EOS stuff
  call set_eos_variables(eos_variables)
  write(*,*) "Rho: ",eos_variables(rhoindex)," g/ccm"
  write(*,*) "T: ",eos_variables(tempindex)," MeV"
  write(*,*) "Ye: ",eos_variables(yeindex)
  write(*,*) "X_n: ",eos_variables(xnindex)
  write(*,*) "X_p: ",eos_variables(xpindex)
  write(*,*) "X_alpha: ",eos_variables(xaindex)
  write(*,*) "X_muons: ",eos_variables(mumuindex)
  
  
  ! ------------ Absorption-----------
  allocate(local_emissivity(3,mypoint_number_groups))
  allocate(local_absopacity(3,mypoint_number_groups))
  allocate(local_scatopacity(2,mypoint_number_groups))
  allocate(local_delta(mypoint_number_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  matter_muhat0 = eos_variables(muhatindex) - delta_np 
  neutron_number_density = max(1.0d-100,eos_variables(xnindex))* &
       eos_variables(rhoindex)/(m_ref*mev_to_gram) !# neutrons/cm^3
  proton_number_density = max(1.0d-100,eos_variables(xpindex))* &
       eos_variables(rhoindex)/(m_ref*mev_to_gram) !# protons/cm^3
       
       
	local_emissivity =0.0d0
	local_absopacity =0.0d0
	
	
  call return_blackbody_spectra(blackbody_spectra,eos_variables)
  
  
  open(unit = 666, file = "muons_emi.txt"  )
  
  
  do i=1,number_groups
  
  
   local_emissivity(2,i) =   numu_absorption_on_n(energies(i),eos_variables)* & !crosssection, cm^2
				             neutron_number_density* &! # neutrons/cm^3
				             max(0.0d0,(proton_number_density/neutron_number_density-1.0d0)/ & !final state proton blocking, =1 if blocking is irrelevant
				             (exp(-matter_muhat0/eos_variables(tempindex))-1.0d0))
				  
  
   local_emissivity(3,i) =  anumu_absorption_on_p(energies(i),eos_variables)* & !crosssection, cm^2
				             proton_number_density* &! # neutrons/cm^3
				             max(0.0d0,(neutron_number_density/proton_number_density-1.0d0)/ & !final state proton blocking, =1 if blocking is irrelevant
				             (exp(-matter_muhat0/eos_variables(tempindex))-1.0d0))
				            
	local_absopacity(2,i) =  anumu_absorption_on_p(energies(i),eos_variables) / blackbody_spectra(3,i)
	local_absopacity(3,i) =  numu_absorption_on_n(energies(i),eos_variables) / blackbody_spectra(3,i)     
				             
    write(*,*) blackbody_spectra(2,i)   
				
!~ 	local_emissivity(1,i) = local_emissivity(1,i) + anue_absorption_on_p(energies(i),eos_variables)* & !crosssection, cm^2
!~ 				             proton_number_density* & ! # protons/cm^3
!~ 				             max(0.0d0,(neutron_number_density/proton_number_density-1.0d0)/ & !final state neutron blocking, =1 if blocking is irrelevant
!~ 				             (exp(matter_muhat0/eos_variables(tempindex))-1.0d0)) 
				             
			             
			             
	local_emissivity(1,i) = local_emissivity(1,i) + nue_absorption_on_n(energies(i),eos_variables)* & !crosssection, cm^2
				             neutron_number_density* &! # neutons/cm^3
				             max(0.0d0,(proton_number_density/neutron_number_density-1.0d0)/ & !final state proton blocking, =1 if blocking is irrelevant
				             (exp(-matter_muhat0/eos_variables(tempindex))-1.0d0))
  
  write(666,*) energies(i),local_emissivity(1,i),local_emissivity(2,i),local_emissivity(2,i) &
			,local_absopacity(1,i),local_absopacity(2,i),local_absopacity(2,i)
  enddo
  
  
  
  close(666)
  
  deallocate(local_emissivity)
  deallocate(local_absopacity)
  deallocate(local_scatopacity)
  deallocate(local_delta)
  
  ! ------------ Scattering
  
  allocate(local_Phi0(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_Phi1(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_Phi0_muons(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_Phi1_muons(mypoint_number_output_species,mypoint_number_groups))
  
  
  
  
  
  open(unit = 404, file = "muons.txt"  )

  Q=0.0d0
  Q=0.0d0


  do i=1,number_groups 
  call single_muonscat_point_return_all(i,eos_variables(mumuindex)/2.0d0, &
                   eos_variables(tempindex),local_Phi0_muons,local_Phi1_muons,mypoint_neutrino_scheme)
  call single_Ipoint_return_all(i,eos_variables(mueindex), &
                   eos_variables(tempindex),local_Phi0,local_Phi1,mypoint_neutrino_scheme)
                   
  Q =  energies(i)**3	 /(2.0d0* pi * hbarc_mevcm)**3&
	  *2.0d0*pi /clight &
	 * local_Phi0(2,i) 
                   
                   
  Q2 = energies(i)**3	 /(2.0d0* pi * hbarc_mevcm)**3&
	  *2.0d0*pi /clight &
	 * local_Phi0_muons(2,i) 
                   
                   
                   
!~    write(*,*) local_Phi0, local_Phi1
   write(404,*) energies(i),energies(i)**3/(2.0d0* pi * hbarc_mevcm)**3&
					*2.0d0*pi /clight*local_Phi0(:,i),&
				energies(i)**3/(2.0d0* pi * hbarc_mevcm)**3&
					*2.0d0*pi /clight*local_Phi0_muons(:,i)
  enddo
  
  close(404)
deallocate(local_Phi0)
deallocate(local_Phi1)
deallocate(local_Phi0_muons)
deallocate(local_Phi1_muons)

end program muons_point_example
  
