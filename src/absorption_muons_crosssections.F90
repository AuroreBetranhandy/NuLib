!-*-f90-*-
!always keep term dimensionless unless it is the units of the final answer...
function numu_absorption_on_n(neutrino_energy,eos_variables) result(crosssection)
  
  use nulib
  implicit none

  real*8, intent(in) :: eos_variables(total_eos_variables)  
  real*8, intent(in) :: neutrino_energy !MeV
  real*8 :: crosssection !final answer in cm^2

  !local variables
  real*8 :: final_muon_energy !MeV
  real*8 :: feminus_exp_log10ed,SA_exp_log10ed,feminus_over_SA_exp_log10ed !dimensionless
  real*8 :: one_plus_feminus_exp_log10ed,one_plus_SA_exp_log10ed !dimensionless
  real*8 :: weak_mag !dimensionless
  real*8 :: mu_nu_eq !MeV
  real*8 :: logterm,expterm !dimensionless

  !function declarations
  real*8 :: weak_mag_correction_absorption
  real*8 :: fermidirac_exptermonly_log10ed

  !based on Todd Thompson PhD Appendix B, Eq. B8, final state proton
  !blocking is done outside of this subroutine (bottom of
  !absorption_crosssections.F90)

  final_muon_energy = neutrino_energy + delta_np
  mu_nu_eq = eos_variables(mumuindex)-eos_variables(muhatindex)
  feminus_exp_log10ed = fermidirac_exptermonly_log10ed(final_muon_energy, &
       eos_variables(mumuindex),eos_variables(tempindex)) !final state muon blocking
  SA_exp_log10ed = fermidirac_exptermonly_log10ed(neutrino_energy,mu_nu_eq,eos_variables(tempindex))
  if (do_weak_mag_corrections) then
     weak_mag = weak_mag_correction_absorption(neutrino_energy,0) !to first order 1.0d0+1.01d0*neutrino_energy/m_ref, Horowitz 2002
  else
     weak_mag = 1.0d0
  endif

  !Note on the exp's.  The Fermi functions and simulated absorption
  !terms can be huge/small.  The best way to deal with this is to
  !play tricks.  We write them all in terms of the log of the exp, not
  !the fermi functions, and then combine appropiately (taking into
  !account the size of the exp to deal with the +1's appropiately

  !Note, we can also combine the feminus_exp/SA_exp term
  !the arguement of feminus_exp is neutrino_energy + delta_np - matter_mue
  !the arguement of SA_exp is neutrino_energy - matter_nue + matter_muhat
  !feplus_exp/SA_exp = exponential with arguement of
  !delta_np-matter_muhat
  feminus_over_SA_exp_log10ed = (delta_np-eos_variables(muhatindex))/eos_variables(tempindex)*log10exp

  !The full term we are lumping together is:
  !(1-f_{e-})/(1-f_{\nu_e}^{eq}) =
  !fexp_{e-}*(1+fexp_{SA})/(1+fexp_{e-})/fexp_{SA} =
  !10.0d0**(feminus_over_SA_exp_log10ed - one_plus_feminus_exp_log10ed
  !+ one_plus_SA_exp_log10ed)

  !deal with fermi functions +1's
  if (feminus_exp_log10ed.gt.30.0d0) then !exp >> 1
     one_plus_feminus_exp_log10ed = feminus_exp_log10ed 
  else if (feminus_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_feminus_exp_log10ed = 0.0d0
  else 
     one_plus_feminus_exp_log10ed = log10(1.0d0+10.0d0**(feminus_exp_log10ed))
  endif

  if (SA_exp_log10ed.gt.30.0d0) then  !exp >> 1
     one_plus_SA_exp_log10ed = SA_exp_log10ed
  else if (SA_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_SA_exp_log10ed = 0.0d0
  else
     one_plus_SA_exp_log10ed = log10(1.0d0+10.0d0**(SA_exp_log10ed))
  endif

  logterm = min(200.0d0,max(-200.0d0,feminus_over_SA_exp_log10ed - one_plus_feminus_exp_log10ed + one_plus_SA_exp_log10ed))
  expterm = 10.0d0**(logterm)

  crosssection = sigma0 * & !cm^2
          (1.0d0+3.0d0*gA**2) / 4.0d0 * & !dimensionless
          (final_muon_energy/m_mu)**2 * & !dimensionless
          (1.0d0-(m_mu/final_muon_energy)**2)**0.5d0 * & !dimensionless
          expterm*weak_mag !dimensionless

end function numu_absorption_on_n

function anumu_absorption_on_p(neutrino_energy,eos_variables) result(crosssection)
  
  use nulib
  implicit none

  real*8, intent(in) :: eos_variables(total_eos_variables)  
  real*8, intent(in) :: neutrino_energy !MeV
  real*8  :: crosssection !final answer in cm^2

  !local variables
  real*8 :: final_muonplus_energy !MeV
  real*8 :: feplus_exp_log10ed,SA_exp_log10ed,feplus_over_SA_exp_log10ed !dimensionless
  real*8 :: one_plus_SA_exp_log10ed,one_plus_feplus_exp_log10ed !dimensionless
  real*8 :: weak_mag !dimensionless
  real*8 :: mu_nu_eq !MeV
  real*8 :: expterm, logterm !dimensionless

  !function declarations
  real*8 :: fermidirac_exptermonly_log10ed
  real*8 :: weak_mag_correction_absorption

  !based on Todd Thompson PhD Appendix B, Eq. B11, final state neutron
  !blocking is done outside of this subroutine (bottom of
  !absorption_crosssections.F90)
  final_muonplus_energy = neutrino_energy - delta_np
  if (final_muonplus_energy.lt.2.0d0*m_mu) then
     !only happens if neutrino energy < ~2.3 MeV
     !we choose 2m_e because GR1D has trouble with zero opacity
     final_muonplus_energy = 2.0d0*m_mu
  endif

  mu_nu_eq = -(eos_variables(mumuindex)-eos_variables(muhatindex))
  feplus_exp_log10ed = fermidirac_exptermonly_log10ed(final_muonplus_energy,-eos_variables(mumuindex),eos_variables(tempindex))
  SA_exp_log10ed = fermidirac_exptermonly_log10ed(neutrino_energy,mu_nu_eq,eos_variables(tempindex))
  if (do_weak_mag_corrections) then
     weak_mag = weak_mag_correction_absorption(neutrino_energy,1) !to first order 1.0d0-7.1d0*neutrino_energy/m_p, Horowitz 2002
  else
     weak_mag = 1.0d0
  endif

  !Note on the exp's.  The Fermi functions and simulated absorbtion
  !terms can be huges/small.  The best way to deal with this is to
  !play tricks.  We write them all in terms of the log of the exp, not
  !the fermi functions, and then combine appropiately (taking into
  !account the size of the exp to deal with the +1 appropiately

  !Note, we can also combine the feplus_exp/SA_exp term
  !the arguement of feplus_exp is neutrino_energy - delta_np + matter_mue
  !the arguement of SA_exp is neutrino_energy + matter_nue - matter_muhat
  !feplus_exp/SA_exp = exponential with arguement of
  !-delta_np+matter_muhat
  feplus_over_SA_exp_log10ed = (-delta_np+eos_variables(muhatindex))/eos_variables(tempindex)*log10exp

  !The full term we are lumping together is:
  !(1-f_{e+})/(1-f_{\bar{nu}_e}^{eq}) =
  !fexp_{e+}*(1+fexp_{SA})/(1+fexp_{e+})/fexp_{SA} =
  !10.0d0**(feplus_over_SA_exp_log10ed - one_plus_feplus_exp_log10ed
  !+ one_plus_SA_exp_log10ed)

  !deal with fermi functions +1's
  if (feplus_exp_log10ed.gt.30.0d0) then !exp >> 1
     one_plus_feplus_exp_log10ed = feplus_exp_log10ed 
  else if (feplus_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_feplus_exp_log10ed = 0.0d0
  else 
     one_plus_feplus_exp_log10ed = log10(1.0d0+10.0d0**(feplus_exp_log10ed))
  endif

  if (SA_exp_log10ed.gt.30.0d0) then  !exp >> 1
     one_plus_SA_exp_log10ed = SA_exp_log10ed
  else if (SA_exp_log10ed.lt.-30.0d0) then !exp << 1
     one_plus_SA_exp_log10ed = 0.0d0
  else
     one_plus_SA_exp_log10ed = log10(1.0d0+10.0d0**(SA_exp_log10ed))
  endif

  !add all log terms up together, and also apply limits for extreme cases!!!
  logterm = min(200.0d0,max(-200.0d0,feplus_over_SA_exp_log10ed - one_plus_feplus_exp_log10ed + one_plus_SA_exp_log10ed))

  expterm = 10.0d0**(logterm)

  crosssection = sigma0 * & !cm^2
       (1.0d0+3.0d0*gA**2) / 4.0d0 * & !dimensionless
       (final_muonplus_energy/m_mu)**2 * & !dimensionless
       (1.0d0-(m_mu/final_muonplus_energy)**2)**0.5d0 * & !dimensionless
       expterm*weak_mag !dimensionless

end function anumu_absorption_on_p

!~ function nue_absorption_on_A(neutrino_energy,eos_variables) result(crosssection)

!~   use nulib
!~   implicit none

!~   real*8, intent(in) :: eos_variables(total_eos_variables)  
!~   real*8, intent(in) :: neutrino_energy !MeV
!~   real*8 :: crosssection !final answer in cm^2

!~   !local variables
!~   real*8 :: N_p_Z, N_n_N !dimensionless
!~   real*8 :: final_muon_energy !MeV, energy of muon in products of interaction
!~   real*8 :: feminus_exp_log10ed,SA_exp_log10ed !dimensionless
!~   real*8 :: one_plus_feminus_exp_log10ed,one_plus_SA_exp_log10ed !dimensionless
!~   real*8 :: exponential_term_log10ed !dimensionless, place holder
!~   real*8 :: mu_nu_eq
!~   real*8 :: logterm, expterm

!~   !function declarations
!~   real*8 :: fermidirac_exptermonly_log10ed

!~   !based on Todd Thompson PhD Appendix B, Eq. B13, and BRT06 Eq. 12,
!~   !BRT06 say better expressions out there...  nu_energy + Qprime, as
!~   !in Thompson, mus have rest mass difference included
!~   final_muon_energy = neutrino_energy  + eos_variables(muhatindex) + 3.0d0 
!~   !At low energies, muhat can be quite negative, therefore
!~   !final_energy_energy < m_e... Rampp & Janka set crosssection to zero
!~   if (final_muon_energy.lt.m_e+1.d-10) then
!~      crosssection = 0.0d0
!~      return
!~   endif

!~   mu_nu_eq = eos_variables(mueindex)-eos_variables(muhatindex)
!~   feminus_exp_log10ed = fermidirac_exptermonly_log10ed(final_muon_energy,eos_variables(mueindex),eos_variables(tempindex))
!~   SA_exp_log10ed =fermidirac_exptermonly_log10ed(neutrino_energy,mu_nu_eq,eos_variables(tempindex))
!~   !Boltzmann term for neutron to be in 3.0MeV exicted state (see Thompson and Bruenn 1985) 
!~   exponential_term_log10ed = -3.0d0/eos_variables(tempindex)*log10exp

!~   !Zbar_Nterm 
!~   if (eos_variables(zbarindex).lt.20.0d0) then
!~      N_p_Z = 0.0d0
!~   else if (eos_variables(zbarindex).lt.28.0d0) then
!~      N_p_Z = eos_variables(zbarindex) - 20.0d0
!~   else 
!~      N_p_Z = 8.0d0
!~   endif
  
!~   !Nbar_Nterm 
!~   if (eos_variables(abarindex)-eos_variables(zbarindex).lt.34.0d0) then
!~      N_n_N = 6.0d0
!~   else if (eos_variables(abarindex)-eos_variables(zbarindex).lt.40.0d0) then
!~      N_n_N = 40.0d0 - (eos_variables(abarindex)-eos_variables(zbarindex))
!~   else 
!~      N_n_N = 0.0d0
!~   endif

!~   if (SA_exp_log10ed.gt.30.0d0) then  !exp >> 1
!~      one_plus_SA_exp_log10ed = SA_exp_log10ed
!~   else if (SA_exp_log10ed.lt.-30.0d0) then !exp << 1
!~      one_plus_SA_exp_log10ed = 0.0d0
!~   else
!~      one_plus_SA_exp_log10ed = log10(1.0d0+10.0d0**(SA_exp_log10ed))
!~   endif

!~   if (feminus_exp_log10ed.gt.30.0d0) then  !exp >> 1
!~      one_plus_feminus_exp_log10ed = feminus_exp_log10ed
!~   else if (feminus_exp_log10ed.lt.-30.0d0) then !exp << 1
!~      one_plus_feminus_exp_log10ed = 0.0d0
!~   else
!~      one_plus_feminus_exp_log10ed = log10(1.0d0+10.0d0**(feminus_exp_log10ed))
!~   endif

!~   !note, in principle we do this
!~   !       feminus_exp/(1.0d0+feminus_exp)*(1.0d0+SA_exp)/SA_exp * & !dimensionless
!~   !       exponential_term !dimensionless
!~   !but in practice, the expontential_term*feminus_exp/SA_exp = 1, therefore
!~   logterm = one_plus_SA_exp_log10ed - one_plus_feminus_exp_log10ed
  
!~   expterm = 10.0d0**(logterm)

!~   crosssection = sigma0 / 14.0d0 * gA**2 * & !cm^2
!~        N_p_Z * N_n_N * & ! dimensionless
!~        (final_muon_energy/m_e)**2 * & !dimensionless
!~        (1.0d0 - (m_e/final_muon_energy)**2)**0.5d0 * & !dimensionless
!~        expterm

!~ end function nue_absorption_on_A


