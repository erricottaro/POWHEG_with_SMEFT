module mod_amplitudes_parameters
  use mod_types; use common_def; use mod_consts_dp
  implicit none
  private

  public :: set_amplitudes_parameters, update_amplitudes_parameters

contains

  subroutine set_amplitudes_parameters(mz_in,gz_in,mw_in,gw_in,sw2_in,e_in,mh_in,gh_in,mt_in, mb_in, nf_in, nup_in, ndn_in, ng_in, decay1_in, decay2_in, proc_in, contr_in, cggh_in, cb_in, ct_in, cz_in)

    real(dp), intent(in) :: mz_in,gz_in,mw_in,gw_in,sw2_in,e_in,mh_in,gh_in,mt_in,mb_in
    integer, intent(in) :: nf_in,nup_in,ndn_in,ng_in, decay1_in, decay2_in
    character, intent(in) :: proc_in*(10), contr_in*(10)
    ! - SMEFT couplings
    real(dp), intent(in) :: cggh_in, cb_in, ct_in, cz_in
    
  
     mh = mh_in
     mhsq = mh**2
     mz=mz_in
     mzsq = mz**2
     mw=mw_in
     mwsq = mw**2
     mt = mt_in
     mtsq = mt**2
     mb = mb_in
     mbsq = mb**2
     GaW=gw_in
     GaZ=gz_in
     GaH=gh_in
     gwsq = e_in**2/sw2_in
     sinW2=sw2_in
     cosW2 = one-sinW2
     sw = sqrt(sinW2)
     cw = sqrt(cosW2)
     vev = two*mw/sw
     ! - SMEFT couplings
     cggh = cggh_in
     cb = cb_in
     ct = ct_in
     cz = cz_in

     Vup = 1.0_dp/2.0_dp - 4.0_dp/3.0_dp*sinW2
     Vdn = -1.0_dp/2.0_dp + 2.0_dp/3.0_dp*sinW2
     Vel = -half + 2.0_dp*sinW2

     Lup = Vup + Aup
     Rup = Vup - Aup
     Ldn = Vdn + Adn
     Rdn = Vdn - Adn

     Lel = Vel + Ael
     Rel = Vel - Ael

     if (proc_in.eq."ZZ") then
        if(mod(abs(decay1_in),2).eq.1) then
           Llep=Lel
           Rlep=Rel
        else
           Llep=Lnu
           Rlep=Rnu
        endif
     endif

     select case (contr_in)
        case ("full")
           contr="full"
        case ("only_h")
           contr="sigl"
        case ("no_h")
           contr="bkgd"
        case ("interf_h")
           contr="intf"
        case default
           stop "contr not allowed"
     end select

! Higgs couples  with bottom  
     if( mb>0)  hwithb = .true.
     
     nf = nf_in
     nflav = nf
     b0=(xn*11.0_dp - 2.0_dp*nflav)/6.0_dp ! eqn 5.11 of paper (02/04/2024)
     b1=(17.0_dp*Ca**2 - 5.0_dp*Ca*nf - 3.0_dp*Cf*nf)/6.0_dp

     CfNf = CF*nf
     CaNf = CA*nf
     NfNf = nf*nf

     nup=nup_in
     ndn=ndn_in
     ng=ng_in

     expmass=mt_in
     HiggsExp = .false.       
     VVExp = .false.
     mlessloops = .true.
     massloops = .true.

  end subroutine set_amplitudes_parameters


  subroutine update_amplitudes_parameters(mu_in)

    real(dp), intent(in) :: mu_in

    mu = mu_in
    musq = mu_in**2

  end subroutine update_amplitudes_parameters

end module mod_amplitudes_parameters
