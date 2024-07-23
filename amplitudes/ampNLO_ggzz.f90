module mod_ampNLO_ggzz
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  
  implicit none
  private
  
  public :: getewampsNLO_ggzz
  
contains

  !-- g(1) g(2) -> [Z->l(3) lb(4)] [Z->l(5) lb(6)]
  !-- returns both the LO and NLO amplitudes, with delta(a,b) asontwopi[asontwopi**2] factored out
  subroutine getewampsNLO_ggzz(za,zb,sprod,ELL,ELR,ewamp_ggzz)
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(in) :: ELL(9),ELR(9)
    complex(dp), intent(out) :: ewamp_ggzz(-1:1,-1:1,-1:1,-1:1)
    real(dp) :: ma2,mb2,ss,tt
    complex(dp) :: E1LL(9),E1LR(9),E2LL(9),E2LR(9)
    complex(dp) :: propz34,propz56,propg34,propg56
    real(dp) :: Czz,Cgg,Czg,Cww,cgl1(-1:1),cgl2(-1:1),cZl1(-1:1),cZl2(-1:1)
    real(dp) :: pref
    integer :: i3,i5

    !--- propagator factors, 
    !-- z exchange
    propz34=1/(sprod(3,4)-mzsq + ci*mz*gaz)
    propz56=1/(sprod(5,6)-mzsq + ci*mz*gaz)
    !-- gamma exchange
    propg34=1/(sprod(3,4))
    propg56=1/(sprod(5,6))

    !--- form factors
    ma2 = sprod(3,4)
    mb2 = sprod(5,6)
    ss = sprod(1,2)
    tt = ma2 - sprod(1,3)-sprod(1,4)

    !--- couplings
    
    !-- gamma-leptons couplings
    cgl1(-1)=-sqrt2*Qel*sw
    cgl1(+1)=-sqrt2*Qel*sw
    cgl2(-1)=-sqrt2*Qel*sw
    cgl2(+1)=-sqrt2*Qel*sw

    !-- Z-leptons couplings
    cZl1(-1)=Lel/sqrt2/cw
    cZl1(+1)=Rel/sqrt2/cw
    cZl2(-1)=Lel/sqrt2/cw
    cZl2(+1)=Rel/sqrt2/cw

    !-- ZZ fermion loop
    Czz = (nup*(Aup**2+Vup**2)+ndn*(Adn**2+Vdn**2))/two/cosW2
    !-- gamma gamma fermion loop
    Cgg = two/9.0_dp*sinW2*(ndn+four*nup)
    !-- Zg fermion loop
    Czg = -sw/cw*(nup*Vup*Qup+ndn*Vdn*Qdn)

    ewamp_ggzz(:,:,-1,-1) = ampspinors(1,2,3,4,5,6,za,zb,ELL,ELR)
    ewamp_ggzz(:,:,-1,+1) = ampspinors(1,2,3,4,6,5,za,zb,ELL,ELR)
    ewamp_ggzz(:,:,+1,-1) = ampspinors(1,2,4,3,5,6,za,zb,ELL,ELR)
    ewamp_ggzz(:,:,+1,+1) = ampspinors(1,2,4,3,6,5,za,zb,ELL,ELR)

    ! now the couplings, Z + gamma exchange!

    pref = (gwsq/two)**2

    do i3=-1,1,2
    do i5=-1,1,2
       ewamp_ggzz(:,:,i3,i5) = pref*propz34*propz56*Czz*cZl1(i3)*cZl2(i5)*ewamp_ggzz(:,:,i3,i5) &
                              + pref*propg34*propg56*Cgg*cgl1(i3)*cgl2(i5)*ewamp_ggzz(:,:,i3,i5) &
                              + pref*propz34*propg56*Czg*czl1(i3)*cgl2(i5)*ewamp_ggzz(:,:,i3,i5) &
                              + pref*propg34*propz56*Czg*cgl1(i3)*czl2(i5)*ewamp_ggzz(:,:,i3,i5)

    enddo
    enddo

    !-- this should be multiplied by ci to get MCFM notation

    return

  end subroutine getewampsNLO_ggzz

!!  subroutine get_eformfactors(ss,tt,ma2,mb2,E1LL,E1LR,E2LL,E2LR)
!!    use, intrinsic :: iso_c_binding
!!    real(dp), intent(in) :: ss,tt,ma2,mb2
!!    complex(dp), intent(out) :: E1LL(9),E1LR(9),E2LL(9),E2LR(9)
!!    integer :: j
!!    interface
!!       subroutine ampggvv(ss,tt,ma2,mb2,E1LLr,E1LRr,E2LLr,E2LRr,E1LLi,E1LRi,E2LLi,E2LRi) bind (c)
!!         use iso_c_binding
!!         real ( c_double ), VALUE :: ss
!!         real ( c_double ), VALUE :: tt
!!         real ( c_double ), VALUE :: ma2
!!         real ( c_double ), VALUE :: mb2
!!         real ( c_double ) :: E1LLr(*)
!!         real ( c_double ) :: E1LRr(*)
!!         real ( c_double ) :: E2LLr(*)
!!         real ( c_double ) :: E2LRr(*)
!!         real ( c_double ) :: E1LLi(*)
!!         real ( c_double ) :: E1LRi(*)
!!         real ( c_double ) :: E2LLi(*)
!!         real ( c_double ) :: E2LRi(*)
!!       end subroutine ampggvv
!!    end interface
!!    integer ( c_int ), parameter :: n=9
!!    real (c_double) E1LLr(n)
!!    real (c_double) E1LRr(n)
!!    real (c_double) E2LLr(n)
!!    real (c_double) E2LRr(n)
!!    real (c_double) E1LLi(n)
!!    real (c_double) E1LRi(n)
!!    real (c_double) E2LLi(n)
!!    real (c_double) E2LRi(n)
!!
!!    call ampggvv(ss,tt,ma2,mb2,E1LLr,E1LRr,E2LLr,E2LRr,E1LLi,E1LRi,E2LLi,E2LRi)
!!
!!    !---factor "-2" would be needed 
!!    !---to agree with overall normalization from arXiv:1503.08759
!!
!!    do j=1,9
!!       E1LL(j) = ( E1LLr(j) + ci*E1LLi(j) )
!!       E1LR(j) = ( E1LRr(j) + ci*E1LRi(j) )
!!       E2LL(j) = ( E2LLr(j) + ci*E2LLi(j) )
!!       E2LR(j) = ( E2LRr(j) + ci*E2LRi(j) )
!!    enddo
!!
!!    return
!!
!!  end subroutine get_eformfactors


  function ampspinors(j1,j2,j3,j4,j5,j6,za,zb,ELL,ELR)
    complex(dp) :: ampspinors(-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp) :: AmpLL,AmpLR,AmpRR,AmpRL
    complex(dp), intent(in) :: za(6,6),zb(6,6),ELL(9),ELR(9)

    ! compute 4 independent configurations for different helicities of incoming gluons

    AmpLL = ( zb(j1,j3)*za(j3,j2) + zb(j1,j4)*za(j4,j2) )*za(j1,j2)/zb(j1,j2) &
            *(( zb(j2,j3)*za(j3,j1) + zb(j2,j4)*za(j4,j1) ) & 
               *( ELL(1)*za(j3,j5)*zb(j4,j6) & 
                + ELL(2)*za(j1,j3)*za(j1,j5)*zb(j1,j4)*zb(j1,j6) & 
                + ELL(3)*za(j1,j3)*za(j2,j5)*zb(j1,j4)*zb(j2,j6) &  
                + ELL(4)*za(j2,j3)*za(j1,j5)*zb(j2,j4)*zb(j1,j6) & 
                + ELL(5)*za(j2,j3)*za(j2,j5)*zb(j2,j4)*zb(j2,j6) ) & 
              + ELL(6)*za(j1,j3)*za(j1,j5)*zb(j1,j4)*zb(j2,j6) &
              + ELL(7)*za(j1,j3)*za(j1,j5)*zb(j2,j4)*zb(j1,j6) &
              + ELL(8)*za(j1,j3)*za(j2,j5)*zb(j2,j4)*zb(j2,j6) &
              + ELL(9)*za(j2,j3)*za(j1,j5)*zb(j2,j4)*zb(j2,j6) )


    AmpLR = ( zb(j2,j3)*za(j3,j1) + zb(j2,j4)*za(j4,j1) ) &
           *(( zb(j2,j3)*za(j3,j1) + zb(j2,j4)*za(j4,j1) ) & 
               *( ELR(1)*za(j3,j5)*zb(j4,j6) & 
                + ELR(2)*za(j1,j3)*za(j1,j5)*zb(j1,j4)*zb(j1,j6) & 
                + ELR(3)*za(j1,j3)*za(j2,j5)*zb(j1,j4)*zb(j2,j6) &  
                + ELR(4)*za(j2,j3)*za(j1,j5)*zb(j2,j4)*zb(j1,j6) & 
                + ELR(5)*za(j2,j3)*za(j2,j5)*zb(j2,j4)*zb(j2,j6) ) & 
              + ELR(6)*za(j1,j3)*za(j1,j5)*zb(j1,j4)*zb(j2,j6) &
              + ELR(7)*za(j1,j3)*za(j1,j5)*zb(j2,j4)*zb(j1,j6) &
              + ELR(8)*za(j1,j3)*za(j2,j5)*zb(j2,j4)*zb(j2,j6) &
              + ELR(9)*za(j2,j3)*za(j1,j5)*zb(j2,j4)*zb(j2,j6) )


    AmpRR = ( conjg(zb(j1,j4)*za(j4,j2) + zb(j1,j3)*za(j3,j2)) )*conjg(za(j1,j2)/zb(j1,j2)) &
            *(( conjg(zb(j2,j4)*za(j4,j1) + zb(j2,j3)*za(j3,j1)) ) & 
               *( ELL(1)*conjg(za(j4,j6)*zb(j3,j5)) & 
                + ELL(2)*conjg(za(j1,j4)*za(j1,j6)*zb(j1,j3)*zb(j1,j5)) & 
                + ELL(3)*conjg(za(j1,j4)*za(j2,j6)*zb(j1,j3)*zb(j2,j5)) &  
                + ELL(4)*conjg(za(j2,j4)*za(j1,j6)*zb(j2,j3)*zb(j1,j5)) & 
                + ELL(5)*conjg(za(j2,j4)*za(j2,j6)*zb(j2,j3)*zb(j2,j5)) ) & 
              + ELL(6)*conjg(za(j1,j4)*za(j1,j6)*zb(j1,j3)*zb(j2,j5)) &
              + ELL(7)*conjg(za(j1,j4)*za(j1,j6)*zb(j2,j3)*zb(j1,j5)) &
              + ELL(8)*conjg(za(j1,j4)*za(j2,j6)*zb(j2,j3)*zb(j2,j5)) &
              + ELL(9)*conjg(za(j2,j4)*za(j1,j6)*zb(j2,j3)*zb(j2,j5)) )


    AmpRL = ( conjg(zb(j2,j4)*za(j4,j1) + zb(j2,j3)*za(j3,j1)) ) &
            *(( conjg(zb(j2,j4)*za(j4,j1) + zb(j2,j3)*za(j3,j1)) ) & 
               *( ELR(1)*conjg(za(j4,j6)*zb(j3,j5)) & 
                + ELR(2)*conjg(za(j1,j4)*za(j1,j6)*zb(j1,j3)*zb(j1,j5)) & 
                + ELR(3)*conjg(za(j1,j4)*za(j2,j6)*zb(j1,j3)*zb(j2,j5)) &  
                + ELR(4)*conjg(za(j2,j4)*za(j1,j6)*zb(j2,j3)*zb(j1,j5)) & 
                + ELR(5)*conjg(za(j2,j4)*za(j2,j6)*zb(j2,j3)*zb(j2,j5)) ) & 
              + ELR(6)*conjg(za(j1,j4)*za(j1,j6)*zb(j1,j3)*zb(j2,j5)) &
              + ELR(7)*conjg(za(j1,j4)*za(j1,j6)*zb(j2,j3)*zb(j1,j5)) &
              + ELR(8)*conjg(za(j1,j4)*za(j2,j6)*zb(j2,j3)*zb(j2,j5)) &
              + ELR(9)*conjg(za(j2,j4)*za(j1,j6)*zb(j2,j3)*zb(j2,j5)) )


    ampspinors(-1,-1) = AmpLL
    ampspinors(-1,+1) = AmpLR
    ampspinors(+1,+1) = AmpRR
    ampspinors(+1,-1) = AmpRL

    return

  end function ampspinors

end module mod_ampNLO_ggzz

