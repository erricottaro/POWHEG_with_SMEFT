module mod_ampNLO_ggww
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  
  implicit none
  private

  public :: getewampsNLO_ggww
  
contains

  !-- g(1) g(2) -> [Wp->l(3) lb(4)] [Wm->l(5) lb(6)]
  !-- returns both the LO and NLO amplitudes, with delta(a,b) asontwopi[asontwopi**2] factored out
  subroutine getewampsNLO_ggww(za,zb,sprod,E2LL,E2LR,ewampNLO_ggww)
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    complex(dp), intent(in) :: E2LL(9),E2LR(9)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: ewampNLO_ggww(-1:1,-1:1,-1:1,-1:1)
    real(dp) :: ma2,mb2,ss,tt
    complex(dp) :: propw34,propw56,pref
    real(dp) :: Cww


    !-- w exchange
    propw34=one/(sprod(3,4)-mwsq + ci*mw*gaw)
    propw56=one/(sprod(5,6)-mwsq + ci*mw*gaw)

    !--- form factors
    ma2 = sprod(3,4)
    mb2 = sprod(5,6)
    ss = sprod(1,2)
    tt = ma2 - sprod(1,3)-sprod(1,4)

    !--- couplings
    Cww = ng*half
    pref = Cww*(gwsq/two)**2*propw34*propw56

    ewampNLO_ggww(:,:,-1,-1) = pref*ampspinors(1,2,3,4,5,6,za,zb,E2LL,E2LR)

    return

  end subroutine getewampsNLO_ggww



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

  ! !-- g(1) g(2) -> [Wp -> l(3) lb(4)] [Wm -> l(5) lb(6)] @ two loops
  ! !-- return coefficient of (as/twopi)**4
  ! subroutine resNLOV_WW(p,res1l,res2l,poles)
  !   real(dp), intent(in) :: p(4,6)
  !   real(dp), intent(out) :: res1l,res2l,poles(-2:-1)
  !   complex(dp) :: ampwwNLO1l(-1:1,-1:1,-1:1,-1:1),ampwwNLO2l(-1:1,-1:1,-1:1,-1:1)
  !   complex(dp) :: za(6,6),zb(6,6)
  !   real(dp) :: pdum(4,6),sprod(6,6),s12
  !   integer :: i1,i2
  !   real(dp), parameter :: colf = 8.0_dp ! Ca**2-one
    
    
  !   res1l = zero
  !   res2l = zero
    
  !   call awayfromzaxis(6,p,pdum)

  !   !print *, 'pdum(1) ', pdum(:,1)
  !   !print *, 'pdum(2) ', pdum(:,2)
  !   !print *, 'pdum(3) ', pdum(:,3)
  !   !print *, 'pdum(4) ', pdum(:,4)
  !   !print *, 'pdum(5) ', pdum(:,5)
  !   !print *, 'pdum(6) ', pdum(:,6)

  !   call spinorur(6,pdum,za,zb,sprod)

  !   call getewampsLONLO_ggww(za,zb,sprod,ampwwNLO1l,ampwwNLO2l)

    
  !   ! summing helicities squared
  !   do i1=-1,1,2
  !   do i2=-1,1,2
  !   res1l = res1l + ampwwNLO1l(i1,i2,-1,-1)*conjg(ampwwNLO1l(i1,i2,-1,-1))
  !   !-- interference 1loop * 2loop
  !   res2l = res2l + ampwwNLO1l(i1,i2,-1,-1)*conjg(ampwwNLO2l(i1,i2,-1,-1)) & 
  !        + ampwwNLO2l(i1,i2,-1,-1)*conjg(ampwwNLO1l(i1,i2,-1,-1))
  !   !print *, 'ampwwNLO1l(',i1,',',i2,',',i3,',',i4,') = ',ampwwNLO1l(i1,i2,i3,i4)
  !   !print *, 'ampwwNLO2l(',i1,',',i2,',',i3,',',i4,') = ',ampwwNLO2l(i1,i2,i3,i4)
  !   enddo
  !   enddo
  !   ! normalization
  !   s12 = sprod(1,2)
  !   res1l = res1l*avegg*colf

  !   ! the finite piece in Catani scheme (selected in wrapper_c.cpp) 
  !   res2l = res2l*avegg*colf 

  !   ! the double poles of the 2loopX1loop can be reconstructed knowing
  !   ! that the 2loop amplitude has a double pole = -CA*1loop -> 
  !   ! note that a factor 2 comes in the interference with the 1 loop!
  !   poles(-2) = -2*CA*res1l 
  !   poles(-1) = (-two * b0 - two * CA * log(musq/s12))*res1l

  !   return

  ! end subroutine resNLOV_WW

end module mod_ampNLO_ggww

