  module mod_DoubleTriangle
     use mod_types; use common_def
     use mod_consts_dp
    
     implicit none
    
     private
     public :: amp_doubletriangle
  contains

!  Amplitudes for double triangles(massless and massive) 

! 1 ccccccccc
!            |\       3
!            | \      |
!            |  \ cccc|
!            |  /     |
!            | /      4
!            |/
!            c
!            c
!            |\       5
!            | \      |
!            |  \ cccc|
!            |  /     |
!            | /      6
!            |/
! 2 ccccccccc

! -- R.R., Nov. 2015

     subroutine amp_doubletriangle(k1,k2,k3,k4,k5,k6,za,zb,sprod,&
          amp_bb,amp_bt,amp_tb,amp_tt)
! -- returns helicity amplitudes for double triangles
! -- indices: h1,h2,h3,h5
! -- the quarks in the two triangle loops can be bottom-bottom (amp_bb), bottom-top (amp_bt), top-bottom (amp_tb), or top-top (amp_tt)
        implicit none
        integer, intent(in)     :: k1,k2,k3,k4,k5,k6
        complex(dp15), intent(in) :: za(6,6),zb(6,6)
        real(dp15), intent(in)    :: sprod(6,6)
        complex(dp15),intent(out) :: amp_bb(-1:1,-1:1,-1:1,-1:1),amp_bt(-1:1,-1:1,-1:1,-1:1)
        complex(dp15),intent(out) :: amp_tb(-1:1,-1:1,-1:1,-1:1),amp_tt(-1:1,-1:1,-1:1,-1:1)
        complex(dp15)             :: T(-1:1,-1:1,-1:1,-1:1),T_old(-1:1,-1:1),F1m1(0:nheft_max),F1m2(0:nheft_max)
        real(dp15)                :: qsq,s34,s56,twok1Dq,twok2Dq,k1Dq,k2Dq,trimass
        complex(dp15)             :: qcdl,qlI2,qlI3,triloop1_mass,triloop1_mless,triloop2_mass,triloop2_mless
        complex(dp15)             :: triloop1_check(2),triloop2_check(2)
        integer                   :: i,j
        logical                 :: QCDlooptest

        QCDlooptest=.false.

        ! -- get helicity structures
        call hel_doubletriangle_LL(k1,k2,k3,k4,k5,k6,za,zb,sprod,T(:,:,-1,-1))   ! LL 
        call hel_doubletriangle_LL(k1,k2,k4,k3,k5,k6,za,zb,sprod,T(:,:,+1,-1))   ! RL
        call hel_doubletriangle_LL(k1,k2,k3,k4,k6,k5,za,zb,sprod,T(:,:,-1,+1))   ! LR
        call hel_doubletriangle_LL(k1,k2,k4,k3,k6,k5,za,zb,sprod,T(:,:,+1,+1))   ! RR


        qsq=sprod(k1,k3)+sprod(k1,k4)+sprod(k3,k4)
        s34=sprod(k3,k4)
        twok1Dq = (s34-qsq)                         ! for first triangle, k1sq=qsq,k3sq=s34
        k1Dq = twok1Dq/two
        s56=sprod(k5,k6)
        twok2Dq =s56-qsq                            ! for second triangle, k1sq=qsq, k3sq=s56
        k2Dq = twok2Dq/two


        !massless triangles
        triloop1_mless =F1mless_new(0d0,qsq,s34)
        triloop2_mless =F1mless_new(0d0,qsq,s56)

        trimass = expmass
        ! first triangle    F_1(q,k1,m1)  with q+k1+k34=0
        triloop1_mass=1.0_dp/twok1Dq* ( &
             2.0d0 + 4d0*trimass**2*qlI3(qsq,0.0_dp,s34,trimass**2,trimass**2,trimass**2,musq,0)+&
             (2d0 + qsq/k1Dq)*( qlI2(s34,trimass**2,trimass**2,musq,0) - qlI2(qsq,trimass**2,trimass**2,musq,0) ) &
             )

        ! second triangle    F_1(q,k2,m2)  with -q+k2+k56=0
        triloop2_mass=1.0_dp/twok2Dq* ( &
             2.0d0 + 4d0*trimass**2*qlI3(qsq,0.0_dp,s56,trimass**2,trimass**2,trimass**2,musq,0)+&
             (2d0 + qsq/k2Dq)*( qlI2(s56,trimass**2,trimass**2,musq,0) - qlI2(qsq,trimass**2,trimass**2,musq,0) ) &
             )

        if (qcdlooptest) then    ! check expression above against QCDlop for massless FF, and against 1/mt expansion for massive FF

           ! QCDLoop values for massless triangles
           trimass = 0d0
           triloop1_check(1)=1.0_dp/twok1Dq* ( &
             2.0d0 + 4d0*trimass**2*qlI3(qsq,0.0_dp,s34,trimass**2,trimass**2,trimass**2,musq,0)+&
             (2d0 + qsq/k1Dq)*( qlI2(s34,trimass**2,trimass**2,musq,0) - qlI2(qsq,trimass**2,trimass**2,musq,0) ) &
             )

           triloop2_check(1) = 1.0_dp/twok2Dq* ( &
                2.0d0 + 4d0*trimass**2*qlI3(qsq,0.0_dp,s56,trimass**2,trimass**2,trimass**2,musq,0)+&
                (2d0 + qsq/k2Dq)*( qlI2(s56,trimass**2,trimass**2,musq,0) - qlI2(qsq,trimass**2,trimass**2,musq,0) ) &
                )

           ! Mass expansion for massive triangles
           call F2mass(0.0_dp,qsq,s34,expmass,F1m1)
           triloop1_check(2)=sum(F1m1)
           call F2mass(0.0_dp,qsq,s56,expmass,F1m2)
           triloop2_check(2)=sum(F1m2)


           print *, "testing massless triangles vs QCDLoop"
           print *, "ratio",triloop1_mless/triloop1_check(1)
           print *, "ratio",triloop2_mless/triloop2_check(1)
           
           print *, "testing massive triangles vs independent expansion"
           print *, "ratio",triloop1_mass/triloop1_check(2)
           print *, "ratio",triloop2_mass/triloop2_check(2)
        endif


        amp_bb =  T * triloop1_mless*triloop2_mless     ! two massless loops
        amp_bt =  T * triloop1_mless*triloop2_mass     ! massless*massive loops
        amp_tb =  T * triloop1_mass*triloop2_mless    ! massive*massless loops
        amp_tt =  T * triloop1_mass*triloop2_mass     ! two massive loops


        


        
      end subroutine amp_doubletriangle


     subroutine hel_doubletriangle_LL(k1,k2,k3,k4,k5,k6,za,zb,sprod,&
         helstr)
! -- returns helicity structures 
! -- indices are helicities of gluons, both leptons are LH
! -- other lepton helicities can be obtained by usual swaps
        implicit none
        integer, intent(in)     :: k1,k2,k3,k4,k5,k6
        complex(dp15), intent(in) :: za(6,6),zb(6,6)
        real(dp15), intent(in)    :: sprod(6,6)
        complex(dp15),intent(out) :: helstr(-1:1,-1:1)
        complex(dp15)             :: iza,izb,zab,zab2
        complex(dp15)             :: LT(-1:1,-1:1)

        helstr(-1,-1) = + za(k1,k2)*za(k1,k3)**2*za(k2,k5)*zb(k1,k3)*zb(k4,k6) &
             + 1.0_dp/2.0_dp*za(k1,k2)*za(k1,k3)**2*za(k2,k5)*zb(k1,k6)*zb(k3,k4)   &
             + za(k1,k2)*za(k1,k3)*za(k1,k4)*za(k2,k5)*zb(k1,k4)*zb(k4,k6)  &
             + za(k1,k2)*za(k1,k3)*za(k2,k5)*za(k3,k4)*zb(k3,k4)*zb(k4,k6)  &
             - 1.0_dp/2.0_dp*za(k1,k3)**2*za(k2,k3)*za(k2,k5)*zb(k3,k4)*zb(k3,k6)   &
             - 1.0_dp/2.0_dp*za(k1,k3)**2*za(k2,k4)*za(k2,k5)*zb(k3,k4)*zb(k4,k6)   
        

        helstr(-1,+1) = + 1.0_dp/2.0_dp*za(k1,k3)**2*za(k1,k5)*zb(k1,k2)*zb(k2,k6)*zb(k3,k4) &
             - za(k1,k3)**2*za(k1,k5)*zb(k1,k3)*zb(k2,k4)*zb(k2,k6)     &
             - 1.0_dp/2.0_dp*za(k1,k3)**2*za(k3,k5)*zb(k2,k3)*zb(k2,k6)*zb(k3,k4)    &
             - 1.0_dp/2.0_dp*za(k1,k3)**2*za(k4,k5)*zb(k2,k4)*zb(k2,k6)*zb(k3,k4)  &
             - za(k1,k3)*za(k1,k4)*za(k1,k5)*zb(k1,k4)*zb(k2,k4)*zb(k2,k6) &
             - za(k1,k3)*za(k1,k5)*za(k3,k4)*zb(k2,k4)*zb(k2,k6)*zb(k3,k4)


        helstr(+1,-1) = - 1.0_dp/2.0_dp*za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(k1,k4)**2*zb(k1,k6) &
             - za(k1,k3)*za(k2,k3)*za(k2,k5)*zb(k1,k3)*zb(k1,k4)*zb(k1,k6)    &
             - za(k1,k4)*za(k2,k3)*za(k2,k5)*zb(k1,k4)**2*zb(k1,k6)           &
             + 1.0_dp/2.0_dp*za(k2,k3)*za(k2,k5)*za(k3,k4)*zb(k1,k4)**2*zb(k3,k6)     &
             - za(k2,k3)*za(k2,k5)*za(k3,k4)*zb(k1,k4)*zb(k1,k6)*zb(k3,k4)    &
             + 1.0_dp/2.0_dp*za(k2,k4)*za(k2,k5)*za(k3,k4)*zb(k1,k4)**2*zb(k4,k6)


        helstr(+1,+1) = + za(k1,k3)*za(k3,k5)*zb(k1,k2)*zb(k1,k3)*zb(k1,k4)*zb(k2,k6) &
             + za(k1,k4)*za(k3,k5)*zb(k1,k2)*zb(k1,k4)**2*zb(k2,k6)           &
             - 1.0_dp/2.0_dp*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)**2*zb(k2,k6)     &
             + za(k3,k4)*za(k3,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6)*zb(k3,k4)    &
             + 1.0_dp/2.0_dp*za(k3,k4)*za(k3,k5)*zb(k1,k4)**2*zb(k2,k3)*zb(k2,k6)     &
             + 1.0_dp/2.0_dp*za(k3,k4)*za(k4,k5)*zb(k1,k4)**2*zb(k2,k4)*zb(k2,k6)

      helstr = helstr/sprod(k3,k4)/sprod(k5,k6)

    end subroutine hel_doubletriangle_LL
    


    function F1mless_new(k1sq,k2sq,k3sq)
      use mod_types; use mod_consts_dp
      !-- Form factor for massless triangle
      !-- Taken from Eq. 2.14 of Hagiwara, Kuruma, Yamada, Nucl. Phys. B358 (1991),80-96
      !-- excl. factor of 1/pi^2
      !-- include factor of 4 relative to CEZ in integrand, cf. HKY Eqs (2.9)-(2.12) vs CEZ Eq. (A.6)
      implicit none
       real(dp15),intent(in)     :: k1sq,k2sq,k3sq
       complex(dp15)             :: F1mless_new
       real(dp15)                :: k1Dk2

       k1Dk2=(k3sq-k1sq-k2sq)/2.0_dp
       F1mless_new=4.0_dp*(1.0_dp/2.0_dp/k1Dk2 * (1.0_dp/2.0_dp + k3sq/4.0_dp/k1Dk2*(log(-k2sq/k3sq)+ci*pi)))

     end function F1mless_new
     


   end module mod_DoubleTriangle
