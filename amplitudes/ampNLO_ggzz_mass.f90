module mod_ampNLO_ggzz_mass
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_ampLO_ggzz_mass
  use mod_DoubleTriangle
  
  implicit none
  private
  
  logical :: polecheck=.false.


!  public :: getewampsNLO_ggzz_heft_eps
  public :: getewampsNLO_ggzz_heft, getewampsNLO_ggzz_rewgt
  
contains
  
  !-- g(1) g(2) -> [Z->l(3) lb(4)] [Z->l(5) lb(6)]
  !-- returns the  NLO amplitudes calling SUBTRACTED 2-loop amplitudes
  !-- first index is the expansion index (0=mt->infty limit)
  !-- other index are h1,h2,h3,h5
  subroutine getewampsNLO_ggzz_heft(za,zb,sprod,mass,ewampNLO_ggzz_heft)
    complex(dp15), intent(in)  :: za(6,6),zb(6,6)
    real(dp15), intent(in)     :: sprod(6,6),mass
    complex(dp15), intent(out) :: ewampNLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    integer                   :: i1,i2,h1,h2,h3,h5,nheft
    complex(dp15)              :: A_nlo_ax(20,0:nheft_max),A_nlo_vec(20,0:nheft_max),T(20,-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: A_nlo_ax_sub(20,0:nheft_max,-2:0),A_nlo_vec_sub(20,0:nheft_max,-2:0)
    complex(dp15)              :: ewampNLO_ggzz_ax(0:nheft_max,-1:1,-1:1,-1:1,-1:1),ewampNLO_ggzz_vec(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: ewampDT_bb(-1:1,-1:1,-1:1,-1:1),ewampDT_bt(-1:1,-1:1,-1:1,-1:1),ewampDT_tb(-1:1,-1:1,-1:1,-1:1),&
         ewampDT_tt(-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: ewampDT_bb3456(-1:1,-1:1,-1:1,-1:1),ewampDT_bt3456(-1:1,-1:1,-1:1,-1:1),ewampDT_tb3456(-1:1,-1:1,-1:1,-1:1),&
         ewampDT_tt3456(-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: ewampDT_bb5634(-1:1,-1:1,-1:1,-1:1),ewampDT_bt5634(-1:1,-1:1,-1:1,-1:1),ewampDT_tb5634(-1:1,-1:1,-1:1,-1:1),&
         ewampDT_tt5634(-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: prop34,prop56
    real(dp15)                 :: twosinwcosw,cl1(-1:1),cl2(-1:1),cvec(1:2),cax(1:2),q1,q2,Qu,Qd
    complex(dp)                :: A_lo_ax(20,0:5,-2:2),A_lo_vec(20,0:5,-2:2)

    ewampNLO_ggzz_heft = czero
    ewampNLO_ggzz_ax = czero
    ewampNLO_ggzz_vec = czero
    ewampDT_bb3456 = czero
    ewampDT_bt3456 = czero
    ewampDT_tb3456 = czero
    ewampDT_tt3456 = czero
    ewampDT_bb5634 = czero
    ewampDT_bt5634 = czero
    ewampDT_tb5634 = czero
    ewampDT_tt5634 = czero
    ewampDT_bb = czero
    ewampDT_bt = czero
    ewampDT_tb = czero
    ewampDT_tt = czero

    if (polecheck) then
       print *, "*********************************************************************"
       print *, "* WARNING: substituting double pole values for two-loop amplitudes! *"
       print *, "*********************************************************************"
    endif

! subtracted 2-loop form factors -- these return finite parts only!
    call ffaxNLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_nlo_ax)
    call ffvecNLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_nlo_vec)

    if (polecheck) then
       call ffaxLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_lo_ax)
       call ffvecLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_lo_vec)
       A_nlo_vec = -CA * A_lo_vec(:,:,0)
       A_nlo_ax = -CA * A_lo_ax(:,:,0)
    endif

    

   !-- populate the helicity array
    call heli_zz_heft(1,2,3,4,5,6,za,zb,sprod,T)

    !-- build the amplitudes
    do nheft=0,nheft_max
    do h1=-1,1,2
    do h2=-1,1,2
    do h3=-1,1,2
    do h5=-1,1,2
    do i1=1,20

       ewampNLO_ggzz_ax(nheft,h1,h2,h3,h5) = ewampNLO_ggzz_ax(nheft,h1,h2,h3,h5) + &
            A_nlo_ax(i1,nheft)*T(i1,h1,h2,h3,h5)
       ewampNLO_ggzz_vec(nheft,h1,h2,h3,h5) = ewampNLO_ggzz_vec(nheft,h1,h2,h3,h5) + &
            A_nlo_vec(i1,nheft)*T(i1,h1,h2,h3,h5)

    enddo
    enddo
    enddo
    enddo
    enddo
    enddo


    ! double triangles
    if (massloops .and. massloops) then
       call amp_doubletriangle(1,2,3,4,5,6,za,zb,sprod,ewampDT_bb3456,ewampDT_bt3456,ewampDT_tb3456,ewampDT_tt3456)    !(1,34), (2,56)
       call amp_doubletriangle(1,2,5,6,3,4,za,zb,sprod,ewampDT_bb5634,ewampDT_bt5634,ewampDT_tb5634,ewampDT_tt5634)    !(1,56), (2,34)
    else                              ! remove double triangles unless we use the full 3rd generation                       
       ewampDT_bb3456 = zero
       ewampDT_tb3456 = zero
       ewampDT_bt3456 = zero
       ewampDT_tt3456 = zero
       ewampDT_bb5634 = zero
       ewampDT_tb5634 = zero
       ewampDT_bt5634 = zero
       ewampDT_tt5634 = zero
    endif

! now include couplings, as done in ampNLO_ggzz
    
 !--- propagator factors
    prop34=sprod(3,4)/(sprod(3,4)-mzsq + ci*mz*gaz)
    prop56=sprod(5,6)/(sprod(5,6)-mzsq + ci*mz*gaz)

    !-- couplings, see src/Need/couplz.f and coupling.f
    !-- xw = sin(tw)^2, sin2w = 2*sin(tw)*cos(tw)

    twosinwcosw=two*sw*cw
    
!-orig    cl1(-1)=Lel/(twosinwcosw)
!-orig    cl1(+1)=Rel/(twosinwcosw)
!-orig    cl2(-1)=Lel/(twosinwcosw)
!-orig    cl2(+1)=Rel/(twosinwcosw)
!-orig
!-orig    cvec(up)=Vup/(twosinwcosw)
!-orig    cvec(dn)=Vdn/(twosinwcosw)
!-orig    cax(up)=Aup/(twosinwcosw)
!-orig    cax(dn)=Adn/(twosinwcosw)

    cl1(-1)=Lel/(2.0_dp*cw)
    cl1(+1)=Rel/(2.0_dp*cw)
    cl2(-1)=Lel/(2.0_dp*cw)
    cl2(+1)=Rel/(2.0_dp*cw)

    cvec(up)=Vup/(2.0_dp*cw)
    cvec(dn)=Vdn/(2.0_dp*cw)
    cax(up)=Aup/(2.0_dp*cw)
    cax(dn)=Adn/(2.0_dp*cw)

    q1=-one
    q2=-one
    Qu=Qup
    Qd=Qdn

    !--- dress vector and axial amplitudes with appropriate couplings
    do h1=-1,1,2
    do h2=-1,1,2
    do h3=-1,1,2
    do h5=-1,1,2
       ewampNLO_ggzz_heft(:,h1,h2,h3,h5)=ewampNLO_ggzz_vec(:,h1,h2,h3,h5)*( &
            (Qu*q1*sinW2+cvec(up)*cl1(h3)*prop34)*(Qu*q2*sinW2+cvec(up)*cl2(h5)*prop56)) &
            +ewampNLO_ggzz_ax(:,h1,h2,h3,h5)*(cax(up)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56)

       ewampDT_bb3456(h1,h2,h3,h5)=ewampDT_bb3456(h1,h2,h3,h5)*(cax(dn)*cl1(h3)*prop34)*(cax(dn)*cl2(h5)*prop56)
       ewampDT_bt3456(h1,h2,h3,h5)=ewampDT_bt3456(h1,h2,h3,h5)*(cax(dn)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56)
       ewampDT_tb3456(h1,h2,h3,h5)=ewampDT_tb3456(h1,h2,h3,h5)*(cax(up)*cl1(h3)*prop34)*(cax(dn)*cl2(h5)*prop56)
       ewampDT_tt3456(h1,h2,h3,h5)=ewampDT_tt3456(h1,h2,h3,h5)*(cax(up)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56)

       ewampDT_bb5634(h1,h2,h3,h5)=ewampDT_bb5634(h1,h2,h3,h5)*(cax(dn)*cl1(h5)*prop34)*(cax(dn)*cl2(h3)*prop56)
       ewampDT_bt5634(h1,h2,h3,h5)=ewampDT_bt5634(h1,h2,h3,h5)*(cax(dn)*cl1(h5)*prop34)*(cax(up)*cl2(h3)*prop56)
       ewampDT_tb5634(h1,h2,h3,h5)=ewampDT_tb5634(h1,h2,h3,h5)*(cax(up)*cl1(h5)*prop34)*(cax(dn)*cl2(h3)*prop56)
       ewampDT_tt5634(h1,h2,h3,h5)=ewampDT_tt5634(h1,h2,h3,h5)*(cax(up)*cl1(h5)*prop34)*(cax(up)*cl2(h3)*prop56)


    enddo
    enddo
    enddo
    enddo

    ewampDT_bb= ewampDT_bb3456  + ewampDT_bb5634
    ewampDT_bt= ewampDT_bt3456  + ewampDT_bt5634
    ewampDT_tb= ewampDT_tb3456  + ewampDT_tb5634
    ewampDT_tt= ewampDT_tt3456  + ewampDT_tt5634

! -- change normalizations from Tr[TaTb]=1/2*delta_[ab] (used by Matthew)
! -- to Tr[TaTb]=delta_[ab] (used in code)
! -- should include factor of sqrt(2) for each Ta
!    ewampNLO_ggzz_heft=2d0*ewampNLO_ggzz_heft


    !CHECK: normalization to agree with MCFM -- ?????
    ewampNLO_ggzz_heft=4d0*ewampNLO_ggzz_heft
! -- now add in phase to agree with MCFM
    ewampNLO_ggzz_heft=ci*ewampNLO_ggzz_heft

! -- double triangle phases as for the others
    ewampDT_bb = ci*ewampDT_bb
    ewampDT_bt = ci*ewampDT_bt
    ewampDT_tb = ci*ewampDT_tb
    ewampDT_tt = ci*ewampDT_tt

! add double triangles
    if (.not. polecheck) then
       ewampNLO_ggzz_heft(0,:,:,:,:) = ewampNLO_ggzz_heft(0,:,:,:,:) + ewampDT_bb(:,:,:,:) + ewampDT_bt(:,:,:,:) + ewampDT_tb(:,:,:,:) + ewampDT_tt(:,:,:,:)
    endif

! --overall  normalization
!-orig    ewampNLO_ggzz_heft = ewampNLO_ggZZ_heft * two * (gwsq * sinW2)**2
    ewampNLO_ggzz_heft = ewampNLO_ggZZ_heft * gwsq**2



    return
  end subroutine getewampsNLO_ggzz_heft



  !-- g(1) g(2) -> [Z->l(3) lb(4)] [Z->l(5) lb(6)]
  !-- returns the massive two-loop amplitudes obtained via reweighting, plus the double-triangle contribution
  subroutine getewampsNLO_ggzz_rewgt(za,zb,sprod,ewampNLO_ggzz_mless,ewampLO_ggzz_mless,ewampLO_ggzz_mass,ewampNLO_ggzz_rewgt)
    complex(dp15), intent(in)  :: za(6,6),zb(6,6)
    real(dp15), intent(in)     :: sprod(6,6)
    complex(dp15), intent(in)  :: ewampLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1),ewampLO_ggzz_mass(-1:1,-1:1,-1:1,-1:1),ewampNLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
    complex(dp15), intent(out) :: ewampNLO_ggzz_rewgt(-1:1,-1:1,-1:1,-1:1)
    integer                   :: i1,i2,h1,h2,h3,h5,nheft
    complex(dp15)              :: A_nlo_ax(20,0:nheft_max),A_nlo_vec(20,0:nheft_max),T(20,-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: A_nlo_ax_sub(20,0:nheft_max,-2:0),A_nlo_vec_sub(20,0:nheft_max,-2:0)
    complex(dp15)              :: ewampNLO_ggzz_ax(0:nheft_max,-1:1,-1:1,-1:1,-1:1),ewampNLO_ggzz_vec(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: ewampDT_bb(-1:1,-1:1,-1:1,-1:1),ewampDT_bt(-1:1,-1:1,-1:1,-1:1),ewampDT_tb(-1:1,-1:1,-1:1,-1:1),&
         ewampDT_tt(-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: ewampDT_bb3456(-1:1,-1:1,-1:1,-1:1),ewampDT_bt3456(-1:1,-1:1,-1:1,-1:1),ewampDT_tb3456(-1:1,-1:1,-1:1,-1:1),&
         ewampDT_tt3456(-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: ewampDT_bb5634(-1:1,-1:1,-1:1,-1:1),ewampDT_bt5634(-1:1,-1:1,-1:1,-1:1),ewampDT_tb5634(-1:1,-1:1,-1:1,-1:1),&
         ewampDT_tt5634(-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: prop34,prop56
    real(dp15)                 :: twosinwcosw,cl1(-1:1),cl2(-1:1),cvec(1:2),cax(1:2),q1,q2,Qu,Qd
    complex(dp)                :: A_lo_ax(20,0:5,-2:2),A_lo_vec(20,0:5,-2:2)


    ewampNLO_ggzz_rewgt = czero
    ewampDT_bb3456 = czero
    ewampDT_bt3456 = czero
    ewampDT_tb3456 = czero
    ewampDT_tt3456 = czero
    ewampDT_bb5634 = czero
    ewampDT_bt5634 = czero
    ewampDT_tb5634 = czero
    ewampDT_tt5634 = czero
    ewampDT_bb = czero
    ewampDT_bt = czero
    ewampDT_tb = czero
    ewampDT_tt = czero

    ! double triangles
    if (massloops .and. massloops) then
       call amp_doubletriangle(1,2,3,4,5,6,za,zb,sprod,ewampDT_bb3456,ewampDT_bt3456,ewampDT_tb3456,ewampDT_tt3456)    !(1,34), (2,56)
       call amp_doubletriangle(1,2,5,6,3,4,za,zb,sprod,ewampDT_bb5634,ewampDT_bt5634,ewampDT_tb5634,ewampDT_tt5634)    !(1,56), (2,34)
    else                              ! remove double triangles unless we use the full 3rd generation                       
       ewampDT_bb3456 = zero
       ewampDT_tb3456 = zero
       ewampDT_bt3456 = zero
       ewampDT_tt3456 = zero
       ewampDT_bb5634 = zero
       ewampDT_tb5634 = zero
       ewampDT_bt5634 = zero
       ewampDT_tt5634 = zero
    endif

! now include couplings, as done in ampNLO_ggzz
    
 !--- propagator factors
    prop34=sprod(3,4)/(sprod(3,4)-mzsq + ci*mz*gaz)
    prop56=sprod(5,6)/(sprod(5,6)-mzsq + ci*mz*gaz)

    twosinwcosw=two*sw*cw
    
    cl1(-1)=Lel/(2.0_dp*cw)
    cl1(+1)=Rel/(2.0_dp*cw)
    cl2(-1)=Lel/(2.0_dp*cw)
    cl2(+1)=Rel/(2.0_dp*cw)

    cvec(up)=Vup/(2.0_dp*cw)
    cvec(dn)=Vdn/(2.0_dp*cw)
    cax(up)=Aup/(2.0_dp*cw)
    cax(dn)=Adn/(2.0_dp*cw)

    q1=-one
    q2=-one
    Qu=Qup
    Qd=Qdn

    !--- dress double triangles vector and axial amplitudes with appropriate couplings,
    !--- and compute reweighted massive amplitude for each helicity configuration
    do h1=-1,1,2
    do h2=-1,1,2
    do h3=-1,1,2
    do h5=-1,1,2

       ewampDT_bb3456(h1,h2,h3,h5)=ewampDT_bb3456(h1,h2,h3,h5)*(cax(dn)*cl1(h3)*prop34)*(cax(dn)*cl2(h5)*prop56)
       ewampDT_bt3456(h1,h2,h3,h5)=ewampDT_bt3456(h1,h2,h3,h5)*(cax(dn)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56)
       ewampDT_tb3456(h1,h2,h3,h5)=ewampDT_tb3456(h1,h2,h3,h5)*(cax(up)*cl1(h3)*prop34)*(cax(dn)*cl2(h5)*prop56)
       ewampDT_tt3456(h1,h2,h3,h5)=ewampDT_tt3456(h1,h2,h3,h5)*(cax(up)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56)

       ewampDT_bb5634(h1,h2,h3,h5)=ewampDT_bb5634(h1,h2,h3,h5)*(cax(dn)*cl1(h5)*prop34)*(cax(dn)*cl2(h3)*prop56)
       ewampDT_bt5634(h1,h2,h3,h5)=ewampDT_bt5634(h1,h2,h3,h5)*(cax(dn)*cl1(h5)*prop34)*(cax(up)*cl2(h3)*prop56)
       ewampDT_tb5634(h1,h2,h3,h5)=ewampDT_tb5634(h1,h2,h3,h5)*(cax(up)*cl1(h5)*prop34)*(cax(dn)*cl2(h3)*prop56)
       ewampDT_tt5634(h1,h2,h3,h5)=ewampDT_tt5634(h1,h2,h3,h5)*(cax(up)*cl1(h5)*prop34)*(cax(up)*cl2(h3)*prop56)

       ! reweight
       if (abs(ewampLO_ggzz_mless(h1,h2,h3,h5)) .gt. zero) then
          ewampNLO_ggzz_rewgt(h1,h2,h3,h5) = ewampNLO_ggzz_mless(h1,h2,h3,h5) * ewampLO_ggzz_mass(h1,h2,h3,h5)/ewampLO_ggzz_mless(h1,h2,h3,h5)
       endif

    enddo
    enddo
    enddo
    enddo

    ewampDT_bb= ewampDT_bb3456  + ewampDT_bb5634
    ewampDT_bt= ewampDT_bt3456  + ewampDT_bt5634
    ewampDT_tb= ewampDT_tb3456  + ewampDT_tb5634
    ewampDT_tt= ewampDT_tt3456  + ewampDT_tt5634

! -- double triangle phases as for the others, and  factors of gwsq
    ewampDT_bb = ci*ewampDT_bb*gwsq**2
    ewampDT_bt = ci*ewampDT_bt*gwsq**2
    ewampDT_tb = ci*ewampDT_tb*gwsq**2
    ewampDT_tt = ci*ewampDT_tt*gwsq**2

! sum reweighted + double triangles
    ewampNLO_ggzz_rewgt = ewampNLO_ggzz_rewgt + ewampDT_bb + ewampDT_bt+ ewampDT_tb+ ewampDT_tt


    


    return
  end subroutine getewampsNLO_ggzz_rewgt


!! -- removed the pole check for powheg implementation -- RR July 2019 -- !!

!!--!!  !-- g(1) g(2) -> [Z->l(3) lb(4)] [Z->l(5) lb(6)]
!!--!!  !-- returns the  NLO amplitudes, 
!!--!!  !-- calling UNSUBTRACTED 2-loop amplitudes and then implementing subtraction numerically
!!--!!  !-- first index is the expansion index (0=mt->infty limit)
!!--!!  !-- other index are h1,h2,h3,h5
!!--!!  !-- final index is eps^j -- should be zero for poles due to subtr
!!--!!  subroutine getewampsNLO_ggzz_heft_eps(za,zb,sprod,mass,ewampNLO_ggzz_eps)
!!--!!    complex(dp), intent(in)   :: za(6,6),zb(6,6)
!!--!!    real(dp), intent(in)      :: sprod(6,6),mass
!!--!!    complex(dp), intent(out)  :: ewampNLO_ggzz_eps(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2)
!!--!!    integer                   :: i1,i2,h1,h2,h3,h5,nheft,neps
!!--!!    complex(dp)               :: A_nlo_ax_unsub(20,0:5,-2:2),A_nlo_vec_unsub(20,0:5,-2:2),A_lo_ax(20,0:5,-2:2),A_lo_vec(20,0:5,-2:2)
!!--!!    complex(dp)               :: T(20,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp)               :: prop34,prop56
!!--!!    real(dp)                  :: twosinwcosw,cl1(-1:1),cl2(-1:1),cvec(1:2),cax(1:2),q1,q2,Qu,Qd
!!--!!    complex(dp)               :: ewampLO_ggzz_ax(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2),ewampLO_ggzz_vec(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2)
!!--!!    complex(dp)               :: ewampNLO_ggzz_ax_unsub(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2),ewampNLO_ggzz_vec_unsub(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2)
!!--!!    complex(dp)               :: ewampNLO_ggzz_ax_sub(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2),ewampNLO_ggzz_vec_sub(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2)
!!--!!    complex(dp)               :: L,Subtr_coeff(-2:0)
!!--!!    real(dp)                  :: delta_qT
!!--!!    
!!--!!    ewampNLO_ggzz_eps = czero
!!--!!    ewampLO_ggzz_ax = czero
!!--!!    ewampLO_ggzz_vec = czero
!!--!!    ewampNLO_ggzz_ax_unsub = czero
!!--!!    ewampNLO_ggzz_vec_unsub = czero
!!--!!    ewampNLO_ggzz_ax_sub = czero
!!--!!    ewampNLO_ggzz_vec_sub = czero
!!--!!
!!--!!    !-- populate the form factors seperately for vector and axial, LO and NLO
!!--!!    call ffaxLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_lo_ax)
!!--!!    call ffvecLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_lo_vec)
!!--!!
!!--!!    call ffaxNLO_ggzz_heft_unsub(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_nlo_ax_unsub)
!!--!!    call ffvecNLO_ggzz_heft_unsub(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_nlo_vec_unsub)
!!--!!
!!--!!    !-- populate the helicity array
!!--!!    call heli_zz_heft(1,2,3,4,5,6,za,zb,sprod,T)
!!--!!
!!--!!    !-- build the amplitudes
!!--!!    do neps=-2,2
!!--!!    do nheft=0,nheft_max
!!--!!    do h1=-1,1,2
!!--!!    do h2=-1,1,2
!!--!!    do h3=-1,1,2
!!--!!    do h5=-1,1,2
!!--!!    do i1=1,20
!!--!!
!!--!!! 1-loop
!!--!!       ewampLO_ggzz_ax(nheft,h1,h2,h3,h5,neps) = ewampLO_ggzz_ax(nheft,h1,h2,h3,h5,neps) + &
!!--!!            A_lo_ax(i1,nheft,neps)*T(i1,h1,h2,h3,h5)
!!--!!       ewampLO_ggzz_vec(nheft,h1,h2,h3,h5,neps) = ewampLO_ggzz_vec(nheft,h1,h2,h3,h5,neps) + &
!!--!!            A_lo_vec(i1,nheft,neps)*T(i1,h1,h2,h3,h5)
!!--!!
!!--!!! 2-loop
!!--!!       ewampNLO_ggzz_ax_unsub(nheft,h1,h2,h3,h5,neps) = ewampNLO_ggzz_ax_unsub(nheft,h1,h2,h3,h5,neps) + &
!!--!!            A_nlo_ax_unsub(i1,nheft,neps)*T(i1,h1,h2,h3,h5)
!!--!!       ewampNLO_ggzz_vec_unsub(nheft,h1,h2,h3,h5,neps) = ewampNLO_ggzz_vec_unsub(nheft,h1,h2,h3,h5,neps) + &
!!--!!            A_nlo_vec_unsub(i1,nheft,neps)*T(i1,h1,h2,h3,h5)
!!--!!
!!--!!
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!
!!--!!
!!--!!! -- change normalizations from Tr[TaTb]=1/2*delta_[ab] (used by Matthew)
!!--!!! -- to Tr[TaTb]=delta_[ab] (used in code)
!!--!!! -- should include factor of sqrt(2) for each Ta
!!--!!    ewampLO_ggzz_ax=4d0*ewampLO_ggzz_ax
!!--!!    ewampLO_ggzz_vec=4d0*ewampLO_ggzz_vec
!!--!!    ewampNLO_ggzz_ax_unsub=8d0*ewampNLO_ggzz_ax_unsub
!!--!!    ewampNLO_ggzz_vec_unsub=8d0*ewampNLO_ggzz_vec_unsub
!!--!!! -- now add in phase to agree with MCFM
!!--!!    ewampLO_ggzz_ax=ci*ewampLO_ggzz_ax
!!--!!    ewampLO_ggzz_vec=ci*ewampLO_ggzz_vec
!!--!!    ewampNLO_ggzz_ax_unsub=ci*ewampNLO_ggzz_ax_unsub
!!--!!    ewampNLO_ggzz_vec_unsub=ci*ewampNLO_ggzz_vec_unsub
!!--!!
!!--!!
!!--!!
!!--!!! subtraction coefficients  -- cf. 1503.08835, Eq. 4.5-4.10 and R.R. notes
!!--!!    L =dcmplx(1.0_dp,0.0_dp)*musq/sprod(1,2)
!!--!!    L = log(L)
!!--!!! these are into qT-scheme
!!--!!    delta_qT=0.0_dp
!!--!!    Subtr_Coeff(-2) = -CA
!!--!!    Subtr_Coeff(-1) = -CA*(L+ci*pi)-b0
!!--!!    Subtr_Coeff(0)  = -CA*(delta_qT+ci*pi*L+half*(L**2-pisq/6.0_dp))-b0*L
!!--!!!    ! Uncomment for catani check
!!--!!!    Subtr_Coeff(-2) = 0._dp
!!--!!!    Subtr_Coeff(-1) = 0._dp
!!--!!!    Subtr_Coeff(0)  = 0._dp
!!--!!
!!--!!! now subtract:
!!--!!    ewampNLO_ggzz_ax_sub(:,:,:,:,:,-2) =ewampNLO_ggzz_ax_unsub(:,:,:,:,:,-2)  &
!!--!!         -ewampLO_ggzz_ax(:,:,:,:,:,0)*Subtr_Coeff(-2)
!!--!!    ewampNLO_ggzz_vec_sub(:,:,:,:,:,-2)=ewampNLO_ggzz_vec_unsub(:,:,:,:,:,-2) &
!!--!!         -ewampLO_ggzz_vec(:,:,:,:,:,0)*Subtr_Coeff(-2)
!!--!!
!!--!!    ewampNLO_ggzz_ax_sub(:,:,:,:,:,-1) =ewampNLO_ggzz_ax_unsub(:,:,:,:,:,-1)  &
!!--!!         -ewampLO_ggzz_ax(:,:,:,:,:,0)*Subtr_Coeff(-1) &
!!--!!         -ewampLO_ggzz_ax(:,:,:,:,:,1)*Subtr_Coeff(-2) 
!!--!!    ewampNLO_ggzz_vec_sub(:,:,:,:,:,-1)=ewampNLO_ggzz_vec_unsub(:,:,:,:,:,-1) &
!!--!!         -ewampLO_ggzz_vec(:,:,:,:,:,0)*Subtr_Coeff(-1) &
!!--!!         -ewampLO_ggzz_vec(:,:,:,:,:,1)*Subtr_Coeff(-2)
!!--!!
!!--!!    ewampNLO_ggzz_ax_sub(:,:,:,:,:,0) =ewampNLO_ggzz_ax_unsub(:,:,:,:,:,0)  &
!!--!!         -ewampLO_ggzz_ax(:,:,:,:,:,0)*Subtr_Coeff(0) &
!!--!!         -ewampLO_ggzz_ax(:,:,:,:,:,1)*Subtr_Coeff(-1) &
!!--!!         -ewampLO_ggzz_ax(:,:,:,:,:,2)*Subtr_Coeff(-2) 
!!--!!    ewampNLO_ggzz_vec_sub(:,:,:,:,:,0)=ewampNLO_ggzz_vec_unsub(:,:,:,:,:,0) &
!!--!!         -ewampLO_ggzz_vec(:,:,:,:,:,0)*Subtr_Coeff(0) &
!!--!!         -ewampLO_ggzz_vec(:,:,:,:,:,1)*Subtr_Coeff(-1) &
!!--!!         -ewampLO_ggzz_vec(:,:,:,:,:,2)*Subtr_Coeff(-2)
!!--!!
!!--!!! now include couplings, as done in ampNLO_ggzz
!!--!!    
!!--!! !--- propagator factors
!!--!!    prop34=sprod(3,4)/(sprod(3,4)-mzsq + ci*mz*gaz)
!!--!!    prop56=sprod(5,6)/(sprod(5,6)-mzsq + ci*mz*gaz)
!!--!!
!!--!!    !-- couplings, see src/Need/couplz.f and coupling.f
!!--!!    !-- xw = sin(tw)^2, sin2w = 2*sin(tw)*cos(tw)
!!--!!
!!--!!    twosinwcosw=two*sw*cw
!!--!!    
!!--!!
!!--!!
!!--!!    cl1(-1)=Lel/(2.0_dp*cw)
!!--!!    cl1(+1)=Rel/(2.0_dp*cw)
!!--!!    cl2(-1)=Lel/(2.0_dp*cw)
!!--!!    cl2(+1)=Rel/(2.0_dp*cw)
!!--!!
!!--!!    cvec(up)=Vup/(2.0_dp*cw)
!!--!!    cvec(dn)=Vdn/(2.0_dp*cw)
!!--!!    cax(up)=Aup/(2.0_dp*cw)
!!--!!    cax(dn)=Adn/(2.0_dp*cw)
!!--!!
!!--!!    q1=-one
!!--!!    q2=-one
!!--!!    Qu=Qup
!!--!!    Qd=Qdn
!!--!!
!!--!!    !--- dress vector and axial amplitudes with appropriate couplings
!!--!!    do h1=-1,1,2
!!--!!    do h2=-1,1,2
!!--!!    do h3=-1,1,2
!!--!!    do h5=-1,1,2
!!--!!       ewampNLO_ggzz_eps(:,h1,h2,h3,h5,:)=ewampNLO_ggzz_vec_sub(:,h1,h2,h3,h5,:)*( &
!!--!!            (Qu*q1*sinW2+cvec(up)*cl1(h3)*prop34)*(Qu*q2*sinW2+cvec(up)*cl2(h5)*prop56)) &
!!--!!            +ewampNLO_ggzz_ax_sub(:,h1,h2,h3,h5,:)*(cax(up)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56)
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!
!!--!!
!!--!!! --overall  normalization
!!--!!!-orig    ewampNLO_ggzz_eps = ewampNLO_ggZZ_eps * two * (gwsq * sinW2)**2
!!--!!    ewampNLO_ggzz_eps = ewampNLO_ggZZ_eps * gwsq**2
!!--!!
!!--!!  end subroutine getewampsNLO_ggzz_heft_eps
!!--!!
!!--!!  subroutine pole_check(ampLO,ampNLO,s12)
!!--!!    implicit none
!!--!!    complex(dp),intent(in) :: ampLO(0:2),ampNLO(-2:0)
!!--!!    real(dp), intent(in)   :: s12
!!--!!    complex(dp)            :: doublepole,singlepole
!!--!!    logical                :: verbose,correctdp,correctsp
!!--!!    
!!--!!    verbose=.false.
!!--!!    doublepole = -CA*ampLO(0)
!!--!!    singlepole = -(b0+CA*(log(musq/s12) + ci*pi))*ampLO(0)-CA*ampLO(1)
!!--!!
!!--!!    correctdp = .false.
!!--!!    correctsp = .false.
!!--!!    if ( (abs(doublepole) .lt. 1d-30)) then
!!--!!       correctdp = .true.
!!--!!    else
!!--!!       if (verbose) print *,"doublepole", abs( (ampNLO(-2) -doublepole)/doublepole )
!!--!!       if ( (abs( (ampNLO(-2) -doublepole)/doublepole ) .lt. 1d-10)) then
!!--!!          correctdp = .true.
!!--!!       endif
!!--!!    endif
!!--!!
!!--!!    if (abs(singlepole) .lt. 1d-30) then
!!--!!       correctsp = .true.
!!--!!    else
!!--!!       if (verbose) print *, "singlepole",abs( (ampNLO(-1) -singlepole)/singlepole )
!!--!!       if (abs( (ampNLO(-1) -singlepole)/singlepole ) .lt. 1d-10) then 
!!--!!          correctsp = .true.
!!--!!       endif
!!--!!       
!!--!!    endif
!!--!!    
!!--!!    if (.not. correctdp) then
!!--!!       print *, "ERROR: double poles are incorrect!"
!!--!!       print *, "double pole",ampNLO(-2),doublepole
!!--!!       stop
!!--!!    endif
!!--!!    
!!--!!    if (.not. correctsp) then
!!--!!       print *, "ERROR: single poles are incorrect!"
!!--!!       print *, "single pole",ampNLO(-1),singlepole
!!--!!       stop
!!--!!    endif
!!--!!
!!--!!    return
!!--!!  end subroutine pole_check

end module mod_ampNLO_ggzz_mass
  

