module mod_ampLO_ggww
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_avec_mcfm
  use mod_func_for_h2vv
  implicit none
  private

  public :: ewampLO_ggww, ewampLO_ggWWmass

contains

  !-- 0 -> g(1) g(2) [Wp->l(3) lb(4)] [Wp->l(5) lb(6)]
  !-- coefficient of asontwopi * delta(a,b)
  !-- phases as in MCFM
  function ewampLO_ggww(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: ewampLO_ggww(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp) :: myavec(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: prop34,prop56
    real(dp) :: twosinwcosw,cl1(-1:1),cl2(-1:1),cvec(1:2),cax(1:2),q1,q2,Qu,Qd
    integer :: h1,h2,h34,h56
    complex(dp) :: Mloop_uptype(-1:1,-1:1,-1:1,-1:1),Mloop_dntype(-1:1,-1:1,-1:1,-1:1)

    ewampLO_ggww = czero

    myavec = avec(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    !--- propagator factors
    prop34=sprod(j3,j4)/(sprod(j3,j4)-mwsq + ci*mw*gaw)
    prop56=sprod(j5,j6)/(sprod(j5,j6)-mwsq + ci*mw*gaw)

    ewampLO_ggww(:,:,-1,-1) = myavec(:,:,-1,-1)*prop34*prop56*(gwsq/two)**2*ng

    return

  end function ewampLO_ggww

  function ewampLO_ggWWmass(p,za,zb,sprod,mass)
  !-- 0 -> g(1) g(2) [Wp->l(3) lb(4)] [Wp->l(5) lb(6)] through a massive loop
  !-- coefficient of asontwopi * delta(a,b)
  !-- routines due to Campbell, Ellis, Williams; taken from MCFM
  !-- phases as in MCFM
  !-- doesnt have the usual momenta swap potential -- can be included if necessary
    use mod_types; use mod_consts_dp
    implicit none
    real(dp)    :: p(4,6),mass,sprod(6,6)
    complex(dp) :: za(6,6),zb(6,6)
    complex(dp) :: ewampLO_ggWWmass(-1:1,-1:1,-1:1,-1:1)
    real(dp)    :: s12,s56,s34,dot1256,afac,bfac,delta,gden,dot1234,dot3456,pflat(4,12),pdum(4,12)
    complex(dp) :: zaflat(12,12),zbflat(12,12)
    real(dp)    :: sprodflat(12,12)
    complex(dp) :: box(2,2,-2:0),triang(2,2,-2:0),bub(2,2,-2:0),d(2,2,6)
    complex(dp) :: prop34, prop56, spinorweight(-1:1,-1:1,-1:1,-1:1)
    integer :: h1,h2,h34,h56


    ewampLO_ggWWmass = czero
    box    = czero
    triang = czero
    bub    = czero

!--create "flattened" momenta
    s12=sprod(1,2)
    s34=sprod(3,4)
    s56=sprod(5,6)
    dot1256=0.5d0*(s34-s12-s56)
    delta=dot1256**2-s12*s56
    gden=dot1256+sqrt(delta)
    afac=s12/gden
    bfac=s56/gden
    pflat(1:4,1:6)=p(1:4,1:6)
    pflat(1:4,1)=-pflat(1:4,1)
    pflat(1:4,2)=-pflat(1:4,2)


    pflat(1:4,7)=one/(one-afac*bfac) &
         *(pflat(1:4,5)+pflat(1:4,6)-bfac*(pflat(1:4,1)+pflat(1:4,2)))
    pflat(1:4,8)=one/(one-afac*bfac) &
         *(pflat(1:4,1)+pflat(1:4,2)-afac*(pflat(1:4,5)+pflat(1:4,6)))

    dot1234=0.5d0*(s56-s12-s34)
    delta=dot1234**2-s12*s34
    gden=dot1234+sqrt(delta)
    afac=s12/gden
    bfac=s34/gden
    pflat(1:4,9)=one/(one-afac*bfac) &
         *(pflat(1:4,1)+pflat(1:4,2)-afac*(pflat(1:4,3)+pflat(1:4,4)))
    pflat(1:4,10)=one/(one-afac*bfac) &
         *(pflat(1:4,3)+pflat(1:4,4)-bfac*(pflat(1:4,1)+pflat(1:4,2)))

    dot3456=0.5d0*(s12-s56-s34)
    delta=dot3456**2-s56*s34
    gden=dot3456+sqrt(delta)
    afac=s56/gden
    bfac=s34/gden
    pflat(1:4,11)=one/(one-afac*bfac) &
         *(pflat(1:4,5)+pflat(1:4,6)-afac*(pflat(1:4,3)+pflat(1:4,4)))
    pflat(1:4,12)=one/(one-afac*bfac) &
         *(pflat(1:4,3)+pflat(1:4,4)-bfac*(pflat(1:4,5)+pflat(1:4,6)))

    call spinorur(12,pflat,zaflat,zbflat,sprodflat)

    ! put in the pTW safety cut at some point...
    call WWmassivebox6(1,2,3,4,5,6,zaflat,zbflat,sprodflat,mass,d,box)
    call WWmassivetri6(1,2,3,4,5,6,zaflat,zbflat,sprodflat,mass,d,triang)
    call WWmassivebub(1,2,3,4,5,6,zaflat,zbflat,sprodflat,mass,bub)

    ewampLO_ggWWmass(-1,-1,-1,-1) = box(1,1,0) + triang(1,1,0) + bub(1,1,0)
    ewampLO_ggWWmass(-1,+1,-1,-1) = box(1,2,0) + triang(1,2,0) + bub(1,2,0)
    ewampLO_ggWWmass(+1,-1,-1,-1) = box(2,1,0) + triang(2,1,0) + bub(2,1,0)
    ewampLO_ggWWmass(+1,+1,-1,-1) = box(2,2,0) + triang(2,2,0) + bub(2,2,0)

  !--- propagator factors
    prop34=sprodflat(3,4)/(sprodflat(3,4)-mwsq + ci*mw*gaw)
    prop56=sprodflat(5,6)/(sprodflat(5,6)-mwsq + ci*mw*gaw)


! RR --check against MCFM pt1, for mt=173 GeV
!    call spinorur(6,p,za,zb,sprod)
!    include "spinorweight.f90"
!    print *, "ratio -1,-1",ewampLO_ggWWmass(-1,-1,-1,-1)/spinorweight(-1,-1,-1,-1)/dcmplx(  1.6465507106420087d-013,  1.6571966858356929d-013)
!    print *, "ratio -1,+1",ewampLO_ggWWmass(-1,+1,-1,-1)/spinorweight(-1,+1,-1,-1)/dcmplx( -7.9176367891982149d-017, -1.7211745126616654d-016)
!    print *, "ratio +1,-1",ewampLO_ggWWmass(+1,-1,-1,-1)/spinorweight(+1,-1,-1,-1)/dcmplx(  3.0436656697831753d-015, -5.5335081755309963d-015)
!    print *, "ratio +1,+1",ewampLO_ggWWmass(+1,+1,-1,-1)/spinorweight(+1,+1,-1,-1)/dcmplx(  1.8394685133137157d-013,  1.6709164060889685d-013)


    ewampLO_ggWWmass = ewampLO_ggWWmass * prop34 * prop56 * gwsq**2/two



    return
  end function ewampLO_ggWWmass

end module mod_ampLO_ggww
