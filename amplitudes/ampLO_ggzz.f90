module mod_ampLO_ggzz
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_avec_mcfm
  implicit none
  private

  public :: ewampLO_ggzz

contains

  !-- from src/ZZ/getggZZamps.f

  !-- 0 -> g(1) g(2) [Z->l(3) lb(4)] [Z->l(5) lb(6)]
  !-- coefficient of asontwopi * delta(a,b)
  !-- phases as in MCFM
  function ewampLO_ggzz(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: ewampLO_ggzz(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp) :: myavec(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: prop34,prop56
    real(dp) :: twosinwcosw,cl1(-1:1),cl2(-1:1),cvec(1:2),cax(1:2),q1,q2,Qu,Qd
    integer :: h1,h2,h34,h56
    complex(dp) :: Mloop_uptype(-1:1,-1:1,-1:1,-1:1),Mloop_dntype(-1:1,-1:1,-1:1,-1:1)

    ewampLO_ggzz = czero

    myavec = avec(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    !--- propagator factors
    prop34=sprod(j3,j4)/(sprod(j3,j4)-mzsq + ci*mz*gaz)
    prop56=sprod(j5,j6)/(sprod(j5,j6)-mzsq + ci*mz*gaz)

    !-- couplings, see src/Need/couplz.f and coupling.f
    !-- xw = sin(tw)^2, sin2w = 2*sin(tw)*cos(tw)

    twosinwcosw=two*sw*cw !maybe I'll have to modify this (date: 23/03/24)

    cl1(-1)=Llep/(twosinwcosw)
    cl1(+1)=Rlep/(twosinwcosw)
    cl2(-1)=Llep/(twosinwcosw)
    cl2(+1)=Rlep/(twosinwcosw)

    cvec(up)=Vup/(twosinwcosw)
    cvec(dn)=Vdn/(twosinwcosw)
    cax(up)=Aup/(twosinwcosw)
    cax(dn)=Adn/(twosinwcosw)

    ! to be changed!!!!
    q1=-one
    q2=-one
    Qu=Qup
    Qd=Qdn

    !--- dress vector and axial amplitudes with appropriate couplings
    do h1=-1,1,2
    do h2=-1,1,2
    do h34=-1,1,2
    do h56=-1,1,2

       !--- internal loops of massless up-type quarks
       Mloop_uptype(h1,h2,h34,h56)=myavec(h1,h2,h34,h56)*( &
            (Qu*q1+cvec(up)*cl1(h34)*prop34)*(Qu*q2+cvec(up)*cl2(h56)*prop56) &
            +(cax(up)*cl1(h34)*prop34)*(cax(up)*cl2(h56)*prop56))

       !--- internal loops of massless down-type quarks
       Mloop_dntype(h1,h2,h34,h56)=myavec(h1,h2,h34,h56)*( &
            (Qd*q1+cvec(dn)*cl1(h34)*prop34)*(Qd*q2+cvec(dn)*cl2(h56)*prop56) &
            +(cax(dn)*cl1(h34)*prop34)*(cax(dn)*cl2(h56)*prop56))

    enddo
    enddo
    enddo
    enddo

    !--
    ewampLO_ggzz(:,:,:,:) = (nup*Mloop_uptype(:,:,:,:)+ndn*Mloop_dntype(:,:,:,:)) * (two * (gwsq * sinW2)**2)

    return

  end function ewampLO_ggzz

end module mod_ampLO_ggzz
