module mod_avec_mcfm
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_func_for_amp_1loop
  implicit none
  private

  public :: avec

contains

  !-- from src/ZZ/getggZZamps.f
  !-- MCFM pre-factor: 4d0*esq*gsq/(16d0*pisq)*esq * delta(a,b)
  !-- 0 -> g(1) g(2) [Z->l(3) lb(4)] [Z->l(5) lb(6)]

  function avec(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: avec(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)

    avec = czero

    !-- from src/ZZ/getggZZamps.f
    !-- 2 --> +1, 1 --> -1
    avec(+1,+1,-1,-1)=a64v('q+qb-g-g-',j3,j4,j1,j2,j6,j5,zb,za,sprod)*(-ci)
    avec(+1,-1,-1,-1)=a64v('q+qb-g-g+',j3,j4,j1,j2,j6,j5,zb,za,sprod)*(-ci)
    avec(-1,+1,-1,-1)=a64v('q+qb-g+g-',j3,j4,j1,j2,j6,j5,zb,za,sprod)*(-ci)
    avec(-1,-1,-1,-1)=a64v('q+qb-g+g+',j3,j4,j1,j2,j6,j5,zb,za,sprod)*(-ci)

    avec(+1,+1,+1,-1)=a64v('q+qb-g+g+',j3,j4,j1,j2,j5,j6,za,zb,sprod)*(-ci)
    avec(+1,-1,+1,-1)=a64v('q+qb-g+g-',j3,j4,j1,j2,j5,j6,za,zb,sprod)*(-ci)
    avec(-1,+1,+1,-1)=a64v('q+qb-g-g+',j3,j4,j1,j2,j5,j6,za,zb,sprod)*(-ci)
    avec(-1,-1,+1,-1)=a64v('q+qb-g-g-',j3,j4,j1,j2,j5,j6,za,zb,sprod)*(-ci)

    avec(+1,+1,-1,+1)=a64v('q+qb-g-g-',j3,j4,j1,j2,j5,j6,zb,za,sprod)*(-ci)
    avec(+1,-1,-1,+1)=a64v('q+qb-g-g+',j3,j4,j1,j2,j5,j6,zb,za,sprod)*(-ci)
    avec(-1,+1,-1,+1)=a64v('q+qb-g+g-',j3,j4,j1,j2,j5,j6,zb,za,sprod)*(-ci)
    avec(-1,-1,-1,+1)=a64v('q+qb-g+g+',j3,j4,j1,j2,j5,j6,zb,za,sprod)*(-ci)

    avec(+1,+1,+1,+1)=a64v('q+qb-g+g+',j3,j4,j1,j2,j6,j5,za,zb,sprod)*(-ci)
    avec(+1,-1,+1,+1)=a64v('q+qb-g+g-',j3,j4,j1,j2,j6,j5,za,zb,sprod)*(-ci)
    avec(-1,+1,+1,+1)=a64v('q+qb-g-g+',j3,j4,j1,j2,j6,j5,za,zb,sprod)*(-ci)
    avec(-1,-1,+1,+1)=a64v('q+qb-g-g-',j3,j4,j1,j2,j6,j5,za,zb,sprod)*(-ci)

    !--- a64v amplitudes have overall color factor delta(A,B);
    !--- put in factor of two to extract delta(A,B)/2 as in massive case
    avec(:,:,:,:)=two*avec(:,:,:,:)

    return

  end function avec


  !-- from src/Zbb/xzqqgg_v.f
  function a64v(st,j1,j4,j2,j3,j5,j6,za,zb,sprod)
    complex(dp) :: a64v
    character, intent(in) :: st*(9)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    !----definition (2.13) of BDK, writes in terms of fvs and fvf

    if     (st.eq.'q+qb-g-g-') then

       a64v=-fvs('q+qb-g+g+',j4,j1,j3,j2,j6,j5,zb,za,sprod) &
            -fvf('q+qb-g+g+',j4,j1,j3,j2,j6,j5,zb,za,sprod)

       !print *, 'fvs(q+qb-g+g+,j4,j1,j3,j2,j6,j5,zb,za)',abs(fvs('q+qb-g+g+',j4,j1,j3,j2,j6,j5,zb,za,sprod))
       !print *, 'fvf(q+qb-g+g+,j4,j1,j3,j2,j6,j5,zb,za)',abs(fvf('q+qb-g+g+',j4,j1,j3,j2,j6,j5,zb,za,sprod))

    elseif (st.eq.'q+qb-g-g+') then

       a64v=-fvs('q+qb-g+g-',j1,j4,j3,j2,j5,j6,za,zb,sprod) &
            -fvf('q+qb-g+g-',j1,j4,j3,j2,j5,j6,za,zb,sprod)

       !print *, 'fvs(q+qb-g+g-,j1,j4,j3,j2,j5,j6,za,zb)',abs(fvs('q+qb-g+g-',j1,j4,j3,j2,j5,j6,za,zb,sprod))
       !print *, 'fvf(q+qb-g+g-,j1,j4,j3,j2,j5,j6,za,zb)',abs(fvf('q+qb-g+g-',j1,j4,j3,j2,j5,j6,za,zb,sprod))

    else

       a64v=-fvs(st,j1,j4,j2,j3,j5,j6,za,zb,sprod) &
            -fvf(st,j1,j4,j2,j3,j5,j6,za,zb,sprod)

       !print *, 'fvs(st,j1,j4,j2,j3,j5,j6,za,zb)',abs(fvs(st,j1,j4,j2,j3,j5,j6,za,zb,sprod))
       !print *, 'fvf(st,j1,j4,j2,j3,j5,j6,za,zb)',abs(fvf(st,j1,j4,j2,j3,j5,j6,za,zb,sprod))

    endif

    return

  end function a64v


  !-- from src/W2jetvirt/fvf.f
  function fvf(st,j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: fvf
    character, intent(in) :: st*(9)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp) :: I3m123456,I3m563412,zab2
    real(dp) :: t

    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
    t(j1,j2,j3) = sprod(j1,j2)+sprod(j1,j3)+sprod(j2,j3)

    I3m123456=I3m(sprod(j1,j2),sprod(j3,j4),sprod(j5,j6))
    I3m563412=I3m123456

    if(st.eq.'q+qb-g+g-') then

       Fvf= &
            -zab2(j5,j2,j4,j1)**2/(za(j5,j6)*zb(j1,j2)*zab2(j3,j1,j2,j4)**2) &
            *Lsm1_2mh(sprod(j3,j4),t(j1,j2,j4),sprod(j1,j2),sprod(j5,j6)) &
            -zab2(j2,j1,j3,j6)**2/(za(j1,j2)*zb(j5,j6)*zab2(j3,j1,j2,j4)**2) &
            *Lsm1_2mh(sprod(j3,j4),t(j1,j2,j3),sprod(j1,j2),sprod(j5,j6)) &
            !
            -I3m123456*za(j4,j5)*zb(j1,j3) &
            *zab2(j2,j1,j3,j6)/(2.0_dp*zab2(j3,j1,j2,j4)*t(j1,j2,j3)) &
            !
            -I3m123456*za(j2,j4)*zb(j3,j6) &
            *zab2(j5,j2,j4,j1)/(2.0_dp*zab2(j3,j1,j2,j4)*t(j1,j2,j4)) &
            !
            +I3m563412*za(j2,j4)*zb(j3,j6) &
            *zab2(j5,j3,j6,j1)/(2.0_dp*zab2(j3,j1,j2,j4)*t(j3,j5,j6)) &
            !
            +I3m563412*za(j4,j5)*zb(j1,j3) &
            *zab2(j2,j4,j5,j6)/(2.0_dp*zab2(j3,j1,j2,j4)*t(j4,j5,j6))

    elseif(st.eq.'q+qb-g+g+') then

       Fvf=-za(j2,j5)**2/(za(j1,j2)*za(j5,j6)*za(j3,j4)**2) &
            *Lsm1_2me(t(j1,j2,j3),t(j1,j2,j4),sprod(j1,j2),sprod(j5,j6))

    endif

    return

  end function fvf

  !-- from src/W2jetvirt/fvf.f
  function fvs(st,j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: fvs
    character, intent(in) :: st*(9)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    !--- This function returns the result of BDK
    !--- Published in Nucl.Phys.B513:3-86,1998.
    !--- e-Print: hep-ph/9708239
    !--- Eqs.(11.1) and Eqs.(11.6) as appropriate dependent on choice of "st"
    real(dp) :: t
    complex(dp) :: Brackppa,Brackpp

    t(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)

    !---This is the square bracket in Eq.(11.1) subject to exchange 16,25
    Brackppa(j1,j2,j3,j4,j5,j6)= &
         za(j1,j2)*(za(j5,j3)*zb(j3,j1))**2/za(j5,j6)/za(j3,j4)**2 &
         *L1(-t(j1,j2,j3),-sprod(j1,j2))/sprod(j1,j2)**2 &
         +za(j5,j2)*za(j5,j3)*zb(j3,j1)/za(j5,j6)/za(j3,j4)**2 &
         *L0(-t(j1,j2,j3),-sprod(j1,j2))/sprod(j1,j2)

    !---This is the whole expression in Eq.(11.1) subject to exchange 34
    Brackpp(j1,j2,j3,j4,j5,j6)= &
         -za(j2,j3)*za(j2,j4)*za(j3,j5)*za(j4,j5) &
         /za(j1,j2)/za(j5,j6)/za(j3,j4)**4 &
         *Lsm1_2me(t(j1,j2,j3),t(j1,j2,j4),sprod(j1,j2),sprod(j5,j6)) &
         -(za(j2,j4)*za(j3,j5)+za(j2,j3)*za(j4,j5)) &
         /za(j1,j2)/za(j5,j6)/za(j3,j4)**3 &
         *(za(j5,j3)*zb(j3,j1)*za(j1,j2) &
         *L0(-t(j1,j2,j3),-sprod(j1,j2))/sprod(j1,j2) &
         +za(j5,j6)*zb(j6,j4)*za(j4,j2) &
         *L0(-t(j1,j2,j3),-sprod(j5,j6))/sprod(j5,j6) &
         +0.5_dp*za(j5,j2)*(Lnrat(-t(j1,j2,j3),-t(j1,j2,j4)))) &
         -Brackppa(j1,j2,j3,j4,j5,j6)-Brackppa(j6,j5,j3,j4,j2,j1) &
         +0.5_dp/za(j3,j4)**2 &
         *(za(j2,j5)**2/za(j1,j2)/za(j5,j6) &
         -zb(j1,j6)**2/zb(j1,j2)/zb(j5,j6))

    if(st.eq.'q+qb-g+g-') then
       !---flip2:( 1<-->2, 3<-->4, 5<-->6, za<-->zb)
       Fvs= &
            +Brackpm(j1,j2,j3,j4,j5,j6,za,zb,sprod) &
            +Brackpm(j2,j1,j4,j3,j6,j5,zb,za,sprod)

       !write(6,*) 'Brackpm(j1,j2,j3,j4,j5,j6,za,zb,sprod)',abs(Brackpm(j1,j2,j3,j4,j5,j6,za,zb,sprod))

    elseif(st.eq.'q+qb-g+g+') then

       !---exch34:(3<-->4)
       Fvs= &
            +Brackpp(j1,j2,j3,j4,j5,j6) &
            +Brackpp(j1,j2,j4,j3,j5,j6)

    endif

    return

  end function fvs

  function Brackpm(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: brackpm
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    !---This is the whole expression in Eq.(11.6) subject to exchange flip_2
    !---ie  flip2:( 1<-->2, 3<-->4, 5<-->6, za<-->zb)
    real(dp) :: t,delta,IDelta
    complex(dp) :: zab2

    t(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
    delta(j1,j2,j3,j4,j5,j6)=sprod(j1,j2)-sprod(j3,j4)-sprod(j5,j6)
    zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

    IDelta=1.0_dp/(sprod(j1,j2)**2+sprod(j3,j4)**2+sprod(j5,j6)**2 &
         -2.0_dp*(+sprod(j1,j2)*sprod(j3,j4)+sprod(j1,j2)*sprod(j5,j6)+sprod(j5,j6)*sprod(j3,j4)))

    Brackpm= &
         -2.0_dp*za(j2,j3)*zb(j4,j6)*zab2(j3,j1,j2,j6) &
         *zab2(j2,j5,j6,j4)*t(j1,j2,j3) &
         /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**4 &
         *Lsm1_2mh(sprod(j3,j4),t(j1,j2,j3),sprod(j1,j2),sprod(j5,j6)) &
         !
         &+(2.0_dp*za(j2,j3)*zb(j4,j6)*(zb(j1,j2)*za(j5,j6)*za(j2,j3)*zb(j4,j6) &
         !
         & -zab2(j5,j2,j4,j1)*zab2(j3,j1,j2,j4))*zab2(j4,j1,j2,j3) &
         & *(t(j1,j2,j3)-t(j1,j2,j4))/zab2(j3,j1,j2,j4)**3*IDelta &
         !
         -3.0_dp*(sprod(j3,j4)*delta(j3,j4,j1,j2,j5,j6)*zab2(j2,j3,j4,j1) &
         *zab2(j5,j3,j4,j6)*IDelta &
         -za(j2,j3)*zb(j3,j1)*za(j5,j4)*zb(j4,j6) &
         -za(j2,j4)*zb(j4,j1)*za(j5,j3)*zb(j3,j6)) &
         !
         *zab2(j4,j1,j2,j3)/zab2(j3,j1,j2,j4)*IDelta &
         -zb(j1,j4)*za(j3,j5)*za(j2,j3)*zb(j4,j6) &
         *zab2(j4,j1,j2,j3)**2/zab2(j3,j1,j2,j4)**2*IDelta &
         -zb(j1,j3)*za(j4,j5)*za(j2,j4)*zb(j3,j6)*IDelta) &
         *I3m(sprod(j1,j2),sprod(j3,j4),sprod(j5,j6)) &
         !
         +(2.0_dp*za(j2,j3)*zb(j4,j6)*t(j1,j2,j3)*zab2(j2,j1,j4,j6) &
         /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**3 &
         !
         -(zab2(j2,j1,j4,j6)**2 &
         +2.0_dp*za(j2,j3)*zb(j3,j6)*za(j2,j4)*zb(j4,j6)) &
         /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**2) &
         *Lnrat(-t(j1,j2,j3),-sprod(j3,j4)) &
         !
         !---ie  flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)
         +Brackpma(j1,j2,j3,j4,j5,j6,za,zb,sprod) &
         +Brackpma(j5,j6,j4,j3,j1,j2,zb,za,sprod) &
         !
         +zb(j1,j6)*zab2(j4,j1,j2,j3) &
         *(zb(j1,j6)*delta(j3,j4,j5,j6,j1,j2) &
         -2.0_dp*zb(j1,j2)*za(j2,j5)*zb(j5,j6)) &
         /zb(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)*IDelta &
         !
         +(zab2(j2,j1,j4,j6)**2+za(j2,j1)*zb(j1,j6)*za(j2,j5)*zb(j5,j6)) &
         /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**2

    return

  end function Brackpm

  function Brackpma(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: Brackpma
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    !---This is the curly braces expression in Eq.(11.6) subject to flip_3
    !---ie  flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)
    real(dp) :: t
    complex(dp) :: delta,zab2,IDelta

    t(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
    delta(j1,j2,j3,j4,j5,j6)=sprod(j1,j2)-sprod(j3,j4)-sprod(j5,j6)
    zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

    IDelta=1.0_dp/(sprod(j1,j2)**2+sprod(j3,j4)**2+sprod(j5,j6)**2 &
         -2.0_dp*(sprod(j1,j2)*sprod(j3,j4)+sprod(j1,j2)*sprod(j5,j6)+sprod(j5,j6)*sprod(j3,j4)))

    Brackpma= &
         +2.0_dp*zb(j1,j2)*za(j2,j3)*zb(j3,j6) &
         /zb(j5,j6)/zab2(j3,j1,j2,j4)**2 &
         *(zab2(j3,j1,j2,j6)*zab2(j2,j1,j3,j4)/zab2(j3,j1,j2,j4) &
         *L0(-t(j1,j2,j3),-sprod(j1,j2))/sprod(j1,j2) &
         -za(j2,j3)*zb(j3,j6) &
         *(L0(-t(j1,j2,j3),-sprod(j1,j2))/sprod(j1,j2) &
         -1/2.0_dp*L1(-t(j1,j2,j3),-sprod(j1,j2))/sprod(j1,j2))) &
         !
         -(3.0_dp*delta(j5,j6,j3,j4,j1,j2)*zab2(j2,j3,j4,j1) &
         *zab2(j5,j3,j4,j6)*zab2(j4,j1,j2,j3) &
         /zab2(j3,j1,j2,j4)*IDelta**2 &
         !
         +2.0_dp*za(j3,j2)*zb(j2,j1)*zb(j4,j6)*(t(j1,j2,j3)-t(j1,j2,j4)) &
         *(za(j2,j5)*t(j1,j2,j3)+za(j2,j1)*zb(j1,j6)*za(j6,j5)) &
         /zab2(j3,j1,j2,j4)**3*IDelta &
         !
         -(2.0_dp*za(j2,j3)*zb(j3,j6)*(t(j1,j2,j3)-t(j1,j2,j4)) &
         +delta(j1,j2,j3,j4,j5,j6) &
         *(za(j2,j5)*t(j1,j2,j3) &
         +za(j2,j1)*zb(j1,j6)*za(j6,j5))/za(j5,j6)) &
         *za(j5,j2)*zb(j2,j1)/zab2(j3,j1,j2,j4)**2*IDelta &
         !
         +0.5_dp*(za(j1,j2)*zb(j1,j6)**2/zb(j5,j6) &
         +zb(j1,j2)*za(j2,j5)**2/za(j5,j6) &
         -2.0_dp*za(j2,j5)*zb(j1,j6))*zab2(j4,j1,j2,j3) &
         /zab2(j3,j1,j2,j4)*IDelta) &
         *(Lnrat(-sprod(j1,j2),-sprod(j3,j4)))

    return

  end function Brackpma

end module mod_avec_mcfm
