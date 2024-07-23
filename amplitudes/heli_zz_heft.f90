!-- 0 -> g(1) g(2) [Z -> l(3) lb(4)] [Z -> l(5) lb(6)]
!-- returns the helicity dependent parts of the amplitude with no couplings included
!-- 1-20 is the number of the form factor that it multiplies
!-- h1,h2,h3,h5 are the helicity indices
subroutine heli_zz_heft(p1,p2,p3,p4,p5,p6,za,zb,s1,T)
  use mod_types; use mod_consts_dp
  implicit none
  integer, intent(in) :: p1,p2,p3,p4,p5,p6
  complex(dp), intent(in) :: za(6,6),zb(6,6)
  real(dp), intent(in) :: s1(6,6)
  complex(dp), intent(out) :: T(20,-1:1,-1:1,-1:1,-1:1)

  T = czero

  !-- T1 =
  T(1,-1,-1,-1,-1) = 1._dp/2._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p4,p6)*za(p2,p1)*za(p5,p3)

  T(1,-1,-1,-1,1) = 1._dp/2._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p4,p5)*za(p2,p1)*za(p6,p3)

  T(1,-1,-1,1,-1) = 1._dp/2._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p6)*za(p2,p1)*za(p5,p4)

  T(1,-1,-1,1,1) = 1._dp/2._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p5)*za(p2,p1)*za(p6,p4)

  T(1,-1,1,-1,-1) =  0

  T(1,-1,1,-1,1) =  0

  T(1,-1,1,1,-1) =  0

  T(1,-1,1,1,1) =  0

  T(1,1,-1,-1,-1) =  0

  T(1,1,-1,-1,1) =  0

  T(1,1,-1,1,-1) =  0

  T(1,1,-1,1,1) =  0

  T(1,1,1,-1,-1) = 1._dp/2._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p2)*zb(p4,p6)*za(p5,p3)

  T(1,1,1,-1,1) = 1._dp/2._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p4,p5)*za(p6,p3)

  T(1,1,1,1,-1) = 1._dp/2._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p3,p6)*za(p5,p4)

  T(1,1,1,1,1) = 1._dp/2._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p2)*zb(p3,p5)*za(p6,p4)

  !-- T2 =
  T(2,-1,-1,-1,-1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p6)*zb(p2,p4)*za(p3,p1)*za(p5,p2)

  T(2,-1,-1,-1,1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p5)*zb(p2,p4)*za(p3,p1)*za(p6,p2)

  T(2,-1,-1,1,-1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p6)*zb(p2,p3)*za(p4,p1)*za(p5,p2)

  T(2,-1,-1,1,1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p5)*zb(p2,p3)*za(p4,p1)*za(p6,p2)

  T(2,-1,1,-1,-1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p4)*zb(p2,p6)*za(p3,p1)*za(p5,p1)

  T(2,-1,1,-1,1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p4)*zb(p2,p5)*za(p3,p1)*za(p6,p1)

  T(2,-1,1,1,-1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p2,p6)*za(p4,p1)*za(p5,p1)

  T(2,-1,1,1,1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p2,p5)*za(p4,p1)*za(p6,p1)

  T(2,1,-1,-1,-1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)*zb(p1,p6)*za(p3,p2)*za(p5,p2)

  T(2,1,-1,-1,1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p4)*zb(p1,p5)*za(p3,p2)*za(p6,p2)

  T(2,1,-1,1,-1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p3)*zb(p1,p6)*za(p4,p2)*za(p5,p2)

  T(2,1,-1,1,1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p3)*zb(p1,p5)*za(p4,p2)*za(p6,p2)

  T(2,1,1,-1,-1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)*zb(p2,p6)*za(p3,p2)*za(p5,p1)

  T(2,1,1,-1,1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p4)*zb(p2,p5)*za(p3,p2)*za(p6,p1)

  T(2,1,1,1,-1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)*zb(p2,p6)*za(p4,p2)*za(p5,p1)

  T(2,1,1,1,1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)*zb(p2,p5)*za(p4,p2)*za(p6,p1)

  !-- T3 =
  T(3,-1,-1,-1,-1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p4)*zb(p2,p6)*za(p3,p2)*za(p5,p1)

  T(3,-1,-1,-1,1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p4)*zb(p2,p5)*za(p3,p2)*za(p6,p1)

  T(3,-1,-1,1,-1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)*zb(p2,p6)*za(p4,p2)*za(p5,p1)

  T(3,-1,-1,1,1) =  - 1._dp/2._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p2,p5)*za(p4,p2)*za(p6,p1)

  T(3,-1,1,-1,-1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p4)*zb(p2,p6)*za(p3,p1)*za(p5,p1)

  T(3,-1,1,-1,1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p4)*zb(p2,p5)*za(p3,p1)*za(p6,p1)

  T(3,-1,1,1,-1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p2,p6)*za(p4,p1)*za(p5,p1)

  T(3,-1,1,1,1) =  - 1._dp/2._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p2,p5)*za(p4,p1)*za(p6,p1)

  T(3,1,-1,-1,-1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)*zb(p1,p6)*za(p3,p2)*za(p5,p2)

  T(3,1,-1,-1,1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p4)*zb(p1,p5)*za(p3,p2)*za(p6,p2)

  T(3,1,-1,1,-1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p3)*zb(p1,p6)*za(p4,p2)*za(p5,p2)

  T(3,1,-1,1,1) = 1._dp/2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p3)*zb(p1,p5)*za(p4,p2)*za(p6,p2)

  T(3,1,1,-1,-1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p6)*zb(p2,p4)*za(p3,p1)*za(p5,p2)

  T(3,1,1,-1,1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p5)*zb(p2,p4)*za(p3,p1)*za(p6,p2)

  T(3,1,1,1,-1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p6)*zb(p2,p3)*za(p4,p1)*za(p5,p2)

  T(3,1,1,1,1) =  - 1._dp/2._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p2,p3)*za(p4,p1)*za(p6,p2)

  !-- T4 =
  T(4,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p4,p1)*zb(p6,p1)*za(p1,p3)*za(p1,p5)*za(p2,p1)

  T(4,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p4,p1)*zb(p5,p1)*za(p1,p3)*za(p1,p6)*za(p2,p1)

  T(4,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p1)*zb(p6,p1)*za(p1,p4)*za(p1,p5)*za(p2,p1)

  T(4,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p1)*zb(p5,p1)*za(p1,p4)*za(p1,p6)*za(p2,p1)

  T(4,-1,1,-1,-1) =  0

  T(4,-1,1,-1,1) =  0

  T(4,-1,1,1,-1) =  0

  T(4,-1,1,1,1) =  0

  T(4,1,-1,-1,-1) =  0

  T(4,1,-1,-1,1) =  0

  T(4,1,-1,1,-1) =  0

  T(4,1,-1,1,1) =  0

  T(4,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p2)*zb(p4,p1)*zb(p6,p1)*za(p1,p3)*za(p1,p5)

  T(4,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p4,p1)*zb(p5,p1)*za(p1,p3)*za(p1,p6)

  T(4,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p3,p1)*zb(p6,p1)*za(p1,p4)*za(p1,p5)

  T(4,1,1,1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p2)*zb(p3,p1)*zb(p5,p1)*za(p1,p4)*za(p1,p6)

  !-- T5 =
  T(5,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p4,p1)*zb(p6,p2)*za(p1,p3)*za(p2,p1)*za(p2,p5)

  T(5,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p4,p1)*zb(p5,p2)*za(p1,p3)*za(p2,p1)*za(p2,p6)

  T(5,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p1)*zb(p6,p2)*za(p1,p4)*za(p2,p1)*za(p2,p5)

  T(5,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p1)*zb(p5,p2)*za(p1,p4)*za(p2,p1)*za(p2,p6)

  T(5,-1,1,-1,-1) =  0

  T(5,-1,1,-1,1) =  0

  T(5,-1,1,1,-1) =  0

  T(5,-1,1,1,1) =  0

  T(5,1,-1,-1,-1) =  0

  T(5,1,-1,-1,1) =  0

  T(5,1,-1,1,-1) =  0

  T(5,1,-1,1,1) =  0

  T(5,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p2)*zb(p4,p1)*zb(p6,p2)*za(p1,p3)*za(p2,p5)

  T(5,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p4,p1)*zb(p5,p2)*za(p1,p3)*za(p2,p6)

  T(5,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p3,p1)*zb(p6,p2)*za(p1,p4)*za(p2,p5)

  T(5,1,1,1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p2)*zb(p3,p1)*zb(p5,p2)*za(p1,p4)*za(p2,p6)

  !-- T6 =
  T(6,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p4,p2)*zb(p6,p1)*za(p1,p5)*za(p2,p1)*za(p2,p3)

  T(6,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p4,p2)*zb(p5,p1)*za(p1,p6)*za(p2,p1)*za(p2,p3)

  T(6,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p2)*zb(p6,p1)*za(p1,p5)*za(p2,p1)*za(p2,p4)

  T(6,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p2)*zb(p5,p1)*za(p1,p6)*za(p2,p1)*za(p2,p4)

  T(6,-1,1,-1,-1) =  0

  T(6,-1,1,-1,1) =  0

  T(6,-1,1,1,-1) =  0

  T(6,-1,1,1,1) =  0

  T(6,1,-1,-1,-1) =  0

  T(6,1,-1,-1,1) =  0

  T(6,1,-1,1,-1) =  0

  T(6,1,-1,1,1) =  0

  T(6,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p2)*zb(p4,p2)*zb(p6,p1)*za(p1,p5)*za(p2,p3)

  T(6,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p4,p2)*zb(p5,p1)*za(p1,p6)*za(p2,p3)

  T(6,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p3,p2)*zb(p6,p1)*za(p1,p5)*za(p2,p4)

  T(6,1,1,1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p2)*zb(p3,p2)*zb(p5,p1)*za(p1,p6)*za(p2,p4)

  !-- T7 =
  T(7,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p4,p2)*zb(p6,p2)*za(p2,p1)*za(p2,p3)*za(p2,p5)

  T(7,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p4,p2)*zb(p5,p2)*za(p2,p1)*za(p2,p3)*za(p2,p6)

  T(7,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p2)*zb(p6,p2)*za(p2,p1)*za(p2,p4)*za(p2,p5)

  T(7,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p3,p2)*zb(p5,p2)*za(p2,p1)*za(p2,p4)*za(p2,p6)

  T(7,-1,1,-1,-1) =  0

  T(7,-1,1,-1,1) =  0

  T(7,-1,1,1,-1) =  0

  T(7,-1,1,1,1) =  0

  T(7,1,-1,-1,-1) =  0

  T(7,1,-1,-1,1) =  0

  T(7,1,-1,1,-1) =  0

  T(7,1,-1,1,1) =  0

  T(7,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p2)*zb(p4,p2)*zb(p6,p2)*za(p2,p3)*za(p2,p5)

  T(7,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p4,p2)*zb(p5,p2)*za(p2,p3)*za(p2,p6)

  T(7,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1&
       ,p2)*zb(p3,p2)*zb(p6,p2)*za(p2,p4)*za(p2,p5)

  T(7,1,1,1,1) = 1._dp/4._dp/(za(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p2)*zb(p3,p2)*zb(p5,p2)*za(p2,p4)*za(p2,p6)

  !-- T8 =
  T(8,-1,-1,-1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
      /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p1)*&
       za(p1,p5)*za(p3,p2)*s1(p1,p3)

  T(8,-1,-1,-1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p1)*za(&
       p1,p6)*za(p3,p2)*s1(p1,p3)

  T(8,-1,-1,1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p3)*zb(p6,p1)*za(p1,p5)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p1)*za(&
       p1,p5)*za(p4,p2)*s1(p1,p4)

  T(8,-1,-1,1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p3)*zb(p5,p1)*za(p1,p6)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p1)*za(&
       p1,p6)*za(p4,p2)*s1(p1,p4)

  T(8,-1,1,-1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)*zb(p2,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p6,p1)*za(&
       p1,p5)*za(p3,p1)*za(p4,p1)

  T(8,-1,1,-1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)*zb(p2,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p5,p1)*za(&
       p1,p6)*za(p3,p1)*za(p4,p1)

  T(8,-1,1,1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p6,p1)*za(p1,p5)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p6,&
       p1)*za(p1,p5)*za(p4,p1)**2

  T(8,-1,1,1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p5,p1)*za(p1,p6)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p5,&
       p1)*za(p1,p6)*za(p4,p1)**2

  T(8,1,-1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p2)**2&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p4)**2*zb(p6,p1)*za(p1,p5)*za(p3,p2)*za(p4,p2)

  T(8,1,-1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p2)**2 + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p4)**2*zb(p5,p1)*za(p1,p6)*za(p3,p2)*za(p4,p2)

  T(8,1,-1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)**2*zb(p6,p1)*za(p1,p5)*za(p3,p2)*za(p4,p2) + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p3)*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p4,p2)**2

  T(8,1,-1,1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)**2*zb(p5,p1)*za(p1,p6)*za(p3,p2)*za(p4,p2) + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p3)*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p4,p2)**2

  T(8,1,1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p1)*za(&
       p1,p5)*za(p3,p2)*s1(p1,p4)

  T(8,1,1,-1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p1)*za(&
       p1,p6)*za(p3,p2)*s1(p1,p4)

  T(8,1,1,1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p1,p3)*zb(p6,p1)*za(p1,p5)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p1)*za(&
       p1,p5)*za(p4,p2)*s1(p1,p3)

  T(8,1,1,1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p3)*zb(p5,p1)*za(p1,p6)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p1)*za(&
       p1,p6)*za(p4,p2)*s1(p1,p3)

  !-- T9 =
  T(9,-1,-1,-1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p2)*&
       za(p2,p5)*za(p3,p2)*s1(p1,p3)

  T(9,-1,-1,-1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p2)*za(&
       p2,p6)*za(p3,p2)*s1(p1,p3)

  T(9,-1,-1,1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p3)*zb(p6,p2)*za(p2,p5)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p2)*za(&
       p2,p5)*za(p4,p2)*s1(p1,p4)

  T(9,-1,-1,1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p3)*zb(p5,p2)*za(p2,p6)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p2)*za(&
       p2,p6)*za(p4,p2)*s1(p1,p4)

  T(9,-1,1,-1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)*zb(p2,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p6,p2)*za(&
       p2,p5)*za(p3,p1)*za(p4,p1)

  T(9,-1,1,-1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)*zb(p2,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p5,p2)*za(&
       p2,p6)*za(p3,p1)*za(p4,p1)

  T(9,-1,1,1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p6,p2)*za(p2,p5)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p6,&
       p2)*za(p2,p5)*za(p4,p1)**2

  T(9,-1,1,1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p5,p2)*za(p2,p6)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p5,&
       p2)*za(p2,p6)*za(p4,p1)**2

  T(9,1,-1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p2)**2&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p4)**2*zb(p6,p2)*za(p2,p5)*za(p3,p2)*za(p4,p2)

  T(9,1,-1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p2)**2 + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p4)**2*zb(p5,p2)*za(p2,p6)*za(p3,p2)*za(p4,p2)

  T(9,1,-1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)**2*zb(p6,p2)*za(p2,p5)*za(p3,p2)*za(p4,p2) + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p3)*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p4,p2)**2

  T(9,1,-1,1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)**2*zb(p5,p2)*za(p2,p6)*za(p3,p2)*za(p4,p2) + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p3)*zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p4,p2)**2

  T(9,1,1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p2)*za(&
       p2,p5)*za(p3,p2)*s1(p1,p4)

  T(9,1,1,-1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p2)*za(&
       p2,p6)*za(p3,p2)*s1(p1,p4)

  T(9,1,1,1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p1,p3)*zb(p6,p2)*za(p2,p5)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p2)*za(&
       p2,p5)*za(p4,p2)*s1(p1,p3)

  T(9,1,1,1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p3)*zb(p5,p2)*za(p2,p6)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p2)*za(&
       p2,p6)*za(p4,p2)*s1(p1,p3)

  !-- T10 =
  T(10,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)**2*zb(p2,p6)*za(p1,p3)*za(p4,p2)*za(p5,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p4,p1)*&
       za(p3,p2)*za(p5,p1)*s1(p1,p3)

  T(10,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p4)**2*zb(p2,p5)*za(p1,p3)*za(p4,p2)*za(p6,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p4,p1)*&
       za(p3,p2)*za(p6,p1)*s1(p1,p3)

  T(10,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)**2*zb(p2,p6)*za(p1,p4)*za(p3,p2)*za(p5,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p3,p1)*&
       za(p4,p2)*za(p5,p1)*s1(p1,p4)

  T(10,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)**2*zb(p2,p5)*za(p1,p4)*za(p3,p2)*za(p6,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p3,p1)*za(&
       p4,p2)*za(p6,p1)*s1(p1,p4)

  T(10,-1,1,-1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p2,p6)*zb(p4,p1)*za(p1,p3)**2*za(p5,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p6)*za(p1,&
       p3)*za(p5,p1)*s1(p4,p1)

  T(10,-1,1,-1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p2,p5)*zb(p4,p1)*za(p1,p3)**2*za(p6,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p5)*za(p1,&
       p3)*za(p6,p1)*s1(p4,p1)

  T(10,-1,1,1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p2,p6)*za(p1,p4)*za(p5,p1)*s1(p3,p1) + 1._dp/4._dp/(&
       s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p6)*zb(p3,p1)&
       *za(p1,p4)**2*za(p5,p1)

  T(10,-1,1,1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p2,p5)*za(p1,p4)*za(p6,p1)*s1(p3,p1) + 1._dp/4._dp/(&
       s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p5)*zb(p3,p1)&
       *za(p1,p4)**2*za(p6,p1)

  T(10,1,-1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p4)**2*zb(p1,p6)*za(p1,p3)*za(p4,p2)*za(p5,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p4,p1)*za(p3,p2)*za(p5,p2)*s1(p1,p3)

  T(10,1,-1,-1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p4)**2*zb(p1,p5)*za(p1,p3)*za(p4,p2)*za(p6,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p4,p1)*za(p3,p2)*za(p6,p2)*s1(p1,p3)

  T(10,1,-1,1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p1,p6)*za(p1,p4)*za(p3,p2)*za(p5,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p3,p1)*za(p4,p2)*za(p5,p2)*s1(p1,p4)

  T(10,1,-1,1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**2*zb(p1,p5)*za(p1,p4)*za(p3,p2)*za(p6,p2)&
        - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p3,p1)*za(p4,p2)*za(p6,p2)*s1(p1,p4)

  T(10,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p2,p3)*zb(p4,p1)*za(p1,p3)**2*za(p5,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p2,p4)*za(&
       p1,p3)*za(p5,p2)*s1(p4,p1)

  T(10,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p2,p3)*zb(p4,p1)*za(p1,p3)**2*za(p6,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p2,p4)*za(&
       p1,p3)*za(p6,p2)*s1(p4,p1)

  T(10,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p2,p3)*za(p1,p4)*za(p5,p2)*s1(p3,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p2,p4)*zb(&
       p3,p1)*za(p1,p4)**2*za(p5,p2)

  T(10,1,1,1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p2,p3)*za(p1,p4)*za(p6,p2)*s1(p3,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p2,p4)*zb(&
       p3,p1)*za(p1,p4)**2*za(p6,p2)

  !-- T11 =
  T(11,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p2,p6)*zb(p4,p2)*za(p2,p3)**2*za(p5,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p6)*&
       za(p2,p3)*za(p5,p1)*s1(p4,p2)

  T(11,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)*zb(p2,p5)*zb(p4,p2)*za(p2,p3)**2*za(p6,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p5)*&
       za(p2,p3)*za(p6,p1)*s1(p4,p2)

  T(11,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)*zb(p2,p6)*za(p2,p4)*za(p5,p1)*s1(p3,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p6)*zb(&
       p3,p2)*za(p2,p4)**2*za(p5,p1)

  T(11,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)*zb(p2,p5)*za(p2,p4)*za(p6,p1)*s1(p3,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p5)*zb(&
       p3,p2)*za(p2,p4)**2*za(p6,p1)

  T(11,-1,1,-1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p4)**2*zb(p2,p6)*za(p2,p3)*za(p4,p1)*za(p5,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p4,p2)*za(p3,&
       p1)*za(p5,p1)*s1(p2,p3)

  T(11,-1,1,-1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p4)**2*zb(p2,p5)*za(p2,p3)*za(p4,p1)*za(p6,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p4,p2)*za(p3,&
       p1)*za(p6,p1)*s1(p2,p3)

  T(11,-1,1,1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p2,p6)*za(p2,p4)*za(p3,p1)*za(p5,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p3,p2)*za(p4,&
       p1)*za(p5,p1)*s1(p2,p4)

  T(11,-1,1,1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p2,p5)*za(p2,p4)*za(p3,p1)*za(p6,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p3,p2)*za(p4,&
       p1)*za(p6,p1)*s1(p2,p4)

  T(11,1,-1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p6)*zb(p4,p2)*za(p2,p3)**2*za(p5,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p1,p6)*za(p2,p3)*za(p5,p2)*s1(p4,p2)

  T(11,1,-1,-1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p5)*zb(p4,p2)*za(p2,p3)**2*za(p6,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p1,p5)*za(p2,p3)*za(p6,p2)*s1(p4,p2)

  T(11,1,-1,1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p6)*za(p2,p4)*za(p5,p2)*s1(p3,p2)&
        - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p1,p6)*zb(p3,p2)*za(p2,p4)**2*za(p5,p2)

  T(11,1,-1,1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)*zb(p1,p5)*za(p2,p4)*za(p6,p2)*s1(p3,p2) - 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p4)*zb(p1,p5)*zb(p3,p2)*za(p2,p4)**2*za(p6,p2)

  T(11,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p2,p4)**2*za(p2,p3)*za(p4,p1)*za(p5,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p4,p2)*za(&
       p3,p1)*za(p5,p2)*s1(p2,p3)

  T(11,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p2,p4)**2*za(p2,p3)*za(p4,p1)*za(p6,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p4,p2)*za(&
       p3,p1)*za(p6,p2)*s1(p2,p3)

  T(11,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p2,p3)**2*za(p2,p4)*za(p3,p1)*za(p5,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p3,p2)*za(&
       p4,p1)*za(p5,p2)*s1(p2,p4)

  T(11,1,1,1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p2,p3)**2*za(p2,p4)*za(p3,p1)*za(p6,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p3,p2)*za(&
       p4,p1)*za(p6,p2)*s1(p2,p4)

  !-- T12 =
  T(12,-1,-1,-1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p1)*&
       za(p1,p5)*za(p3,p2)*s1(p1,p4)

  T(12,-1,-1,-1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p1)*&
       za(p1,p6)*za(p3,p2)*s1(p1,p4)

  T(12,-1,-1,1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p3)*zb(p6,p1)*za(p1,p5)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p1)*&
       za(p1,p5)*za(p4,p2)*s1(p1,p3)

  T(12,-1,-1,1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p3)*zb(p5,p1)*za(p1,p6)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p1)*za(&
       p1,p6)*za(p4,p2)*s1(p1,p3)

  T(12,-1,1,-1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p2,p3)*zb(p2,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p6,p1)*&
       za(p1,p5)*za(p3,p1)*za(p4,p1)

  T(12,-1,1,-1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)*zb(p2,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p5,p1)*za(&
       p1,p6)*za(p3,p1)*za(p4,p1)

  T(12,-1,1,1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p6,p1)*za(p1,p5)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p6,&
       p1)*za(p1,p5)*za(p4,p1)**2

  T(12,-1,1,1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p5,p1)*za(p1,p6)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p5,&
       p1)*za(p1,p6)*za(p4,p1)**2

  T(12,1,-1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p2)**2&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p4)**2*zb(p6,p1)*za(p1,p5)*za(p3,p2)*za(p4,p2)

  T(12,1,-1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p2)**2&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p4)**2*zb(p5,p1)*za(p1,p6)*za(p3,p2)*za(p4,p2)

  T(12,1,-1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)**2*zb(p6,p1)*za(p1,p5)*za(p3,p2)*za(p4,p2)&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p3)*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p4,p2)**2

  T(12,1,-1,1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)**2*zb(p5,p1)*za(p1,p6)*za(p3,p2)*za(p4,p2) + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p3)*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p4,p2)**2

  T(12,1,1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p1)*za(&
       p1,p5)*za(p3,p2)*s1(p1,p3)

  T(12,1,1,-1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p1)*za(&
       p1,p6)*za(p3,p2)*s1(p1,p3)

  T(12,1,1,1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p3)*zb(p6,p1)*za(p1,p5)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p1)*za(&
       p1,p5)*za(p4,p2)*s1(p1,p4)

  T(12,1,1,1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p1,p3)*zb(p5,p1)*za(p1,p6)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p1)*za(&
       p1,p6)*za(p4,p2)*s1(p1,p4)

  !-- T13 =
  T(13,-1,-1,-1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p2)*&
       za(p2,p5)*za(p3,p2)*s1(p1,p4)

  T(13,-1,-1,-1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p2)*&
       za(p2,p6)*za(p3,p2)*s1(p1,p4)

  T(13,-1,-1,1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,&
       p4))*zb(p1,p3)*zb(p6,p2)*za(p2,p5)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p2)*&
       za(p2,p5)*za(p4,p2)*s1(p1,p3)

  T(13,-1,-1,1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p3)*zb(p5,p2)*za(p2,p6)*za(p4,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p2)*za(&
       p2,p6)*za(p4,p2)*s1(p1,p3)

  T(13,-1,1,-1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p2,p3)*zb(p2,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p6,p2)*&
       za(p2,p5)*za(p3,p1)*za(p4,p1)

  T(13,-1,1,-1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)*zb(p2,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)**2*zb(p5,p2)*za(&
       p2,p6)*za(p3,p1)*za(p4,p1)

  T(13,-1,1,1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p6,p2)*za(p2,p5)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p6,&
       p2)*za(p2,p5)*za(p4,p1)**2

  T(13,-1,1,1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p2,p3)**2*zb(p5,p2)*za(p2,p6)*za(p3,p1)*za(p4,p1) - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p2,p4)*zb(p5,&
       p2)*za(p2,p6)*za(p4,p1)**2

  T(13,1,-1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p2)**2&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p4)**2*zb(p6,p2)*za(p2,p5)*za(p3,p2)*za(p4,p2)

  T(13,1,-1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)*zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p2)**2&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p4)**2*zb(p5,p2)*za(p2,p6)*za(p3,p2)*za(p4,p2)

  T(13,1,-1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(&
       s1(p3,p4))*zb(p1,p3)**2*zb(p6,p2)*za(p2,p5)*za(p3,p2)*za(p4,p2)&
        + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*&
       zb(p1,p3)*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p4,p2)**2

  T(13,1,-1,1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(&
       p3,p4))*zb(p1,p3)**2*zb(p5,p2)*za(p2,p6)*za(p3,p2)*za(p4,p2) + 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p5,p6))/(s1(p3,p4))*zb(p1,&
       p3)*zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p4,p2)**2

  T(13,1,1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4&
       ))*zb(p1,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p6,p2)*za(&
       p2,p5)*za(p3,p2)*s1(p1,p3)

  T(13,1,1,-1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)*s1(p2,p4) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p4)*zb(p5,p2)*za(&
       p2,p6)*za(p3,p2)*s1(p1,p3)

  T(13,1,1,1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4)&
       )*zb(p1,p3)*zb(p6,p2)*za(p2,p5)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p6,p2)*za(&
       p2,p5)*za(p4,p2)*s1(p1,p4)

  T(13,1,1,1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))&
       *zb(p1,p3)*zb(p5,p2)*za(p2,p6)*za(p4,p1)*s1(p2,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p5,p6))/(s1(p3,p4))*zb(p2,p3)*zb(p5,p2)*za(&
       p2,p6)*za(p4,p2)*s1(p1,p4)

  !-- T14 =
  T(14,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p6)*zb(p2,p3)*zb(p4,p1)*za(p1,p3)**2*za(p5,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p2,p4)*&
       za(p1,p3)*za(p5,p2)*s1(p4,p1)

  T(14,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p5)*zb(p2,p3)*zb(p4,p1)*za(p1,p3)**2*za(p6,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p2,p4)*&
       za(p1,p3)*za(p6,p2)*s1(p4,p1)

  T(14,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p6)*zb(p2,p3)*za(p1,p4)*za(p5,p2)*s1(p3,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p2,p4)*zb(&
       p3,p1)*za(p1,p4)**2*za(p5,p2)

  T(14,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p2,p3)*za(p1,p4)*za(p6,p2)*s1(p3,p1) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p2,p4)*zb(&
       p3,p1)*za(p1,p4)**2*za(p6,p2)

  T(14,-1,1,-1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p2,p6)*zb(p4,p1)*za(p1,p3)**2*za(p5,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p6)*za(p1,&
       p3)*za(p5,p1)*s1(p4,p1)

  T(14,-1,1,-1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p2,p5)*zb(p4,p1)*za(p1,p3)**2*za(p6,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p5)*za(p1,&
       p3)*za(p6,p1)*s1(p4,p1)

  T(14,-1,1,1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p2,p6)*za(p1,p4)*za(p5,p1)*s1(p3,p1) + 1._dp/4._dp/(&
       s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p6)*zb(p3,p1)&
       *za(p1,p4)**2*za(p5,p1)

  T(14,-1,1,1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p2,p5)*za(p1,p4)*za(p6,p1)*s1(p3,p1) + 1._dp/4._dp/(&
       s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p2,p5)*zb(p3,p1)&
       *za(p1,p4)**2*za(p6,p1)

  T(14,1,-1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p4)**2*zb(p1,p6)*za(p1,p3)*za(p4,p2)*za(p5,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p4,p1)*za(p3,p2)*za(p5,p2)*s1(p1,p3)

  T(14,1,-1,-1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p4)**2*zb(p1,p5)*za(p1,p3)*za(p4,p2)*za(p6,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p4,p1)*za(p3,p2)*za(p6,p2)*s1(p1,p3)

  T(14,1,-1,1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p1,p6)*za(p1,p4)*za(p3,p2)*za(p5,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p6)*zb(p3,p1)*za(p4,p2)*za(p5,p2)*s1(p1,p4)

  T(14,1,-1,1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**2*zb(p1,p5)*za(p1,p4)*za(p3,p2)*za(p6,p2)&
        - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p3,p1)*za(p4,p2)*za(p6,p2)*s1(p1,p4)

  T(14,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)**2*zb(p2,p6)*za(p1,p3)*za(p4,p2)*za(p5,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p4,p1)*za(&
       p3,p2)*za(p5,p1)*s1(p1,p3)

  T(14,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)**2*zb(p2,p5)*za(p1,p3)*za(p4,p2)*za(p6,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p4,p1)*za(&
       p3,p2)*za(p6,p1)*s1(p1,p3)

  T(14,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)**2*zb(p2,p6)*za(p1,p4)*za(p3,p2)*za(p5,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p3,p1)*za(&
       p4,p2)*za(p5,p1)*s1(p1,p4)

  T(14,1,1,1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)**2*zb(p2,p5)*za(p1,p4)*za(p3,p2)*za(p6,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p3,p1)*za(&
       p4,p2)*za(p6,p1)*s1(p1,p4)

  !-- T15 =
  T(15,-1,-1,-1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p6)*zb(p2,p4)**2*za(p2,p3)*za(p4,p1)*za(p5,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p4,p2)*&
       za(p3,p1)*za(p5,p2)*s1(p2,p3)

  T(15,-1,-1,-1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p5)*zb(p2,p4)**2*za(p2,p3)*za(p4,p1)*za(p6,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p4,p2)*&
       za(p3,p1)*za(p6,p2)*s1(p2,p3)

  T(15,-1,-1,1,-1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p6)*zb(p2,p3)**2*za(p2,p4)*za(p3,p1)*za(p5,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p6)*zb(p3,p2)*&
       za(p4,p1)*za(p5,p2)*s1(p2,p4)

  T(15,-1,-1,1,1) = 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p5)*zb(p2,p3)**2*za(p2,p4)*za(p3,p1)*za(p6,p2) + 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p5)*zb(p3,p2)*za(&
       p4,p1)*za(p6,p2)*s1(p2,p4)

  T(15,-1,1,-1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p4)**2*zb(p2,p6)*za(p2,p3)*za(p4,p1)*za(p5,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p4,p2)*za(p3,&
       p1)*za(p5,p1)*s1(p2,p3)

  T(15,-1,1,-1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p4)**2*zb(p2,p5)*za(p2,p3)*za(p4,p1)*za(p6,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p4,p2)*za(p3,&
       p1)*za(p6,p1)*s1(p2,p3)

  T(15,-1,1,1,-1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p2,p6)*za(p2,p4)*za(p3,p1)*za(p5,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p6)*zb(p3,p2)*za(p4,&
       p1)*za(p5,p1)*s1(p2,p4)

  T(15,-1,1,1,1) = 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p2,p5)*za(p2,p4)*za(p3,p1)*za(p6,p1) + 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p5)*zb(p3,p2)*za(p4,&
       p1)*za(p6,p1)*s1(p2,p4)

  T(15,1,-1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p6)*zb(p4,p2)*za(p2,p3)**2*za(p5,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p1,p6)*za(p2,p3)*za(p5,p2)*s1(p4,p2)

  T(15,1,-1,-1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p5)*zb(p4,p2)*za(p2,p3)**2*za(p6,p2&
       ) - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p1,p5)*za(p2,p3)*za(p6,p2)*s1(p4,p2)

  T(15,1,-1,1,-1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p6)*za(p2,p4)*za(p5,p2)*s1(p3,p2)&
        - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p1,p6)*zb(p3,p2)*za(p2,p4)**2*za(p5,p2)

  T(15,1,-1,1,1) =  - 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)*zb(p1,p5)*za(p2,p4)*za(p6,p2)*s1(p3,p2) - 1._dp&
       /4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p4)*zb(p1,p5)*zb(p3,p2)*za(p2,p4)**2*za(p6,p2)

  T(15,1,1,-1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)*zb(p2,p6)*zb(p4,p2)*za(p2,p3)**2*za(p5,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p6)*za(&
       p2,p3)*za(p5,p1)*s1(p4,p2)

  T(15,1,1,-1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)*zb(p2,p5)*zb(p4,p2)*za(p2,p3)**2*za(p6,p1) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p5)*za(&
       p2,p3)*za(p6,p1)*s1(p4,p2)

  T(15,1,1,1,-1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)*zb(p2,p6)*za(p2,p4)*za(p5,p1)*s1(p3,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p6)*zb(&
       p3,p2)*za(p2,p4)**2*za(p5,p1)

  T(15,1,1,1,1) = 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p3)*zb(p2,p5)*za(p2,p4)*za(p6,p1)*s1(p3,p2) + 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2,p5)*zb(&
       p3,p2)*za(p2,p4)**2*za(p6,p1)

  !-- T16 =
  T(16,-1,-1,-1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p2,p4)*zb(p4,p6)*za(p3,p2)*za(p4,p1)*za(p5,p3)&
        - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)*zb(p4,p6)*za(p3,p1)*za(p4,p2)*za(p5,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p6)*za(p5,p3)*s1(&
       p1,p3)*s1(p2,p3) - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5&
       ,p6))*zb(p4,p6)*za(p5,p3)*s1(p1,p4)*s1(p2,p4)

  T(16,-1,-1,-1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p2,p4)*zb(p4,p5)*za(p3,p2)*za(p4,p1)*za(p6,p3)&
        - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)*zb(p4,p5)*za(p3,p1)*za(p4,p2)*za(p6,p3) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p5)*za(p6,p3)*s1(&
       p1,p3)*s1(p2,p3) - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5&
       ,p6))*zb(p4,p5)*za(p6,p3)*s1(p1,p4)*s1(p2,p4)

  T(16,-1,-1,1,-1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p2,p4)*zb(p3,p6)*za(p3,p2)*za(p4,p1)*za(p5,p4)&
        - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)*zb(p3,p6)*za(p3,p1)*za(p4,p2)*za(p5,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p6)*za(p5,p4)*s1(&
       p1,p3)*s1(p2,p3) - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5&
       ,p6))*zb(p3,p6)*za(p5,p4)*s1(p1,p4)*s1(p2,p4)

  T(16,-1,-1,1,1) =  - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)*zb(p2,p4)*zb(p3,p5)*za(p3,p2)*za(p4,p1)*za(p6,p4)&
        - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)*zb(p3,p5)*za(p3,p1)*za(p4,p2)*za(p6,p4) - 1._dp/4._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p5)*za(p6,p4)*s1(&
       p1,p3)*s1(p2,p3) - 1._dp/4._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5&
       ,p6))*zb(p3,p5)*za(p6,p4)*s1(p1,p4)*s1(p2,p4)

  T(16,-1,1,-1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p2,p3)**2*zb(p4,p6)*za(p3,p1)**2*za(p5,p3) - 1._dp/2._dp/(&
       s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)*zb(p4,p6)&
       *za(p3,p1)*za(p4,p1)*za(p5,p3) - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,&
       p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p4,p6)*za(p4,p1)**2*za(p5,p3)

  T(16,-1,1,-1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)**2*zb(p4,p5)*za(p3,p1)**2*za(p6,p3) - 1._dp/2._dp/(s1(&
       p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)*zb(p4,p5)*&
       za(p3,p1)*za(p4,p1)*za(p6,p3) - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,&
       p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p4,p5)*za(p4,p1)**2*za(p6,p3)

  T(16,-1,1,1,-1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)**2*zb(p3,p6)*za(p3,p1)**2*za(p5,p4) - 1._dp/2._dp/(s1(&
       p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)*zb(p3,p6)*&
       za(p3,p1)*za(p4,p1)*za(p5,p4) - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,&
       p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p3,p6)*za(p4,p1)**2*za(p5,p4)

  T(16,-1,1,1,1) =  - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)**2*zb(p3,p5)*za(p3,p1)**2*za(p6,p4) - 1._dp/2._dp/(s1(&
       p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)*zb(p3,p5)*&
       za(p3,p1)*za(p4,p1)*za(p6,p4) - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,&
       p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p3,p5)*za(p4,p1)**2*za(p6,p4)

  T(16,1,-1,-1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**2*zb(p4,p6)*za(p3,p2)**2*za(p5,p3) + 1._dp/&
       2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p3)*&
       zb(p1,p4)*zb(p4,p6)*za(p3,p2)*za(p4,p2)*za(p5,p3) + 1._dp/4._dp&
       /(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)**2*&
       zb(p4,p6)*za(p4,p2)**2*za(p5,p3)

  T(16,1,-1,-1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**2*zb(p4,p5)*za(p3,p2)**2*za(p6,p3) + 1._dp/&
       2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p3)*&
       zb(p1,p4)*zb(p4,p5)*za(p3,p2)*za(p4,p2)*za(p6,p3) + 1._dp/4._dp&
       /(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)**2*&
       zb(p4,p5)*za(p4,p2)**2*za(p6,p3)

  T(16,1,-1,1,-1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**2*zb(p3,p6)*za(p3,p2)**2*za(p5,p4) + 1._dp/&
       2._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p3)*&
       zb(p1,p4)*zb(p3,p6)*za(p3,p2)*za(p4,p2)*za(p5,p4) + 1._dp/4._dp&
       /(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)**2*&
       zb(p3,p6)*za(p4,p2)**2*za(p5,p4)

  T(16,1,-1,1,1) = 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p3)**2*zb(p3,p5)*za(p3,p2)**2*za(p6,p4) + 1._dp/2._dp&
       /(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p3)*zb(&
       p1,p4)*zb(p3,p5)*za(p3,p2)*za(p4,p2)*za(p6,p4) + 1._dp/4._dp/(&
       za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)**2*zb(&
       p3,p5)*za(p4,p2)**2*za(p6,p4)

  T(16,1,1,-1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)*zb(p2,p4)*zb(p4,p6)*za(p3,p2)*za(p4,p1)*za(p5,p3)&
        - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)*zb(p4,p6)*za(p3,p1)*za(p4,p2)*za(p5,p3) - 1._dp/4._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p6)*za(p5,p3)*s1(&
       p1,p3)*s1(p2,p3) - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5&
       ,p6))*zb(p4,p6)*za(p5,p3)*s1(p1,p4)*s1(p2,p4)

  T(16,1,1,-1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p2,p4)*zb(p4,p5)*za(p3,p2)*za(p4,p1)*za(p6,p3) - &
       1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(&
       p2,p3)*zb(p4,p5)*za(p3,p1)*za(p4,p2)*za(p6,p3) - 1._dp/4._dp/(&
       za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p5)*za(p6,p3)*s1(p1,&
       p3)*s1(p2,p3) - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p4,p5)*za(p6,p3)*s1(p1,p4)*s1(p2,p4)

  T(16,1,1,1,-1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p2,p4)*zb(p3,p6)*za(p3,p2)*za(p4,p1)*za(p5,p4) - &
       1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(&
       p2,p3)*zb(p3,p6)*za(p3,p1)*za(p4,p2)*za(p5,p4) - 1._dp/4._dp/(&
       za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p6)*za(p5,p4)*s1(p1,&
       p3)*s1(p2,p3) - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p3,p6)*za(p5,p4)*s1(p1,p4)*s1(p2,p4)

  T(16,1,1,1,1) =  - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)*zb(p2,p4)*zb(p3,p5)*za(p3,p2)*za(p4,p1)*za(p6,p4) - 1._dp&
       /4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2&
       ,p3)*zb(p3,p5)*za(p3,p1)*za(p4,p2)*za(p6,p4) - 1._dp/4._dp/(za(&
       p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p5)*za(p6,p4)*s1(p1,p3)&
       *s1(p2,p3) - 1._dp/4._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p3,p5)*za(p6,p4)*s1(p1,p4)*s1(p2,p4)

  !-- T17 =
  T(17,-1,-1,-1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p4)**2*zb(p2,p3)*zb(p6,p1)*za(p1,p3)**2*za(p1,p5)*za(&
       p4,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p1)&
       *za(p1,p3)*za(p1,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p1)*za(p1,p3)*za(&
       p1,p5)*s1(p1,p4)*s1(p2,p4)

  T(17,-1,-1,-1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p4)**2*zb(p2,p3)*zb(p5,p1)*za(p1,p3)**2*za(p1,p6)*za(&
       p4,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p1)&
       *za(p1,p3)*za(p1,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p1)*za(p1,p3)*za(&
       p1,p6)*s1(p1,p4)*s1(p2,p4)

  T(17,-1,-1,1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)**2*zb(p2,p4)*zb(p6,p1)*za(p1,p4)**2*za(p1,p5)*za(&
       p3,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p6,p1)*za(p1,p5)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p1)&
       *za(p1,p4)*za(p1,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p1)*za(p1,p4)*za(&
       p1,p5)*s1(p1,p4)*s1(p2,p4)

  T(17,-1,-1,1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)**2*zb(p2,p4)*zb(p5,p1)*za(p1,p4)**2*za(p1,p6)*za(p3&
       ,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p3)*zb(p5,p1)*za(p1,p6)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/8._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p1)*&
       za(p1,p4)*za(p1,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p1)*za(p1,p4)*za(p1&
       ,p6)*s1(p1,p4)*s1(p2,p4)

  T(17,-1,1,-1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)**2*zb(p4,p1)*zb(p6,p1)*za(p1,p3)**2*za(p1,p5)*za(p3,p1&
       ) - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*&
       zb(p2,p4)*zb(p6,p1)*za(p1,p3)**2*za(p1,p5)*s1(p4,p1) + 1._dp/8._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p6,p1)*za(&
       p1,p3)*za(p1,p5)*za(p4,p1)*s1(p4,p1)

  T(17,-1,1,-1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p4,p1)*zb(p5,p1)*za(p1,p3)**2*za(p1,p6)*za(p3,p1)&
        - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(&
       p2,p4)*zb(p5,p1)*za(p1,p3)**2*za(p1,p6)*s1(p4,p1) + 1._dp/8._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p5,p1)*za(&
       p1,p3)*za(p1,p6)*za(p4,p1)*s1(p4,p1)

  T(17,-1,1,1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p6,p1)*za(p1,p4)*za(p1,p5)*za(p3,p1)*s1(p3,p1) - 1._dp&
       /4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)&
       *zb(p6,p1)*za(p1,p4)**2*za(p1,p5)*s1(p3,p1) + 1._dp/8._dp/(s1(p1&
       ,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p3,p1)*zb(p6,p1)*&
       za(p1,p4)**2*za(p1,p5)*za(p4,p1)

  T(17,-1,1,1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p5,p1)*za(p1,p4)*za(p1,p6)*za(p3,p1)*s1(p3,p1) - 1._dp&
       /4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)&
       *zb(p5,p1)*za(p1,p4)**2*za(p1,p6)*s1(p3,p1) + 1._dp/8._dp/(s1(p1&
       ,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p3,p1)*zb(p5,p1)*&
       za(p1,p4)**2*za(p1,p6)*za(p4,p1)

  T(17,1,-1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p4,p1)*zb(p6,p1)*za(p1,p5)*za(p3,p2)**&
       2*s1(p1,p3) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**3*zb(p6,p1)*za(p1,p3)*za(p1,p5)*za(p4,p2)&
       **2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)**2*zb(p6,p1)*za(p1,p5)*za(p3,p2)*za(p4,p2)*s1(p1,p3)

  T(17,1,-1,-1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p4,p1)*zb(p5,p1)*za(p1,p6)*za(p3,p2)**&
       2*s1(p1,p3) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**3*zb(p5,p1)*za(p1,p3)*za(p1,p6)*za(p4,p2)&
       **2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)**2*zb(p5,p1)*za(p1,p6)*za(p3,p2)*za(p4,p2)*s1(p1,p3)

  T(17,1,-1,1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**3*zb(p6,p1)*za(p1,p4)*za(p1,p5)*za(p3,p2&
       )**2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)**2*zb(p6,p1)*za(p1,p5)*za(p3,p2)*za(p4,p2)*s1(p1,p4&
       ) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p3,p1)*zb(p6,p1)*za(p1,p5)*za(p4,p2)**2*s1(p1,p4)

  T(17,1,-1,1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**3*zb(p5,p1)*za(p1,p4)*za(p1,p6)*za(p3,p2)&
       **2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)**2*zb(p5,p1)*za(p1,p6)*za(p3,p2)*za(p4,p2)*s1(p1,p4)&
        - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p3,p1)*zb(p5,p1)*za(p1,p6)*za(p4,p2)**2*s1(p1,p4)

  T(17,1,1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p4)**2*zb(p2,p3)*zb(p6,p1)*za(p1,p3)**2*za(p1,p5)*za(p4&
       ,p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p4)*zb(p6,p1)*za(p1,p5)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p1)*&
       za(p1,p3)*za(p1,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p1)*za(p1,p3)*za(p1&
       ,p5)*s1(p1,p4)*s1(p2,p4)

  T(17,1,1,-1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)**2*zb(p2,p3)*zb(p5,p1)*za(p1,p3)**2*za(p1,p6)*za(p4,&
       p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p4)*zb(p5,p1)*za(p1,p6)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p1)*&
       za(p1,p3)*za(p1,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p1)*za(p1,p3)*za(p1&
       ,p6)*s1(p1,p4)*s1(p2,p4)

  T(17,1,1,1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)**2*zb(p2,p4)*zb(p6,p1)*za(p1,p4)**2*za(p1,p5)*za(p3,&
       p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p3)*zb(p6,p1)*za(p1,p5)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p1)*&
       za(p1,p4)*za(p1,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p1)*za(p1,p4)*za(p1&
       ,p5)*s1(p1,p4)*s1(p2,p4)

  T(17,1,1,1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)**2*zb(p2,p4)*zb(p5,p1)*za(p1,p4)**2*za(p1,p6)*za(p3,&
       p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p3)*zb(p5,p1)*za(p1,p6)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p1)*&
       za(p1,p4)*za(p1,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p1)*za(p1,p4)*za(p1&
       ,p6)*s1(p1,p4)*s1(p2,p4)

  !-- T18 =
  T(18,-1,-1,-1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p4)**2*zb(p2,p3)*zb(p6,p2)*za(p1,p3)**2*za(p2,p5)*za(&
       p4,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p2)&
       *za(p1,p3)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p2)*za(p1,p3)*za(&
       p2,p5)*s1(p1,p4)*s1(p2,p4)

  T(18,-1,-1,-1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p4)**2*zb(p2,p3)*zb(p5,p2)*za(p1,p3)**2*za(p2,p6)*za(&
       p4,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p2)&
       *za(p1,p3)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p2)*za(p1,p3)*za(&
       p2,p6)*s1(p1,p4)*s1(p2,p4)

  T(18,-1,-1,1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)**2*zb(p2,p4)*zb(p6,p2)*za(p1,p4)**2*za(p2,p5)*za(&
       p3,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p6,p2)*za(p2,p5)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p2)&
       *za(p1,p4)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p2)*za(p1,p4)*za(&
       p2,p5)*s1(p1,p4)*s1(p2,p4)

  T(18,-1,-1,1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)**2*zb(p2,p4)*zb(p5,p2)*za(p1,p4)**2*za(p2,p6)*za(p3&
       ,p2) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p3)*zb(p5,p2)*za(p2,p6)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/8._dp&
       /(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p2)*&
       za(p1,p4)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p2)*za(p1,p4)*za(p2&
       ,p6)*s1(p1,p4)*s1(p2,p4)

  T(18,-1,1,-1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)**2*zb(p4,p1)*zb(p6,p2)*za(p1,p3)**2*za(p2,p5)*za(p3,p1&
       ) - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*&
       zb(p2,p4)*zb(p6,p2)*za(p1,p3)**2*za(p2,p5)*s1(p4,p1) + 1._dp/8._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p6,p2)*za(&
       p1,p3)*za(p2,p5)*za(p4,p1)*s1(p4,p1)

  T(18,-1,1,-1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p4,p1)*zb(p5,p2)*za(p1,p3)**2*za(p2,p6)*za(p3,p1)&
        - 1._dp/4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(&
       p2,p4)*zb(p5,p2)*za(p1,p3)**2*za(p2,p6)*s1(p4,p1) + 1._dp/8._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p5,p2)*za(&
       p1,p3)*za(p2,p6)*za(p4,p1)*s1(p4,p1)

  T(18,-1,1,1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p6,p2)*za(p1,p4)*za(p2,p5)*za(p3,p1)*s1(p3,p1) - 1._dp&
       /4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)&
       *zb(p6,p2)*za(p1,p4)**2*za(p2,p5)*s1(p3,p1) + 1._dp/8._dp/(s1(p1&
       ,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p3,p1)*zb(p6,p2)*&
       za(p1,p4)**2*za(p2,p5)*za(p4,p1)

  T(18,-1,1,1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**2*zb(p5,p2)*za(p1,p4)*za(p2,p6)*za(p3,p1)*s1(p3,p1) - 1._dp&
       /4._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)*zb(p2,p4)&
       *zb(p5,p2)*za(p1,p4)**2*za(p2,p6)*s1(p3,p1) + 1._dp/8._dp/(s1(p1&
       ,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p3,p1)*zb(p5,p2)*&
       za(p1,p4)**2*za(p2,p6)*za(p4,p1)

  T(18,1,-1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p4,p1)*zb(p6,p2)*za(p2,p5)*za(p3,p2)**&
       2*s1(p1,p3) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**3*zb(p6,p2)*za(p1,p3)*za(p2,p5)*za(p4,p2)&
       **2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)**2*zb(p6,p2)*za(p2,p5)*za(p3,p2)*za(p4,p2)*s1(p1,p3)

  T(18,1,-1,-1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p4,p1)*zb(p5,p2)*za(p2,p6)*za(p3,p2)**&
       2*s1(p1,p3) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**3*zb(p5,p2)*za(p1,p3)*za(p2,p6)*za(p4,p2)&
       **2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)**2*zb(p5,p2)*za(p2,p6)*za(p3,p2)*za(p4,p2)*s1(p1,p3)

  T(18,1,-1,1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**3*zb(p6,p2)*za(p1,p4)*za(p2,p5)*za(p3,p2&
       )**2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)**2*zb(p6,p2)*za(p2,p5)*za(p3,p2)*za(p4,p2)*s1(p1,p4&
       ) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p3,p1)*zb(p6,p2)*za(p2,p5)*za(p4,p2)**2*s1(p1,p4)

  T(18,1,-1,1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**3*zb(p5,p2)*za(p1,p4)*za(p2,p6)*za(p3,p2)&
       **2 + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)**2*zb(p5,p2)*za(p2,p6)*za(p3,p2)*za(p4,p2)*s1(p1,p4)&
        - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p1,p4)*zb(p3,p1)*zb(p5,p2)*za(p2,p6)*za(p4,p2)**2*s1(p1,p4)

  T(18,1,1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p4)**2*zb(p2,p3)*zb(p6,p2)*za(p1,p3)**2*za(p2,p5)*za(p4&
       ,p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p4)*zb(p6,p2)*za(p2,p5)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p2)*&
       za(p1,p3)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p6,p2)*za(p1,p3)*za(p2&
       ,p5)*s1(p1,p4)*s1(p2,p4)

  T(18,1,1,-1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p4)**2*zb(p2,p3)*zb(p5,p2)*za(p1,p3)**2*za(p2,p6)*za(p4,&
       p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p4)*zb(p5,p2)*za(p2,p6)*za(p3,p2)*s1(p1,p3)*s1(p4,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p2)*&
       za(p1,p3)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p1)*zb(p5,p2)*za(p1,p3)*za(p2&
       ,p6)*s1(p1,p4)*s1(p2,p4)

  T(18,1,1,1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)**2*zb(p2,p4)*zb(p6,p2)*za(p1,p4)**2*za(p2,p5)*za(p3,&
       p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p3)*zb(p6,p2)*za(p2,p5)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p2)*&
       za(p1,p4)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p6,p2)*za(p1,p4)*za(p2&
       ,p5)*s1(p1,p4)*s1(p2,p4)

  T(18,1,1,1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)**2*zb(p2,p4)*zb(p5,p2)*za(p1,p4)**2*za(p2,p6)*za(p3,&
       p2) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p2,&
       p3)*zb(p5,p2)*za(p2,p6)*za(p4,p2)*s1(p1,p4)*s1(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p2)*&
       za(p1,p4)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p1)*zb(p5,p2)*za(p1,p4)*za(p2&
       ,p6)*s1(p1,p4)*s1(p2,p4)

  !-- T19 =
  T(19,-1,-1,-1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p2,p4)**2*zb(p6,p1)*za(p1,p5)*za(p2,p3)**2*za(&
       p4,p1) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p1)&
       *za(p1,p5)*za(p2,p3)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p1)*za(p1,p5)*za(&
       p2,p3)*s1(p1,p4)*s1(p2,p4)

  T(19,-1,-1,-1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p2,p4)**2*zb(p5,p1)*za(p1,p6)*za(p2,p3)**2*za(&
       p4,p1) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p1)&
       *za(p1,p6)*za(p2,p3)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p1)*za(p1,p6)*za(&
       p2,p3)*s1(p1,p4)*s1(p2,p4)

  T(19,-1,-1,1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p6,p1)*za(p1,p5)*za(p4,p1)*s1(p2,p4)*s1(p3,p2)&
        - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)**2*zb(p6,p1)*za(p1,p5)*za(p2,p4)**2*za(p3,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p1)&
       *za(p1,p5)*za(p2,p4)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p1)*za(p1,p5)*za(&
       p2,p4)*s1(p1,p4)*s1(p2,p4)

  T(19,-1,-1,1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)*zb(p5,p1)*za(p1,p6)*za(p4,p1)*s1(p2,p4)*s1(p3,p2)&
        - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)**2*zb(p5,p1)*za(p1,p6)*za(p2,p4)**2*za(p3,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p1)&
       *za(p1,p6)*za(p2,p4)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p1)*za(p1,p6)*za(&
       p2,p4)*s1(p1,p4)*s1(p2,p4)

  T(19,-1,1,-1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p4,p2)*zb(p6,p1)*za(p1,p5)*za(p3,p1)**2*s1(p2,p3)&
        + 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**3*&
       zb(p6,p1)*za(p1,p5)*za(p2,p3)*za(p4,p1)**2 - 1._dp/4._dp/(s1(p1,&
       p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p6,p1)*za(p1,p5)*&
       za(p3,p1)*za(p4,p1)*s1(p2,p3)

  T(19,-1,1,-1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p4,p2)*zb(p5,p1)*za(p1,p6)*za(p3,p1)**2*s1(p2,p3) + 1._dp&
       /8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**3*zb(p5,&
       p1)*za(p1,p6)*za(p2,p3)*za(p4,p1)**2 - 1._dp/4._dp/(s1(p1,p2))/(&
       s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p5,p1)*za(p1,p6)*za(p3,p1&
       )*za(p4,p1)*s1(p2,p3)

  T(19,-1,1,1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**3*zb(p6,p1)*za(p1,p5)*za(p2,p4)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)**2*zb(p6,p1)*za(&
       p1,p5)*za(p3,p1)*za(p4,p1)*s1(p2,p4) + 1._dp/8._dp/(s1(p1,p2))/(&
       s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p3,p2)*zb(p6,p1)*za(p1,p5)*&
       za(p4,p1)**2*s1(p2,p4)

  T(19,-1,1,1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**3*zb(p5,p1)*za(p1,p6)*za(p2,p4)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)**2*zb(p5,p1)*za(&
       p1,p6)*za(p3,p1)*za(p4,p1)*s1(p2,p4) + 1._dp/8._dp/(s1(p1,p2))/(&
       s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p3,p2)*zb(p5,p1)*za(p1,p6)*&
       za(p4,p1)**2*s1(p2,p4)

  T(19,1,-1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p4,p2)*zb(p6,p1)*za(p1,p5)*za(p2,p3&
       )**2*za(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p2,p3)**&
       2*s1(p4,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**2*zb(p6,p1)*za(p1,p5)*za(p2,p3)*za(p4,p2)*&
       s1(p4,p2)

  T(19,1,-1,-1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p4,p2)*zb(p5,p1)*za(p1,p6)*za(p2,p3&
       )**2*za(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p2,p3)**&
       2*s1(p4,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**2*zb(p5,p1)*za(p1,p6)*za(p2,p3)*za(p4,p2)*&
       s1(p4,p2)

  T(19,1,-1,1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p6,p1)*za(p1,p5)*za(p2,p4)*za(p3,p2&
       )*s1(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p6,p1)*za(p1,p5)*za(p2,p4)**2*&
       s1(p3,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p4)**2*zb(p3,p2)*zb(p6,p1)*za(p1,p5)*za(p2,p4)**2*&
       za(p4,p2)

  T(19,1,-1,1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**2*zb(p5,p1)*za(p1,p6)*za(p2,p4)*za(p3,p2)*&
       s1(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p5,p1)*za(p1,p6)*za(p2,p4)**2*s1(&
       p3,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p4)**2*zb(p3,p2)*zb(p5,p1)*za(p1,p6)*za(p2,p4)**2*za(&
       p4,p2)

  T(19,1,1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)*zb(p2,p4)**2*zb(p6,p1)*za(p1,p5)*za(p2,p3)**2*za(p4&
       ,p1) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p4)*zb(p6,p1)*za(p1,p5)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p1)*&
       za(p1,p5)*za(p2,p3)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p1)*za(p1,p5)*za(p2&
       ,p3)*s1(p1,p4)*s1(p2,p4)

  T(19,1,1,-1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p2,p4)**2*zb(p5,p1)*za(p1,p6)*za(p2,p3)**2*za(p4,&
       p1) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p4)*zb(p5,p1)*za(p1,p6)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p1)*&
       za(p1,p6)*za(p2,p3)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p1)*za(p1,p6)*za(p2&
       ,p3)*s1(p1,p4)*s1(p2,p4)

  T(19,1,1,1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p6,p1)*za(p1,p5)*za(p4,p1)*s1(p2,p4)*s1(p3,p2) - &
       1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(&
       p2,p3)**2*zb(p6,p1)*za(p1,p5)*za(p2,p4)**2*za(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p1)*za(&
       p1,p5)*za(p2,p4)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)**2&
       )/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p1)*za(p1,p5)*za(p2,p4&
       )*s1(p1,p4)*s1(p2,p4)

  T(19,1,1,1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)*zb(p5,p1)*za(p1,p6)*za(p4,p1)*s1(p2,p4)*s1(p3,p2) - 1._dp&
       /8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2&
       ,p3)**2*zb(p5,p1)*za(p1,p6)*za(p2,p4)**2*za(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p1)*za(&
       p1,p6)*za(p2,p4)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)**2&
       )/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p1)*za(p1,p6)*za(p2,p4&
       )*s1(p1,p4)*s1(p2,p4)

  !-- T20 =
  T(20,-1,-1,-1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p2,p4)**2*zb(p6,p2)*za(p2,p3)**2*za(p2,p5)*za(&
       p4,p1) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p2)&
       *za(p2,p3)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p2)*za(p2,p3)*za(&
       p2,p5)*s1(p1,p4)*s1(p2,p4)

  T(20,-1,-1,-1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p2,p4)**2*zb(p5,p2)*za(p2,p3)**2*za(p2,p6)*za(&
       p4,p1) - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p1,p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p2)&
       *za(p2,p3)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p2)*za(p2,p3)*za(&
       p2,p6)*s1(p1,p4)*s1(p2,p4)

  T(20,-1,-1,1,-1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p3)*zb(p6,p2)*za(p2,p5)*za(p4,p1)*s1(p2,p4)*s1(p3,p2)&
        - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)**2*zb(p6,p2)*za(p2,p4)**2*za(p2,p5)*za(p3,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p2)&
       *za(p2,p4)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p2)*za(p2,p4)*za(&
       p2,p5)*s1(p1,p4)*s1(p2,p4)

  T(20,-1,-1,1,1) =  - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)*zb(p5,p2)*za(p2,p6)*za(p4,p1)*s1(p2,p4)*s1(p3,p2)&
        - 1._dp/8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*&
       zb(p2,p3)**2*zb(p5,p2)*za(p2,p4)**2*za(p2,p6)*za(p3,p1) - 1._dp/&
       8._dp/(zb(p1,p2)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p2)&
       *za(p2,p4)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(zb(p1,p2&
       )**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p2)*za(p2,p4)*za(&
       p2,p6)*s1(p1,p4)*s1(p2,p4)

  T(20,-1,1,-1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*&
       zb(p2,p3)*zb(p4,p2)*zb(p6,p2)*za(p2,p5)*za(p3,p1)**2*s1(p2,p3)&
        + 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**3*&
       zb(p6,p2)*za(p2,p3)*za(p2,p5)*za(p4,p1)**2 - 1._dp/4._dp/(s1(p1,&
       p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p6,p2)*za(p2,p5)*&
       za(p3,p1)*za(p4,p1)*s1(p2,p3)

  T(20,-1,1,-1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)*zb(p4,p2)*zb(p5,p2)*za(p2,p6)*za(p3,p1)**2*s1(p2,p3) + 1._dp&
       /8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**3*zb(p5,&
       p2)*za(p2,p3)*za(p2,p6)*za(p4,p1)**2 - 1._dp/4._dp/(s1(p1,p2))/(&
       s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)**2*zb(p5,p2)*za(p2,p6)*za(p3,p1&
       )*za(p4,p1)*s1(p2,p3)

  T(20,-1,1,1,-1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**3*zb(p6,p2)*za(p2,p4)*za(p2,p5)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)**2*zb(p6,p2)*za(&
       p2,p5)*za(p3,p1)*za(p4,p1)*s1(p2,p4) + 1._dp/8._dp/(s1(p1,p2))/(&
       s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p3,p2)*zb(p6,p2)*za(p2,p5)*&
       za(p4,p1)**2*s1(p2,p4)

  T(20,-1,1,1,1) = 1._dp/8._dp/(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(&
       p2,p3)**3*zb(p5,p2)*za(p2,p4)*za(p2,p6)*za(p3,p1)**2 - 1._dp/4._dp&
       /(s1(p1,p2))/(s1(p3,p4))/(s1(p5,p6))*zb(p2,p3)**2*zb(p5,p2)*za(&
       p2,p6)*za(p3,p1)*za(p4,p1)*s1(p2,p4) + 1._dp/8._dp/(s1(p1,p2))/(&
       s1(p3,p4))/(s1(p5,p6))*zb(p2,p4)*zb(p3,p2)*zb(p5,p2)*za(p2,p6)*&
       za(p4,p1)**2*s1(p2,p4)

  T(20,1,-1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p4,p2)*zb(p6,p2)*za(p2,p3)**2*za(p2&
       ,p5)*za(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p6,p2)*za(p2,p3)**2*za(p2,p5&
       )*s1(p4,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**2*zb(p6,p2)*za(p2,p3)*za(p2,p5)*za(p4,p2)*&
       s1(p4,p2)

  T(20,1,-1,-1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p4,p2)*zb(p5,p2)*za(p2,p3)**2*za(p2&
       ,p6)*za(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p5,p2)*za(p2,p3)**2*za(p2,p6&
       )*s1(p4,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p4)**2*zb(p5,p2)*za(p2,p3)*za(p2,p6)*za(p4,p2)*&
       s1(p4,p2)

  T(20,1,-1,1,-1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))&
       /(s1(p5,p6))*zb(p1,p3)**2*zb(p6,p2)*za(p2,p4)*za(p2,p5)*za(p3,p2&
       )*s1(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p6,p2)*za(p2,p4)**2*za(p2,p5)*&
       s1(p3,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p4)**2*zb(p3,p2)*zb(p6,p2)*za(p2,p4)**2*za(p2,p5)*&
       za(p4,p2)

  T(20,1,-1,1,1) =  - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(&
       s1(p5,p6))*zb(p1,p3)**2*zb(p5,p2)*za(p2,p4)*za(p2,p6)*za(p3,p2)*&
       s1(p3,p2) + 1._dp/4._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(&
       p5,p6))*zb(p1,p3)*zb(p1,p4)*zb(p5,p2)*za(p2,p4)**2*za(p2,p6)*s1(&
       p3,p2) - 1._dp/8._dp/(za(p2,p1))/(zb(p2,p1))/(s1(p3,p4))/(s1(p5,&
       p6))*zb(p1,p4)**2*zb(p3,p2)*zb(p5,p2)*za(p2,p4)**2*za(p2,p6)*za(&
       p4,p2)

  T(20,1,1,-1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6&
       ))*zb(p1,p3)*zb(p2,p4)**2*zb(p6,p2)*za(p2,p3)**2*za(p2,p5)*za(p4&
       ,p1) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p4)*zb(p6,p2)*za(p2,p5)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p2)*&
       za(p2,p3)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p6,p2)*za(p2,p3)*za(p2&
       ,p5)*s1(p1,p4)*s1(p2,p4)

  T(20,1,1,-1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p2,p4)**2*zb(p5,p2)*za(p2,p3)**2*za(p2,p6)*za(p4,&
       p1) - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,&
       p4)*zb(p5,p2)*za(p2,p6)*za(p3,p1)*s1(p2,p3)*s1(p4,p2) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p2)*&
       za(p2,p3)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)&
       **2)/(s1(p3,p4))/(s1(p5,p6))*zb(p4,p2)*zb(p5,p2)*za(p2,p3)*za(p2&
       ,p6)*s1(p1,p4)*s1(p2,p4)

  T(20,1,1,1,-1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6)&
       )*zb(p1,p3)*zb(p6,p2)*za(p2,p5)*za(p4,p1)*s1(p2,p4)*s1(p3,p2) - &
       1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(&
       p2,p3)**2*zb(p6,p2)*za(p2,p4)**2*za(p2,p5)*za(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p2)*za(&
       p2,p4)*za(p2,p5)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)**2&
       )/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p6,p2)*za(p2,p4)*za(p2,p5&
       )*s1(p1,p4)*s1(p2,p4)

  T(20,1,1,1,1) =  - 1._dp/8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))&
       *zb(p1,p3)*zb(p5,p2)*za(p2,p6)*za(p4,p1)*s1(p2,p4)*s1(p3,p2) - 1._dp&
       /8._dp/(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p1,p4)*zb(p2&
       ,p3)**2*zb(p5,p2)*za(p2,p4)**2*za(p2,p6)*za(p3,p1) - 1._dp/8._dp&
       /(za(p2,p1)**2)/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p2)*za(&
       p2,p4)*za(p2,p6)*s1(p1,p3)*s1(p2,p3) - 1._dp/8._dp/(za(p2,p1)**2&
       )/(s1(p3,p4))/(s1(p5,p6))*zb(p3,p2)*zb(p5,p2)*za(p2,p4)*za(p2,p6&
       )*s1(p1,p4)*s1(p2,p4)

end subroutine heli_zz_heft
