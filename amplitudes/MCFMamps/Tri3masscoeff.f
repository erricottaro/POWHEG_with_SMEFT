      function Tri3masscoeff(j1,j2,j3,j4,j5,j6,za,zb,sprod) 
      use mod_types; use mod_consts_dp
C-----Tri3massmass calculates the coefficient of the 3 mass triangle,
C-----in the massless limit using the formula of BDK (11.6) and (11.7)
      implicit none
      integer j1,j2,j3,j4,j5,j6
      complex(dp)  :: Tri3masscoeff
      complex(dp) :: Fvfbit1,Fvsbit1,Fvf,Fvs
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)     :: sprod(6,6)

c---flip2:( j4<-->j3, j2<-->j1, j5<-->j6, za<-->zb)
      Fvf=
     & +Fvfbit1(j1,j2,j3,j4,j5,j6,za,zb,sprod) 
     & +Fvfbit1(j2,j1,j4,j3,j6,j5,zb,za,sprod) 
      
      Fvs=
     & +Fvsbit1(j1,j2,j3,j4,j5,j6,za,zb,sprod) 
     & +Fvsbit1(j2,j1,j4,j3,j6,j5,zb,za,sprod) 


      Tri3masscoeff=(Fvf+Fvs)


      return
      end 

      function Fvfbit1(j1,j2,j3,j4,j5,j6,za,zb,sprod) 
      use mod_types; use mod_consts_dp
      implicit none
      integer j1,j2,j3,j4,j5,j6
      complex(dp) ::  Fvfbit0,zab2
      complex(dp) :: Fvfbit1
      real(dp) :: t,ss,tt,m1sq,m2sq,Lsm1_2mhbit
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)

C---This performs the flip3 on the explicit I3m piece.
C---and adds in the piece from the Lsm1_2mh
c---flip3:( j4<-->j5, j3<-->j6, j2<-->j1, za<-->zb)

C---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      t(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
      Lsm1_2mhbit(ss,tt,m1sq,m2sq)=0.5d0*(ss-m1sq-m2sq)+m1sq*m2sq/tt
C---end statement functions

      Fvfbit1=
     & +Fvfbit0(j1,j2,j3,j4,j5,j6,za,zb,sprod) 
     & +Fvfbit0(j2,j1,j6,j5,j4,j3,zb,za,sprod) 
     & +zab2(j3,j4,j2,j6)**2/(za(j4,j3)*zb(j5,j6)*zab2(j2,j3,j4,j1)**2)
     & *Lsm1_2mhbit(sprod(j1,j2),t(j2,j3,j4),sprod(j3,j4),sprod(j5,j6))
      return
      end
      


      function Fvsbit1(j1,j2,j3,j4,j5,j6,za,zb,sprod) 
C---This is the whole expression in Eq.(11.6) subject to exchange flip_2
c---ie  flip2:( j1<-->j2, j3<-->j4, j5<-->j6, za<-->zb)
      use mod_types; use mod_consts_dp
      implicit none
      integer j1,j2,j3,j4,j5,j6
      complex(dp) :: zab2
      complex(dp) :: Fvsbit1
      real(dp)    ::  t,delta,IDelta,ss,tt,m1sq,m2sq,Lsm1_2mhbit
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)

C---statement functions
      t(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
      delta(j1,j2,j3,j4,j5,j6)=sprod(j1,j2)-sprod(j3,j4)-sprod(j5,j6)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      Lsm1_2mhbit(ss,tt,m1sq,m2sq)=0.5d0*(ss-m1sq-m2sq)+m1sq*m2sq/tt
C---end statement functions

      IDelta=1d0/(sprod(j3,j4)**2+sprod(j1,j2)**2+sprod(j5,j6)**2
     & -2d0*(+sprod(j3,j4)*sprod(j1,j2)+sprod(j3,j4)*sprod(j5,j6)
     & +sprod(j5,j6)*sprod(j1,j2)))

      Fvsbit1=
     & +2d0*za(j3,j2)*zb(j1,j6)*zab2(j2,j4,j3,j6)
     & *zab2(j3,j5,j6,j1)*t(j2,j3,j4)
     & /(za(j4,j3)*zb(j5,j6)*zab2(j2,j3,j4,j1)**4)
     & *Lsm1_2mhbit(sprod(j1,j2),t(j2,j3,j4),sprod(j3,j4),sprod(j5,j6))

     &-(2d0*za(j3,j2)*zb(j1,j6)*(zb(j4,j3)*za(j5,j6)*za(j3,j2)*zb(j1,j6)
     & -zab2(j5,j3,j1,j4)*zab2(j2,j3,j4,j1))*zab2(j1,j3,j4,j2)
     & *(t(j2,j3,j4)-t(j1,j3,j4))/zab2(j2,j3,j4,j1)**3*IDelta

     & -3d0*(sprod(j1,j2)*delta(j2,j1,j4,j3,j5,j6)*zab2(j3,j2,j1,j4)
     &  *zab2(j5,j2,j1,j6)*IDelta
     & -za(j3,j2)*zb(j2,j4)*za(j5,j1)*zb(j1,j6)
     & -za(j3,j1)*zb(j1,j4)*za(j5,j2)*zb(j2,j6))
     & *zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)*IDelta

     & -zb(j4,j1)*za(j2,j5)*za(j3,j2)*zb(j1,j6)
     & *zab2(j1,j3,j4,j2)**2/zab2(j2,j3,j4,j1)**2*IDelta
     & -zb(j4,j2)*za(j1,j5)*za(j3,j1)*zb(j2,j6)*IDelta)

       return
       end
       
       
      function Fvfbit0(j1,j2,j3,j4,j5,j6,za,zb,sprod) 
      use mod_types; use mod_consts_dp
      implicit none
      integer j1,j2,j3,j4,j5,j6
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)
      complex(dp) :: zab2
      complex(dp) :: Fvfbit0

C---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
C---end statement functions
c
      Fvfbit0=0.5d0*zb(j4,j2)*za(j1,j5)
     & *zab2(j3,j2,j4,j6)/zab2(j2,j3,j4,j1)
     & /(sprod(j3,j4)+sprod(j3,j2)+sprod(j4,j2))
      return
      end
      


