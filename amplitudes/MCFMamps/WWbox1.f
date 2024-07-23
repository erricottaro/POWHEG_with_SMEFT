      subroutine WWbox1(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,
     &     app,apm,amp,amm)
      use mod_types; use mod_consts_dp
      implicit none
      real(dp):: masssq,mass
      integer k1,k2,k3,k4,k5,k6
      complex(dp):: app,apm,amp,amm,boxamp,zab2
      real(dp):: s134
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
c--- end statement functions

      masssq=mass**2

c--- Box 1
      app=czero
      amm=czero
c      apm=czero
c      amp=czero

c--- Improved expressions
      s134=sprod(k1,k3)+sprod(k1,k4)+sprod(k3,k4)
      boxamp =sprod(k1,k2)*zb(k2,k6)*za(k1,k3)
     & *(masssq-s134)**3/zab2(k1,k3,k4,k2)**4*
     & (masssq*zb(k4,k2)*za(k1,k5)
     & +za(k1,k3)*zb(k3,k4)*za(k5,k6)*zb(k2,k6))
      apm=boxamp*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))

      boxamp=(masssq-s134)*sprod(k1,k2)/zab2(k2,k3,k4,k1)**4
     & *(masssq*zb(k4,k1)*za(k2,k5)+zab2(k2,k1,k3,k4)*zab2(k5,k3,k4,k1))
     & *(masssq*zb(k1,k6)-zb(k5,k6)*zab2(k5,k2,k6,k1))
     & *(masssq*za(k2,k3)+za(k3,k4)*zab2(k2,k1,k3,k4))

      amp=boxamp*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))

      return
      end
      
