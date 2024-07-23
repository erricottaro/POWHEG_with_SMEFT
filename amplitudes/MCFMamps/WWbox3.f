      subroutine WWbox3(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,
     &     app,apm,amp,amm)
      use mod_types; use mod_consts_dp
      implicit none
      real(dp):: mass,masssq,s12,s34,s56,s134,s156
      integer k1,k2,k3,k4,k5,k6
      complex(dp):: app,apm,amp,amm,boxamp,zab2
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
c--- statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
c--- end statement functions

      apm=czero
      amp=czero

      masssq=mass**2
      s12=sprod(k1,k2)
      s34=sprod(k3,k4)
      s56=sprod(k5,k6)
      s134=sprod(k1,k3)+sprod(k1,k4)+sprod(k3,k4)
      s156=sprod(k1,k5)+sprod(k1,k6)+sprod(k5,k6)

      boxamp=-(s134*s156-s56*s34+masssq*s12)*za(k2,k3)/s12**2
     & /za(k1,k2)**2
     & *(za(k1,k5)*zb(k5,k6)*zb(k1,k2)
     & +zb(k1,k6)*masssq*s12/zab2(k2,k5,k6,k1))
     & *(za(k1,k5)*za(k2,k3)*zb(k3,k4)*zb(k1,k2)
     &  -za(k2,k5)*zb(k1,k4)*masssq*s12/zab2(k2,k5,k6,k1))
      app=boxamp*(-half*ci)/(s34*s56)
c      write(6,*)
c      write(6,*) 'box3:app',app
      boxamp=-(s134*s156-s56*s34+masssq*s12)*zb(k2,k6)/s12**2
     &     /zb(k1,k2)**2
     &  *(za(k3,k4)*za(k1,k2)*zb(k1,k4)
     &  +za(k1,k3)*masssq*s12/zab2(k1,k5,k6,k2))
     & *(za(k5,k6)*za(k1,k2)*zb(k2,k6)*zb(k1,k4)
     & -za(k1,k5)*zb(k2,k4)*masssq*s12/zab2(k1,k5,k6,k2))
      amm=boxamp*(-half*ci)/(s34*s56)
c      write(6,*) 'box3:amm',amm


      boxamp=-(s134*s156-s56*s34+masssq*s12)*zb(k2,k6)
     & *masssq/zab2(k1,k5,k6,k2)**2/s12**2
     & *(za(k1,k2)*za(k4,k3)*zb(k1,k4)
     & -za(k1,k3)*masssq*s12/zab2(k1,k5,k6,k2))
     & *(za(k1,k2)*zb(k1,k4)*zab2(k5,k3,k4,k2)
     &  -za(k1,k5)*zb(k2,k4)*masssq*s12/zab2(k1,k5,k6,k2))
      apm=boxamp*(-half*ci)/(s34*s56)
c      write(6,*) 'box3:apm',apm
      boxamp=-(s134*s156-s56*s34+masssq*s12)*za(k2,k3)
     & *masssq/zab2(k2,k5,k6,k1)**2/s12**2
     & *(za(k1,k5)*zb(k1,k2)*zb(k5,k6)
     &  +zb(k1,k6)*masssq*s12/zab2(k2,k5,k6,k1))
     & *(za(k1,k5)*zb(k1,k2)*zab2(k2,k5,k6,k4)
     & +za(k2,k5)*zb(k1,k4)*masssq*s12/zab2(k2,k5,k6,k1))
      amp=boxamp*(-half*ci)/(s34*s56)
c      write(6,*) 'box3:amp',amp


      return
      end
