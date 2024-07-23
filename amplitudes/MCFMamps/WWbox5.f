      subroutine WWbox5(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,
     &     app,apm,amp,amm)
      use mod_types; use mod_consts_dp
      implicit none
      real(dp):: mass,masssq,s134,rat,masssqps134
      integer k1,k2,k3,k4,k5,k6,i1,i2
      complex(dp):: app,apm,amp,amm,zab2,izab2,iza,izb
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
c--- statement functions
      iza(i1,i2)=cone/za(i1,i2)
      izb(i1,i2)=cone/zb(i1,i2)
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      izab2(k1,k2,k3,k4)=cone/(za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4))
c--- end statement functions

      masssq=mass**2
      s134=sprod(k1,k3)+sprod(k1,k4)+sprod(k3,k4)
      rat=masssq/sprod(k1,k2)
      masssqps134=masssq+s134
      app =  + iza(k1,k2)*masssq * (  - 2.D0*za(k1,k3)*za(k1,k5)
     &     *zb(k1,k2
     &    )*zb(k1,k4)*zb(k1,k6) + 2.D0*za(k1,k3)*za(k2,k5)*zb(k1,k2)*
     &    zb(k1,k4)*zb(k2,k6)*rat + 2.D0*za(k1,k5)*za(k2,k3)*zb(k1,k2)*
     &    zb(k1,k6)*zb(k2,k4)*rat + 2.D0*za(k1,k5)*za(k3,k4)*zb(k1,k2)*
     &    zb(k1,k4)*zb(k4,k6) )
      app = app + iza(k1,k2)*izab2(k1,k3,k4,k2)*masssq*masssqps134 * 
     & ( za(
     &    k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) + za(k1,k3)*
     &    za(k1,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4) - za(k1,k5)*za(k3,k4)
     &    *zb(k1,k2)*zb(k2,k4)*zb(k4,k6) )
      app = app + iza(k1,k2)*izab2(k1,k3,k4,k2)**2*masssq*masssqps134**2
     &  * (  - za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6) )
      app = app + iza(k1,k2)*izab2(k1,k3,k4,k2)*zab2(k2,k3,k4,k1)*masssq
     &  * (  - 2.D0*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*
     &    rat )
      app = app + iza(k1,k2)*izab2(k2,k3,k4,k1)*masssq*masssqps134 * 
     &  ( za(
     &    k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k1,k5)*
     &    za(k2,k3)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k2,k5)*za(k3,k4)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      app = app + iza(k1,k2)*izab2(k2,k3,k4,k1)**2*masssq*masssqps134**2
     &  * (  - za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) )
      app = app + iza(k1,k2)*izab2(k2,k3,k4,k1)*zab2(k1,k3,k4,k2)*masssq
     &  * (  - 2.D0*za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*
     &    rat )
      amm =  + izb(k1,k2)*masssq * ( 2.D0*za(k1,k2)*za(k1,k3)*za(k2,k5)*
     &    zb(k1,k4)*zb(k2,k6)*rat + 2.D0*za(k1,k2)*za(k1,k3)*za(k2,k5)*
     &    zb(k1,k6)*zb(k2,k4) + 2.D0*za(k1,k2)*za(k1,k5)*za(k2,k3)*zb(
     &    k1,k6)*zb(k2,k4)*rat - 2.D0*za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(
     &    k2,k4)*zb(k4,k6) )
      amm = amm + izb(k1,k2)*izab2(k1,k3,k4,k2)*masssq*masssqps134 * 
     &    ( za(
     &    k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4) - za(k1,k2)*
     &    za(k1,k3)*za(k2,k5)*zb(k2,k4)*zb(k2,k6) - za(k1,k2)*za(k1,k5)
     &    *za(k3,k4)*zb(k2,k4)*zb(k4,k6) )
      amm = amm + izb(k1,k2)*izab2(k1,k3,k4,k2)**2*masssq*masssqps134**2
     &  * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6) )
      amm = amm + izb(k1,k2)*izab2(k1,k3,k4,k2)*zab2(k2,k3,k4,k1)*masssq
     &  * (  - 2.D0*za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)*
     &    rat )
      amm = amm + izb(k1,k2)*izab2(k2,k3,k4,k1)*masssq*masssqps134 * 
     &    ( za(
     &    k1,k2)*za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*
     &    za(k2,k3)*za(k2,k5)*zb(k1,k6)*zb(k2,k4) - za(k1,k2)*za(k2,k5)
     &    *za(k3,k4)*zb(k1,k4)*zb(k4,k6) )
      amm = amm + izb(k1,k2)*izab2(k2,k3,k4,k1)**2*masssq*masssqps134**2
     &  * (  - za(k1,k2)*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6) )
      amm = amm + izb(k1,k2)*izab2(k2,k3,k4,k1)*zab2(k1,k3,k4,k2)*masssq
     &  * (  - 2.D0*za(k1,k2)*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*
     &    rat )
      apm =  + masssq * (  - 2.D0*za(k2,k3)*za(k2,k5)*zb(k1,k4)*
     &   zb(k1,k6)
     & - 2.D0*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*rat )
      apm = apm + izab2(k1,k3,k4,k2)*masssqps134 * 
     &     ( za(k1,k2)*za(k1,k3)*
     &    za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*za(k2,k5)
     &    *za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)*masssq*masssqps134 * 
     &    (  - za(k1,k3)*
     &    za(k2,k5)*zb(k1,k4)*zb(k1,k6) - za(k1,k5)*za(k2,k3)*zb(k1,k4)
     &    *zb(k1,k6) + za(k2,k3)*za(k2,k5)*zb(k1,k6)*zb(k2,k4) + za(k2,
     &    k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**2*masssqps134**2 * ( za(k1,k2)*za(
     &    k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*
     &    za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) - za(k1,k2)
     &    *za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4) - za(k1,k2
     &    )*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6) + za(k1,
     &    k2)*za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**2*masssq*masssqps134**2 * 
     &   ( za(k1,k3)
     &    *za(k2,k5)*zb(k1,k4)*zb(k2,k6) + za(k1,k5)*za(k2,k3)*zb(k1,k6
     &    )*zb(k2,k4) )
      apm = apm + izab2(k1,k3,k4,k2)**3*masssqps134**3 * (  - za(k1,k2)*
     &    za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) - za(k1,k2)
     &    *za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4) + za(k1,k2
     &    )*za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6) + za(k1,
     &    k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**4*masssqps134**4 * ( za(k1,k2)*za(
     &    k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**3*zab2(k2,k3,k4,k1)*masssq*
     & masssqps134**2 * (  - 4.D0*za(k1,k3)*za(k1,k5)*zb(k2,k4)
     & *zb(k2,k6))
      apm = apm + izab2(k1,k3,k4,k2)**2*zab2(k2,k3,k4,k1)*masssq*
     & masssqps134 * ( 3.D0*za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k2,k6) 
     & + 3.D0
     &    *za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4) - 3.D0*za(k1,k3)*za(
     &    k2,k5)*zb(k2,k4)*zb(k2,k6) - 3.D0*za(k1,k5)*za(k3,k4)*zb(k2,
     &    k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**2*zab2(k2,k3,k4,k1)**2*masssq * ( 
     & - 2.D0*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)*rat )
      apm = apm + izab2(k1,k3,k4,k2)*zab2(k2,k3,k4,k1)*masssq * 
     &    (  - 2.D0
     &    *za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6) + 2.D0*za(k1,k3)*za(
     &    k2,k5)*zb(k1,k4)*zb(k2,k6) + 2.D0*za(k1,k3)*za(k2,k5)*zb(k1,
     &    k4)*zb(k2,k6)*rat + 2.D0*za(k1,k3)*za(k2,k5)*zb(k1,k6)*zb(k2,
     &    k4) + 2.D0*za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*rat + 2.D0
     &    *za(k1,k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6) - 2.D0*za(k2,k5)*za(
     &    k3,k4)*zb(k2,k4)*zb(k4,k6) )
      amp =  + masssq * (  - 2.D0*za(k1,k3)*za(k1,k5)*zb(k2,k4)*
     &    zb(k2,k6)
     &    *rat )
      amp = amp + izab2(k2,k3,k4,k1)*masssq*masssqps134 * 
     &    (  - za(k1,k3)*
     &    za(k1,k5)*zb(k1,k6)*zb(k2,k4) + za(k1,k5)*za(k3,k4)*zb(k2,k4)
     &    *zb(k4,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**2*masssq*masssqps134**2 * 
     &    ( za(k1,k3)
     &    *za(k2,k5)*zb(k1,k4)*zb(k2,k6) + za(k1,k5)*za(k2,k3)*zb(k1,k6
     &    )*zb(k2,k4) )
      amp = amp + izab2(k2,k3,k4,k1)**3*masssqps134**3 * (  - za(k1,k2)*
     &    za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k1,k2)
     &    *za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**4*masssqps134**4 * ( za(k1,k2)*za(
     &    k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**3*zab2(k1,k3,k4,k2)*masssq*
     & masssqps134**2 * (  - 4.D0*za(k2,k3)*za(k2,k5)*zb(k1,k4)
     &  *zb(k1,k6))
      amp = amp + izab2(k2,k3,k4,k1)**2*zab2(k1,k3,k4,k2)*masssq*
     & masssqps134 * ( 3.D0*za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6) 
     &  - 3.D0
     &    *za(k2,k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**2*zab2(k1,k3,k4,k2)**2*masssq * ( 
     &     - 2.D0*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*rat )
      amp = amp + izab2(k2,k3,k4,k1)*zab2(k1,k3,k4,k2)*masssq * ( 2.D0*
     &    za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)*rat + 2.D0*za(k1,k5)*
     &    za(k2,k3)*zb(k1,k6)*zb(k2,k4)*rat )
      app=app*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))
      amm=amm*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))
      apm=apm*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))
      amp=amp*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))
      return
      end
