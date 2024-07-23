      subroutine WWtriangle12(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,
     &     app,apm)
      use mod_types; use mod_consts_dp; use common_def
      implicit none
      real(dp) ::  masssq,mass
      integer k1,k2,k3,k4,k5,k6,i1,i2
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
      complex(dp) ::  app,apm,triamp,iza,A134,A134c,s134h,s234h
c--- statement functions
      iza(i1,i2)=cone/za(i1,i2)
c--- end statement functions

      masssq=mass**2
      apm=czero

c--- Triangle 12
      A134=-(za(k1,k3)*zb(k3,k2)+za(k1,k4)*zb(k4,k2))
      A134c=dconjg(A134)
      s134h=sprod(k1,k3)+sprod(k1,k4)+sprod(k3,k4)+masssq
      s234h=sprod(k2,k3)+sprod(k2,k4)+sprod(k3,k4)+masssq

      triamp =  + iza(k1,k2)*masssq*A134**(-2) * ( za(k1,k3)*za(k1,k5)*
     &    zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*s234h + za(k1,k3)*za(k1,k5)*zb(
     &    k1,k2)*zb(k2,k4)*zb(k2,k6)*s134h )
      triamp = triamp + iza(k1,k2)*masssq*A134**(-1) * ( za(k1,k3)*
     &     za(k1,
     &    k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) + za(k1,k3)*za(k1,k5)*zb(k1
     &    ,k2)*zb(k1,k6)*zb(k2,k4) + za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(
     &    k2,k4)*zb(k2,k6) + za(k1,k5)*za(k2,k3)*zb(k1,k2)*zb(k2,k4)*
     &    zb(k2,k6) - 2.D0*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(
     &    k4,k6) )
      triamp = triamp + iza(k1,k2)*masssq*A134c**(-2) * 
     &     ( za(k2,k3)*za(k2
     &    ,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*s234h + za(k2,k3)*za(k2,k5
     &    )*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*s134h )
      triamp = triamp + iza(k1,k2)*masssq*A134c**(-1) * ( za(k1,k3)*
     &     za(k2
     &    ,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k1,k5)*za(k2,k3)*zb(
     &    k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k2,k3)*za(k2,k5)*zb(k1,k2)*
     &    zb(k1,k4)*zb(k2,k6) + za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k6)
     &    *zb(k2,k4) - 2.D0*za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(
     &    k4,k6) )

      app=triamp*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))

      triamp =  + A134**(-4) * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1
     &    ,k2)*zb(k2,k4)*zb(k2,k6)*s234h**3 - za(k1,k2)*za(k1,k3)*za(k1
     &    ,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*s134h**3 )
      triamp = triamp + A134**(-3) * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k2,k6)*s134h**2 - za(k1,k2)*za(k1,k3)
     &    *za(k1,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4)*s134h**2 + za(k1,k2)
     &    *za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*s134h**2
     &     - za(k1,k2)*za(k1,k5)*za(k2,k3)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6
     &    )*s234h**2 + za(k1,k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4
     &    )*zb(k4,k6)*s234h**2 + za(k1,k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2
     &    )*zb(k2,k4)*zb(k4,k6)*s134h**2 )
      triamp = triamp + A134**(-2) * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*s134h + za(k1,k2)*za(k1,k3)*
     &    za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6)*s134h + za(k1,k2)*za(
     &    k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4)*s134h + za(k1,
     &    k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6)*s134h
     &     - za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(k4,k6
     &    )*s134h )
      triamp = triamp + A134**(-1) * ( za(k1,k2)*za(k1,k3)*za(k2,k5)*
     &    zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*za(k2,k5)*za(k3,k4)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      triamp = triamp + masssq*A134**(-3)*A134c * ( 2.D0*za(k1,k3)*
     &     za(k1,
     &    k5)*zb(k2,k4)*zb(k2,k6)*s234h + 2.D0*za(k1,k3)*za(k1,k5)*zb(
     &    k2,k4)*zb(k2,k6)*s134h )
      triamp = triamp + masssq*A134**(-2) * (  - za(k1,k3)*za(k2,k5)*zb(
     &    k1,k4)*zb(k2,k6)*s234h - za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,
     &    k6)*s134h - za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*s234h - 
     &    za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*s134h )
      triamp = triamp + masssq*A134**(-2)*A134c * ( za(k1,k3)*za(k1,k5)*
     &    zb(k1,k4)*zb(k2,k6) + za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4)
     &     - za(k1,k3)*za(k2,k5)*zb(k2,k4)*zb(k2,k6) + za(k1,k5)*za(k2,
     &    k3)*zb(k2,k4)*zb(k2,k6) - 2.D0*za(k1,k5)*za(k3,k4)*zb(k2,k4)*
     &    zb(k4,k6) )
      triamp = triamp + masssq*A134**(-1) * (  - za(k1,k3)*za(k2,k5)*zb(
     &    k1,k4)*zb(k1,k6) - za(k1,k5)*za(k2,k3)*zb(k1,k4)*zb(k1,k6) - 
     &    za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6) + za(k2,k3)*za(k2,k5)
     &    *zb(k1,k6)*zb(k2,k4) + 2.D0*za(k2,k5)*za(k3,k4)*zb(k1,k4)*zb(
     &    k4,k6) )

      apm=triamp*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))

      return
      end
