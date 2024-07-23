      subroutine WWtriangle6(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,app,apm)
      use mod_types; use mod_consts_dp; use common_def
      implicit none
      real(dp) ::  masssq,mass
      integer k1,k2,k3,k4,k5,k6,i1,i2,k12h,k34h
      complex(dp) ::  app,apm,triamp,iza,izb,
     & s12,s34,dot1234,delta,ga,x,y
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
      parameter(k12h=9,k34h=10)

c--- statement functions
      iza(i1,i2)=cone/za(i1,i2)
      izb(i1,i2)=cone/zb(i1,i2)
c--- end statement functions

      masssq=mass**2
      apm=czero

c--- Triangle 6
      s12=za(k1,k2)*zb(k2,k1)
      s34=za(k3,k4)*zb(k4,k3)
      dot1234=(za(k1,k3)*zb(k3,k1)+za(k1,k4)*zb(k4,k1)
     &        +za(k2,k3)*zb(k3,k2)+za(k2,k4)*zb(k4,k2))/2d0
      delta=sqrt(dot1234**2-s12*s34)
      ga=dot1234+delta

      x=(s12*s34+ga*(s34-masssq))/(s12*s34-ga**2)
      y=s12*(masssq-s34-ga)/(s12*s34-ga**2)

      triamp =  + y * ( za(k1,k34h)*za(k2,k12h)*za(k3,k5)*za(k3,k12h)*
     &    zb(k2,k6)*zb(k2,k34h)*zb(k3,k4)*iza(k1,k2)*iza(k1,k12h)**2 + 
     &    za(k1,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k3
     &    ,k4)*iza(k1,k2)*iza(k2,k34h) + za(k1,k34h)*za(k3,k5)*za(k3,
     &    k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k3,k4)*iza(k1,k2)*iza(k2,k12h)
     &     - za(k1,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k2)*zb(k2,k6)*
     &    zb(k4,k34h)*iza(k1,k12h)**2 - za(k1,k12h)*za(k2,k34h)*za(k3,
     &    k5)*za(k3,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k3,k4)*iza(k1,k2)*
     &    iza(k2,k12h)**2 + za(k1,k12h)*za(k3,k5)*za(k3,k34h)*zb(k1,k6)
     &    *zb(k1,k34h)*zb(k3,k4)*iza(k1,k2)*iza(k2,k12h) - za(k2,k34h)*
     &    za(k3,k5)*za(k3,k34h)*zb(k2,k6)*zb(k2,k34h)*zb(k3,k4)*iza(k1,
     &    k2)*iza(k1,k34h) - za(k2,k34h)*za(k3,k5)*za(k3,k12h)*zb(k2,k6
     &    )*zb(k2,k34h)*zb(k3,k4)*iza(k1,k2)*iza(k1,k12h) + za(k2,k34h)
     &    *za(k3,k12h)*za(k5,k12h)*zb(k1,k2)*zb(k1,k6)*zb(k4,k34h)*iza(
     &    k2,k12h)**2 - za(k2,k12h)*za(k3,k5)*za(k3,k34h)*zb(k2,k6)*zb(
     &    k2,k34h)*zb(k3,k4)*iza(k1,k2)*iza(k1,k12h) )
      triamp = triamp + y * ( 4.D0*za(k3,k5)*za(k3,k34h)*zb(k1,k2)*zb(
     &    k3,k4)*zb(k6,k34h)*iza(k1,k2) - za(k3,k34h)*za(k5,k34h)*zb(k1
     &    ,k2)*zb(k1,k6)*zb(k4,k34h)*iza(k2,k34h) + za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k2)*zb(k2,k6)*zb(k4,k34h)*iza(k1,k34h) - za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k2)*zb(k1,k6)*zb(k4,k34h)*iza(k2,k12h
     &    ) + za(k3,k34h)*za(k5,k12h)*zb(k1,k2)*zb(k2,k6)*zb(k4,k34h)*
     &    iza(k1,k12h) - za(k3,k12h)*za(k5,k34h)*zb(k1,k2)*zb(k1,k6)*
     &    zb(k4,k34h)*iza(k2,k12h) + za(k3,k12h)*za(k5,k34h)*zb(k1,k2)*
     &    zb(k2,k6)*zb(k4,k34h)*iza(k1,k12h) )
      triamp = triamp + y**2 * (  - za(k1,k34h)**2*za(k2,k12h)*za(k3,
     &    k12h)*za(k5,k12h)*zb(k2,k6)*zb(k2,k34h)*zb(k4,k34h)*iza(k1,k2
     &    )*iza(k1,k12h)**3 - za(k1,k34h)*za(k2,k34h)*za(k3,k12h)*za(k5
     &    ,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,
     &    k12h)**2 + za(k1,k34h)*za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*
     &    zb(k2,k6)*zb(k2,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k12h)**2
     &     + za(k1,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k2,k6)*
     &    zb(k2,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k12h)**2 + za(k1,
     &    k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k2,k6)*zb(k2,
     &    k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k12h)**2 - 2.D0*za(k1,
     &    k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k34h)
     &    *iza(k1,k2)**2 + za(k1,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k6
     &    )*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,k34h) + za(k1,
     &    k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,
     &    k34h)*iza(k1,k2)*iza(k2,k12h) + za(k1,k34h)*za(k3,k12h)*za(k5
     &    ,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,
     &    k12h) )
      triamp = triamp + y**2 * ( za(k1,k12h)*za(k2,k34h)**2*za(k3,k12h)
     &    *za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k2)*
     &    iza(k2,k12h)**3 - za(k1,k12h)*za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,
     &    k12h)**2 - za(k1,k12h)*za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*
     &    zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,k12h)**2
     &     + za(k1,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*
     &    zb(k4,k34h)*iza(k1,k2)*iza(k2,k12h) - 2.D0*za(k2,k34h)*za(k3,
     &    k5)*za(k3,k34h)*zb(k2,k34h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)
     &    **2 - za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k6)*zb(k2,
     &    k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k34h) - za(k2,k34h)*za(k3
     &    ,k34h)*za(k5,k12h)*zb(k2,k6)*zb(k2,k34h)*zb(k4,k34h)*iza(k1,
     &    k2)*iza(k1,k12h) - za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k2,
     &    k6)*zb(k2,k34h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k12h) - za(k2,
     &    k12h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k6)*zb(k2,k34h)*zb(k4,
     &    k34h)*iza(k1,k2)*iza(k1,k12h) )
      triamp = triamp + y**2 * ( 4.D0*za(k3,k34h)*za(k5,k34h)*zb(k1,k2)
     &    *zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2) )
      triamp = triamp + y**3 * (  - 2.D0*za(k1,k34h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0
     &    *za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k2)**2 )
      triamp = triamp + x * (  - za(k1,k34h)*za(k2,k12h)*za(k3,k5)*za(
     &    k3,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k3,k4)*iza(k1,k2)*iza(k2,
     &    k34h)**2 + za(k1,k34h)*za(k3,k5)*za(k3,k12h)*zb(k1,k6)*zb(k1,
     &    k12h)*zb(k3,k4)*iza(k1,k2)*iza(k2,k34h) + za(k1,k12h)*za(k2,
     &    k34h)*za(k3,k5)*za(k3,k34h)*zb(k2,k6)*zb(k2,k12h)*zb(k3,k4)*
     &    iza(k1,k2)*iza(k1,k34h)**2 + za(k1,k12h)*za(k3,k5)*za(k3,k34h
     &    )*zb(k1,k6)*zb(k1,k12h)*zb(k3,k4)*iza(k1,k2)*iza(k2,k34h) + 
     &    za(k1,k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k3
     &    ,k4)*iza(k1,k2)*iza(k2,k12h) - za(k1,k12h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k2)*zb(k2,k6)*zb(k4,k12h)*iza(k1,k34h)**2 - za(k2
     &    ,k34h)*za(k3,k5)*za(k3,k12h)*zb(k2,k6)*zb(k2,k12h)*zb(k3,k4)*
     &    iza(k1,k2)*iza(k1,k34h) - za(k2,k12h)*za(k3,k5)*za(k3,k34h)*
     &    zb(k2,k6)*zb(k2,k12h)*zb(k3,k4)*iza(k1,k2)*iza(k1,k34h) - za(
     &    k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(k2,k6)*zb(k2,k12h)*zb(k3,k4
     &    )*iza(k1,k2)*iza(k1,k12h) + za(k2,k12h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k2)*zb(k1,k6)*zb(k4,k12h)*iza(k2,k34h)**2 )
      triamp = triamp + x * ( 4.D0*za(k3,k5)*za(k3,k12h)*zb(k1,k2)*zb(
     &    k3,k4)*zb(k6,k12h)*iza(k1,k2) - za(k3,k34h)*za(k5,k12h)*zb(k1
     &    ,k2)*zb(k1,k6)*zb(k4,k12h)*iza(k2,k34h) + za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k2)*zb(k2,k6)*zb(k4,k12h)*iza(k1,k34h) - za(k3,
     &    k12h)*za(k5,k34h)*zb(k1,k2)*zb(k1,k6)*zb(k4,k12h)*iza(k2,k34h
     &    ) + za(k3,k12h)*za(k5,k34h)*zb(k1,k2)*zb(k2,k6)*zb(k4,k12h)*
     &    iza(k1,k34h) - za(k3,k12h)*za(k5,k12h)*zb(k1,k2)*zb(k1,k6)*
     &    zb(k4,k12h)*iza(k2,k12h) + za(k3,k12h)*za(k5,k12h)*zb(k1,k2)*
     &    zb(k2,k6)*zb(k4,k12h)*iza(k1,k12h) )
      triamp = triamp + x*y * (  - za(k1,k34h)*za(k2,k12h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(
     &    k2,k34h)**2 - za(k1,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,k34h)**2
     &     + za(k1,k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k2,k6)*
     &    zb(k2,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(k1,k12h)**2 + za(k1,
     &    k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k2,k6)*zb(k2,
     &    k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k12h)**2 - 2.D0*za(k1,
     &    k34h)*za(k3,k5)*za(k3,k12h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k12h)
     &    *iza(k1,k2)**2 - 2.D0*za(k1,k34h)*za(k3,k5)*za(k3,k12h)*zb(k1
     &    ,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)**2 + za(k1,k34h)*za(
     &    k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1
     &    ,k2)*iza(k2,k34h) + za(k1,k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1
     &    ,k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,k34h) + za(k1,
     &    k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,
     &    k12h)*iza(k1,k2)*iza(k2,k34h) )
      triamp = triamp + x*y * ( za(k1,k34h)*za(k3,k12h)*za(k5,k34h)*zb(
     &    k1,k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,k34h) + za(
     &    k1,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,
     &    k12h)*iza(k1,k2)*iza(k2,k12h) + za(k1,k34h)*za(k3,k12h)*za(k5
     &    ,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,
     &    k12h) + za(k1,k12h)*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k2
     &    ,k6)*zb(k2,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(k1,k34h)**2 + za(
     &    k1,k12h)*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k6)*zb(k2,
     &    k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k34h)**2 - za(k1,k12h)*
     &    za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(
     &    k4,k12h)*iza(k1,k2)*iza(k2,k12h)**2 - za(k1,k12h)*za(k2,k34h)
     &    *za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k34h)*
     &    iza(k1,k2)*iza(k2,k12h)**2 - 2.D0*za(k1,k12h)*za(k3,k5)*za(k3
     &    ,k34h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0
     &    *za(k1,k12h)*za(k3,k5)*za(k3,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(
     &    k6,k34h)*iza(k1,k2)**2 )
      triamp = triamp + x*y * ( za(k1,k12h)*za(k3,k34h)*za(k5,k34h)*zb(
     &    k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(k2,k34h) + za(
     &    k1,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,
     &    k34h)*iza(k1,k2)*iza(k2,k34h) + za(k1,k12h)*za(k3,k34h)*za(k5
     &    ,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(k2,
     &    k12h) + za(k1,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,
     &    k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,k12h) + za(k1,k12h)*za(k3
     &    ,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,
     &    k2)*iza(k2,k12h) + za(k1,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,
     &    k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k2,k12h) - 2.D0*
     &    za(k2,k34h)*za(k3,k5)*za(k3,k12h)*zb(k2,k34h)*zb(k3,k4)*zb(k6
     &    ,k12h)*iza(k1,k2)**2 - 2.D0*za(k2,k34h)*za(k3,k5)*za(k3,k12h)
     &    *zb(k2,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)**2 - za(k2,k34h
     &    )*za(k3,k34h)*za(k5,k12h)*zb(k2,k6)*zb(k2,k34h)*zb(k4,k12h)*
     &    iza(k1,k2)*iza(k1,k34h) - za(k2,k34h)*za(k3,k34h)*za(k5,k12h)
     &    *zb(k2,k6)*zb(k2,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k34h) )
      triamp = triamp + x*y * (  - za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*
     &    zb(k2,k6)*zb(k2,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(k1,k34h) - 
     &    za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k2,k6)*zb(k2,k12h)*zb(
     &    k4,k34h)*iza(k1,k2)*iza(k1,k34h) - za(k2,k34h)*za(k3,k12h)*
     &    za(k5,k12h)*zb(k2,k6)*zb(k2,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(
     &    k1,k12h) - za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k2,k6)*zb(
     &    k2,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k12h) - 2.D0*za(k2,
     &    k12h)*za(k3,k5)*za(k3,k34h)*zb(k2,k34h)*zb(k3,k4)*zb(k6,k12h)
     &    *iza(k1,k2)**2 - 2.D0*za(k2,k12h)*za(k3,k5)*za(k3,k34h)*zb(k2
     &    ,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)**2 - za(k2,k12h)*za(
     &    k3,k34h)*za(k5,k34h)*zb(k2,k6)*zb(k2,k34h)*zb(k4,k12h)*iza(k1
     &    ,k2)*iza(k1,k34h) - za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k2
     &    ,k6)*zb(k2,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,k34h) - za(k2,
     &    k12h)*za(k3,k34h)*za(k5,k12h)*zb(k2,k6)*zb(k2,k34h)*zb(k4,
     &    k12h)*iza(k1,k2)*iza(k1,k12h) - za(k2,k12h)*za(k3,k34h)*za(k5
     &    ,k12h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,k34h)*iza(k1,k2)*iza(k1,
     &    k12h) )
      triamp = triamp + x*y * (  - za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*
     &    zb(k2,k6)*zb(k2,k34h)*zb(k4,k12h)*iza(k1,k2)*iza(k1,k12h) - 
     &    za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k2,k6)*zb(k2,k12h)*zb(
     &    k4,k34h)*iza(k1,k2)*iza(k1,k12h) + 4.D0*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k2)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2) + 4.D0*za(
     &    k3,k34h)*za(k5,k12h)*zb(k1,k2)*zb(k4,k12h)*zb(k6,k34h)*iza(k1
     &    ,k2) + 4.D0*za(k3,k12h)*za(k5,k34h)*zb(k1,k2)*zb(k4,k34h)*zb(
     &    k6,k12h)*iza(k1,k2) + 4.D0*za(k3,k12h)*za(k5,k34h)*zb(k1,k2)*
     &    zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2) )
      triamp = triamp + x*y**2 * (  - 2.D0*za(k1,k34h)*za(k3,k34h)*za(
     &    k5,k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)**2 - 
     &    2.D0*za(k1,k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,
     &    k12h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0*za(k1,k34h)*za(k3,k34h
     &    )*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)
     &    **2 - 2.D0*za(k1,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*
     &    zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0*za(k1,k34h)*za(
     &    k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(
     &    k1,k2)**2 - 2.D0*za(k1,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,
     &    k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0*za(k1,k12h
     &    )*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)
     &    *iza(k1,k2)**2 - 2.D0*za(k1,k12h)*za(k3,k34h)*za(k5,k34h)*zb(
     &    k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0*za(k1,
     &    k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)**2 - 2.D0*za(k2,k34h)*za(k3,k34h)*za(k5,k12h
     &    )*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)**2 )
      triamp = triamp + x*y**2 * (  - 2.D0*za(k2,k34h)*za(k3,k34h)*za(
     &    k5,k12h)*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)**2 - 
     &    2.D0*za(k2,k34h)*za(k3,k34h)*za(k5,k12h)*zb(k2,k12h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0*za(k2,k34h)*za(k3,k12h
     &    )*za(k5,k34h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)
     &    **2 - 2.D0*za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k2,k34h)*
     &    zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0*za(k2,k34h)*za(
     &    k3,k12h)*za(k5,k34h)*zb(k2,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(
     &    k1,k2)**2 - 2.D0*za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k2,
     &    k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0*za(k2,k12h
     &    )*za(k3,k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k34h)
     &    *iza(k1,k2)**2 - 2.D0*za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(
     &    k2,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2 )
      triamp = triamp + x**2 * ( za(k1,k34h)*za(k2,k12h)**2*za(k3,k34h)
     &    *za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,k2)*
     &    iza(k2,k34h)**3 - za(k1,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,k2)*iza(k2,
     &    k34h)**2 - za(k1,k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*
     &    zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,k2)*iza(k2,k34h)**2
     &     + za(k1,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)*
     &    zb(k4,k12h)*iza(k1,k2)*iza(k2,k34h) - za(k1,k12h)**2*za(k2,
     &    k34h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,
     &    k12h)*iza(k1,k2)*iza(k1,k34h)**3 + za(k1,k12h)*za(k2,k34h)*
     &    za(k3,k34h)*za(k5,k12h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,k12h)*
     &    iza(k1,k2)*iza(k1,k34h)**2 + za(k1,k12h)*za(k2,k34h)*za(k3,
     &    k12h)*za(k5,k34h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,k12h)*iza(k1,k2
     &    )*iza(k1,k34h)**2 - za(k1,k12h)*za(k2,k12h)*za(k3,k34h)*za(k5
     &    ,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,k2)*iza(k2,
     &    k34h)**2 )
      triamp = triamp + x**2 * ( za(k1,k12h)*za(k2,k12h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,k12h)*iza(k1,k2)*iza(
     &    k1,k34h)**2 - 2.D0*za(k1,k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,
     &    k12h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)**2 + za(k1,k12h)*za(k3
     &    ,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,
     &    k2)*iza(k2,k34h) + za(k1,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,
     &    k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,k2)*iza(k2,k34h) + za(k1,
     &    k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,
     &    k12h)*iza(k1,k2)*iza(k2,k12h) - za(k2,k34h)*za(k3,k12h)*za(k5
     &    ,k12h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,k12h)*iza(k1,k2)*iza(k1,
     &    k34h) - 2.D0*za(k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(k2,k12h)*
     &    zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)**2 - za(k2,k12h)*za(k3,k34h)
     &    *za(k5,k12h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,k12h)*iza(k1,k2)*
     &    iza(k1,k34h) - za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k2,k6)*
     &    zb(k2,k12h)*zb(k4,k12h)*iza(k1,k2)*iza(k1,k34h) - za(k2,k12h)
     &    *za(k3,k12h)*za(k5,k12h)*zb(k2,k6)*zb(k2,k12h)*zb(k4,k12h)*
     &    iza(k1,k2)*iza(k1,k12h) )
      triamp = triamp + x**2 * ( 4.D0*za(k3,k12h)*za(k5,k12h)*zb(k1,k2)
     &    *zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2) )
      triamp = triamp + x**2*y * (  - 2.D0*za(k1,k34h)*za(k3,k12h)*za(
     &    k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)**2 - 
     &    2.D0*za(k1,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,
     &    k34h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0*za(k1,k34h)*za(k3,k12h
     &    )*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)
     &    **2 - 2.D0*za(k1,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*
     &    zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0*za(k1,k12h)*za(
     &    k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(
     &    k1,k2)**2 - 2.D0*za(k1,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,
     &    k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0*za(k1,k12h
     &    )*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k12h)
     &    *iza(k1,k2)**2 - 2.D0*za(k1,k12h)*za(k3,k12h)*za(k5,k34h)*zb(
     &    k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0*za(k1,
     &    k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,
     &    k34h)*iza(k1,k2)**2 - 2.D0*za(k2,k34h)*za(k3,k12h)*za(k5,k12h
     &    )*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)**2 )
      triamp = triamp + x**2*y * (  - 2.D0*za(k2,k34h)*za(k3,k12h)*za(
     &    k5,k12h)*zb(k2,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)**2 - 
     &    2.D0*za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k2,k12h)*zb(k4,
     &    k12h)*zb(k6,k34h)*iza(k1,k2)**2 - 2.D0*za(k2,k12h)*za(k3,k34h
     &    )*za(k5,k12h)*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)
     &    **2 - 2.D0*za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k2,k12h)*
     &    zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0*za(k2,k12h)*za(
     &    k3,k34h)*za(k5,k12h)*zb(k2,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(
     &    k1,k2)**2 - 2.D0*za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k2,
     &    k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0*za(k2,k12h
     &    )*za(k3,k12h)*za(k5,k34h)*zb(k2,k12h)*zb(k4,k34h)*zb(k6,k12h)
     &    *iza(k1,k2)**2 - 2.D0*za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(
     &    k2,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)**2 )
      triamp = triamp + x**3 * (  - 2.D0*za(k1,k12h)*za(k3,k12h)*za(k5,
     &    k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)**2 - 2.D0
     &    *za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k2,k12h)*zb(k4,k12h)*
     &    zb(k6,k12h)*iza(k1,k2)**2 )
      triamp = triamp - za(k3,k5)*za(k3,k34h)*zb(k1,k2)*zb(k1,k6)*zb(k3
     & ,k4)*iza(k2,k34h) + za(k3,k5)*za(k3,k34h)*zb(k1,k2)*zb(k2,k6)*
     &    zb(k3,k4)*iza(k1,k34h) - za(k3,k5)*za(k3,k12h)*zb(k1,k2)*zb(
     &    k1,k6)*zb(k3,k4)*iza(k2,k12h) + za(k3,k5)*za(k3,k12h)*zb(k1,
     &    k2)*zb(k2,k6)*zb(k3,k4)*iza(k1,k12h)

      app=triamp*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))

      triamp =  + y * ( za(k1,k34h)*za(k2,k12h)*za(k3,k5)*za(k3,k12h)*
     &    zb(k1,k6)*zb(k1,k34h)*zb(k3,k4)*iza(k1,k12h)**2*izb(k1,k2) - 
     &    za(k2,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k3
     &    ,k4)*iza(k1,k34h)*izb(k1,k2) - za(k2,k34h)*za(k3,k5)*za(k3,
     &    k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k3,k4)*iza(k1,k12h)*izb(k1,k2)
     &     - za(k2,k12h)*za(k3,k5)*za(k3,k34h)*zb(k1,k6)*zb(k1,k34h)*
     &    zb(k3,k4)*iza(k1,k12h)*izb(k1,k2) - za(k3,k5)*za(k3,k34h)*zb(
     &    k1,k6)*zb(k1,k34h)**2*zb(k3,k4)*izb(k1,k2)*izb(k2,k34h) - 2.D0
     &    *za(k3,k5)*za(k3,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k1,k12h)*zb(
     &    k3,k4)*izb(k1,k2)*izb(k2,k12h) + za(k3,k5)*za(k3,k34h)*zb(k1,
     &    k6)*zb(k1,k12h)**2*zb(k2,k34h)*zb(k3,k4)*izb(k1,k2)*izb(k2,
     &    k12h)**2 )
      triamp = triamp + y**2 * ( za(k1,k34h)**2*za(k2,k12h)**2*za(k3,k5
     &    )*za(k3,k12h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*
     &    iza(k1,k12h)**3*izb(k1,k2) - za(k1,k34h)**2*za(k2,k12h)*za(k3
     &    ,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,
     &    k12h)**3*izb(k1,k2) - 2.D0*za(k1,k34h)*za(k2,k34h)*za(k2,k12h
     &    )*za(k3,k5)*za(k3,k12h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) + za(k1,k34h)*za(k2,
     &    k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,
     &    k34h)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*za(k2,k12h)**2
     &    *za(k3,k5)*za(k3,k34h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k34h)*iza(
     &    k1,k2)*iza(k1,k12h)**2*izb(k1,k2) + za(k1,k34h)*za(k2,k12h)*
     &    za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*
     &    iza(k1,k12h)**2*izb(k1,k2) + za(k1,k34h)*za(k2,k12h)*za(k3,
     &    k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,
     &    k12h)**2*izb(k1,k2) + za(k2,k34h)**2*za(k3,k5)*za(k3,k34h)*
     &    zb(k1,k34h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)*
     &    izb(k1,k2) )
      triamp = triamp + y**2 * ( za(k2,k34h)**2*za(k3,k5)*za(k3,k12h)*
     &    zb(k1,k34h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*
     &    izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,k5)*za(k3,
     &    k34h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k12h)*izb(k1,k2) + za(k2,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,
     &    k34h)**2*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k34h) + za(k2,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,k34h)**2*zb(
     &    k3,k4)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) - 2.D0*
     &    za(k2,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(
     &    k2,k34h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k12h)**2 + 2.D0*za(k2,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,k34h)
     &    *zb(k1,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(
     &    k2,k12h) + za(k2,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,k12h)**2*
     &    zb(k2,k34h)**2*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k12h)**3 - za(k2,k34h)*za(k3,k5)*za(k3,k34h)*zb(k1,
     &    k12h)**2*zb(k2,k34h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,
     &    k2)*izb(k2,k12h)**2 )
      triamp = triamp + y**2 * (  - za(k2,k34h)*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k34h)*izb(k1,k2) - 
     &    za(k2,k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(
     &    k4,k34h)*iza(k1,k12h)*izb(k1,k2) - za(k2,k34h)*za(k3,k12h)*
     &    za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k34h)*iza(k1,k12h)*
     &    izb(k1,k2) - za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k6)*
     &    zb(k1,k34h)*zb(k4,k34h)*iza(k1,k12h)*izb(k1,k2) - za(k3,k34h)
     &    *za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)**2*zb(k4,k34h)*izb(k1,k2)*
     &    izb(k2,k34h) - za(k3,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)
     &    **2*zb(k4,k12h)*izb(k1,k2)*izb(k2,k12h) + 2.D0*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k34h)*zb(
     &    k4,k12h)*izb(k1,k2)*izb(k2,k12h)**2 - 2.D0*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*izb(k1,k2
     &    )*izb(k2,k12h) - za(k3,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h
     &    )**2*zb(k2,k34h)**2*zb(k4,k12h)*izb(k1,k2)*izb(k2,k12h)**3 + 
     &    za(k3,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)**2*zb(k2,k34h)*
     &    zb(k4,k34h)*izb(k1,k2)*izb(k2,k12h)**2 )
      triamp = triamp + y**3 * (  - za(k1,k34h)**3*za(k2,k12h)**2*za(k3
     &    ,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1
     &    ,k2)*iza(k1,k12h)**4*izb(k1,k2) + 2.D0*za(k1,k34h)**2*za(k2,
     &    k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)**3*izb(k1,k2) + za(
     &    k1,k34h)**2*za(k2,k12h)**2*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h
     &    )*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)**3*izb(k1,
     &    k2) + za(k1,k34h)**2*za(k2,k12h)**2*za(k3,k12h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)**
     &    3*izb(k1,k2) - za(k1,k34h)*za(k2,k34h)**2*za(k3,k12h)*za(k5,
     &    k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k12h)**2*izb(k1,k2) - 2.D0*za(k1,k34h)*za(k2,k34h)*za(k2,k12h
     &    )*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)
     &    *iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) - 2.D0*za(k1,k34h)*za(
     &    k2,k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) )
      triamp = triamp + y**3 * (  - za(k1,k34h)*za(k2,k12h)**2*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)*iza(k1,k12h)**2*izb(k1,k2) + za(k2,k34h)**2*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*
     &    iza(k1,k34h)*izb(k1,k2) + za(k2,k34h)**2*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k12h)*izb(k1,k2) + za(k2,k34h)**2*za(k3,k12h)*za(k5,k34h)*zb(
     &    k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*izb(
     &    k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*
     &    izb(k1,k2) - za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)
     &    **2*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)
     &    *izb(k2,k12h)**2 + za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2
     &    ,k34h) + za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)**2*
     &    zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) )
      triamp = triamp + y**3 * ( za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)**2*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k12h) + 2.D0*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1
     &    ,k34h)*zb(k1,k12h)*zb(k2,k34h)**2*zb(k4,k12h)*zb(k6,k12h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**3 - 2.D0*za(k2,k34h)*za(
     &    k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k34h)*zb(
     &    k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 - 
     &    2.D0*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,
     &    k12h)*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,
     &    k2)*izb(k2,k12h)**2 + 2.D0*za(k2,k34h)*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*
     &    izb(k1,k2)*izb(k2,k12h) - za(k2,k34h)*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k12h)**2*zb(k2,k34h)**3*zb(k4,k12h)*zb(k6,k12h)*iza(k1
     &    ,k2)*izb(k1,k2)*izb(k2,k12h)**4 + za(k2,k34h)*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k12h)**2*zb(k2,k34h)**2*zb(k4,k34h)*zb(k6,k12h
     &    )*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**3 )
      triamp = triamp + y**3 * ( za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k12h)**2*zb(k2,k34h)**2*zb(k4,k12h)*zb(k6,k34h)*iza(k1,
     &    k2)*izb(k1,k2)*izb(k2,k12h)**3 - za(k2,k34h)*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k12h)**2*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 )
      triamp = triamp + x * ( za(k1,k12h)*za(k2,k34h)*za(k3,k5)*za(k3,
     &    k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k3,k4)*iza(k1,k34h)**2*izb(k1,
     &    k2) - za(k2,k34h)*za(k3,k5)*za(k3,k12h)*zb(k1,k6)*zb(k1,k12h)
     &    *zb(k3,k4)*iza(k1,k34h)*izb(k1,k2) - za(k2,k12h)*za(k3,k5)*
     &    za(k3,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k3,k4)*iza(k1,k34h)*izb(
     &    k1,k2) - za(k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,k6)*zb(k1,
     &    k12h)*zb(k3,k4)*iza(k1,k12h)*izb(k1,k2) + za(k3,k5)*za(k3,
     &    k12h)*zb(k1,k6)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k3,k4)*izb(k1,
     &    k2)*izb(k2,k34h)**2 - 2.D0*za(k3,k5)*za(k3,k12h)*zb(k1,k6)*
     &    zb(k1,k34h)*zb(k1,k12h)*zb(k3,k4)*izb(k1,k2)*izb(k2,k34h) - 
     &    za(k3,k5)*za(k3,k12h)*zb(k1,k6)*zb(k1,k12h)**2*zb(k3,k4)*izb(
     &    k1,k2)*izb(k2,k12h) )
      triamp = triamp + x*y * (  - za(k1,k34h)*za(k2,k12h)**2*za(k3,k5)
     &    *za(k3,k12h)*zb(k1,k34h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*
     &    iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*za(k2,k12h)**2*za(k3
     &    ,k5)*za(k3,k12h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)
     &    *iza(k1,k12h)**2*izb(k1,k2) + za(k1,k34h)*za(k2,k12h)*za(k3,
     &    k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,
     &    k12h)**2*izb(k1,k2) + za(k1,k34h)*za(k2,k12h)*za(k3,k12h)*za(
     &    k5,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,k12h)**2*
     &    izb(k1,k2) - za(k1,k12h)*za(k2,k34h)**2*za(k3,k5)*za(k3,k34h)
     &    *zb(k1,k34h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)**2
     &    *izb(k1,k2) - za(k1,k12h)*za(k2,k34h)**2*za(k3,k5)*za(k3,k34h
     &    )*zb(k1,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)**
     &    2*izb(k1,k2) + za(k1,k12h)*za(k2,k34h)*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k34h)**2*izb(k1,k2
     &    ) + za(k1,k12h)*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k6)
     &    *zb(k1,k12h)*zb(k4,k34h)*iza(k1,k34h)**2*izb(k1,k2) )
      triamp = triamp + x*y * ( za(k2,k34h)**2*za(k3,k5)*za(k3,k12h)*
     &    zb(k1,k34h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)*
     &    izb(k1,k2) + za(k2,k34h)**2*za(k3,k5)*za(k3,k12h)*zb(k1,k12h)
     &    *zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + 2.
     &    D0*za(k2,k34h)*za(k2,k12h)*za(k3,k5)*za(k3,k34h)*zb(k1,k34h)*
     &    zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + 2.D
     &    0*za(k2,k34h)*za(k2,k12h)*za(k3,k5)*za(k3,k34h)*zb(k1,k12h)*
     &    zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + 2.D
     &    0*za(k2,k34h)*za(k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,k34h)*
     &    zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + 2.D
     &    0*za(k2,k34h)*za(k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,k12h)*
     &    zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) - 
     &    za(k2,k34h)*za(k3,k5)*za(k3,k12h)*zb(k1,k34h)**2*zb(k2,k12h)*
     &    zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2
     &     + za(k2,k34h)*za(k3,k5)*za(k3,k12h)*zb(k1,k34h)**2*zb(k3,k4)
     &    *zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) )
      triamp = triamp + x*y * ( 2.D0*za(k2,k34h)*za(k3,k5)*za(k3,k12h)*
     &    zb(k1,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(
     &    k1,k2)*izb(k2,k34h) + 2.D0*za(k2,k34h)*za(k3,k5)*za(k3,k12h)*
     &    zb(k1,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*izb(
     &    k1,k2)*izb(k2,k12h) - za(k2,k34h)*za(k3,k5)*za(k3,k12h)*zb(k1
     &    ,k12h)**2*zb(k2,k34h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*izb(k1
     &    ,k2)*izb(k2,k12h)**2 + za(k2,k34h)*za(k3,k5)*za(k3,k12h)*zb(
     &    k1,k12h)**2*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(
     &    k2,k12h) - za(k2,k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(
     &    k1,k34h)*zb(k4,k12h)*iza(k1,k34h)*izb(k1,k2) - za(k2,k34h)*
     &    za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k34h)*
     &    iza(k1,k34h)*izb(k1,k2) - za(k2,k34h)*za(k3,k12h)*za(k5,k34h)
     &    *zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k34h)*izb(k1,k2) - 
     &    za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(
     &    k4,k34h)*iza(k1,k34h)*izb(k1,k2) - za(k2,k34h)*za(k3,k12h)*
     &    za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k12h)*
     &    izb(k1,k2) )
      triamp = triamp + x*y * (  - za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*
     &    zb(k1,k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,k12h)*izb(k1,k2) + 
     &    za(k2,k12h)**2*za(k3,k5)*za(k3,k34h)*zb(k1,k34h)*zb(k3,k4)*
     &    zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + za(k2,k12h)
     &    **2*za(k3,k5)*za(k3,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) - za(k2,k12h)*za(k3,k5)*
     &    za(k3,k34h)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k3,k4)*zb(k6,k34h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + za(k2,k12h)*za(k3,k5)
     &    *za(k3,k34h)*zb(k1,k34h)**2*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*
     &    izb(k1,k2)*izb(k2,k34h) + 2.D0*za(k2,k12h)*za(k3,k5)*za(k3,
     &    k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2
     &    )*izb(k1,k2)*izb(k2,k34h) + 2.D0*za(k2,k12h)*za(k3,k5)*za(k3,
     &    k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2
     &    )*izb(k1,k2)*izb(k2,k12h) - za(k2,k12h)*za(k3,k5)*za(k3,k34h)
     &    *zb(k1,k12h)**2*zb(k2,k34h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*
     &    izb(k1,k2)*izb(k2,k12h)**2 )
      triamp = triamp + x*y * ( za(k2,k12h)*za(k3,k5)*za(k3,k34h)*zb(k1
     &    ,k12h)**2*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k12h) - za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,
     &    k34h)*zb(k4,k12h)*iza(k1,k34h)*izb(k1,k2) - za(k2,k12h)*za(k3
     &    ,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k34h)*iza(k1,
     &    k34h)*izb(k1,k2) - za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,
     &    k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k12h)*izb(k1,k2) - za(k2,
     &    k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,
     &    k34h)*iza(k1,k12h)*izb(k1,k2) - za(k2,k12h)*za(k3,k12h)*za(k5
     &    ,k34h)*zb(k1,k6)*zb(k1,k34h)*zb(k4,k12h)*iza(k1,k12h)*izb(k1,
     &    k2) - za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,
     &    k12h)*zb(k4,k34h)*iza(k1,k12h)*izb(k1,k2) + za(k3,k34h)*za(k5
     &    ,k12h)*zb(k1,k6)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k4,k34h)*izb(
     &    k1,k2)*izb(k2,k34h)**2 - za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*
     &    zb(k1,k34h)**2*zb(k4,k12h)*izb(k1,k2)*izb(k2,k34h) - 2.D0*za(
     &    k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,
     &    k34h)*izb(k1,k2)*izb(k2,k34h) )
      triamp = triamp + x*y * (  - 2.D0*za(k3,k34h)*za(k5,k12h)*zb(k1,
     &    k6)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k12h)*izb(k1,k2)*izb(k2,
     &    k12h) + za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)**2*zb(
     &    k2,k34h)*zb(k4,k12h)*izb(k1,k2)*izb(k2,k12h)**2 - za(k3,k34h)
     &    *za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)**2*zb(k4,k34h)*izb(k1,k2)*
     &    izb(k2,k12h) + za(k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)
     &    **2*zb(k2,k12h)*zb(k4,k34h)*izb(k1,k2)*izb(k2,k34h)**2 - za(
     &    k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)**2*zb(k4,k12h)*
     &    izb(k1,k2)*izb(k2,k34h) - 2.D0*za(k3,k12h)*za(k5,k34h)*zb(k1,
     &    k6)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*izb(k1,k2)*izb(k2,
     &    k34h) - 2.D0*za(k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k34h)*
     &    zb(k1,k12h)*zb(k4,k12h)*izb(k1,k2)*izb(k2,k12h) + za(k3,k12h)
     &    *za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)**2*zb(k2,k34h)*zb(k4,k12h)
     &    *izb(k1,k2)*izb(k2,k12h)**2 - za(k3,k12h)*za(k5,k34h)*zb(k1,
     &    k6)*zb(k1,k12h)**2*zb(k4,k34h)*izb(k1,k2)*izb(k2,k12h) )
      triamp = triamp + x*y**2 * ( za(k1,k34h)**2*za(k2,k12h)**2*za(k3,
     &    k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k12h)**3*izb(k1,k2) + za(k1,k34h)**2*za(k2,k12h)**
     &    2*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)
     &    *iza(k1,k2)*iza(k1,k12h)**3*izb(k1,k2) + za(k1,k34h)**2*za(k2
     &    ,k12h)**2*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(
     &    k6,k34h)*iza(k1,k2)*iza(k1,k12h)**3*izb(k1,k2) - 2.D0*za(k1,
     &    k34h)*za(k2,k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,
     &    k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)**2*izb(
     &    k1,k2) - 2.D0*za(k1,k34h)*za(k2,k34h)*za(k2,k12h)*za(k3,k12h)
     &    *za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*
     &    iza(k1,k12h)**2*izb(k1,k2) - 2.D0*za(k1,k34h)*za(k2,k34h)*za(
     &    k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(
     &    k6,k34h)*iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*
     &    za(k2,k12h)**2*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k34h
     &    )*zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) )
      triamp = triamp + x*y**2 * (  - za(k1,k34h)*za(k2,k12h)**2*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,
     &    k2)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*za(k2,k12h)**2*
     &    za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*za(k2,
     &    k12h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(
     &    k6,k12h)*iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*
     &    za(k2,k12h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h
     &    )*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,
     &    k34h)*za(k2,k12h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) - 
     &    za(k1,k12h)*za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h
     &    )*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)**2*izb(k1,
     &    k2) - za(k1,k12h)*za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*zb(
     &    k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)**2*
     &    izb(k1,k2) )
      triamp = triamp + x*y**2 * (  - za(k1,k12h)*za(k2,k34h)**2*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)*iza(k1,k34h)**2*izb(k1,k2) + za(k2,k34h)**2*za(k3,k34h)*
     &    za(k5,k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*
     &    iza(k1,k34h)*izb(k1,k2) + za(k2,k34h)**2*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k34h)*izb(k1,k2) + za(k2,k34h)**2*za(k3,k34h)*za(k5,k12h)*zb(
     &    k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)*izb(
     &    k1,k2) + za(k2,k34h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*
     &    zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + 
     &    za(k2,k34h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h
     &    )*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k34h
     &    )**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k34h)**2*za(
     &    k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(
     &    k1,k2)*iza(k1,k12h)*izb(k1,k2) )
      triamp = triamp + x*y**2 * ( za(k2,k34h)**2*za(k3,k12h)*za(k5,
     &    k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k12h)*izb(k1,k2) + za(k2,k34h)**2*za(k3,k12h)*za(k5,k12h)*zb(
     &    k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*izb(
     &    k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)*
     &    izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*
     &    iza(k1,k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k12h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*
     &    za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,
     &    k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) )
      triamp = triamp + x*y**2 * ( 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,
     &    k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k12h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*
     &    za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,
     &    k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) - za(k2,k34h)*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k4,k34h)*zb(
     &    k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + za(k2,k34h)*
     &    za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)**2*zb(k4,k34h)*zb(k6,k12h
     &    )*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) + za(k2,k34h)*za(k3,k34h
     &    )*za(k5,k12h)*zb(k1,k34h)**2*zb(k4,k12h)*zb(k6,k34h)*iza(k1,
     &    k2)*izb(k1,k2)*izb(k2,k34h) + za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k34h)**2*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(
     &    k1,k2)*izb(k2,k12h) - 2.D0*za(k2,k34h)*za(k3,k34h)*za(k5,k12h
     &    )*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)
     &    *iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 )
      triamp = triamp + x*y**2 * ( 2.D0*za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)*izb(k1,k2)*izb(k2,k34h) + 2.D0*za(k2,k34h)*za(k3,k34h)*
     &    za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) + 2.D0*za(k2,k34h)*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,
     &    k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) + za(k2,k34h)*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k12h)**2*zb(k2,k34h)**2*zb(k4,k12h)*
     &    zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**3 - za(k2,
     &    k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)**2*zb(k2,k34h)*zb(
     &    k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 - 
     &    za(k2,k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)**2*zb(k2,k34h
     &    )*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)
     &    **2 + za(k2,k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)**2*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) - za(
     &    k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k2,k12h)*
     &    zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2
     &     )
      triamp = triamp + x*y**2 * ( za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*
     &    zb(k1,k34h)**2*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k34h) + za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h
     &    )**2*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k34h) + za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)**2*
     &    zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) - 
     &    2.D0*za(k2,k34h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,
     &    k12h)*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,
     &    k2)*izb(k2,k12h)**2 + 2.D0*za(k2,k34h)*za(k3,k12h)*za(k5,k34h
     &    )*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*
     &    izb(k1,k2)*izb(k2,k34h) + 2.D0*za(k2,k34h)*za(k3,k12h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,
     &    k2)*izb(k1,k2)*izb(k2,k12h) + 2.D0*za(k2,k34h)*za(k3,k12h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) + za(k2,k34h)*za(k3,k12h)*
     &    za(k5,k34h)*zb(k1,k12h)**2*zb(k2,k34h)**2*zb(k4,k12h)*zb(k6,
     &    k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**3 )
      triamp = triamp + x*y**2 * (  - za(k2,k34h)*za(k3,k12h)*za(k5,
     &    k34h)*zb(k1,k12h)**2*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(
     &    k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 - za(k2,k34h)*za(k3,k12h)*
     &    za(k5,k34h)*zb(k1,k12h)**2*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k34h
     &    )*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 + za(k2,k34h)*za(k3,
     &    k12h)*za(k5,k34h)*zb(k1,k12h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(
     &    k1,k2)*izb(k1,k2)*izb(k2,k12h) + za(k2,k12h)**2*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*
     &    iza(k1,k12h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k12h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,k34h)*za(k5,k34h)*zb(
     &    k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*izb(
     &    k1,k2) - za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)**2*
     &    zb(k2,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k34h)**2 + za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)**2*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2
     &    ,k34h) )
      triamp = triamp + x*y**2 * ( za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)**2*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k34h) + za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h
     &    )**2*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k12h) - 2.D0*za(k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*
     &    zb(k1,k12h)*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*
     &    izb(k1,k2)*izb(k2,k12h)**2 + 2.D0*za(k2,k12h)*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k34h)*iza(
     &    k1,k2)*izb(k1,k2)*izb(k2,k34h) + 2.D0*za(k2,k12h)*za(k3,k34h)
     &    *za(k5,k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) + 2.D0*za(k2,k12h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,
     &    k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) + za(k2,k12h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k12h)**2*zb(k2,k34h)**2*zb(k4,k12h)*
     &    zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**3 - za(k2,
     &    k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k12h)**2*zb(k2,k34h)*zb(
     &    k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 )
      triamp = triamp + x*y**2 * (  - za(k2,k12h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k12h)**2*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k34h)*iza(
     &    k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 + za(k2,k12h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k12h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)
     &    *izb(k1,k2)*izb(k2,k12h) )
      triamp = triamp + x**2 * ( za(k1,k12h)**2*za(k2,k34h)**2*za(k3,k5
     &    )*za(k3,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*
     &    iza(k1,k34h)**3*izb(k1,k2) - za(k1,k12h)**2*za(k2,k34h)*za(k3
     &    ,k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,
     &    k34h)**3*izb(k1,k2) - za(k1,k12h)*za(k2,k34h)**2*za(k3,k5)*
     &    za(k3,k12h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*iza(
     &    k1,k34h)**2*izb(k1,k2) - 2.D0*za(k1,k12h)*za(k2,k34h)*za(k2,
     &    k12h)*za(k3,k5)*za(k3,k34h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k12h)
     &    *iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) + za(k1,k12h)*za(k2,
     &    k34h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,
     &    k12h)*iza(k1,k34h)**2*izb(k1,k2) + za(k1,k12h)*za(k2,k34h)*
     &    za(k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*
     &    iza(k1,k34h)**2*izb(k1,k2) + za(k1,k12h)*za(k2,k12h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,
     &    k34h)**2*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,k5)*
     &    za(k3,k12h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*iza(
     &    k1,k34h)*izb(k1,k2) )
      triamp = triamp + x**2 * (  - za(k2,k34h)*za(k3,k12h)*za(k5,k12h)
     &    *zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,k34h)*izb(k1,k2) + 
     &    za(k2,k12h)**2*za(k3,k5)*za(k3,k34h)*zb(k1,k12h)*zb(k3,k4)*
     &    zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k12h)
     &    **2*za(k3,k5)*za(k3,k12h)*zb(k1,k12h)*zb(k3,k4)*zb(k6,k12h)*
     &    iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + za(k2,k12h)*za(k3,k5)*
     &    za(k3,k12h)*zb(k1,k34h)**2*zb(k2,k12h)**2*zb(k3,k4)*zb(k6,
     &    k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**3 - za(k2,k12h)*za(
     &    k3,k5)*za(k3,k12h)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k3,k4)*zb(k6
     &    ,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 - 2.D0*za(k2,
     &    k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,
     &    k12h)*zb(k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h
     &    )**2 + 2.D0*za(k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,k34h)*zb(
     &    k1,k12h)*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k34h) + za(k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(k1,k12h)**2*zb(
     &    k3,k4)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) )
      triamp = triamp + x**2 * ( za(k2,k12h)*za(k3,k5)*za(k3,k12h)*zb(
     &    k1,k12h)**2*zb(k3,k4)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(
     &    k2,k12h) - za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k6)*zb(
     &    k1,k12h)*zb(k4,k12h)*iza(k1,k34h)*izb(k1,k2) - za(k2,k12h)*
     &    za(k3,k12h)*za(k5,k34h)*zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*
     &    iza(k1,k34h)*izb(k1,k2) - za(k2,k12h)*za(k3,k12h)*za(k5,k12h)
     &    *zb(k1,k6)*zb(k1,k12h)*zb(k4,k12h)*iza(k1,k12h)*izb(k1,k2) - 
     &    za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k34h)**2*zb(k2,k12h)
     &    **2*zb(k4,k34h)*izb(k1,k2)*izb(k2,k34h)**3 + za(k3,k12h)*za(
     &    k5,k12h)*zb(k1,k6)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k4,k12h)*
     &    izb(k1,k2)*izb(k2,k34h)**2 + 2.D0*za(k3,k12h)*za(k5,k12h)*zb(
     &    k1,k6)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k12h)*zb(k4,k34h)*izb(k1
     &    ,k2)*izb(k2,k34h)**2 - 2.D0*za(k3,k12h)*za(k5,k12h)*zb(k1,k6)
     &    *zb(k1,k34h)*zb(k1,k12h)*zb(k4,k12h)*izb(k1,k2)*izb(k2,k34h)
     &     - za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*zb(k1,k12h)**2*zb(k4,
     &    k34h)*izb(k1,k2)*izb(k2,k34h) )
      triamp = triamp + x**2 * (  - za(k3,k12h)*za(k5,k12h)*zb(k1,k6)*
     &    zb(k1,k12h)**2*zb(k4,k12h)*izb(k1,k2)*izb(k2,k12h) )
      triamp = triamp + x**2*y * (  - za(k1,k34h)*za(k2,k12h)**2*za(k3,
     &    k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*za(k2,k12h)**2*
     &    za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*
     &    iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) - za(k1,k34h)*za(k2,
     &    k12h)**2*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(
     &    k6,k34h)*iza(k1,k2)*iza(k1,k12h)**2*izb(k1,k2) + za(k1,k12h)
     &    **2*za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,
     &    k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)**3*izb(k1,k2) + za(
     &    k1,k12h)**2*za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k12h
     &    )*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)**3*izb(k1,
     &    k2) + za(k1,k12h)**2*za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)**
     &    3*izb(k1,k2) - za(k1,k12h)*za(k2,k34h)**2*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,
     &    k34h)**2*izb(k1,k2) )
      triamp = triamp + x**2*y * (  - za(k1,k12h)*za(k2,k34h)**2*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k34h)**2*izb(k1,k2) - za(k1,k12h)*za(k2,k34h)**2*
     &    za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) - za(k1,k12h)*za(k2,
     &    k34h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h)*zb(
     &    k6,k12h)*iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) - za(k1,k12h)*
     &    za(k2,k34h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h
     &    )*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) - za(k1,
     &    k12h)*za(k2,k34h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(
     &    k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) - 
     &    2.D0*za(k1,k12h)*za(k2,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,
     &    k34h)**2*izb(k1,k2) - 2.D0*za(k1,k12h)*za(k2,k34h)*za(k2,k12h
     &    )*za(k3,k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)
     &    *iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) )
      triamp = triamp + x**2*y * (  - 2.D0*za(k1,k12h)*za(k2,k34h)*za(
     &    k2,k12h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(
     &    k6,k34h)*iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) + za(k2,k34h)
     &    **2*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,
     &    k12h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k34h)**2*za(
     &    k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(
     &    k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k34h)**2*za(k3,k12h)*
     &    za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*
     &    iza(k1,k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*
     &    za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*
     &    iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,
     &    k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,
     &    k34h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*
     &    za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h)*
     &    zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) )
      triamp = triamp + x**2*y * ( 2.D0*za(k2,k34h)*za(k2,k12h)*za(k3,
     &    k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)*
     &    za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,
     &    k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,
     &    k12h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + 2.D0*za(k2,k34h)*
     &    za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k34h)*
     &    zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + 2.D0*za(k2,
     &    k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,
     &    k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + za(k2,
     &    k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)**2*zb(k2,k12h)**2*
     &    zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**3
     &     - za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)**2*zb(k2,
     &    k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k34h)**2 )
      triamp = triamp + x**2*y * (  - za(k2,k34h)*za(k3,k12h)*za(k5,
     &    k12h)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(
     &    k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + za(k2,k34h)*za(k3,k12h)*
     &    za(k5,k12h)*zb(k1,k34h)**2*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)
     &    *izb(k1,k2)*izb(k2,k34h) - 2.D0*za(k2,k34h)*za(k3,k12h)*za(k5
     &    ,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k12h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + 2.D0*za(k2,k34h
     &    )*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k34h)
     &    *zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) + 2.D0*za(k2,
     &    k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,
     &    k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) + 2.D0*
     &    za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*
     &    zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) - 
     &    za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)**2*zb(k2,k34h
     &    )*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h)
     &    **2 )
      triamp = triamp + x**2*y * ( za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*
     &    zb(k1,k12h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k34h) + za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h
     &    )**2*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k12h) + za(k2,k34h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)**2*
     &    zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) + 
     &    za(k2,k12h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k12h
     &    )*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k12h
     &    )**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,
     &    k12h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k12h)**2*za(
     &    k3,k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(
     &    k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,k34h)*
     &    za(k5,k12h)*zb(k1,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*
     &    iza(k1,k12h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,
     &    k12h)*izb(k1,k2) )
      triamp = triamp + x**2*y * ( za(k2,k12h)**2*za(k3,k34h)*za(k5,
     &    k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k12h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,k12h)*za(k5,k34h)*zb(
     &    k1,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)*izb(
     &    k1,k2) + za(k2,k12h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*
     &    zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + 
     &    za(k2,k12h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h
     &    )*zb(k6,k34h)*iza(k1,k2)*iza(k1,k12h)*izb(k1,k2) + za(k2,k12h
     &    )*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)**2*zb(k2,k12h)**2*zb(k4
     &    ,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**3 - 
     &    za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)**2*zb(k2,k12h
     &    )*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)
     &    **2 - za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)**2*zb(
     &    k2,k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2
     &    ,k34h)**2 + za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)**
     &    2*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)
     &     )
      triamp = triamp + x**2*y * (  - 2.D0*za(k2,k12h)*za(k3,k34h)*za(
     &    k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k12h)*zb(k4,k34h)*zb(
     &    k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + 2.D0*za(k2,
     &    k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,
     &    k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) + 2.D0*
     &    za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*
     &    zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) + 
     &    2.D0*za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,
     &    k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k12h) - za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)**2*
     &    zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k12h)**2 + za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,
     &    k12h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2
     &    ,k34h) + za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)**2*
     &    zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) + 
     &    za(k2,k12h)*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)**2*zb(k4,k12h
     &    )*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k12h) )
      triamp = triamp + x**2*y * ( za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*
     &    zb(k1,k34h)**2*zb(k2,k12h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)*izb(k1,k2)*izb(k2,k34h)**3 - za(k2,k12h)*za(k3,k12h)*za(
     &    k5,k34h)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k4,k34h)*zb(k6,k12h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 - za(k2,k12h)*za(k3,
     &    k12h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k2,k12h)*zb(k4,k12h)*zb(
     &    k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + za(k2,k12h)*
     &    za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k4,k12h)*zb(k6,k12h
     &    )*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) - 2.D0*za(k2,k12h)*za(k3
     &    ,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k12h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + 2.D0
     &    *za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k12h)*
     &    zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) + 
     &    2.D0*za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,
     &    k12h)*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k34h) )
      triamp = triamp + x**2*y * ( 2.D0*za(k2,k12h)*za(k3,k12h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,
     &    k2)*izb(k1,k2)*izb(k2,k12h) - za(k2,k12h)*za(k3,k12h)*za(k5,
     &    k34h)*zb(k1,k12h)**2*zb(k2,k34h)*zb(k4,k12h)*zb(k6,k12h)*iza(
     &    k1,k2)*izb(k1,k2)*izb(k2,k12h)**2 + za(k2,k12h)*za(k3,k12h)*
     &    za(k5,k34h)*zb(k1,k12h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)
     &    *izb(k1,k2)*izb(k2,k34h) + za(k2,k12h)*za(k3,k12h)*za(k5,k34h
     &    )*zb(k1,k12h)**2*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2
     &    )*izb(k2,k12h) + za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,
     &    k12h)**2*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2
     &    ,k12h) )
      triamp = triamp + x**3 * (  - za(k1,k12h)**3*za(k2,k34h)**2*za(k3
     &    ,k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1
     &    ,k2)*iza(k1,k34h)**4*izb(k1,k2) + za(k1,k12h)**2*za(k2,k34h)
     &    **2*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,
     &    k12h)*iza(k1,k2)*iza(k1,k34h)**3*izb(k1,k2) + za(k1,k12h)**2*
     &    za(k2,k34h)**2*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h
     &    )*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)**3*izb(k1,k2) + 2.D0*
     &    za(k1,k12h)**2*za(k2,k34h)*za(k2,k12h)*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)
     &    **3*izb(k1,k2) - za(k1,k12h)*za(k2,k34h)**2*za(k3,k12h)*za(k5
     &    ,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,
     &    k34h)**2*izb(k1,k2) - 2.D0*za(k1,k12h)*za(k2,k34h)*za(k2,k12h
     &    )*za(k3,k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)
     &    *iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) - 2.D0*za(k1,k12h)*za(
     &    k2,k34h)*za(k2,k12h)*za(k3,k12h)*za(k5,k34h)*zb(k1,k12h)*zb(
     &    k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k2) )
      triamp = triamp + x**3 * (  - za(k1,k12h)*za(k2,k12h)**2*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k34h)**2*izb(k1,k2) + 2.D0*za(k2,k34h)*za(k2,k12h)
     &    *za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*
     &    iza(k1,k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,
     &    k34h)*za(k5,k12h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,
     &    k2)*iza(k1,k34h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,k12h)*za(
     &    k5,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(
     &    k1,k34h)*izb(k1,k2) + za(k2,k12h)**2*za(k3,k12h)*za(k5,k12h)*
     &    zb(k1,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*iza(k1,k12h)*
     &    izb(k1,k2) - za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)
     &    **2*zb(k2,k12h)**3*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,
     &    k2)*izb(k2,k34h)**4 + za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(
     &    k1,k34h)**2*zb(k2,k12h)**2*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)
     &    *izb(k1,k2)*izb(k2,k34h)**3 + za(k2,k12h)*za(k3,k12h)*za(k5,
     &    k12h)*zb(k1,k34h)**2*zb(k2,k12h)**2*zb(k4,k12h)*zb(k6,k34h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**3 )
      triamp = triamp + x**3 * (  - za(k2,k12h)*za(k3,k12h)*za(k5,k12h)
     &    *zb(k1,k34h)**2*zb(k2,k12h)*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2
     &    )*izb(k1,k2)*izb(k2,k34h)**2 + 2.D0*za(k2,k12h)*za(k3,k12h)*
     &    za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k2,k12h)**2*zb(k4,k34h
     &    )*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**3 - 2.D0*
     &    za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*
     &    zb(k2,k12h)*zb(k4,k34h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k34h)**2 - 2.D0*za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*
     &    zb(k1,k34h)*zb(k1,k12h)*zb(k2,k12h)*zb(k4,k12h)*zb(k6,k34h)*
     &    iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + 2.D0*za(k2,k12h)*za(
     &    k3,k12h)*za(k5,k12h)*zb(k1,k34h)*zb(k1,k12h)*zb(k4,k12h)*zb(
     &    k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) - za(k2,k12h)*za(
     &    k3,k12h)*za(k5,k12h)*zb(k1,k12h)**2*zb(k2,k12h)*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h)**2 + za(k2,
     &    k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h)**2*zb(k4,k34h)*zb(
     &    k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,k34h) )
      triamp = triamp + x**3 * ( za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*
     &    zb(k1,k12h)**2*zb(k4,k12h)*zb(k6,k34h)*iza(k1,k2)*izb(k1,k2)*
     &    izb(k2,k34h) + za(k2,k12h)*za(k3,k12h)*za(k5,k12h)*zb(k1,k12h
     &    )**2*zb(k4,k12h)*zb(k6,k12h)*iza(k1,k2)*izb(k1,k2)*izb(k2,
     &    k12h) )

      apm=triamp*(-half*ci)/(sprod(k3,k4)*sprod(k5,k6))

      return
      end
