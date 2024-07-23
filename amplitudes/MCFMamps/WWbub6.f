      function WWbub6(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass)
      use mod_types; use mod_consts_dp; use common_def
!----- C. Williams June 2011
!----- This function returns the coefficient of the log(s_gg) with j1=p and j2=m 
      implicit none 
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12),mass
      integer j1,j2,j3,j4,j5,j6
      integer del1,del2
      complex(dp)::  sym,antisym,WWbub6
      parameter(del1=7,del2=8)
      

      WWbub6= 
     &+sym(j1,j2,j3,j4,j5,j6,del1,del2,za,zb,mass)
     &+sym(j1,j2,j3,j4,j5,j6,del2,del1,za,zb,mass)
     &+antisym(j1,j2,j3,j4,j5,j6,del1,del2,za,zb,mass)
     &-antisym(j1,j2,j3,j4,j5,j6,del2,del1,za,zb,mass)
      
!---- W propagators 
    
      WWbub6=WWbub6/(sprod(j3,j4)*sprod(j5,j6))
      return 
      end


      function sym(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,mass) 
      use mod_types; use mod_consts_dp; use common_def
!---- symmetric piece of coefficient 
      implicit none 
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: mass
      integer j1,j2,j3,j4,j5,j6,j7,j8 
      complex(dp) :: masssq , X, sym

      masssq=mass**2

      X=-masssq+za(j5,j6)*zb(j6,j5) 

      sym=  -((X*za(j2,j3)*za(j2,j5)*za(j2,j7)*za(j7,j8)**2*zb(j2,j1)*
     &     zb(j7,j4)*zb(j7,j6)*
     &      (-(za(j1,j2)*zb(j7,j2)) + (za(j1,j2)*zb(j2,j1)*(-(za(j1,j5)*
     &     zb(j7,j5)) - za(j1,j6)*zb(j7,j6)))/X)*zb(j8,j1)*zb(j8,j7))/
     &     (za(j1,j2)*zb(j7,j2)*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*
     &     zb(j7,j2))*(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) 
     &     - za(j7,j8)**2*zb(j8,j7)**2)**2)) + 
     &(X*za(j2,j3)*za(j2,j5)*za(j2,j7)*za(j7,j8)**2*zb(j2,j1)*zb(j7,j6)*
     &     (za(j1,j2)*zb(j7,j1) + (za(j1,j2)*zb(j2,j1)*(-(za(j2,j5)*
     &    zb(j7,j5)) - za(j2,j6)*zb(j7,j6)))/X)*zb(j8,j4)*zb(j8,j7))/
     &   (za(j1,j2)*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*zb(j7,j2))*
     & (za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) 
     &    - za(j7,j8)**2*zb(j8,j7)**2)**2) + 
     &(X*za(j1,j3)*za(j2,j5)*za(j2,j7)*za(j7,j8)**2*zb(j2,j1)*zb(j7,j1)*
     &     (za(j1,j2)*zb(j7,j1) + (za(j1,j2)*zb(j2,j1)*(-(za(j2,j5)*
     &   zb(j7,j5)) - za(j2,j6)*zb(j7,j6)))/X)*(zb(j7,j6)*zb(j8,j4) 
     &   + zb(j7,j4)*zb(j8,j6))*
     &     zb(j8,j7))/(za(j1,j2)*zb(j7,j2)*(za(j1,j7)*zb(j7,j1) 
     &   + za(j2,j7)*zb(j7,j2))*
     &     (za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5)
     &     - za(j7,j8)**2*zb(j8,j7)**2)**2) + 
     &     (X*za(j1,j5)*za(j2,j7)*za(j7,j8)**2*zb(j2,j1)*zb(j7,j1)*
     &     (za(j1,j3)*zb(j7,j1) + za(j2,j3)*zb(j7,j2))*
     &     (za(j1,j2)*zb(j7,j1) + (za(j1,j2)*zb(j2,j1)*
     & (-(za(j2,j5)*zb(j7,j5)) - za(j2,j6)*zb(j7,j6)))/X)*
     &(zb(j7,j6)*zb(j8,j4) + zb(j7,j4)*zb(j8,j6))*
     &     zb(j8,j7))/(za(j1,j2)*zb(j7,j2)**2*(za(j1,j7)*zb(j7,j1) 
     &     + za(j2,j7)*zb(j7,j2))*
     &     (za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) 
     &- za(j7,j8)**2*zb(j8,j7)**2)**2) - 
     &(X*za(j2,j5)*za(j2,j7)*za(j7,j8)**2*zb(j2,j1)*(za(j1,j3)*zb(j2,j1)
     &+ (za(j1,j2)*zb(j2,j1)*(za(j3,j5)*zb(j5,j2)
     & + za(j3,j6)*zb(j6,j2)))/X)*zb(j7,j1)*
     &     zb(j7,j4)*zb(j7,j6)*zb(j8,j7)**2)/
     &   (zb(j7,j2)**2*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*zb(j7,j2))*
     &(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) - za(j7,j8)**2*
     & zb(j8,j7)**2)**2)  
     
      sym=sym
     & -(X*za(j1,j5)*za(j2,j7)*za(j7,j8)**2*zb(j2,j1)*(za(j1,j2)*
     & zb(j2,j1) + (za(j1,j2)*zb(j2,j1)*(za(j2,j5)*zb(j5,j2)
     & + za(j2,j6)*zb(j6,j2)))/X)*zb(j7,j1)*
     &     (za(j1,j3)*zb(j7,j1) + za(j2,j3)*zb(j7,j2))*
     &  zb(j7,j4)*zb(j7,j6)*zb(j8,j7)**2)/
     &   (za(j1,j2)*zb(j7,j2)**3*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*
     & zb(j7,j2))*(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) 
     & - za(j7,j8)**2*zb(j8,j7)**2)**2) + 
     &  (X*za(j2,j3)*za(j2,j5)*za(j2,j7)*za(j7,j8)**2*
     & zb(j2,j1)*zb(j6,j2)*zb(j7,j4)*
     &     (za(j1,j2)*zb(j7,j1) + (za(j1,j2)*zb(j2,j1)*
     &(-(za(j2,j5)*zb(j7,j5)) - za(j2,j6)*zb(j7,j6)))/X)*zb(j8,j7)**2)/
     &   (za(j1,j2)*zb(j7,j2)*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*
     &     zb(j7,j2))*(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) 
     & - za(j7,j8)**2*zb(j8,j7)**2)**2) - 
     &  (X*za(j1,j5)*za(j2,j7)*za(j7,j8)**2*zb(j2,j1)**2*(za(j1,j3)*
     & zb(j7,j1) + za(j2,j3)*zb(j7,j2))*zb(j7,j4)*zb(j7,j6)*
     &     (za(j1,j2)*zb(j7,j1) + (za(j1,j2)*zb(j2,j1)*(-(za(j2,j5)*
     & zb(j7,j5)) - za(j2,j6)*zb(j7,j6)))/X)*zb(j8,j7)**2)/
     &   (za(j1,j2)*zb(j7,j2)**3*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*
     & zb(j7,j2))*(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5)
     & - za(j7,j8)**2*zb(j8,j7)**2)**2) - 
     &  (X*za(j2,j7)*za(j3,j7)*za(j7,j8)**2*zb(j2,j1)**2*(za(j1,j5)*
     & zb(j7,j1) + za(j2,j5)*zb(j7,j2))*zb(j7,j4)*zb(j7,j6)*
     &     (za(j1,j2)*zb(j7,j1) + (za(j1,j2)*zb(j2,j1)*(-(za(j2,j5)*
     & zb(j7,j5)) - za(j2,j6)*zb(j7,j6)))/X)*zb(j8,j7)**2)/
     &   (zb(j7,j2)**2*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*zb(j7,j2))**2*
     & (-(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5)) + za(j7,j8)**2*
     &zb(j8,j7)**2)**2)

      return 
      end
      

      double complex function antisym(j1,j2,j3,j4,j5,j6,j7,j8,
     &     za,zb,mass) 
      use mod_types; use mod_consts_dp; use common_def
!---- antisymmetric piece of coefficient 
      implicit none 
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: mass
      integer j1,j2,j3,j4,j5,j6,j7,j8 

      antisym=
     &-((za(j2,j7)*za(j3,j5)*za(j7,j8)*zb(j2,j1)*zb(j6,j5)*
     & zb(j7,j1)*(za(j1,j5)*zb(j7,j1) + za(j2,j5)*zb(j7,j2))*
     & zb(j7,j4)*zb(j8,j7))/
     &     (zb(j7,j2)**2*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*zb(j7,j2))*
     &     (za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) - 
     & za(j7,j8)**2*zb(j8,j7)**2))) - 
     &  (za(j2,j7)*za(j7,j8)*zb(j2,j1)*zb(j4,j1)*(za(j1,j3)*zb(j7,j1) + 
     & za(j2,j3)*zb(j7,j2))*(za(j1,j5)*zb(j7,j1) + za(j2,j5)*zb(j7,j2))*
     & zb(j7,j6)*zb(j8,j7))/
     &   (zb(j7,j2)**2*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*zb(j7,j2))*
     & (za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) - 
     & za(j7,j8)**2*zb(j8,j7)**2)) + 
     &  (za(j2,j7)**2*za(j7,j8)*zb(j2,j1)**2*(za(j1,j3)*zb(j7,j1) 
     &+ za(j2,j3)*zb(j7,j2))*(za(j1,j5)*zb(j7,j1) + za(j2,j5)*zb(j7,j2))
     & *zb(j7,j4)*zb(j7,j6)*
     &     zb(j8,j7))/(2d0*zb(j7,j2)**2*(za(j1,j7)*zb(j7,j1)
     &+ za(j2,j7)*zb(j7,j2))**2*(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) 
     &     - za(j7,j8)**2*zb(j8,j7)**2))
     &   + (za(j1,j3)*za(j1,j5)*za(j2,j7)*za(j7,j8)*zb(j2,j1)**2*
     & zb(j7,j1)**2*zb(j7,j4)*zb(j7,j6)*zb(j8,j7))/
     &   (zb(j7,j2)**3*(za(j1,j7)*zb(j7,j1) + za(j2,j7)*zb(j7,j2))*
     &(za(j1,j2)*za(j5,j6)*zb(j2,j1)*zb(j6,j5) 
     &- za(j7,j8)**2*zb(j8,j7)**2))

      return 
      end
