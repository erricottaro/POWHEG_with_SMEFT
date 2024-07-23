      subroutine WWmbc(MCFMst,j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,bcoeff)
      use mod_types; use mod_consts_dp; use common_def
      implicit none
c----Bubble coefficients extracted from BDK 11.5, 11.8
c----after performing the transformation 
c    (1-->4)
c    (2-->3)
c    (3-->1)
c    (4-->2)
c    (5-->5)
c    (6-->6)
      integer j1,j2,j3,j4,j5,j6,j,k
      include 'blabels.f'
      character*9 MCFMst
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
      complex(dp) ::  bcoeff(8),zab2,izab2,WWbub6
      complex(dp) ::  iza(6,6),izb(6,6)
      real(dp) ::  isprod(6,6),t,t134ms34,t234ms34,t134ms56,t234ms56
      real(dp) ::  IDelta,mass
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      izab2(j1,j2,j3,j4)=cone/zab2(j1,j2,j3,j4)
      t(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j1,j3)
      do j=1,6
      do k=j+1,6
      isprod(j,k)=1d0/sprod(j,k)
      iza(j,k)=cone/za(j,k)
      izb(j,k)=cone/zb(j,k)
      isprod(k,j)=isprod(j,k)
      iza(k,j)=-iza(j,k)
      izb(k,j)=-izb(j,k)
      enddo
      enddo

      IDelta=sprod(j1,j2)**2+sprod(j3,j4)**2+sprod(j5,j6)**2
     & -2d0*sprod(j1,j2)*sprod(j3,j4)
     & -2d0*sprod(j3,j4)*sprod(j5,j6)
     & -2d0*sprod(j1,j2)*sprod(j5,j6)
      IDelta=1d0/IDelta

      t134ms34=t(j1,j3,j4)-sprod(j3,j4)
      t234ms34=t(j2,j3,j4)-sprod(j3,j4)
      t134ms56=t(j1,j3,j4)-sprod(j5,j6)
      t234ms56=t(j2,j3,j4)-sprod(j5,j6)

      do j=1,8
      bcoeff(j)=czero
      enddo
      
      if     (MCFMst.eq.'q+qb-g+g+') then

      bcoeff(b34)= + t134ms34**(-2) * (  - za(j4,j3)*za(j1,j5)**2*zb(j4
     &    ,j1)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + t134ms34**(-1) * ( za(j3,j1)*za(j1,j5
     &    )*za(j2,j5)*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j2)*
     &    za(j1,j5)**2*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j5)*
     &    za(j1,j5)*zb(j4,j1)*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-2) * (  - za(j4,j3)*za(j2
     &    ,j5)**2*zb(j4,j2)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-1) * (  - za(j3,j1)*za(j2
     &    ,j5)**2*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j2)*za(j1,
     &    j5)*za(j2,j5)*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j5)*
     &    za(j2,j5)*zb(j4,j2)*iza(j1,j2)**2*iza(j5,j6) )

      bcoeff(b56)= + t234ms56**(-2) * (  - za(j3,j1)**2*za(j5,j6)*zb(j1
     &    ,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b56) = bcoeff(b56) + t234ms56**(-1) * (  - za(j3,j1)**2*
     &    za(j2,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,j1)*za(
     &    j3,j2)*za(j1,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,
     &    j1)*za(j3,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-2) * (  - za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-1) * ( za(j3,j1)*za(j3,j2
     &    )*za(j2,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,j2)**2
     &    *za(j1,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,j2)*za(
     &    j3,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**2 )

      bcoeff(b134)= + t134ms34**(-2) * ( za(j4,j3)*za(j1,j5)**2*zb(j4,
     &    j1)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b134) = bcoeff(b134) + t134ms34**(-1) * (  - za(j3,j1)*za(
     &    j1,j5)*za(j2,j5)*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) - za(j3,
     &    j2)*za(j1,j5)**2*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) + za(j3,
     &    j5)*za(j1,j5)*zb(j4,j1)*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-2) * ( za(j3,j2)**2*za(
     &    j5,j6)*zb(j2,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-1) * (  - za(j3,j1)*za(
     &    j3,j2)*za(j2,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,
     &    j2)**2*za(j1,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,
     &    j2)*za(j3,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b134) = bcoeff(b134) - za(j3,j1)*za(j3,j5)*za(j2,j5)*iza(
     & j4,j3)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j2)*za(j3,j5)*za(j1,j5)*
     &    iza(j4,j3)*iza(j1,j2)**3*iza(j5,j6)

      bcoeff(b234)= + t234ms34**(-2) * ( za(j4,j3)*za(j2,j5)**2*zb(j4,
     &    j2)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b234) = bcoeff(b234) + t234ms34**(-1) * ( za(j3,j1)*za(j2,
     &    j5)**2*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j2)*za(j1,
     &    j5)*za(j2,j5)*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j5)*
     &    za(j2,j5)*zb(j4,j2)*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-2) * ( za(j3,j1)**2*za(
     &    j5,j6)*zb(j1,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-1) * ( za(j3,j1)**2*za(
     &    j2,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,j1)*za(j3,
     &    j2)*za(j1,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,j1)*
     &    za(j3,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b234) = bcoeff(b234) + za(j3,j1)*za(j3,j5)*za(j2,j5)*iza(
     & j4,j3)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j2)*za(j3,j5)*za(j1,j5)*
     &    iza(j4,j3)*iza(j1,j2)**3*iza(j5,j6)

      bcoeff(rat)= + t134ms34**(-1) * ( za(j4,j3)*za(j1,j5)**2*zb(j4,j1
     &    )**2*iza(j1,j2)**2*iza(j5,j6)*isprod(j4,j3) )
      bcoeff(rat) = bcoeff(rat) + t234ms34**(-1) * ( za(j4,j3)*za(j2,j5
     &    )**2*zb(j4,j2)**2*iza(j1,j2)**2*iza(j5,j6)*isprod(j4,j3) )
      bcoeff(rat) = bcoeff(rat) + t234ms56**(-1) * ( za(j3,j1)**2*za(j5
     &    ,j6)*zb(j1,j6)**2*iza(j4,j3)*iza(j1,j2)**2*isprod(j5,j6) )
      bcoeff(rat) = bcoeff(rat) + t134ms56**(-1) * ( za(j3,j2)**2*za(j5
     &    ,j6)*zb(j2,j6)**2*iza(j4,j3)*iza(j1,j2)**2*isprod(j5,j6) )
      bcoeff(rat) = bcoeff(rat) + za(j3,j5)**2*iza(j4,j3)*iza(j1,j2)**2
     & *iza(j5,j6) - zb(j4,j6)**2*iza(j1,j2)**2*izb(j4,j3)*izb(j5,j6)

c      if     ((case .eq. 'HWWint') .or. (case .eq. 'WWqqbr')
c     &   .or. (case .eq. 'HWWH+i') .or. (case .eq. 'ggWW4l')
c     &   .or. (case .eq. 'ggWWbx')
c     &   .or. (case .eq. 'ggVV4l') .or. (case .eq. 'ggVVbx')) then
      bcoeff(b12)=czero
      bcoeff(b12zm)=czero
      bcoeff(b12m)=czero
c      elseif ((case .eq. 'HZZint') .or. (case .eq. 'ZZlept')
c     &   .or. (case .eq. 'HZZH+i') .or. (case .eq. 'ggZZ4l')
c     &   .or. (case .eq. 'ggZZbx')) then
c      bcoeff(b12zm)=czero
c      bcoeff(b12)=bcoeff(b12zm)
c      bcoeff(b1)=-bcoeff(b12)-bcoeff(b34)-bcoeff(b56)
c     & -bcoeff(b134)-bcoeff(b234)
c      else
c      write(6,*) 'Unimplemented case in mbc.f'
c      stop
c      endif

      elseif (MCFMst.eq.'q+qb-g+g-') then

      bcoeff(b34)= + t134ms34**(-2) * ( za(j3,j1)**2*zb(j4,j3)*zb(j1,j6
     &    )**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2*sprod(j4,j3) )
      bcoeff(b34) = bcoeff(b34) + t134ms34**(-1) * ( 2.D0*za(j3,j1)**2*
     &    zb(j4,j3)*zb(j1,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 - 2.D0
     &    *za(j3,j1)*zb(j4,j3)*zb(j1,j6)*izb(j5,j6)*zab2(j3,j4,j1,j2)*
     &    zab2(j1,j4,j3,j6)*izab2(j1,j4,j3,j2)**3 )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-2) * ( za(j4,j3)*za(j2,j5
     &    )**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2*
     &     sprod(j4,j3) )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-1) * ( 2.D0*za(j4,j3)*za(
     &    j2,j5)**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2 - 2.D0
     &    *za(j4,j3)*za(j2,j5)*zb(j4,j2)*iza(j5,j6)*zab2(j1,j3,j2,j4)*
     &    zab2(j5,j3,j4,j2)*izab2(j1,j3,j4,j2)**3 )
      bcoeff(b34) = bcoeff(b34) + IDelta**2 * (  - 3.D0*zab2(j3,j1,j2,
     &    j4)*zab2(j2,j4,j3,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)
     &     *sprod(
     &    j4,j3) - 3.D0*zab2(j3,j1,j2,j4)*zab2(j2,j4,j3,j1)*zab2(j5,j1,
     &    j2,j6)*izab2(j1,j4,j3,j2)*sprod(j1,j2) + 3.D0*
     &     zab2(j3,j1,j2,j4)*
     &    zab2(j2,j4,j3,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)*
     &     sprod(j5,
     &    j6) - 3.D0*zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j3,j4,j2)*sprod(j4,j3) - 3.D0*zab2(j3,j2,j1,j4)*
     &    zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,j6)*izab2(j1,j3,j4,j2)*
     &     sprod(j1,
     &    j2) + 3.D0*zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j3,j4,j2)*sprod(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + IDelta * ( 2.D0*za(j4,j3)*za(j3,j1)*
     &    za(j5,j6)*zb(j4,j3)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3
     &    *t(j4,j3,j1) - 2.D0*za(j4,j3)*za(j3,j1)*za(j5,j6)*zb(j4,j3)*
     &    zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j2) + 2.D0*
     &    za(j4,j3)*za(j3,j5)*za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(j5,j6)*
     &    izab2(j1,j3,j4,j2)**3*t(j3,j4,j1) - 2.D0*za(j4,j3)*za(j3,j5)*
     &    za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j3,j4,j2)**3
     &    *t(j3,j4,j2) - za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(
     &    j1,j4,j3,j2)**2*sprod(j4,j3) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*
     &     zb(
     &    j4,j6)*izab2(j1,j4,j3,j2)**2*sprod(j1,j2) + za(j4,j3)*
     &     za(j3,j5)*
     &    zb(j4,j3)*zb(j4,j6)*izab2(j1,j4,j3,j2)**2*sprod(j5,j6) - 
     &     za(j4,j3
     &    )*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2*
     &     sprod(j4,j3
     &    ) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j3,j4,j2
     &    )**2*sprod(j1,j2) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*
     &    izab2(j1,j3,j4,j2)**2*sprod(j5,j6) + 2.D0*za(j4,j3)*za(j1,j5)
     &     *zb(
     &    j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j1)*t(j3,j4,j2
     &    ) )
      bcoeff(b34) = bcoeff(b34) + IDelta * (  - 2.D0*za(j4,j3)*za(j1,j5
     &    )*zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j2)**2 - 
     &    2.D0*za(j4,j3)*za(j2,j5)*zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,
     &    j2)**2*t(j3,j4,j1) + 2.D0*za(j4,j3)*za(j2,j5)*zb(j4,j2)*zb(j4
     &    ,j6)*izab2(j1,j3,j4,j2)**2*t(j3,j4,j2) + 1.D0/2.D0*za(j4,j3)*
     &    zb(j4,j6)**2*izb(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2)
     &     + 1.D0/2.D0*za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*zab2(j2,j3,j4,
     &    j1)*izab2(j1,j3,j4,j2) - za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*
     &    izab2(j1,j3,j4,j2)**2*sprod(j4,j3)*t(j3,j4,j2) + za(j4,j3)*
     &     zb(j4,
     &    j6)**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*sprod(j1,j2)*
     &     t(j3,j4,j2)
     &     + za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*
     &     sprod(
     &    j5,j6)*t(j3,j4,j2) - 2.D0*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j1
     &    ,j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j1) + 2.D0*za(j3,j1)*za(j3
     &    ,j5)*zb(j4,j3)*zb(j1,j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j2) + 
     &    2.D0*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j2,j6)*izab2(j1,j4,j3,
     &    j2)**3*t(j4,j3,j1)**2 )
      bcoeff(b34) = bcoeff(b34) + IDelta * (  - 2.D0*za(j3,j1)*za(j3,j5
     &    )*zb(j4,j3)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j1)*t(j4,
     &    j3,j2) + 1.D0/2.D0*za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*zab2(j2,
     &    j4,j3,j1)*izab2(j1,j4,j3,j2) + 1.D0/2.D0*za(j3,j5)**2*zb(j4,
     &    j3)*iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) - za(j3,
     &    j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,j3,j2)**2*sprod(j4,j3)
     &     *t(
     &    j4,j3,j1) + za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,j3,
     &    j2)**2*sprod(j1,j2)*t(j4,j3,j1) + za(j3,j5)**2*zb(j4,j3)*
     &     iza(j5,
     &    j6)*izab2(j1,j4,j3,j2)**2*sprod(j5,j6)*t(j4,j3,j1) - 
     &     za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2) - za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) )

      bcoeff(b12zm)= + IDelta**2 * ( 3.D0*zab2(j3,j1,j2,j4)*zab2(j2,j4,
     &    j3,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)*sprod(j4,j3) + 
     &     3.D0*
     &    zab2(j3,j1,j2,j4)*zab2(j2,j4,j3,j1)*zab2(j5,j1,j2,j6)*izab2(
     &    j1,j4,j3,j2)*sprod(j1,j2) - 3.D0*zab2(j3,j1,j2,j4)*
     &     zab2(j2,j4,j3,
     &    j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)*sprod(j5,j6) - 3.D0*
     &    zab2(j3,j1,j2,j4)*zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(
     &    j1,j6,j5,j2)*sprod(j4,j3) + 3.D0*zab2(j3,j1,j2,j4)*
     &     zab2(j2,j6,j5,
     &    j1)*zab2(j5,j1,j2,j6)*izab2(j1,j6,j5,j2)*sprod(j1,j2) + 3.D0*
     &    zab2(j3,j1,j2,j4)*zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(
     &    j1,j6,j5,j2)*sprod(j5,j6) + 3.D0*zab2(j3,j2,j1,j4)*
     &     zab2(j2,j3,j4,
     &    j1)*zab2(j5,j2,j1,j6)*izab2(j1,j3,j4,j2)*sprod(j4,j3) + 3.D0*
     &    zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,j6)*izab2(
     &    j1,j3,j4,j2)*sprod(j1,j2) - 3.D0*zab2(j3,j2,j1,j4)*
     &     zab2(j2,j3,j4,
     &    j1)*zab2(j5,j2,j1,j6)*izab2(j1,j3,j4,j2)*sprod(j5,j6) - 3.D0*
     &    zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,j6)*izab2(
     &    j1,j5,j6,j2)*sprod(j4,j3) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta**2 * ( 3.D0*zab2(j3,j2,j1,
     &    j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,j6)*izab2(j1,j5,j6,j2)*
     &     sprod(
     &    j1,j2) + 3.D0*zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,
     &    j1,j6)*izab2(j1,j5,j6,j2)*sprod(j5,j6) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * (  - 2.D0*za(j4,j3)*za(
     &    j3,j1)*za(j5,j6)*zb(j4,j3)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3
     &    ,j2)**3*t(j4,j3,j1) + 2.D0*za(j4,j3)*za(j3,j1)*za(j5,j6)*zb(
     &    j4,j3)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j2)
     &     - 2.D0*za(j4,j3)*za(j3,j5)*za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(
     &    j5,j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j1) + 2.D0*za(j4,j3)*za(
     &    j3,j5)*za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j3,j4
     &    ,j2)**3*t(j3,j4,j2) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)
     &    *izab2(j1,j4,j3,j2)**2*sprod(j4,j3) - za(j4,j3)*za(j3,j5)*
     &     zb(j4,
     &    j3)*zb(j4,j6)*izab2(j1,j4,j3,j2)**2*sprod(j1,j2) - za(j4,j3)*
     &     za(
     &    j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j4,j3,j2)**2*
     &     sprod(j5,j6) + 
     &    za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2
     &    *sprod(j4,j3) - za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*
     &     izab2(j1,
     &    j3,j4,j2)**2*sprod(j1,j2) - za(j4,j3)*za(j3,j5)*zb(j4,j3)*
     &     zb(j4,
     &    j6)*izab2(j1,j3,j4,j2)**2*sprod(j5,j6) - 2.D0*za(j4,j3)*
     &     za(j1,j5)
     &    *za(j5,j6)*zb(j4,j2)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**
     &    3*t(j6,j5,j1) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( 2.D0*za(j4,j3)*za(j1,
     &    j5)*za(j5,j6)*zb(j4,j2)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2
     &    )**3*t(j6,j5,j2) - 2.D0*za(j4,j3)*za(j1,j5)*zb(j4,j2)*zb(j4,
     &    j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j1)*t(j3,j4,j2) + 2.D0*za(
     &    j4,j3)*za(j1,j5)*zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**3*t(
     &    j3,j4,j2)**2 + 2.D0*za(j4,j3)*za(j2,j5)*zb(j4,j2)*zb(j4,j6)*
     &    izab2(j1,j3,j4,j2)**2*t(j3,j4,j1) - 2.D0*za(j4,j3)*za(j2,j5)*
     &    zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2*t(j3,j4,j2) - 1.D0/
     &    2.D0*za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*zab2(j2,j4,j3,j1)*
     &    izab2(j1,j4,j3,j2) - 1.D0/2.D0*za(j4,j3)*zb(j4,j6)**2*izb(j5,
     &    j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) + za(j4,j3)*zb(j4,j6
     &    )**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*sprod(j4,j3)*
     &     t(j3,j4,j2) - 
     &    za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*
     &     sprod(j1,
     &    j2)*t(j3,j4,j2) - za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*izab2(j1,
     &    j3,j4,j2)**2*sprod(j5,j6)*t(j3,j4,j2) - 2.D0*za(j3,j1)*
     &     za(j3,j5)*
     &    za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**3
     &    *t(j5,j6,j1) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( 2.D0*za(j3,j1)*za(j3,
     &    j5)*za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2
     &    )**3*t(j5,j6,j2) + 2.D0*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j1,
     &    j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j1) - 2.D0*za(j3,j1)*za(j3,
     &    j5)*zb(j4,j3)*zb(j1,j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j2) - 2.
     &    D0*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j2,j6)*izab2(j1,j4,j3,j2)
     &    **3*t(j4,j3,j1)**2 + 2.D0*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j2
     &    ,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j1)*t(j4,j3,j2) - 2.D0*za(
     &    j3,j1)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**3*t(
     &    j5,j6,j1)*t(j5,j6,j2) + 2.D0*za(j3,j1)*za(j5,j6)*zb(j4,j6)*
     &    zb(j2,j6)*izab2(j1,j5,j6,j2)**3*t(j5,j6,j2)**2 + 2.D0*za(j3,
     &    j2)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**2*t(j5,
     &    j6,j1) - 2.D0*za(j3,j2)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(
     &    j1,j5,j6,j2)**2*t(j5,j6,j2) - 1.D0/2.D0*za(j3,j5)**2*zb(j4,j3
     &    )*iza(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2) - 1.D0/2.D0
     &    *za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1
     &    ,j3,j4,j2) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( za(j3,j5)**2*zb(j4,j3)
     &    *iza(j5,j6)*izab2(j1,j4,j3,j2)**2*sprod(j4,j3)*t(j4,j3,j1) 
     &     - za(
     &    j3,j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,j3,j2)**2*
     &     sprod(j1,j2)
     &    *t(j4,j3,j1) - za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,
     &    j3,j2)**2*sprod(j5,j6)*t(j4,j3,j1) - 1.D0/2.D0*za(j3,j5)**2*
     &     zb(j5
     &    ,j6)*iza(j4,j3)*zab2(j2,j5,j6,j1)*izab2(j1,j5,j6,j2) - 1.D0/2.
     &    D0*za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*zab2(j2,j6,j5,j1)*izab2(
     &    j1,j6,j5,j2) - za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,
     &    j5,j2)**2*sprod(j4,j3)*t(j6,j5,j1) - za(j3,j5)**2*zb(j5,j6)*
     &     iza(
     &    j4,j3)*izab2(j1,j6,j5,j2)**2*sprod(j1,j2)*t(j6,j5,j1) + 
     &     za(j3,j5)
     &    **2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,j5,j2)**2*sprod(j5,j6)*
     &     t(j6,
     &    j5,j1) + 2.D0*za(j3,j5)*za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izab2(
     &    j1,j6,j5,j2)**2*t(j6,j5,j1) - 2.D0*za(j3,j5)*za(j1,j5)*zb(j4,
     &    j1)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*t(j6,j5,j2) - 2.D0*za(j3,
     &    j5)*za(j1,j5)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,
     &    j5,j1)**2 )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( 2.D0*za(j3,j5)*za(j1,
     &    j5)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,j5,j1)*t(
     &    j6,j5,j2) - za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,
     &    j5,j6,j2)**2*sprod(j4,j3) - za(j3,j5)*za(j5,j6)*zb(j4,j6)*
     &     zb(j5,
     &    j6)*izab2(j1,j5,j6,j2)**2*sprod(j1,j2) + za(j3,j5)*za(j5,j6)*
     &     zb(
     &    j4,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**2*sprod(j5,j6) - 
     &     za(j3,j5)*
     &    za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*
     &     sprod(j4,j3)
     &     - za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)
     &    **2*sprod(j1,j2) + za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*
     &     izab2(
     &    j1,j6,j5,j2)**2*sprod(j5,j6) + za(j3,j5)*zb(j4,j6)*
     &     zab2(j2,j4,j3,
     &    j1)*izab2(j1,j4,j3,j2) + za(j3,j5)*zb(j4,j6)*zab2(j2,j3,j4,j1
     &    )*izab2(j1,j3,j4,j2) + za(j3,j5)*zb(j4,j6)*zab2(j2,j5,j6,j1)*
     &    izab2(j1,j5,j6,j2) + za(j3,j5)*zb(j4,j6)*zab2(j2,j6,j5,j1)*
     &    izab2(j1,j6,j5,j2) - 1.D0/2.D0*za(j5,j6)*zb(j4,j6)**2*izb(j4,
     &    j3)*zab2(j2,j5,j6,j1)*izab2(j1,j5,j6,j2) - 1.D0/2.D0*za(j5,j6
     &    )*zb(j4,j6)**2*izb(j4,j3)*zab2(j2,j6,j5,j1)*izab2(j1,j6,j5,j2
     &    ) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * (  - za(j5,j6)*zb(j4,j6)
     &    **2*izb(j4,j3)*izab2(j1,j5,j6,j2)**2*sprod(j4,j3)*t(j5,j6,j2)
     &     - 
     &    za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*izab2(j1,j5,j6,j2)**2*
     &     sprod(j1,
     &    j2)*t(j5,j6,j2) + za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*izab2(j1,
     &    j5,j6,j2)**2*sprod(j5,j6)*t(j5,j6,j2) )
      bcoeff(b12zm) = bcoeff(b12zm) - 2.D0*za(j3,j1)*za(j3,j2)*zb(j1,j6
     & )*zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 + 2.D0*
     &    za(j3,j1)*zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*zab2(j3,j4,j2,j6)*
     &    izab2(j1,j4,j3,j2)**3*t(j4,j3,j1) - 2.D0*za(j1,j5)*za(j2,j5)*
     &    zb(j4,j1)*zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*izab2(j1,j3,j4,j2)
     &    **2 + 2.D0*za(j1,j5)*zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*zab2(j5,
     &    j3,j1,j4)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j2) - iza(j4,j3)*izb(
     &    j5,j6)*zab2(j3,j4,j2,j6)**2*izab2(j1,j4,j3,j2)**2 - iza(j5,j6
     &    )*izb(j4,j3)*zab2(j5,j3,j1,j4)**2*izab2(j1,j3,j4,j2)**2

      bcoeff(b56)= + t234ms56**(-2) * ( za(j1,j5)**2*zb(j4,j1)**2*zb(j5
     &    ,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2*sprod(j5,j6) )
      bcoeff(b56) = bcoeff(b56) + t234ms56**(-1) * ( 2.D0*za(j1,j5)**2*
     &    zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2 - 2.D0
     &    *za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izb(j4,j3)*zab2(j1,j6,j5,j4)*
     &    zab2(j5,j6,j1,j2)*izab2(j1,j6,j5,j2)**3 )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-2) * ( za(j3,j2)**2*za(j5
     &    ,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2*
     &     sprod(j5,j6) )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-1) * ( 2.D0*za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2 - 2.D0
     &    *za(j3,j2)*za(j5,j6)*zb(j2,j6)*iza(j4,j3)*zab2(j3,j5,j6,j2)*
     &    zab2(j1,j5,j2,j6)*izab2(j1,j5,j6,j2)**3 )
      bcoeff(b56) = bcoeff(b56) + IDelta**2 * ( 3.D0*zab2(j3,j1,j2,j4)*
     &    zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j6,j5,j2)*
     &     sprod(j4,
     &    j3) - 3.D0*zab2(j3,j1,j2,j4)*zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,
     &    j6)*izab2(j1,j6,j5,j2)*sprod(j1,j2) - 3.D0*zab2(j3,j1,j2,j4)*
     &    zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j6,j5,j2)*
     &     sprod(j5,
     &    j6) + 3.D0*zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j5,j6,j2)*sprod(j4,j3) - 3.D0*zab2(j3,j2,j1,j4)*
     &    zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,j6)*izab2(j1,j5,j6,j2)*
     &     sprod(j1,
     &    j2) - 3.D0*zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j5,j6,j2)*sprod(j5,j6) )
      bcoeff(b56) = bcoeff(b56) + IDelta * ( 2.D0*za(j4,j3)*za(j1,j5)*
     &    za(j5,j6)*zb(j4,j2)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3
     &    *t(j6,j5,j1) - 2.D0*za(j4,j3)*za(j1,j5)*za(j5,j6)*zb(j4,j2)*
     &    zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,j5,j2) + 2.D0*
     &    za(j3,j1)*za(j3,j5)*za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*
     &    izab2(j1,j5,j6,j2)**3*t(j5,j6,j1) - 2.D0*za(j3,j1)*za(j3,j5)*
     &    za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**3
     &    *t(j5,j6,j2) + 2.D0*za(j3,j1)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*
     &    izab2(j1,j5,j6,j2)**3*t(j5,j6,j1)*t(j5,j6,j2) - 2.D0*za(j3,j1
     &    )*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**3*t(j5,j6
     &    ,j2)**2 - 2.D0*za(j3,j2)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(
     &    j1,j5,j6,j2)**2*t(j5,j6,j1) + 2.D0*za(j3,j2)*za(j5,j6)*zb(j4,
     &    j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**2*t(j5,j6,j2) + 1.D0/2.D0*
     &    za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*zab2(j2,j5,j6,j1)*izab2(j1,
     &    j5,j6,j2) + 1.D0/2.D0*za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*zab2(
     &    j2,j6,j5,j1)*izab2(j1,j6,j5,j2) )
      bcoeff(b56) = bcoeff(b56) + IDelta * ( za(j3,j5)**2*zb(j5,j6)*
     &    iza(j4,j3)*izab2(j1,j6,j5,j2)**2*sprod(j4,j3)*t(j6,j5,j1) + 
     &     za(j3
     &    ,j5)**2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,j5,j2)**2*
     &     sprod(j1,j2)*
     &    t(j6,j5,j1) - za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,
     &    j5,j2)**2*sprod(j5,j6)*t(j6,j5,j1) - 2.D0*za(j3,j5)*za(j1,j5)*
     &     zb(
     &    j4,j1)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*t(j6,j5,j1) + 2.D0*za(
     &    j3,j5)*za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*t(
     &    j6,j5,j2) + 2.D0*za(j3,j5)*za(j1,j5)*zb(j4,j2)*zb(j5,j6)*
     &    izab2(j1,j6,j5,j2)**3*t(j6,j5,j1)**2 - 2.D0*za(j3,j5)*za(j1,
     &    j5)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,j5,j1)*t(
     &    j6,j5,j2) + za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,
     &    j5,j6,j2)**2*sprod(j4,j3) + za(j3,j5)*za(j5,j6)*zb(j4,j6)*
     &     zb(j5,
     &    j6)*izab2(j1,j5,j6,j2)**2*sprod(j1,j2) - za(j3,j5)*za(j5,j6)*
     &     zb(
     &    j4,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**2*sprod(j5,j6) + 
     &     za(j3,j5)*
     &    za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*
     &     sprod(j4,j3)
     &     + za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)
     &    **2*sprod(j1,j2) )
      bcoeff(b56) = bcoeff(b56) + IDelta * (  - za(j3,j5)*za(j5,j6)*zb(
     &    j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*sprod(j5,j6) - 
     &     za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j5,j6,j1)*izab2(j1,j5,j6,j2) - za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j6,j5,j1)*izab2(j1,j6,j5,j2) + 1.D0/2.D0*
     &    za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*zab2(j2,j5,j6,j1)*izab2(j1,
     &    j5,j6,j2) + 1.D0/2.D0*za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*zab2(
     &    j2,j6,j5,j1)*izab2(j1,j6,j5,j2) + za(j5,j6)*zb(j4,j6)**2*izb(
     &    j4,j3)*izab2(j1,j5,j6,j2)**2*sprod(j4,j3)*t(j5,j6,j2) + 
     &     za(j5,j6)
     &    *zb(j4,j6)**2*izb(j4,j3)*izab2(j1,j5,j6,j2)**2*sprod(j1,j2)*
     &     t(j5,
     &    j6,j2) - za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*izab2(j1,j5,j6,j2)
     &    **2*sprod(j5,j6)*t(j5,j6,j2) )

      bcoeff(b134)= + t134ms34**(-2) * (  - za(j3,j1)**2*zb(j4,j3)*zb(
     &    j1,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2*sprod(j4,j3) )
      bcoeff(b134) = bcoeff(b134) + t134ms34**(-1) * (  - 2.D0*za(j3,j1
     &    )**2*zb(j4,j3)*zb(j1,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2
     &     + 2.D0*za(j3,j1)*zb(j4,j3)*zb(j1,j6)*izb(j5,j6)*zab2(j3,j4,
     &    j1,j2)*zab2(j1,j4,j3,j6)*izab2(j1,j4,j3,j2)**3 )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-2) * (  - za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2*
     &     sprod(j5,
     &    j6) )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-1) * (  - 2.D0*za(j3,j2
     &    )**2*za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2
     &     + 2.D0*za(j3,j2)*za(j5,j6)*zb(j2,j6)*iza(j4,j3)*zab2(j3,j5,
     &    j6,j2)*zab2(j1,j5,j2,j6)*izab2(j1,j5,j6,j2)**3 )
      bcoeff(b134) = bcoeff(b134) + 2.D0*za(j3,j1)*za(j3,j2)*zb(j1,j6)*
     & zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 - 2.D0*za(
     &    j3,j1)*zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*zab2(j3,j4,j2,j6)*
     &    izab2(j1,j4,j3,j2)**3*t(j4,j3,j1) + iza(j4,j3)*izb(j5,j6)*
     &    zab2(j3,j4,j2,j6)**2*izab2(j1,j4,j3,j2)**2

      bcoeff(b234)= + t234ms34**(-2) * (  - za(j4,j3)*za(j2,j5)**2*zb(
     &    j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2*sprod(j4,j3) )
      bcoeff(b234) = bcoeff(b234) + t234ms34**(-1) * (  - 2.D0*za(j4,j3
     &    )*za(j2,j5)**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2
     &     + 2.D0*za(j4,j3)*za(j2,j5)*zb(j4,j2)*iza(j5,j6)*zab2(j1,j3,
     &    j2,j4)*zab2(j5,j3,j4,j2)*izab2(j1,j3,j4,j2)**3 )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-2) * (  - za(j1,j5)**2*
     &    zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2*
     &     sprod(j5,
     &    j6) )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-1) * (  - 2.D0*za(j1,j5
     &    )**2*zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2
     &     + 2.D0*za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izb(j4,j3)*zab2(j1,j6,
     &    j5,j4)*zab2(j5,j6,j1,j2)*izab2(j1,j6,j5,j2)**3 )
      bcoeff(b234) = bcoeff(b234) + 2.D0*za(j1,j5)*za(j2,j5)*zb(j4,j1)*
     & zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*izab2(j1,j3,j4,j2)**2 - 2.D0*za(
     &    j1,j5)*zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*zab2(j5,j3,j1,j4)*
     &    izab2(j1,j3,j4,j2)**3*t(j3,j4,j2) + iza(j5,j6)*izb(j4,j3)*
     &    zab2(j5,j3,j1,j4)**2*izab2(j1,j3,j4,j2)**2

      bcoeff(rat)= + t134ms34**(-1) * (  - za(j3,j1)**2*zb(j4,j3)*zb(j1
     &    ,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + t234ms34**(-1) * (  - za(j4,j3)*za(j2
     &    ,j5)**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + t234ms56**(-1) * (  - za(j1,j5)**2*
     &    zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + t134ms56**(-1) * (  - za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + IDelta * (  - za(j3,j5)**2*iza(j4,j3)
     &    *iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2)*sprod(j4,j3)
     &     + 
     &    za(j3,j5)**2*iza(j4,j3)*iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1
     &    ,j3,j4,j2)*sprod(j1,j2) - za(j3,j5)**2*iza(j4,j3)*iza(j5,j6)*
     &    zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2)*sprod(j5,j6) -
     &     2.D0*za(j3,j5
     &    )*zb(j4,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2) - 2.D0*za(j3
     &    ,j5)*zb(j4,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) - zb(j4,
     &    j6)**2*izb(j4,j3)*izb(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3
     &    ,j2)*sprod(j4,j3) + zb(j4,j6)**2*izb(j4,j3)*izb(j5,j6)*
     &     zab2(j2,j4
     &    ,j3,j1)*izab2(j1,j4,j3,j2)*sprod(j1,j2) - zb(j4,j6)**2*
     &     izb(j4,j3)
     &    *izb(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2)*
     &     sprod(j5,j6) )
      bcoeff(rat) = bcoeff(rat) - za(j3,j5)*zb(j4,j6)*izab2(j1,j4,j3,j2
     & )**2 - za(j3,j5)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2 + iza(j4,j3)*
     &    izb(j5,j6)*zab2(j3,j4,j2,j6)**2*izab2(j1,j4,j3,j2)**2 + iza(
     &    j5,j6)*izb(j4,j3)*zab2(j5,j3,j1,j4)**2*izab2(j1,j3,j4,j2)**2

c      if     ((case .eq. 'HWWint') .or. (case .eq. 'WWqqbr')
c     &   .or. (case .eq. 'HWWH+i') .or. (case .eq. 'ggWW4l')
c     &   .or. (case .eq. 'ggWWbx')
c     &   .or. (case .eq. 'ggVV4l') .or. (case .eq. 'ggVVbx')) then
      bcoeff(b12)=WWbub6(j2,j1,j4,j3,j6,j5,zb,za,sprod,mass)
      bcoeff(b12m)=bcoeff(b12zm)-bcoeff(b12)
c      elseif ((case .eq. 'HZZint') .or. (case .eq. 'ZZlept')
c     &   .or. (case .eq. 'HZZH+i') .or. (case .eq. 'ggZZ4l')
c     &   .or. (case .eq. 'ggZZbx')) then
c      bcoeff(b12)=bcoeff(b12zm)
c      bcoeff(b1)=-bcoeff(b12)-bcoeff(b34)-bcoeff(b56)
c     & -bcoeff(b134)-bcoeff(b234)
c      else
c      write(6,*) 'Unimplemented case in mbc.f'
c      stop
c      endif

      endif
      
      do j=1,8
      bcoeff(j)=ci*bcoeff(j)
      enddo

      return
      end




