      subroutine ZZC01x2LLmp(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,Xmp,Xpm)
      use mod_types; use mod_consts_dp
      implicit none
c--- Author: J. M. Campbell, October 2013
      integer k1,k2,k3,k4,k5,k6
      integer h3,h5,j1,j2,j3,j4,j5,j6,itot,irat
      real(dp)  :: mass,masssq,t,s134,s234
      complex(dp)  :: zab2,amp2,Xmp(2,2,2),Xpm(2,2,2)
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)
      parameter(itot=1,irat=2)
      
C---statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      t(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j1,j3)
C---end statement functions

      s134=t(j1,j3,j4)
      s234=t(j2,j3,j4)
      k1=j1
      k2=j2
      masssq=mass**2

      do h3=1,2
      do h5=1,2
      if (h3 .eq. 1) then
        k3=j3
        k4=j4
      elseif (h3 .eq. 2) then
        k3=j4
        k4=j3
      endif
      if (h5 .eq. 1) then
        k5=j5
        k6=j6
      elseif (h5 .eq. 2) then
        k5=j6
        k6=j5
      endif
      amp2=(
     & -2d0*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)**3
     &     *(s134+s234)*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)   
     & +1d0/zab2(k2,k3,k4,k1)**2*( 
     &    +(s134+s234)*(za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)
     &                 +za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4))
     &    +za(k2,k3)*zb(k1,k6)*zab2(k5,k1,k6,k4)*zab2(k1,k3,k4,k2)
     &    -za(k2,k5)*zb(k1,k4)*zab2(k3,k1,k4,k6)*zab2(k1,k3,k4,k2))  
     & -za(k3,k5)*zb(k4,k6)*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)
     &      )/sprod(k3,k4)/sprod(k5,k6)

      Xmp(h3,h5,itot)=amp2*masssq
      Xmp(h3,h5,irat)=amp2*(+0.5d0)

c      if (docheck) call ggZZcapture('1x2',h3,h5,j1,j2,j3,j4,j5,j6,
c     &                              czero,amp2,czero)

      enddo
      enddo

c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5,:)=dconjg(Xmp(3-h3,3-h5,:))
      enddo
      enddo

      return
      end
