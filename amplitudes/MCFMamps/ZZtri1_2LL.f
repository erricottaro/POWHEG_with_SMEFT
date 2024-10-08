      subroutine ZZtri1_2LL(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,
     &     Xpp,Xmp,Xpm,Xmm,Xrat)
      use mod_types; use mod_consts_dp
C-----Author: R.K. Ellis (September 2013)
C-----Trianglecoefficient for LL coupling 
C-----Triangle C0(p1,p2,mass,mass,mass)
C-----Xpp and Xmp refer to the initial state gluon polarizations
      implicit none
      include 'ggZZcomputemp.f'
      real(dp)  :: mass,masssq,t
      integer j1,j2,j3,j4,j5,j6,k1,k2,k3,k4,k5,k6,h3,h5
      complex(dp)  :: Xpp(2,2),Xmp(2,2),Xmm(2,2),Xpm(2,2),zab2,
     & Funcpp_2,Funcmp_0,Funcmp_2,amp0,amp2,app2,Xrat(2,2,2,2)
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)

c--- statement functions
      t(k1,k2,k3)=sprod(k1,k2)+sprod(k2,k3)+sprod(k1,k3)
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)

C----Functions for the mp amplitude
      Funcmp_0(k1,k2,k3,k4,k5,k6)=
     & 0.5d0*sprod(k1,k2)*((zab2(k2,k1,k3,k4)*zab2(k5,k3,k4,k1))**2
     & +(t(k1,k3,k4)*za(k2,k5)*zb(k4,k1))**2)
     & /(za(k5,k6)*zb(k4,k3)*zab2(k2,k3,k4,k1)**4)

      Funcmp_2(k1,k2,k3,k4,k5,k6)=0.5d0*(
     & (-2*(t(k1,k3,k4)+t(k2,k3,k4))*za(k2,k3)*za(k2,k5)
     & *zab2(k1,k3,k4,k2)*zb(k4,k1)*zb(k6,k1)/zab2(k2,k3,k4,k1)**3)
     & -(za(k2,k3)*zab2(k1,k3,k4,k2)
     & *(za(k1,k5)*zb(k4,k1)-za(k2,k5)*zb(k4,k2))*zb(k6,k1))
     & /zab2(k2,k3,k4,k1)**2
     & +((t(k1,k3,k4)+t(k2,k3,k4))
     & *(za(k1,k5)*za(k2,k3)*zb(k4,k2)*zb(k6,k1)
     & +za(k1,k3)*za(k2,k5)*zb(k4,k1)*zb(k6,k2)))/zab2(k2,k3,k4,k1)**2
     & +(za(k1,k3)*za(k1,k5)*zb(k2,k1)*zb(k6,k4))/zab2(k2,k3,k4,k1)
     & -(za(k1,k5)*zb(k4,k2)*(za(k2,k3)*zb(k6,k2)-za(k3,k4)*zb(k6,k4)))
     & /zab2(k2,k3,k4,k1)
     & -(zb(k4,k2)
     & *(za(k1,k3)*za(k2,k5)*zb(k6,k2)-za(k1,k5)*za(k3,k4)*zb(k6,k4)))
     & /zab2(k2,k3,k4,k1)
     & -(za(k2,k5)*zab2(k1,k3,k4,k2)*zb(k4,k1)
     & *(za(k3,k4)*zb(k6,k4)-za(k3,k5)*zb(k6,k5)))
     & /zab2(k2,k3,k4,k1)**2)/(sprod(k3,k4)*sprod(k5,k6))

C----Functions for the pp amplitude
      Funcpp_2(k1,k2,k3,k4,k5,k6)=
     & -zb(k1,k2)/(za(k1,k2)*sprod(k3,k4)*sprod(k5,k6))
     &  *(za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)/zab2(k1,k3,k4,k2)**2
     &  *(sprod(k1,k3)+sprod(k1,k4)+sprod(k2,k3)+sprod(k2,k4)+
     & 2d0*sprod(k3,k4))
     &  +(za(k1,k3)*za(k5,k6)*zb(k2,k6)+za(k1,k5)*za(k3,k4)*zb(k2,k4))
     & *zb(k4,k6)/zab2(k1,k3,k4,k2))

c--old fully symmetrized version
C      Funcpp_2(k1,k2,k3,k4,k5,k6)=
c     & -zb(k1,k2)/(za(k1,k2)*sprod(k3,k4)*sprod(k5,k6))
c     &  *(izab2(k1,k3,k4,k2)**2*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)
c     &      *(sprod(k1,k3)+sprod(k1,k4)+sprod(k2,k3)+sprod(k2,k4)+2d0*sprod(k3,k4))
c     &  +izab2(k1,k3,k4,k2)
c     &  *(za(k1,k3)*za(k5,k6)*zb(k2,k6)
c     &   +za(k1,k5)*za(k3,k4)*zb(k2,k4))*zb(k4,k6)
c     & +izab2(k2,k3,k4,k1)**2*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)
c     &    *(sprod(k1,k3)+sprod(k1,k4)+sprod(k2,k3)+sprod(k2,k4)+2d0*sprod(k3,k4))
c     &  +izab2(k2,k3,k4,k1)
c     & *(za(k2,k3)*za(k5,k6)*zb(k1,k6)
c     &  +za(k2,k5)*za(k3,k4)*zb(k1,k4))*zb(k4,k6))

c--- end statement functions

      masssq=mass**2

      k1=j1
      k2=j2

      do h3=1,2
         if (h3 .eq.1) then
             k3=j3
             k4=j4
         elseif (h3 .eq.2) then
             k3=j4
             k4=j3
         endif 
      do h5=1,2
         if (h5 .eq.1) then
             k5=j5
             k6=j6
         elseif (h5 .eq.2) then
             k5=j6
             k6=j5
         endif 

      app2=Funcpp_2(k1,k2,k3,k4,k5,k6)+Funcpp_2(k2,k1,k5,k6,k3,k4)

       Xpp(h3,h5)=masssq*app2
c--- contribution to rational part (coefficient of -mass^2)
      Xrat(2,2,h3,h5)=-app2     
      
      if (computemp) then
        amp0=Funcmp_0(k1,k2,k3,k4,k5,k6)+Funcmp_0(k1,k2,k5,k6,k3,k4)
        amp2=Funcmp_2(k1,k2,k3,k4,k5,k6)+Funcmp_2(k1,k2,k5,k6,k3,k4)
        Xmp(h3,h5)=amp0+masssq*amp2
c--- contribution to rational part (coefficient of -mass^2)
        Xrat(1,2,h3,h5)=-amp2     
      else
        Xmp(h3,h5)=czero
        Xrat(1,2,h3,h5)=czero
      endif

c      if (docheck) call ggZZcapture('1x2pp',h3,h5,j1,j2,j3,j4,j5,j6,
c     &                              czero,app2,czero)

      enddo
      enddo

c      write(6,*) 'Xpp(1,1)',Xpp(1,1)
c      write(6,*) 'Xmp(1,1)',Xmp(1,1)

c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5)=dconjg(Xmp(3-h3,3-h5))
      Xmm(h3,h5)=dconjg(Xpp(3-h3,3-h5))
      Xrat(2,1,h3,h5)=dconjg(Xrat(1,2,3-h3,3-h5))
      Xrat(1,1,h3,h5)=dconjg(Xrat(2,2,3-h3,3-h5))
      enddo
      enddo

c      if (docheck) then
c        call ggZZparsecheck('LL',-1,+1,Xmp,
c     &   '( s12,   0,   0;masssq,masssq,masssq)')
c        call ggZZparsecheck('LL',+1,+1,Xpp,
c     &   '( s12,   0,   0;masssq,masssq,masssq)')
c        call ggZZparsecheck('LL',+1,-1,Xpm,
c     &   '( s12,   0,   0;masssq,masssq,masssq)')
c        call ggZZparsecheck('LL',-1,-1,Xmm,
c     &   '( s12,   0,   0;masssq,masssq,masssq)')
c      endif

      return
      end

