      subroutine ZZtri12_34LL(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,
     & Xpp,Xmp,Xpm,Xmm)
      use mod_types; use mod_consts_dp
C-----Author: R.K. Ellis (September 2013)
C-----Triangle coefficient for LL coupling 
C-----Triangle C0(p12,p34,0,0,0)
C-----Tri3masscoeff is taken from BDK (11.6) and (11.7)
C-----The full mass^2 dependence is calculated in ZZmassivetri
C-----Xpp and Xmp refer to the initial state gluon polarizations
      implicit none
      include 'ggZZcomputemp.f'
      real(dp)  :: mass
      integer j1,j2,j3,j4,j5,j6,k1,k2,k3,k4,k5,k6,h3,h5
      complex(dp)  :: Xpp(2,2),Xmp(2,2),Tri3masscoeff,
     & Xpm(2,2),Xmm(2,2)
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)

c     & ,Tri3masscoeff_new,Tri3masscoeffmasssq,

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

C----Helicity 1^+ 2^+
      Xpp(:,:)=czero


      if (computemp) then
C----Helicity 1^- 2^+
        Xmp(h3,h5)=Tri3masscoeff(k1,k2,k3,k4,k5,k6,za,zb,sprod) 
c        write(6,*) 'Explicit',h3,h5,
c     &   Tri3masscoeffmasssq(k1,k2,k3,k4,k5,k6,za,zb)\
      else
        Xmp(h3,h5)=czero
      endif
      enddo
      enddo
      
c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5)=dconjg(Xmp(3-h3,3-h5))
      Xmm(h3,h5)=dconjg(Xpp(3-h3,3-h5))
      enddo
      enddo

      
      return
      end
      
