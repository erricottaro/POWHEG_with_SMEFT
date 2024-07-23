      subroutine WWmassivetri6(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,
     &     d,triang)
      use mod_types; use mod_consts_dp; use common_def
      implicit none
      complex(dp) ::  c(2,2,12),d(2,2,6),Cint(12,-2:0),triang(2,2,-2:0),
     & qlI3,tmp
      real(dp) ::  s12,s34,s56,s134,s156,masssq,Delta,shift,
     & cred13,cred23,mass
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
      integer j,k1,k2,k3,k4,k5,k6,e,h1,h2


      masssq=mass**2

      s134=sprod(k1,k3)+sprod(k1,k4)+sprod(k3,k4)
      s156=sprod(k2,k3)+sprod(k2,k4)+sprod(k3,k4)
      s12=sprod(k1,k2)
      s34=sprod(k3,k4)
      s56=sprod(k5,k6)

c--- Initialize all triangle coefficients to zero
      do h1=1,2
      do h2=1,2
      do j=1,12
      c(h1,h2,j)=czero
      enddo     
      enddo     
      enddo     

c--- Triangle 2  ! NO LONGER NEEDED
c      call triangle2(1,2,3,4,5,6,za,zb,c(2,2,2),c(2,1,2),
c     &                                 c(1,2,2),c(1,1,2))
     
c--- Triangle 6
      call WWtriangle6(1,2,3,4,5,6,za,zb,sprod,mass,c(2,2,6),c(2,1,6))
      call WWtriangle6(2,1,3,4,5,6,za,zb,sprod,mass,c(1,1,6),c(1,2,6))

c--- Triangles 7 and 9
c      call triangle7(1,2,3,4,5,6,za,zb,c(2,2,7),c(2,1,7))
c      call triangle7(1,2,6,5,4,3,zb,za,c(1,1,9),c(1,2,9))
      call WWtriangle7new(1,2,3,4,5,6,za,zb,sprod,mass,
     &     c(2,2,7),c(2,1,7))
      call WWtriangle7new(1,2,6,5,4,3,zb,za,sprod,mass,
     &     c(1,1,9),c(1,2,9))
      
c      call triangle9(1,2,3,4,5,6,za,zb,c(2,2,9),c(2,1,9))
c      call triangle9(1,2,6,5,4,3,zb,za,c(1,1,7),c(1,2,7))
      call WWtriangle9new(1,2,3,4,5,6,za,zb,sprod,mass,
     &     c(2,2,9),c(2,1,9))
      call WWtriangle9new(1,2,6,5,4,3,zb,za,sprod,mass,
     &     c(1,1,7),c(1,2,7))

c--- Triangles 8 and 10
c      call triangle9(2,1,3,4,5,6,za,zb,c(2,2,8),c(1,2,8))
c      call triangle9(2,1,6,5,4,3,zb,za,c(1,1,10),c(2,1,10))
      call WWtriangle9new(2,1,3,4,5,6,za,zb,sprod,mass,
     &     c(2,2,8),c(1,2,8))
      call WWtriangle9new(2,1,6,5,4,3,zb,za,sprod,mass,
     &     c(1,1,10),c(2,1,10))
      
c      call triangle7(2,1,3,4,5,6,za,zb,c(2,2,10),c(1,2,10))
c      call triangle7(2,1,6,5,4,3,zb,za,c(1,1,8),c(2,1,8))
      call WWtriangle7new(2,1,3,4,5,6,za,zb,sprod,mass,
     &     c(2,2,10),c(1,2,10))
      call WWtriangle7new(2,1,6,5,4,3,zb,za,sprod,mass,
     &     c(1,1,8),c(2,1,8))

c--- Triangle 11
c      call triangle11(1,2,3,4,5,6,za,zb,c(2,2,11),c(2,1,11))
c--- This works for mp
c      call triangle11(2,1,3,4,5,6,za,zb,tmp,c(1,2,11))
c--- This works for mm
c      call triangle11(2,1,4,3,6,5,zb,za,c(1,1,11),tmp)
      call WWtriangle11new(1,2,3,4,5,6,za,zb,sprod,mass,
     &     c(2,2,11),c(2,1,11))
c--- This works for mp
      call WWtriangle11new(2,1,3,4,5,6,za,zb,sprod,mass,tmp,c(1,2,11))
c--- This works for mm
      call WWtriangle11new(2,1,4,3,6,5,zb,za,sprod,mass,c(1,1,11),tmp)

c--- Triangle 12
      call WWtriangle12(1,2,3,4,5,6,za,zb,sprod,mass,
     &     c(2,2,12),c(2,1,12))
      call WWtriangle12(1,2,6,5,4,3,zb,za,sprod,mass,
     &     c(1,1,12),c(1,2,12))
c      call triangle12new(1,2,3,4,5,6,za,zb,c(2,2,12),c(2,1,12))
c      call triangle12new(1,2,6,5,4,3,zb,za,c(1,1,12),c(1,2,12))


      

c--- Perform shifting of triangle coefficients 
      Delta=(masssq-s134)*(masssq-s156)-(masssq-s34)*(masssq-s56)
      shift=(2d0*s12*masssq-Delta)/Delta**2
c      cred31=2d0*(s134-s34)*shift
c      cred34=2d0*(s156-s56)*shift
c      cred42=2d0*(s156-s34)*shift
c      cred43=2d0*(s134-s56)*shift
       cred13=
     &  2d0*((s134+masssq)*(s56+s34-s12)-2d0*s134*masssq-2d0*s34*s56)
     & /(s12*(s134-masssq)**2)
       cred23=
     &  2d0*((s156+masssq)*(s56+s34-s12)-2d0*s156*masssq-2d0*s34*s56)
     & /(s12*(s156-masssq)**2)


      do h1=1,2
      do h2=1,2
      c(h1,h2,6)=c(h1,h2,6)
     & -0.5d0*cred13*d(h1,h2,1)-0.5d0*cred23*d(h1,h2,2)
      c(h1,h2,7) = c(h1,h2,7)-(s134-s34)*shift*d(h1,h2,3)
      c(h1,h2,8) = c(h1,h2,8)-(s134-s56)*shift*d(h1,h2,4)
      c(h1,h2,9) = c(h1,h2,9)-(s156-s56)*shift*d(h1,h2,3)
      c(h1,h2,10)=c(h1,h2,10)-(s156-s34)*shift*d(h1,h2,4)
      enddo
      enddo
      
c      do e=-2,0
c      Cint(1,e)=qlI3(0d0,0d0,s12,0d0,0d0,0d0,musq,e)
c      enddo
c      do e=-2,0
c      Cint(2,e)=qlI3(0d0,s134,s34,0d0,0d0,masssq,musq,e)
c      enddo
c      do e=-2,0
c      Cint(3,e)=qlI3(0d0,s56,s134,0d0,0d0,masssq,musq,e)
c      enddo
c      do e=-2,0
c      Cint(4,e)=qlI3(0d0,s156,s56,0d0,0d0,masssq,musq,e)
c      enddo
c      do e=-2,0
c      Cint(5,e)=qlI3(0d0,s34,s156,0d0,0d0,masssq,musq,e)
c      enddo
      do e=-2,0
      Cint(6,e)=qlI3(s12,s56,s34,0d0,0d0,masssq,musq,e)
      enddo
      do e=-2,0
      Cint(7,e)=qlI3(s134,0d0,s34,0d0,masssq,masssq,musq,e)
      enddo
      do e=-2,0
      Cint(8,e)=qlI3(s56,0d0,s134,0d0,masssq,masssq,musq,e)
      enddo
      do e=-2,0
      Cint(9,e)=qlI3(s156,0d0,s56,0d0,masssq,masssq,musq,e)
      enddo
      do e=-2,0
      Cint(10,e)=qlI3(s34,0d0,s156,0d0,masssq,masssq,musq,e)
      enddo
      do e=-2,0
         Cint(11,e)=qlI3(s56,s12,s34,0d0,masssq,masssq,musq,e)
      enddo
      do e=-2,0
      Cint(12,e)=qlI3(s12,0d0,0d0,masssq,masssq,masssq,musq,e)
      enddo

      do h1=1,2
      do h2=1,2
      do e=-2,0
      triang(h1,h2,e)=czero
      do j=6,12
      triang(h1,h2,e)=triang(h1,h2,e)+c(h1,h2,j)*Cint(j,e)
      enddo
      enddo
      enddo
      enddo

      return
      end

