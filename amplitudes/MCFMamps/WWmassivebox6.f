      subroutine WWmassivebox6(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,d,box)
      use mod_types; use mod_consts_dp; use common_def
      implicit none
      complex(dp):: d(2,2,6),box(2,2,-2:0),
     & qlI4,Dint(6),WWD1six,WWD2six,WWD3six,WWD4six
      real(dp):: s12,s34,s56,s134,s156,masssq,mass
      integer j,k1,k2,k3,k4,k5,k6,h1,h2
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)


      masssq=mass**2

      s134=sprod(k1,k3)+sprod(k1,k4)+sprod(k3,k4)
      s156=sprod(k1,k5)+sprod(k1,k6)+sprod(k5,k6)
      s12=sprod(k1,k2)
      s34=sprod(k3,k4)
      s56=sprod(k5,k6)

c--- Initialize all box coefficients to zero
      do h1=1,2
      do h2=1,2
      do j=1,6
      d(h1,h2,j)=czero
      enddo
      enddo
      enddo

c--- Boxes 1 and 2      
      call WWbox1(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,d(2,2,1),d(2,1,1),
     &                                  d(1,2,1),d(1,1,1))
      call WWbox1(k1,k2,k6,k5,k4,k3,zb,za,sprod,mass,d(1,1,2),d(1,2,2),
     &                                  d(2,1,2),d(2,2,2))

c--- Boxes 5 and 6      
      call WWbox5(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,d(2,2,5),d(2,1,5),
     &                                  d(1,2,5),d(1,1,5))
      call WWbox5(k1,k2,k6,k5,k4,k3,zb,za,sprod,mass,d(1,1,6),d(1,2,6),
     &                                  d(2,1,6),d(2,2,6))

c--- Box 3 and 4  
      call WWbox3(k1,k2,k3,k4,k5,k6,za,zb,sprod,mass,d(2,2,3),d(2,1,3),
     &                                  d(1,2,3),d(1,1,3))

      call WWbox3(k2,k1,k6,k5,k4,k3,zb,za,sprod,mass,d(1,1,4),d(2,1,4),
     &                                  d(1,2,4),d(2,2,4))

c--- Box 4   
c      call box4(k1,k2,k3,k4,k5,k6,za,zb,d(2,2,4),d(2,1,4),
c     &                                  d(1,2,4),d(1,1,4))
c      write(6,*) 'box4 d(2,2,4)',d(2,2,4)
c      write(6,*) 'box4 d(2,1,4)',d(2,1,4)
c      write(6,*) 'box4 d(1,2,4)',d(1,2,4)
c      write(6,*) 'box4 d(1,1,4)',d(1,1,4)
c      pause
      

C----for boxes 1-4 we convert to sixdim box
      Dint(1)=WWD1six(k1,k2,k3,k4,k5,k6,sprod,mass)
      Dint(2)=WWD2six(k1,k2,k3,k4,k5,k6,sprod,mass)
      Dint(3)=WWD3six(k1,k2,k3,k4,k5,k6,sprod,mass)
      Dint(4)=WWD4six(k1,k2,k3,k4,k5,k6,sprod,mass)
      Dint(5)=qlI4(s34,0d0,0d0,s56,s134,s12,0d0,
     &     masssq,masssq,masssq,musq,0)
      Dint(6)=qlI4(s56,0d0,0d0,s34,s156,s12,0d0,
     &     masssq,masssq,masssq,musq,0)
      do h1=1,2
      do h2=1,2
      box(h1,h2,-2)=czero
      box(h1,h2,-1)=czero
      box(h1,h2,0)=czero
      do j=1,6
      box(h1,h2,0)=box(h1,h2,0)+d(h1,h2,j)*Dint(j)
      enddo
      enddo
      enddo

      return
      end
