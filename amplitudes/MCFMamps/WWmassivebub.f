      subroutine WWmassivebub(k1,k4,k2,k3,k5,k6,za,zb,sprod,mass,bub)
      use mod_types; use mod_consts_dp; use common_def
      implicit none
c----Bubble coefficients extracted from BDK 11.5, 11.8
      include 'blabels.f'
      character*9 MCFMst,MCFMst1
      integer j,k1,k2,k3,k4,k5,k6,h1,h2,e
      complex(dp) ::  b(2,2,7),bcoeff(8),bub(2,2,-2:0),Bint(7,-2:0),
     & qlI2,tmp
      real(dp) ::  masssq,s12,s34,s56,s134,s156,mass
      complex(dp)  :: za(12,12),zb(12,12)
      real(dp)  :: sprod(12,12)
      masssq=mass**2

c--- note: these look unnatural due to the
c---       (k1,k4,k2,k3,k5,k6) ordering in the call
      s134=sprod(k1,k2)+sprod(k1,k3)+sprod(k2,k3)
      s156=sprod(k1,k5)+sprod(k1,k6)+sprod(k5,k6)
      s12=sprod(k1,k4)
      s34=sprod(k2,k3)
      s56=sprod(k5,k6)

      do h1=1,2
      do h2=1,2

      if ((h1.eq.1) .and. (h2.eq.1)) MCFMst='q+qb-g-g-' 
      if ((h1.eq.1) .and. (h2.eq.2)) MCFMst='q+qb-g-g+' 
      if ((h1.eq.2) .and. (h2.eq.1)) MCFMst='q+qb-g+g-' 
      if ((h1.eq.2) .and. (h2.eq.2)) MCFMst='q+qb-g+g+' 
      if (MCFMst .eq. 'q+qb-g-g-') then
          MCFMst1='q+qb-g+g+'
          call WWmbc(MCFMst1,k4,k1,k3,k2,k6,k5,zb,za,sprod,mass,bcoeff)
c--- this call interchanges roles of 1,2 so change coeffs. accordingly
          tmp=bcoeff(b134)
        bcoeff(b134)=bcoeff(b234)
        bcoeff(b234)=tmp
c--- (-,-) amplitudes
c      write(6,*) '--'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),cdabs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),cdabs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),cdabs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),cdabs(bcoeff(b56))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),cdabs(bcoeff(b134))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),cdabs(bcoeff(b34))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),cdabs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),cdabs(bcoeff(rat))
c      write(6,*)
        
      elseif (MCFMst.eq.'q+qb-g-g+') then
          MCFMst1='q+qb-g+g-'
          call WWmbc(MCFMst1,k4,k1,k2,k3,k5,k6,za,zb,sprod,mass,bcoeff)
c--- this call interchanges roles of 1,2 so change coeffs. accordingly
          tmp=bcoeff(b134)
        bcoeff(b134)=bcoeff(b234)
        bcoeff(b234)=tmp
c--- (-,+) amplitudes
c      write(6,*) '-+'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),cdabs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),cdabs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),cdabs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),cdabs(bcoeff(b56))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),cdabs(bcoeff(b134))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),cdabs(bcoeff(b34))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),cdabs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),cdabs(bcoeff(rat))
c      write(6,*)
        
      elseif (MCFMst.eq.'q+qb-g+g+') then
          call WWmbc(MCFMst,k1,k4,k2,k3,k5,k6,za,zb,sprod,mass,bcoeff)
c--- (+,+) amplitudes
c      write(6,*) '++'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),cdabs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),cdabs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),cdabs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),cdabs(bcoeff(b56))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),cdabs(bcoeff(b134))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),cdabs(bcoeff(b34))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),cdabs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),cdabs(bcoeff(rat))
c      write(6,*)
        
      elseif (MCFMst.eq.'q+qb-g+g-') then
          call WWmbc(MCFMst,k1,k4,k2,k3,k5,k6,za,zb,sprod,mass,bcoeff)
c--- (+,-) amplitudes
c      write(6,*) '+-'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),cdabs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),cdabs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),cdabs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),cdabs(bcoeff(b56))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),cdabs(bcoeff(b34))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),cdabs(bcoeff(b134))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),cdabs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),cdabs(bcoeff(rat))
c      write(6,*)
        
      endif

      do j=1,7
      b(h1,h2,j)=bcoeff(j)          
      enddo
      
      enddo
      enddo


      do h1=1,2
      do h2=1,2
      bub(h1,h2,-2)=czero
      enddo
      enddo
      

      do e=-1,0
      Bint(1,e)=qlI2(s12,0d0,0d0,musq,e)
      enddo
      do e=-1,0
      Bint(2,e)=qlI2(s34,0d0,masssq,musq,e)
      enddo
      do e=-1,0
      Bint(3,e)=qlI2(s56,0d0,masssq,musq,e)
      enddo
      do e=-1,0
      Bint(4,e)=qlI2(s134,0d0,masssq,musq,e)
      enddo
      do e=-1,0
      Bint(5,e)=qlI2(s156,0d0,masssq,musq,e)
      enddo
      do e=-1,0
      Bint(6,e)=qlI2(s12,masssq,masssq,musq,e)
      enddo
C-----rational part
      Bint(7,-1)=czero
      Bint(7, 0)=cone

c--- Note: Bint(8) is the massless result, should not be included in the sum
      do h1=1,2
      do h2=1,2
      do e=-2,0
      bub(h1,h2,e)=czero
      do j=1,7
      bub(h1,h2,e)=bub(h1,h2,e)+b(h1,h2,j)*Bint(j,e)
      enddo
      enddo
      enddo
      enddo

            
      return
      end

