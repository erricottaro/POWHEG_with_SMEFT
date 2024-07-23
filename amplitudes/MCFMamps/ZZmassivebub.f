      subroutine ZZmassivebub(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,
     &     bub,totrat)
      use mod_types; use mod_consts_dp
      implicit none
c----Bubble coefficients extracted from BDK 11.5, 11.8
      include 'ggZZintegrals.f'
      include 'blabels.f'
c      include 'docheck.f'
c      logical ggZZunstable
      character*9 MCFMst,MCFMst1
      integer j,j1,j3,j4,j2,j5,j6,k1,k3,k4,k2,k5,k6,h1,h2,h3,h4
      complex(dp)  :: b(2,2,2,2,7),bcoeff(8),bub(2,2,2,2,-2:0),
     & totrat(2,2,2,2),tmp
      real(dp)  :: mass
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)


c      common/ggZZunstable/ggZZunstable


c--- QCDLoop already initialized from call to massivebox
c      call qlinit
c      if (first) then
c      first=.false. 
c      write(6,*) 'masssq',masssq
c      write(6,*) 'musq',musq
c      pause
c      endif

c      s134=sprod(j1,j3)+sprod(j1,j4)+sprod(j3,j4)
c      s156=sprod(j1,j5)+sprod(j1,j6)+sprod(j5,j6)
c      s12=sprod(j1,j2)
c      s34=sprod(j3,j4)
c      s56=sprod(j5,j6)



      k1=j1
      k2=j2

      do h1=1,2
c      do h2=1,2
c--- only compute h2=1: will obtain others by c.c.
      h2=2
      do h3=1,2
      do h4=1,2
      
      if (h3 .eq.1) then
          k3=j3
          k4=j4
      elseif (h3 .eq.2) then
          k3=j4
          k4=j3
      endif 
      if (h4 .eq.1) then
          k5=j5
          k6=j6
      elseif (h4 .eq.2) then
          k5=j6
          k6=j5
      endif 

      if ((h1.eq.1) .and. (h2.eq.1)) MCFMst='q+qb-g-g-' 
      if ((h1.eq.1) .and. (h2.eq.2)) MCFMst='q+qb-g-g+' 
      if ((h1.eq.2) .and. (h2.eq.1)) MCFMst='q+qb-g+g-' 
      if ((h1.eq.2) .and. (h2.eq.2)) MCFMst='q+qb-g+g+' 

      if (MCFMst .eq. 'q+qb-g-g-') then
          MCFMst1='q+qb-g+g+'
          call ZZmbc(MCFMst1,k2,k1,k4,k3,k6,k5,zb,za,sprod,bcoeff)
c--- this call interchanges roles of 1,2 so change coeffs. accordingly
          tmp=bcoeff(b134)
          bcoeff(b134)=bcoeff(b234)
          bcoeff(b234)=tmp
        
      elseif (MCFMst.eq.'q+qb-g-g+') then
          call ZZmbc(MCFMst,k1,k2,k3,k4,k5,k6,za,zb,sprod,bcoeff)
          
      elseif (MCFMst.eq.'q+qb-g+g+') then
          call ZZmbc(MCFMst,k1,k2,k3,k4,k5,k6,za,zb,sprod,bcoeff)
        
      elseif (MCFMst.eq.'q+qb-g+g-') then
          MCFMst1='q+qb-g+g-'
          call ZZmbc(MCFMst1,k1,k2,k4,k3,k6,k5,zb,za,sprod,bcoeff)

      endif

      do j=1,7
      b(h1,h2,h3,h4,j)=bcoeff(j)
      enddo
      
      enddo

      enddo
c      enddo
      enddo

c--- obtain h2=1 results by c.c.
      do h3=1,2
      do h4=1,2
      b(2,1,h3,h4,:)=dconjg(b(1,2,3-h3,3-h4,:))
      b(1,1,h3,h4,:)=dconjg(b(2,2,3-h3,3-h4,:))
      enddo
      enddo

c--- copy rational parts to another array to return separately
      totrat(:,:,:,:)=b(:,:,:,:,rat)
     
c      if (docheck) then
c      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b34),'( s34;masssq,masssq)')
c      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b34),'( s34;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b34),'( s34;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b34),'( s34;masssq,masssq)')
c      
c      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b134),'(s134;masssq,masssq)')
c      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b134),'(s134;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b134),'(s134;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b134),'(s134;masssq,masssq)')
c      
c      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b56),'( s56;masssq,masssq)')
c      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b56),'( s56;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b56),'( s56;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b56),'( s56;masssq,masssq)')
c      
c      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b234),'(s156;masssq,masssq)')
c      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b234),'(s156;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b234),'(s156;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b234),'(s156;masssq,masssq)')
c      
c      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b12),'( s12;masssq,masssq)')
c      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b12),'( s12;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b12),'( s12;masssq,masssq)')
c      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b12),'( s12;masssq,masssq)')
c      
c      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,rat),'xe,res/lo  1')
c      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,rat),'xe,res/lo  1')
c      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,rat),'xe,res/lo  1')
c      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,rat),'xe,res/lo  1')
c      endif
      
C--- Check that coefficients sum to zero to test numerical stability
c--- note: only do for (1,2,1,1) and (2,2,1,1) for simplicity
c--- and multiply by s12 to remove dimensions
c      testmm=(b(1,1,1,1,1)+b(1,1,1,1,2)+b(1,1,1,1,3)
c     &       +b(1,1,1,1,4)+b(1,1,1,1,5))*s12
c      testmp=(b(1,2,1,1,1)+b(1,2,1,1,2)+b(1,2,1,1,3)
c     &       +b(1,2,1,1,4)+b(1,2,1,1,5))*s12
c      if ((abs(testmm) .gt. 1d-4) .or. (abs(testmp) .gt. 1d-4)) then
c        write(6,*) 'testmm=',testmm
c        write(6,*) 'testmp=',testmp
c        ggZZunstable=.true.
c      else
c        ggZZunstable=.false.
c      endif

      bub(:,:,:,:,:)=czero

c--- Note: Bint(6) does not appear in this calculation
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do j=1,7
      if (j .ne. 6) then
      bub(h1,h2,h3,h4,0)=bub(h1,h2,h3,h4,0)+b(h1,h2,h3,h4,j)*intB0(j)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      
c      if (docheck) then
c        do j=1,7
c        do h3=1,2
c        do h4=1,2
c        call ggZZcapture('bubmp'//char(ichar('0')+j),h3,h4,1,2,3,4,5,6,
c     &   b(1,2,h3,h4,j),
c     &   czero,czero)
c        call ggZZcapture('bubpp'//char(ichar('0')+j),h3,h4,1,2,3,4,5,6,
c     &   b(2,2,h3,h4,j),
c     &   czero,czero)
c        enddo        
c        enddo        
c        enddo        
c      endif

c      if (docheck) then
c        do j=1,7
c        write(6,*) 'j,Bint(j)',j,Bint(j,0)*b(2,1,j)
c        enddo
c      endif
      
      return
      end

