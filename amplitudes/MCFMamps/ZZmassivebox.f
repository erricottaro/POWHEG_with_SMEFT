      subroutine ZZmassivebox(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,box
     &     ,drat)
      use mod_types; use mod_consts_dp
      implicit none
      include 'ggZZintegrals.f'
      include 'ZZdlabels.f'
      integer j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,j
      real(dp)  :: mass
      complex(dp)  :: Xpp(2,2),Xmp(2,2),Xpm(2,2),Xmm(2,2),
     & box(2,2,2,2,-2:0),d(2,2,2,2,3),
     & Xrat(2,2,2,2),drat(2,2,2,2,3)
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)
      

      call ZZbox2LL(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,
     &     Xpp,Xmp,Xpm,Xmm,Xrat)
      d(2,2,:,:,d2_1_34)=Xpp(:,:)
      d(1,2,:,:,d2_1_34)=Xmp(:,:)
      d(2,1,:,:,d2_1_34)=Xpm(:,:)
      d(1,1,:,:,d2_1_34)=Xmm(:,:)
      drat(:,:,:,:,d2_1_34)=Xrat(:,:,:,:)
c      write(6,*) 'Xmp(1,1)',Xmp(1,1)
c      write(6,*) 'Xmp(1,2)',Xmp(1,2)
c      write(6,*) 'Xmp(2,1)',Xmp(2,1)
c      write(6,*) 'Xmp(2,2)',Xmp(2,2)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)
      
      call ZZbox2LL(j1,j2,j5,j6,j3,j4,za,zb,sprod,mass,
     &     Xpp,Xmp,Xpm,Xmm,Xrat)
      d(2,2,:,:,d1_2_34)=Xpp(:,:)
      d(1,2,:,:,d1_2_34)=Xmp(:,:)
      d(2,1,:,:,d1_2_34)=Xpm(:,:)
      d(1,1,:,:,d1_2_34)=Xmm(:,:)
      drat(:,:,:,:,d1_2_34)=Xrat(:,:,:,:)
c      write(6,*) 'Xmp(1,1)',Xmp(1,1)
c      write(6,*) 'Xmp(1,2)',Xmp(1,2)
c      write(6,*) 'Xmp(2,1)',Xmp(2,1)
c      write(6,*) 'Xmp(2,2)',Xmp(2,2)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)

      call ZZbox1LL(j1,j2,j3,j4,j5,j6,za,zb,sprod,mass,
     &     Xpp,Xmp,Xpm,Xmm,Xrat)
      d(2,2,:,:,d1_34_2)=Xpp(:,:)
      d(1,2,:,:,d1_34_2)=Xmp(:,:)
      d(2,1,:,:,d1_34_2)=Xpm(:,:)
      d(1,1,:,:,d1_34_2)=Xmm(:,:)
      drat(:,:,:,:,d1_34_2)=Xrat(:,:,:,:)
c      write(6,*) 'Xpp(1,1)',Xpp(1,1)
c      write(6,*) 'Xpp(1,2)',Xpp(1,2)
c      write(6,*) 'Xpp(2,1)',Xpp(2,1)
c      write(6,*) 'Xpp(2,2)',Xpp(2,2)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)

c--- These integrals are now initialized in ZZintegraleval      
c      do e=-2,0
c      Dint(1,e)
c     & =qlI4(s34,0d0,0d0,s56,s134,s12,masssq,masssq,masssq,masssq,musq,e)
c      enddo
c      do e=-2,0
c      Dint(2,e)
c     & =qlI4(s56,0d0,0d0,s34,s156,s12,masssq,masssq,masssq,masssq,musq,e)
c      enddo
c      do e=-2,0
c      Dint(3,e)
c     & =qlI4(s56,0d0,s34,0d0,s156,s134,masssq,masssq,masssq,masssq,musq,e)
c      enddo

c      write(6,*) 'check',intD0(d2_1_34)/Dint(1,0)
c      write(6,*) 'check',intD0(d1_2_34)/Dint(2,0)
c      write(6,*) 'check',intD0(d1_34_2)/Dint(3,0)
c      pause

c--- DEBUG
c      write(6,*) 'd(1,1,1,1,d1_2_34)/s12/s234,intD0(d1_2_34)*s12*s234',
c     & d(1,1,1,1,d1_2_34)/s12/s156,intD0(d1_2_34)*s12*s156
c      write(6,*) 'd(1,1,1,1,d2_1_34)/s12/s134,intD0(d2_1_34)*s12*s134',
c     & d(1,1,1,1,d2_1_34)/s12/s134,intD0(d2_1_34)*s12*s134
c--- DEBUG

      box(:,:,:,:,:)=czero

c--- only compute finite part here
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do j=1,3

      box(h1,h2,h3,h4,0)=box(h1,h2,h3,h4,0)+d(h1,h2,h3,h4,j)*intD0(j)

      enddo
      enddo
      enddo
      enddo
      enddo
      
c      write(6,*) 'box(1,1,1,1,0)',box(1,1,1,1,0)

      return
      end
      
