      subroutine ggZZmassamp_new(za,zb,sprod,mass,AmassLL,AmassLR)
      use mod_types; use common_def; use mod_consts_dp
      use mod_auxfunctions
      implicit none
c--- Author: J. Campbell, September 2013
c---
c--- Given momentum p and spinor products za,zb, calculate contribution
c--- of gg->ZZ through a loop of fermions with mass "mass".
c--- Returns amplitudes AmassLL, AmassLR indexed by
c--- helicities (h1=gluon,h2=gluon,h34=Z(34),h56=Z(p56))
c--- with LL and LR couplings in loop
      include 'ggZZcomputemp.f'
      real(dp)  :: mass,p(6,4)
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)
      complex(dp)  :: bub(2,2,2,2,-2:0),box(2,2,2,2,-2:0),
     & tri(2,2,2,2,-2:0),drat(2,2,2,2,3),totrat(2,2,2,2),
     & A(6),Xmp(2,2),Xpp(2,2),Xpm(2,2),Xmm(2,2),
     & AmassLL(2,2,2,2),AmassLR(2,2,2,2)
      
c--- initialize all scalar integrals
!      call ZZintegraleval(p,mass)   

c--- First calculate ALR
      call Acalc(1,2,3,4,5,6,sprod,mass,A)
      call LRcalc(1,2,3,4,5,6,za,zb,sprod,A,xpp,xmp,xpm,xmm)
      AmassLR(1,2,:,:)=Xmp(:,:)
      AmassLR(2,2,:,:)=Xpp(:,:)
      AmassLR(1,1,:,:)=Xmm(:,:)
      AmassLR(2,1,:,:)=Xpm(:,:)

c--- boxes      
      call ZZmassivebox(1,2,3,4,5,6,za,zb,sprod,mass,box,drat)
c--- bubbles
      call ZZmassivebub(1,2,3,4,5,6,za,zb,sprod,mass,bub,totrat)

c--- triangles -  note: calculated last to allow easy calculation
c--- of three-mass triangle from rational and drat (passed in)
      call ZZmassivetri(1,2,3,4,5,6,za,zb,sprod,mass,totrat,drat,tri)
      AmassLL(:,:,:,:)=box(:,:,:,:,0)+tri(:,:,:,:,0)+bub(:,:,:,:,0)

c--- new implementation using 6-d boxes for (-,+) and (+,-) amplitudes
      call ZZmassiveboxtri(1,2,3,4,5,6,za,zb,sprod,mass,totrat,box,tri)

      AmassLL(1,2,:,:)=box(1,2,:,:,0)+tri(1,2,:,:,0)+bub(1,2,:,:,0)
      AmassLL(2,1,:,:)=box(2,1,:,:,0)+tri(2,1,:,:,0)+bub(2,1,:,:,0)

c      if (docheck) then
cc-- check final result for LL
c        call ggZZparsecheck('LL',-1,-1,AmassLL(1,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('LL',-1,+1,AmassLL(1,2,:,:),'0+1,res/lo')
c        call ggZZparsecheck('LL',+1,-1,AmassLL(2,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('LL',+1,+1,AmassLL(2,2,:,:),'0+1,res/lo')
cc-- check that final result for RR is same as for LL
c        call ggZZparsecheck('RR',-1,-1,AmassLL(1,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('RR',-1,+1,AmassLL(1,2,:,:),'0+1,res/lo')
c        call ggZZparsecheck('RR',+1,-1,AmassLL(2,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('RR',+1,+1,AmassLL(2,2,:,:),'0+1,res/lo')
cc-- check final result for LR
c        call ggZZparsecheck('LR',-1,-1,AmassLR(1,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('LR',-1,+1,AmassLR(1,2,:,:),'0+1,res/lo')
c        call ggZZparsecheck('LR',+1,-1,AmassLR(2,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('LR',+1,+1,AmassLR(2,2,:,:),'0+1,res/lo')
cc-- check that final result for RL is same as for LR
c        call ggZZparsecheck('RL',-1,-1,AmassLR(1,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('RL',-1,+1,AmassLR(1,2,:,:),'0+1,res/lo')
c        call ggZZparsecheck('RL',+1,-1,AmassLR(2,1,:,:),'0+1,res/lo')
c        call ggZZparsecheck('RL',+1,+1,AmassLR(2,2,:,:),'0+1,res/lo')
c      endif

      return
      end
      
