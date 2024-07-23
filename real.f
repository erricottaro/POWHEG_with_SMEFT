      subroutine setreal(pin,rflav,amp2)
!      use mod_checkgram, only: grambox,gramtri
      use openloops_powheg, only: ol_loop2real => openloops_loop2real
      use mod_ampsqR_zz, only: resR_ZZ_heft
      implicit none
      real*8 dotp
      external dotp
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_flst.h'
      real *8 pin(0:3,nlegbornexternal+1),amp2,amp2_smeft(-2:0),amp2_ol,acc, accmax !-amp2_ol is the value of the sq amplitude given by openloops
      integer rflav(nlegbornexternal+1),mu,j
      real * 8 plan_cut
      parameter (plan_cut=1d-7)
      logical ini
      data ini/.true./
      save ini
      real * 8 ptVVcut
      save ptVVcut
      real *8 powheginput,coplanarity,coplan
      external powheginput
      logical pwhg_isfinite
      external pwhg_isfinite
      common /OL_stability/acc, accmax
      logical, save :: noexternalvqq = .false.

c     Set OL accuracy to very high value: if OL manages to calculate a stable point, acc is set to a small value
      acc=1d8

      if (ini) then
         ptVVcut=powheginput("#ptVVcut")
         if (ptVVcut.lt.0d0) then
             write(*,*) "For stable results ptVVcut has been set to 0.1"
             ptVVcut=0.1d0
         endif
         if(powheginput("#ol_noexternalvqq").eq.1d0) then
           noexternalvqq = .true.
         endif
         ini=.false.
      endif

      if(flg_btildepart.ne.'R') then
         if (.not.pwhg_isfinite(kn_jacborn)
     $        .or.(kn_jacborn.eq.0d0)) then
            call increasecnt('real discarded no valid jacborn')
            amp2=0d0
            return
         endif
      endif
c      if (0.5d0*sqrt(kn_sreal*(1d0-kn_y**2))*kn_csi.lt.
c     $        ptVVcut) then
      if (sqrt(kn_cmpreal(1,flst_ireallength)**2+kn_cmpreal(2
     $     ,flst_ireallength)**2).lt.ptVVcut) then

            call increasecnt('real discarded ptVV below cut')
            amp2 = 0d0
            return
         endif


c     if the gluon and the two Vs lie on a plane the gram is singular
      coplan = coplanarity(pin,3,4,9)
      if (coplan.le.plan_cut) then
c      too singular even for quad precision to fix it, throw it away
         call increasecnt('real discarded coplan cut')
         amp2=0d0
         return
      endif

c     here we introduce the SMEFT amplitudes (can be used only with only_h)
	  call resR_ZZ_heft(pin, ph_tmass, amp2_smeft)
c     i do not know which component of the output I should use
	  amp2=amp2_smeft(0)
c      if (noexternalvqq .and. (rflav(1) /= 0 .or. rflav(2) /= 0) .and.
c     &      (flg_approx == "" .or. flg_approx == "interf" .or. flg_approx == "noh")
c     &      ) then
c        if (flg_approx == "interf" .or. flg_approx == "noh") then
c          call ol_loop2real(pin,rflav,amp2_ol,acc,trim(flg_approx) // "_noexternalvqq")
c        else
c          call ol_loop2real(pin,rflav,amp2_ol,acc,trim(flg_approx) // "noexternalvqq")
c        end if
c      else
c        call ol_loop2real(pin,rflav,amp2_ol,acc,trim(flg_approx))
c      end if
c      call openloops_real(p,rflavnores,amp2_ol)
c      amp2=amp2_ol
      end


      function coplanarity(p,iz1,iz2,ig)
c     rotates until Z1 has zero x component
c     evaluate x component of gluon
c     normalize to its 3-momentum
       real *8 coplanarity
       real *8 p(0:3,*)
       real *8 pT1,pmod2
       integer i1,i2
       pTz1=sqrt(p(1,iz1)**2+p(2,iz1)**2)
       pMz1=sqrt(p(1,iz1)**2+p(2,iz1)**2+p(3,iz1)**2)
       coplanarity=abs((p(2,iz1)*p(1,iz2)-p(1,iz1)*p(2,iz2))/pTz1/pMz1)
      end

      subroutine regularcolour_lh
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
     
      icolup=0
      
c     colored particles
      if(idup(1).gt.0) then
c     q qbar
         icolup(1,1)=501
         icolup(2,1)=0
         icolup(1,2)=0
         icolup(2,2)=502
      else
c     qbar q
         icolup(1,2)=501
         icolup(2,2)=0
         icolup(1,1)=0
         icolup(2,1)=502
      endif
c     Last particle is the emitted gluon
      icolup(1,nup)=501
      icolup(2,nup)=502
      end
