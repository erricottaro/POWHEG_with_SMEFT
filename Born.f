!#define USE_OPENLOOPS
c     computation of the Born amplitude
      subroutine setborn(pin,bflav,born,bornjk,bmunu)
      use mod_ampsqLO_zz, only: resLO_ZZ_heft, resLO_ZZ_heft_smeft, resLO_sping_ZZ_heft
      use mod_ampsqLO_ww, only: resLO_WW,resLO_sping_WW
      use mod_auxfunctions, only: getPolVectors
      use openloops_powheg, only: ol_loop2 => openloops_loop2
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 pin(0:3,nlegbornexternal),p(0:3,nlegbornexternal)
      integer bflav(nlegbornexternal)
      real * 8 born, OLborn
      real *8 bornjk(nlegbornexternal,nlegbornexternal),bmunu(0:3,0:3
     $     ,nlegbornexternal)
      integer idvecbos,vdecaymode,Vdecmod
      common/cvecbos/idvecbos,vdecaymode,Vdecmod
      integer mu,nu,j,k,l,lp
      complex * 16 resHel(-1:1,-1:1)
      complex * 16 eps1(0:3,-1:1), eps2(0:3,-1:1)
      logical ini
      data ini/.true./
      logical pwhg_isfinite
      external pwhg_isfinite
      common /OLborn/OLborn
c
      if((.not.pwhg_isfinite(kn_jacborn)).or.(kn_jacborn.eq.0d0)) then
         born=0d0
         bornjk=0
         bmunu=0
         return
      endif

	  !-this is just for testing the scaling of virtual cross section with respect to SMEFT coupling parameters (they scale correctly)
      ! call  setvirtual(pin,bflav,born)
      ! born = born * st_alpha/(2d0*pi)
      ! return

#ifdef USE_OPENLOOPS
      call ol_loop2(pin,bflav,OLborn,bornjk,bmunu,approx=trim(flg_approx))
      born=OLborn
      return
#else
      if (ini) then
         call qlinit()
         ini=.false.
      endif

c    change lepton ordering
      p = pin
      if(bflav(3)<0) then
         p(:,3) = pin(:,4)
         p(:,4) = pin(:,3)
      endif
      if(bflav(5)<0) then
         p(:,5) = pin(:,6)
         p(:,6) = pin(:,5)
      endif

      if (flg_bornonly) then
        if (flg_proc.eq."ZZ") then
           call resLO_ZZ_heft_smeft(p,ph_tmass,born) ! done with SMEFT parameters
           !call resLO_ZZ_heft(p,ph_tmass,born)
        elseif (flg_proc.eq."WW")  then
           call resLO_WW(p,ph_tmass,born)
        endif
        born= born*(st_alpha/2d0/pi)**2
      else
        if (flg_proc.eq."ZZ") then
           call resLO_sping_ZZ_heft(p,ph_tmass,resHel)
        elseif (flg_proc.eq."WW")  then
           call resLO_sping_WW(p,ph_tmass,resHel)
        endif
        born=resHel(-1,-1)+resHel(1,1)
        born= born*(st_alpha/2d0/pi)**2

c       initialization of bornjk
        bornjk=0d0
        bornjk(1,2)= ca*born
        bornjk(2,1)=bornjk(1,2)

c       spin correlated born amplitude
        call getPolVectors(p,eps1,eps2)
        bmunu=0d0
        do mu=0,3
           do nu=0,3
              do l=-1,1,2
                 do lp=-1,1,2
                    bmunu(mu,nu,1) = bmunu(mu,nu,1) + resHel(l,lp)
     $                 *conjg(eps1(mu,l))*eps1(nu,lp)
                    bmunu(mu,nu,2) = bmunu(mu,nu,2) + resHel(l,lp)
     $                 *conjg(eps2(mu,l))*eps2(nu,lp)
                 enddo
              enddo
           enddo
        enddo
        bmunu=bmunu*(st_alpha/2d0/pi)**2
      end if

#endif

!      print*, born, OLborn, born/OLborn-1

      end

      subroutine borncolour_lh
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      icolup=0
c     colored particles
      icolup(1,1)=501
      icolup(2,1)=502
      icolup(1,2)=502
      icolup(2,2)=501
      end



      subroutine finalize_lh
      implicit none
      include 'LesHouches.h'
      integer i1,i2
      real * 8 v(5),v1
      integer itmp2(2)
      real * 8 random
      external random
c     give masses to final-state light particles
      call lhefinitemasses
      end


