module mod_ampsqNLOV_zz
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_ampLO_ggzz
  use mod_ampNLO_ggzz
  use mod_ampLO_ggzz_mass
  use mod_ampNLO_ggzz_mass
  use mod_ampLO_ggh2zz
  use mod_ampNLO_ggh2zz


  implicit none
  private

  logical :: polecheck=.false.

  public :: resNLOV_ZZ_heft!,resNLOV_ZZ_heft_swap

  contains 


    ! full amp. sq. of NLO gg->ZZ in heavy top expansion   
    ! g g -> Z[l lb] Z[l' lb']
    ! coefficient of (as/twopi)**2, finite remainder in qT scheme
    ! debug --- original file name: sqNLO_ggzz.f90
    ! the mass is needed for the heavy top expansion
  subroutine resNLOV_ZZ_heft(p,mass,E1LL,E1LR,E2LL,E2LR,f_exp,resLO,resNLO,poles)
    
    implicit none
    real(dp), intent(in) :: p(4,6),mass,f_exp
    complex(dp), intent(in) :: E1LL(9),E1LR(9),E2LL(9),E2LR(9)
    real(dp), intent(out) :: resLO,resNLO,poles(-2:-1)
    real(dp) :: resHLO,resHNLO,resZZLO,resZZNLO,resZZHLO,resZZHNLO
    real(dp) :: resHLO_mt(0:nheft_max),resHNLO_mt(0:nheft_max)
    real(dp) :: resZZLO_mt(0:nheft_max),resZZNLO_mt(0:nheft_max)
    real(dp) :: resZZHLO_mt(0:nheft_max),resZZHNLO_mt(0:nheft_max)

    complex(dp) :: ampLO_zzh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampNLO_zzh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampLO_ggzz_all(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampNLO_ggzz_all(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampNLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampNLO_ggzz_rewgt(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampdum(-1:1,-1:1,-1:1,-1:1),ampdumheft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
 
    complex(dp15) :: ewampLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp15) :: ewampLO_ggzz_fullmass(-1:1,-1:1,-1:1,-1:1)
    complex(dp15) :: ewampNLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp15) :: ewampNLO_ggzz_eps(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2)

    complex(dp) :: ampLO_ggh2zz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampNLO_ggh2zz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_ggh2zz(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampNLO_ggh2zz(-1:1,-1:1,-1:1,-1:1)

    complex(dp) :: za(6,6),zb(6,6)
    complex(dp15) :: zadp(6,6),zbdp(6,6)
    real(dp)    :: pdum(4,6),sprod(6,6),pMCFM(6,4),pTZ,m4l
    real(dp15)  :: pMCFMdp(6,4),sproddp(6,6),massdp
    integer     :: i1,i2,i3,i5,j,nheft,neps
    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one
    complex(dp) :: diff
    complex(dp) :: ampDT(-1:1,-1:1,-1:1,-1:1)
    logical, save :: first=.true.

    resLO     = zero
    resNLO    = zero
    resHLO    = zero
    resHNLO   = zero
    resZZLO   = zero
    resZZNLO  = zero
    resZZHLO  = zero
    resZZHNLO = zero


    resHLO_mt    = zero
    resHNLO_mt   = zero
    resZZLO_mt   = zero
    resZZNLO_mt  = zero
    resZZHLO_mt  = zero
    resZZHNLO_mt = zero

    poles = zero

    ewampLO_ggzz_heft     = czero
    ewampNLO_ggzz_heft    = czero
    ewampNLO_ggzz_eps     = czero
    ewampNLO_ggzz_all     = czero
    ewampLO_ggzz_all      = czero
    ewampLO_ggzz_fullmass = czero
    ewampLO_ggzz_mless    = czero
    ewampNLO_ggzz_mless   = czero
    ampLO_zzh2zz          = czero
    ampNLO_zzh2zz         = czero
    ampLO_ggh2zz          = czero
    ampNLO_ggh2zz         = czero
    ampLO_ggh2zz_heft     = czero
    ampNLO_ggh2zz_heft    = czero
    
! set up spinors
    call awayfromzaxis(6,p,pdum)

    call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)


! calculate gg -> ZZ amps
    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then

       if (mlessloops) then
          if (first) print *, "including massless loops"
          ewampLO_ggzz_mless = ewampLO_ggzz(1,2,3,4,5,6,za,zb,sprod)
          ! massless implemented for gg->VV not 0->ggVV, so recompute za,zb,sprod
          call spinorur(6,pdum,za,zb,sprod)
          call getewampsNLO_ggzz(za,zb,sprod,E2LL,E2LR,ewampNLO_ggzz_mless)

          ! Raoul TODO: do we need this phase shift now?
          ! phase shift to have agreement with LO results (checked against MCFM)
          ewampNLO_ggzz_mless = ci * ewampNLO_ggzz_mless


          ! get "normal" spinors again
          call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)

       endif

       if (massloops) then
          if (first .and. VVExp) print *, "including LO massive loop in HEFT"
          if (first .and. .not. VVExp) print *, "including LO massive loop exactly"
          if (first) print *, "including NLO massive loops via reweighting and/or inverse mt expansion"
          call converttoMCFMmom(p,pMCFM)
          if (kind(p) .ne. 8) then
             print *, "converting to dp for massive amps", kind(p)
             pMCFMdp(1:6,1:4) = real(pMCFM(1:6,1:4),kind=dp15)
             massdp = real(mass,kind=dp15)
             zadp(1:6,1:6) = cmplx(za(1:6,1:6),kind=dp15)
             zbdp(1:6,1:6) = cmplx(zb(1:6,1:6),kind=dp15)
             sproddp(1:6,1:6)=real(sprod(1:6,1:6),kind=dp15)
          else
             zadp=za
             zbdp=zb
             sproddp=sprod
             pMCFMdp=pMCFM
             massdp=mass
          endif
          call ZZintegraleval(pMCFMdp,massdp)
          ! LO
          call getewampsLO_ggzz_heft(zadp,zbdp,sproddp,massdp,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft)

          ! 2-loop amplitudes -- with analytic and numerical subtr
          call getewampsNLO_ggzz_heft(zadp,zbdp,sproddp,massdp,ewampNLO_ggzz_heft)
          call getewampsNLO_ggzz_rewgt(zadp,zbdp,sproddp,ewampNLO_ggzz_mless,ewampLO_ggzz_mless,ewampLO_ggzz_fullmass,ewampNLO_ggzz_rewgt)

              


!          ! Only used for catani check
!          call getewampsNLO_ggzz_heft_eps(za,zb,sprod,mass,ewampNLO_ggzz_eps)
          
! if pTZ < 0.1 GeV, omit top loops
! this should be done at ph space generation level, but it's fine for now...
          pTZ = sqrt( (p(2,3)+p(2,4))**2 + (p(3,3)+p(3,4))**2 )
          if ( pTZ .le. 0.1d0) then
             !print *, "skipping point with pTZ=",pTZ
             ewampLO_ggzz_fullmass = czero
          endif
       endif

! only LO has exact mass dependence
       if (VVExp) then
          ewampLO_ggzz_all = ewampLO_ggzz_heft
          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_all(0,:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
       else
          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_fullmass(:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
       endif

       ! interpolate between 1/mt expansion and reweighting
       ewampNLO_ggzz_all = ewampNLO_ggzz_heft * f_exp
       ewampNLO_ggzz_all(0,:,:,:,:) = ewampNLO_ggzz_all(0,:,:,:,:) +  ewampNLO_ggzz_rewgt(:,:,:,:) *(one-f_exp)

       ! sum massless+(approx) massive 2 loop amps
       ewampNLO_ggzz_all(0,:,:,:,:) = ewampNLO_ggzz_all(0,:,:,:,:)+ewampNLO_ggzz_mless(:,:,:,:)

!       ! Only used for catani check
!       ! We don't have poles for massless so we assume -CA*LO gives correct ep^-2
!       ewampNLO_ggzz_all = ewampNLO_ggzz_eps(:,:,:,:,:,-2)
!       ewampNLO_ggzz_all(0,:,:,:,:) = ewampNLO_ggzz_all(0,:,:,:,:)-CA*ewampLO_ggzz_mless(:,:,:,:)
    endif
          

! calculate gg->H->ZZ amplitudes
    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if (HiggsExp) then     
          ! Higgs expanded
          if (first) print *, "Using heavy top expansion for Higgs amplitudes"
          ampLO_ggh2zz_heft =  ewampLO_ggh2zz_heft(1,2,3,4,5,6,za,zb,sprod)
          call getewamplsLONLO_ggh2zz_heft(1,2,3,4,5,6,za,zb,sprod,ampdumheft,ampNLO_ggh2zz_heft)
       else
          if (first) print *, "Using full mass dependence for Higgs amplitudes"
          ampLO_ggh2zz =  ewampLO_ggh2zz(1,2,3,4,5,6,za,zb,sprod)
          call getewamplsLONLO_ggh2zz(1,2,3,4,5,6,za,zb,sprod,ampdum,ampNLO_ggh2zz(:,:,:,:))
       endif
    endif

! add the Higgs and ZZ amplitudes
    if (contr .eq. "full" .or. contr .eq. "intf") then
       if (HiggsExp) then
          ampLO_zzh2zz = ewampLO_ggzz_all + ampLO_ggh2zz_heft
          ampNLO_zzh2zz = ewampNLO_ggzz_all + ampNLO_ggh2zz_heft
       else
          ampLO_zzh2zz = ewampLO_ggzz_all
          ampLO_zzh2zz(0,:,:,:,:) = ampLO_zzh2zz(0,:,:,:,:) + ampLO_ggh2zz
          ampNLO_zzh2zz = ewampNLO_ggzz_all
          ampNLO_zzh2zz(0,:,:,:,:) = ampNLO_zzh2zz(0,:,:,:,:) + ampNLO_ggh2zz
!       ! Only used for catani check
!       ampNLO_zzh2zz = ewampNLO_ggzz_all
!       ampNLO_zzh2zz(0,:,:,:,:) = ampNLO_zzh2zz(0,:,:,:,:) - CA*ampLO_ggh2zz
       endif
    endif

! ! CHECK THAT THE ANALYTIC SUBTR AND THE NUMERICAL SUBTR GIVE THE SAME RESULTS !!!!
! ! IMPORTANT: to check poles you must remove the Double Triagle contribution in
! !            getewampsNLO_ggzz_heft currently line 502

! ! check poles in numerical subtr implementation
!     do neps=-2,-1
!        print *, "checking pole neps=",neps
!     do i1=-1,+1,2
!     do i2=-1,+1,2
!     do i3=-1,+1,2
!     do i5=-1,+1,2
!     do nheft=0,nheft_max
!        if (abs(ewampNLO_ggzz_eps(nheft,i1,i2,i3,i5,neps)/ewampNLO_ggzz_eps(nheft,i1,i2,i3,i5,0)) .gt. 1d-10) then
!           print *, "WARNING: POLE POSSIBLY NON-ZERO!"
!           print *, nheft,i1,i2,i3,i5,ewampNLO_ggzz_eps(nheft,i1,i2,i3,i5,neps)
!           pause
!        endif
!     enddo
!     enddo
!     enddo
!     enddo
!     enddo
!     print *, "check done!"
!     enddo

! ! check finite piece for numerical vs analytic subtr
!     print *, "checking finite piece"
!     do i1=-1,+1,2
!     do i2=-1,+1,2
!     do i3=-1,+1,2
!     do i5=-1,+1,2
!     do nheft=0,nheft_max
!        if ( abs(ewampNLO_ggzz_heft(nheft,i1,i2,i3,i5)) .lt.1d-14 ) then
!           diff = ewampNLO_ggzz_eps(nheft,i1,i2,i3,i5,0)-ewampNLO_ggzz_heft(nheft,i1,i2,i3,i5)
!           if ( abs(diff) .gt. 1d-14 ) then
!              print *, "WARNING: POSSIBLE DISAGREEMENT IN FINITE PART! A"
!              print *,nheft,i1,i2,i3,i5,ewampNLO_ggzz_eps(nheft,i1,i2,i3,i5,0),ewampNLO_ggzz_heft(nheft,i1,i2,i3,i5)
!           endif
!        else
!           diff = (ewampNLO_ggzz_eps(nheft,i1,i2,i3,i5,0)-ewampNLO_ggzz_heft(nheft,i1,i2,i3,i5))/ewampNLO_ggzz_heft(nheft,i1,i2,i3,i5)
!           if ( abs(diff) .gt.1d-10 ) then
!              print *, "WARNING: POSSIBLE DISAGREEMENT IN FINITE PART! B"
!              print *, "diff=",diff
!              print *,nheft,i1,i2,i3,i5,ewampNLO_ggzz_eps(nheft,i1,i2,i3,i5,0),ewampNLO_ggzz_heft(nheft,i1,i2,i3,i5)
!              pause
!           endif
!        endif
!     enddo
!     enddo
!     enddo
!     enddo
!     enddo
!     pause
!     print *, "check done!"
! !!!!! END OF CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !    pause

! square the amps, keeping the order of the 1/mt expansion
!M the factors of 1/mt already exist in these amplitudes
    do i1=-1,1,2
    do i2=-1,1,2
    do i3=-1,1,2
    do i5=-1,1,2

! |H|^2
       if (HiggsExp) then
          ! LO amp-sq
          resHLO_mt(0) = resHLO_mt(0) + abs(ampLO_ggh2zz_heft(0,i1,i2,i3,i5))**2
          resHLO_mt(1) = resHLO_mt(1) + &
               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(1,i1,i2,i3,i5)))
          resHLO_mt(2) = resHLO_mt(2) + &
               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(2,i1,i2,i3,i5))) + &
               abs(ampLO_ggh2zz_heft(1,i1,i2,i3,i5))**2
          resHLO_mt(3) = resHLO_mt(3) + &
               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(3,i1,i2,i3,i5))) + &
               2.0_dp*real(ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(2,i1,i2,i3,i5)))
          resHLO_mt(4) = resHLO_mt(4) + &
               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(4,i1,i2,i3,i5))) + &
               2.0_dp*real(ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(3,i1,i2,i3,i5))) + &
               abs(ampLO_ggh2zz_heft(2,i1,i2,i3,i5))**2

          ! LO*NLO
          resHNLO_mt(0) = resHNLO_mt(0) + &
               2.0_dp*real( ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5)))
          resHNLO_mt(1) = resHNLO_mt(1) + &
               2.0_dp*real( & 
               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5)) &
               )
          resHNLO_mt(2) = resHNLO_mt(2) + &
               2.0_dp*real( &
               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(2,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(2,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5))   &
               )
          resHNLO_mt(3) = resHNLO_mt(3) + &
               2.0_dp*real( &
               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(3,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(2,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(2,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(3,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5))   &
               )
          resHNLO_mt(4) = resHNLO_mt(4) + &
               2.0_dp*real( &
               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(4,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(3,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(2,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(2,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(3,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
               ampLO_ggh2zz_heft(4,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5))   &
               )
       else
       
          ! LO amp-sq
          ! print *, "ampLO_ggh2zz",i1, i2, i3, i5,  ampLO_ggh2zz(i1,i2,i3,i5) 
          resHLO = resHLO + abs(ampLO_ggh2zz(i1,i2,i3,i5))**2
          ! print *, "Partial resHLO", i1, i2, i3, i5, resHLO
          ! LO*NLO
          resHNLO = resHNLO + 2.0_dp*real(ampLO_ggh2zz(i1,i2,i3,i5)*conjg(ampNLO_ggh2zz(i1,i2,i3,i5)))
       endif
! |ZZ|^2
       ! LO amp-sq
       resZZLO_mt(0) = resZZLO_mt(0) + abs(ewampLO_ggzz_all(0,i1,i2,i3,i5))**2
       resZZLO_mt(1) = resZZLO_mt(1) + &
             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,i1,i2,i3,i5)))
       resZZLO_mt(2) = resZZLO_mt(2) + &
             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5))) + &
             abs(ewampLO_ggzz_all(1,i1,i2,i3,i5))**2
       resZZLO_mt(3) = resZZLO_mt(3) + &
            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5)))
       resZZLO_mt(4) = resZZLO_mt(4) + &
            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(4,i1,i2,i3,i5))) + &
            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
            abs(ewampLO_ggzz_all(2,i1,i2,i3,i5))**2

       ! LO*NLO
       resZZNLO_mt(0) = resZZNLO_mt(0) + &
            2.0_dp*real( ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5)))
       resZZNLO_mt(1) = resZZNLO_mt(1) + &
            2.0_dp*real( & 
            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5)) &
            )
       resZZNLO_mt(2) = resZZNLO_mt(2) + &
            2.0_dp*real( &
            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(2,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5))   &
            )
       resZZNLO_mt(3) = resZZNLO_mt(3) + &
            2.0_dp*real( &
            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(3,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(2,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(3,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5))   &
            )
       
       resZZNLO_mt(4) = resZZNLO_mt(4) + &
            2.0_dp*real( &
            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(4,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(3,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(2,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(3,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
            ewampLO_ggzz_all(4,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5))   &
            )

! | ZZ + HZZ |^2

       ! LO amp-sq
       resZZHLO_mt(0) = resZZHLO_mt(0) + abs(ampLO_zzh2zz(0,i1,i2,i3,i5))**2
       resZZHLO_mt(1) = resZZHLO_mt(1) + &
            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,i1,i2,i3,i5)))
       resZZHLO_mt(2) = resZZHLO_mt(2) + &
            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5))) + &
            abs(ampLO_zzh2zz(1,i1,i2,i3,i5))**2
       resZZHLO_mt(3) = resZZHLO_mt(3) + &
            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5)))
       resZZHLO_mt(4) = resZZHLO_mt(4) + &
            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(4,i1,i2,i3,i5))) + &
            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
            abs(ampLO_zzh2zz(2,i1,i2,i3,i5))**2


       ! LO*NLO
       resZZHNLO_mt(0) = resZZHNLO_mt(0) + &
            2.0_dp*real( ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5)))
       resZZHNLO_mt(1) = resZZHNLO_mt(1) + &
            2.0_dp*real( & 
            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5)) &
            )
       resZZHNLO_mt(2) = resZZHNLO_mt(2) + &
            2.0_dp*real( &
            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(2,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5))   &
            )
       resZZHNLO_mt(3) = resZZHNLO_mt(3) + &
            2.0_dp*real( &
            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(3,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(2,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(3,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5))   &
            )
       
       resZZHNLO_mt(4) = resZZHNLO_mt(4) + &
            2.0_dp*real( &
            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(4,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(3,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(2,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(3,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
            ampLO_zzh2zz(4,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5))   &
            )

    enddo
    enddo
    enddo
    enddo

    ! add the various pieces
    ! if-statements are redundant, just put in to be extra safe

    do nheft=0,expord
          resZZLO = resZZLO + resZZLO_mt(nheft)
          resZZNLO = resZZNLO + resZZNLO_mt(nheft)
       
       if(HiggsExp) then
       	  !print *, "we shouldn't be here"
          resHLO = resHLO + resHLO_mt(nheft)
          resHNLO = resHNLO + resHNLO_mt(nheft)
       endif
          resZZHLO = resZZHLO + resZZHLO_mt(nheft)
          resZZHNLO = resZZHNLO + resZZHNLO_mt(nheft)
    enddo
    
    first=.false.

! define res according to the 'contr' input
    if (contr .eq. "bkgd") then                ! |gg->ZZ|^2
       if (first) print *, "Doing background |gg->ZZ|^2"
       resLO = resZZLO
       resNLO = resZZNLO
    elseif (contr .eq. "sigl") then            ! |gg->H->ZZ|^2
       if (first) print *, "Doing signal |gg->H->ZZ|^2"
       resLO = resHLO
       resNLO = resHNLO
    elseif (contr .eq. "full") then            ! |gg->H->ZZ  +  gg->ZZ|^2
       if (first) print *, "Doing full  |gg->H->ZZ + gg->ZZ|^2"
       resLO = resZZHLO
       resNLO = resZZHNLO
    elseif (contr .eq. "intf") then            ! |gg->H->ZZ + gg->ZZ|^2 - |gg->ZZ|^2 - |gg->H->ZZ|^2 = 2*Real( (gg->H->ZZ)*conj(gg->ZZ) )
       resLO = resZZHLO - resZZLO - resHLO
       resNLO = resZZHNLO - resZZNLO - resHNLO
       if (first) print *, "Doing interference 2Real( (gg->H->ZZ)*conj(gg->ZZ) )"
    endif

    !-- fix pre-factor
    resLO = resLO * avegg * colf 
    resNLO = resNLO * avegg * colf 

!    print *, "resLO",resLO
!    print *, "resNLO",resNLO
!    stop

! poles, as in massless piece
! might want to check on this order-by-order in 1/mt exp.
    poles(-2) = -2*CA*resLO
    poles(-1) = (-two * b0 - two * CA * log(musq/sprod(1,2)))*resLO
!    print *, "DP/1l Catani",poles(-2)/resLO
!    print *, "SP/1l Catani",poles(-1)/resLO
!    stop

    return
  end subroutine resNLOV_ZZ_heft




!!--!!    ! full amp. sq. of NLO gg->ZZ in heavy top expansion   
!!--!!    ! g g -> Z[l lb] Z[l lb]
!!--!!    ! same decays into same flavour leptons (e-e+e-e+ etc)
!!--!!    ! coefficient of (as/twopi)**2, finite remainder in qT scheme
!!--!!    ! identical to above routine but with loop over swap of momenta
!!--!!    ! debug --- original file name: sqNLO_ggzz.f90
!!--!!  subroutine resNLOV_ZZ_heft_swap(p,mass,E1LL,E1LR,E2LL,E2LR,E1LL_swap,E1LR_swap,E2LL_swap,E2LR_swap,resLO,resNLO,poles)
!!--!!    
!!--!!    implicit none
!!--!!    real(dp), intent(in) :: p(4,6),mass
!!--!!    complex(dp), intent(in) :: E1LL(9),E1LR(9),E2LL(9),E2LR(9)
!!--!!    real(dp), intent(out) :: resLO,resNLO,poles(-2:-1)
!!--!!    real(dp) :: resHLO,resHNLO,resZZLO,resZZNLO,resZZHLO,resZZHNLO
!!--!!    real(dp) :: resHLO_mt(0:nheft_max),resHNLO_mt(0:nheft_max)
!!--!!    real(dp) :: resZZLO_mt(0:nheft_max),resZZNLO_mt(0:nheft_max)
!!--!!    real(dp) :: resZZHLO_mt(0:nheft_max),resZZHNLO_mt(0:nheft_max)
!!--!!
!!--!!    complex(dp) :: ampLO_zzh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ampNLO_zzh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ewampLO_ggzz_all(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ewampNLO_ggzz_all(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ewampLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ewampNLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ewampNLO_ggzz_rewgt(-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ampdum(-1:1,-1:1,-1:1,-1:1),ampdumheft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!! 
!!--!!    complex(dp15) :: ewampLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp15) :: ewampLO_ggzz_fullmass(-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp15) :: ewampNLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp15) :: ewampNLO_ggzz_eps(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-2:2)
!!--!!
!!--!!    complex(dp) :: ampLO_ggh2zz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ampNLO_ggh2zz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ampLO_ggh2zz(-1:1,-1:1,-1:1,-1:1)
!!--!!    complex(dp) :: ampNLO_ggh2zz(-1:1,-1:1,-1:1,-1:1)
!!--!!
!!--!!    complex(dp) :: za(6,6),zb(6,6)
!!--!!    complex(dp15) :: zadp(6,6),zbdp(6,6)
!!--!!    real(dp)    :: pdum(4,6),sprod(6,6),pMCFM(6,4),pTZ,m4l
!!--!!    real(dp15)  :: pMCFMdp(6,4),sproddp(6,6),massdp
!!--!!    integer     :: i1,i2,i3,i5,j,nheft,neps,iflav
!!--!!    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one
!!--!!    complex(dp) :: diff
!!--!!    complex(dp) :: ampDT(-1:1,-1:1,-1:1,-1:1)
!!--!!    logical       :: rewgt_2loop,expand_2loop
!!--!!    logical, save :: first=.true.
!!--!!
!!--!!    resLO     = zero
!!--!!    resNLO    = zero
!!--!!    resHLO    = zero
!!--!!    resHNLO   = zero
!!--!!    resZZLO   = zero
!!--!!    resZZNLO  = zero
!!--!!    resZZHLO  = zero
!!--!!    resZZHNLO = zero
!!--!!
!!--!!
!!--!!    resHLO_mt    = zero
!!--!!    resHNLO_mt   = zero
!!--!!    resZZLO_mt   = zero
!!--!!    resZZNLO_mt  = zero
!!--!!    resZZHLO_mt  = zero
!!--!!    resZZHNLO_mt = zero
!!--!!
!!--!!    poles = zero
!!--!!
!!--!!    ewampLO_ggzz_heft     = czero
!!--!!    ewampNLO_ggzz_heft    = czero
!!--!!    ewampNLO_ggzz_eps     = czero
!!--!!    ewampNLO_ggzz_all     = czero
!!--!!    ewampLO_ggzz_all      = czero
!!--!!    ewampLO_ggzz_fullmass = czero
!!--!!    ewampLO_ggzz_mless    = czero
!!--!!    ewampNLO_ggzz_mless   = czero
!!--!!    ampLO_zzh2zz          = czero
!!--!!    ampNLO_zzh2zz         = czero
!!--!!    ampLO_ggh2zz          = czero
!!--!!    ampNLO_ggh2zz         = czero
!!--!!    ampLO_ggh2zz_heft     = czero
!!--!!    ampNLO_ggh2zz_heft    = czero
!!--!!    
!!--!!! set up spinors
!!--!!
!!--!!    do iflav = 1,2
!!--!!
!!--!!       if (iflav .eq. 2) then    ! swap 4 < -- > 6
!!--!!          pswap(:,4) = p(:,6)
!!--!!          pswap(:,6) = p(:,4)
!!--!!          p(:,4) = pswap(:,4)
!!--!!          p(:,6) = pswap(:,6)
!!--!!       endif
!!--!!
!!--!!       call awayfromzaxis(6,p,pdum)
!!--!!
!!--!!
!!--!!       call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)
!!--!!
!!--!!
!!--!!    m4l = sqrt(sprod(1,2))
!!--!!    if (m4l .gt. 2*mt) then
!!--!!       rewgt_2loop = .true.
!!--!!       expand_2loop = .false.
!!--!!    else
!!--!!       rewgt_2loop = .false.
!!--!!       expand_2loop = .true.
!!--!!    endif
!!--!!
!!--!!
!!--!!! calculate gg -> ZZ amps
!!--!!    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then
!!--!!
!!--!!       if (mlessloops) then
!!--!!          if (first) print *, "including massless loops"
!!--!!          ewampLO_ggzz_mless = ewampLO_ggzz(1,2,3,4,5,6,za,zb,sprod)
!!--!!          ! massless implemented for gg->VV not 0->ggVV, so recompute za,zb,sprod
!!--!!          call spinorur(6,pdum,za,zb,sprod)
!!--!!          if (iflav .eq. 1) then
!!--!!             call getewampsNLO_ggzz(za,zb,sprod,E2LL,E2LR,ewampNLO_ggzz_mless)
!!--!!          elseif (iflav .eq. 2 ) then             
!!--!!             call getewampsNLO_ggzz(za,zb,sprod,E2LL_swap,E2LR_swap,ewampNLO_ggzz_mless)
!!--!!          endif
!!--!!
!!--!!          ! Raoul TODO: do we need this phase shift now?
!!--!!          ! phase shift to have agreement with LO results (checked against MCFM)
!!--!!          ewampNLO_ggzz_mless = ci * ewampNLO_ggzz_mless
!!--!!
!!--!!
!!--!!          ! get "normal" spinors again
!!--!!          call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)
!!--!!
!!--!!       endif
!!--!!
!!--!!       if (massloops) then
!!--!!          if (first .and. VVExp) print *, "including LO massive loop in HEFT"
!!--!!          if (first .and. .not. VVExp) print *, "including LO massive loop exactly"
!!--!!          if (first .and. rewgt_2loop) print *, "including NLO massive loops via reweighting"
!!--!!          if (first .and. expand_2loop) print *, "including NLO massive loops via expansion"
!!--!!          call converttoMCFMmom(p,pMCFM)
!!--!!          if (kind(p) .ne. 8) then
!!--!!             print *, "converting to dp for massive amps", kind(p)
!!--!!             pMCFMdp(1:6,1:4) = real(pMCFM(1:6,1:4),kind=dp15)
!!--!!             massdp = real(mass,kind=dp15)
!!--!!             zadp(1:6,1:6) = cmplx(za(1:6,1:6),kind=dp15)
!!--!!             zbdp(1:6,1:6) = cmplx(zb(1:6,1:6),kind=dp15)
!!--!!             sproddp(1:6,1:6)=real(sprod(1:6,1:6),kind=dp15)
!!--!!          else
!!--!!             zadp=za
!!--!!             zbdp=zb
!!--!!             sproddp=sprod
!!--!!             pMCFMdp=pMCFM
!!--!!             massdp=mass
!!--!!          endif
!!--!!          call ZZintegraleval(pMCFMdp,massdp)
!!--!!          ! LO
!!--!!          call getewampsLO_ggzz_heft(zadp,zbdp,sproddp,massdp,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft)
!!--!!
!!--!!          ! 2-loop amplitudes -- with analytic and numerical subtr
!!--!!          if (expand_2loop) then
!!--!!             call getewampsNLO_ggzz_heft(zadp,zbdp,sproddp,massdp,ewampNLO_ggzz_heft)
!!--!!           endif
!!--!!           if (rewgt_2loop) then
!!--!!              call getewampsNLO_ggzz_rewgt(zadp,zbdp,sproddp,ewampNLO_ggzz_mless,ewampLO_ggzz_mless,ewampLO_ggzz_fullmass,ewampNLO_ggzz_rewgt)
!!--!!           endif
!!--!!
!!--!!
!!--!!!          ! Only used for catani check
!!--!!!          call getewampsNLO_ggzz_heft_eps(za,zb,sprod,mass,ewampNLO_ggzz_eps)
!!--!!          
!!--!!! if pTZ < 0.1 GeV, omit top loops
!!--!!! this should be done at ph space generation level, but it's fine for now...
!!--!!          pTZ = sqrt( (p(2,3)+p(2,4))**2 + (p(3,3)+p(3,4))**2 )
!!--!!          if ( pTZ .le. 0.1d0) then
!!--!!             !print *, "skipping point with pTZ=",pTZ
!!--!!             ewampLO_ggzz_fullmass = czero
!!--!!          endif
!!--!!       endif
!!--!!
!!--!!! only LO has exact mass dependence
!!--!!       if (VVExp) then
!!--!!          ewampLO_ggzz_all = ewampLO_ggzz_heft
!!--!!          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_all(0,:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
!!--!!       else
!!--!!          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_fullmass(:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
!!--!!       endif
!!--!!
!!--!!       ! either reweight or use expansion for 2loop amps
!!--!!       if (expand_2loop) then
!!--!!          ewampNLO_ggzz_all = ewampNLO_ggzz_heft
!!--!!       endif
!!--!!       if (rewgt_2loop) then
!!--!!          ewampNLO_ggzz_all(0,:,:,:,:) = ewampNLO_ggzz_rewgt(:,:,:,:)
!!--!!       endif
!!--!!       ! sum massless+(approx) massive 2 loop amps
!!--!!       ewampNLO_ggzz_all(0,:,:,:,:) = ewampNLO_ggzz_all(0,:,:,:,:)+ewampNLO_ggzz_mless(:,:,:,:)
!!--!!
!!--!!    endif
!!--!!          
!!--!!
!!--!!! calculate gg->H->ZZ amplitudes
!!--!!    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
!!--!!       if (HiggsExp) then     
!!--!!          ! Higgs expanded
!!--!!          if (first) print *, "Using heavy top expansion for Higgs amplitudes"
!!--!!          ampLO_ggh2zz_heft =  ewampLO_ggh2zz_heft(1,2,3,4,5,6,za,zb,sprod)
!!--!!          call getewamplsLONLO_ggh2zz_heft(1,2,3,4,5,6,za,zb,sprod,ampdumheft,ampNLO_ggh2zz_heft)
!!--!!       else
!!--!!          if (first) print *, "Using full mass dependence for Higgs amplitudes"
!!--!!          ampLO_ggh2zz =  ewampLO_ggh2zz(1,2,3,4,5,6,za,zb,sprod)
!!--!!          call getewamplsLONLO_ggh2zz(1,2,3,4,5,6,za,zb,sprod,ampdum,ampNLO_ggh2zz(:,:,:,:))
!!--!!       endif
!!--!!    endif
!!--!!
!!--!!! add the Higgs and ZZ amplitudes
!!--!!    if (contr .eq. "full" .or. contr .eq. "intf") then
!!--!!       if (HiggsExp) then
!!--!!          ampLO_zzh2zz = ewampLO_ggzz_all + ampLO_ggh2zz_heft
!!--!!          ampNLO_zzh2zz = ewampNLO_ggzz_all + ampNLO_ggh2zz_heft
!!--!!       else
!!--!!          ampLO_zzh2zz = ewampLO_ggzz_all
!!--!!          ampLO_zzh2zz(0,:,:,:,:) = ampLO_zzh2zz(0,:,:,:,:) + ampLO_ggh2zz
!!--!!          ampNLO_zzh2zz = ewampNLO_ggzz_all
!!--!!          ampNLO_zzh2zz(0,:,:,:,:) = ampNLO_zzh2zz(0,:,:,:,:) + ampNLO_ggh2zz
!!--!!       endif
!!--!!    endif
!!--!!
!!--!!
!!--!!    
!!--!!
!!--!! enddo
!!--!!
!!--!!! square the amps, keeping the order of the 1/mt expansion
!!--!!!M the factors of 1/mt already exist in these amplitudes
!!--!!    do i1=-1,1,2
!!--!!    do i2=-1,1,2
!!--!!    do i3=-1,1,2
!!--!!    do i5=-1,1,2
!!--!!
!!--!!! |H|^2
!!--!!       if (HiggsExp) then
!!--!!          ! LO amp-sq
!!--!!          resHLO_mt(0) = resHLO_mt(0) + abs(ampLO_ggh2zz_heft(0,i1,i2,i3,i5))**2
!!--!!          resHLO_mt(1) = resHLO_mt(1) + &
!!--!!               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(1,i1,i2,i3,i5)))
!!--!!          resHLO_mt(2) = resHLO_mt(2) + &
!!--!!               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(2,i1,i2,i3,i5))) + &
!!--!!               abs(ampLO_ggh2zz_heft(1,i1,i2,i3,i5))**2
!!--!!          resHLO_mt(3) = resHLO_mt(3) + &
!!--!!               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(3,i1,i2,i3,i5))) + &
!!--!!               2.0_dp*real(ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(2,i1,i2,i3,i5)))
!!--!!          resHLO_mt(4) = resHLO_mt(4) + &
!!--!!               2.0_dp*real(ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(4,i1,i2,i3,i5))) + &
!!--!!               2.0_dp*real(ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz_heft(3,i1,i2,i3,i5))) + &
!!--!!               abs(ampLO_ggh2zz_heft(2,i1,i2,i3,i5))**2
!!--!!
!!--!!          ! LO*NLO
!!--!!          resHNLO_mt(0) = resHNLO_mt(0) + &
!!--!!               2.0_dp*real( ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5)))
!!--!!          resHNLO_mt(1) = resHNLO_mt(1) + &
!!--!!               2.0_dp*real( & 
!!--!!               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5)) &
!!--!!               )
!!--!!          resHNLO_mt(2) = resHNLO_mt(2) + &
!!--!!               2.0_dp*real( &
!!--!!               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(2,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(2,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5))   &
!!--!!               )
!!--!!          resHNLO_mt(3) = resHNLO_mt(3) + &
!!--!!               2.0_dp*real( &
!!--!!               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(3,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(2,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(2,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(3,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5))   &
!!--!!               )
!!--!!          resHNLO_mt(4) = resHNLO_mt(4) + &
!!--!!               2.0_dp*real( &
!!--!!               ampLO_ggh2zz_heft(0,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(4,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(1,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(3,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(2,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(2,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(3,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(1,i1,i2,i3,i5)) + &
!!--!!               ampLO_ggh2zz_heft(4,i1,i2,i3,i5)*conjg(ampNLO_ggh2zz_heft(0,i1,i2,i3,i5))   &
!!--!!               )
!!--!!       else
!!--!!       
!!--!!          ! LO amp-sq
!!--!!          resHLO = resHLO + abs(ampLO_ggh2zz(i1,i2,i3,i5))**2
!!--!!          ! LO*NLO
!!--!!          resHNLO = resHNLO + 2.0_dp*real(ampLO_ggh2zz(i1,i2,i3,i5)*conjg(ampNLO_ggh2zz(i1,i2,i3,i5)))
!!--!!       endif
!!--!!! |ZZ|^2
!!--!!       ! LO amp-sq
!!--!!       resZZLO_mt(0) = resZZLO_mt(0) + abs(ewampLO_ggzz_all(0,i1,i2,i3,i5))**2
!!--!!       resZZLO_mt(1) = resZZLO_mt(1) + &
!!--!!             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,i1,i2,i3,i5)))
!!--!!       resZZLO_mt(2) = resZZLO_mt(2) + &
!!--!!             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5))) + &
!!--!!             abs(ewampLO_ggzz_all(1,i1,i2,i3,i5))**2
!!--!!       resZZLO_mt(3) = resZZLO_mt(3) + &
!!--!!            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
!!--!!            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5)))
!!--!!       resZZLO_mt(4) = resZZLO_mt(4) + &
!!--!!            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(4,i1,i2,i3,i5))) + &
!!--!!            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
!!--!!            abs(ewampLO_ggzz_all(2,i1,i2,i3,i5))**2
!!--!!
!!--!!       ! LO*NLO
!!--!!       resZZNLO_mt(0) = resZZNLO_mt(0) + &
!!--!!            2.0_dp*real( ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5)))
!!--!!       resZZNLO_mt(1) = resZZNLO_mt(1) + &
!!--!!            2.0_dp*real( & 
!!--!!            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5)) &
!!--!!            )
!!--!!       resZZNLO_mt(2) = resZZNLO_mt(2) + &
!!--!!            2.0_dp*real( &
!!--!!            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(2,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5))   &
!!--!!            )
!!--!!       resZZNLO_mt(3) = resZZNLO_mt(3) + &
!!--!!            2.0_dp*real( &
!!--!!            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(3,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(2,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(3,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5))   &
!!--!!            )
!!--!!       
!!--!!       resZZNLO_mt(4) = resZZNLO_mt(4) + &
!!--!!            2.0_dp*real( &
!!--!!            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(4,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(3,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(2,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(3,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(1,i1,i2,i3,i5)) + &
!!--!!            ewampLO_ggzz_all(4,i1,i2,i3,i5)*conjg(ewampNLO_ggzz_all(0,i1,i2,i3,i5))   &
!!--!!            )
!!--!!
!!--!!! | ZZ + HZZ |^2
!!--!!
!!--!!       ! LO amp-sq
!!--!!       resZZHLO_mt(0) = resZZHLO_mt(0) + abs(ampLO_zzh2zz(0,i1,i2,i3,i5))**2
!!--!!       resZZHLO_mt(1) = resZZHLO_mt(1) + &
!!--!!            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,i1,i2,i3,i5)))
!!--!!       resZZHLO_mt(2) = resZZHLO_mt(2) + &
!!--!!            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5))) + &
!!--!!            abs(ampLO_zzh2zz(1,i1,i2,i3,i5))**2
!!--!!       resZZHLO_mt(3) = resZZHLO_mt(3) + &
!!--!!            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
!!--!!            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5)))
!!--!!       resZZHLO_mt(4) = resZZHLO_mt(4) + &
!!--!!            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(4,i1,i2,i3,i5))) + &
!!--!!            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
!!--!!            abs(ampLO_zzh2zz(2,i1,i2,i3,i5))**2
!!--!!
!!--!!
!!--!!       ! LO*NLO
!!--!!       resZZHNLO_mt(0) = resZZHNLO_mt(0) + &
!!--!!            2.0_dp*real( ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5)))
!!--!!       resZZHNLO_mt(1) = resZZHNLO_mt(1) + &
!!--!!            2.0_dp*real( & 
!!--!!            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5)) &
!!--!!            )
!!--!!       resZZHNLO_mt(2) = resZZHNLO_mt(2) + &
!!--!!            2.0_dp*real( &
!!--!!            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(2,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5))   &
!!--!!            )
!!--!!       resZZHNLO_mt(3) = resZZHNLO_mt(3) + &
!!--!!            2.0_dp*real( &
!!--!!            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(3,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(2,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(3,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5))   &
!!--!!            )
!!--!!       
!!--!!       resZZHNLO_mt(4) = resZZHNLO_mt(4) + &
!!--!!            2.0_dp*real( &
!!--!!            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(4,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(3,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(2,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(3,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(1,i1,i2,i3,i5)) + &
!!--!!            ampLO_zzh2zz(4,i1,i2,i3,i5)*conjg(ampNLO_zzh2zz(0,i1,i2,i3,i5))   &
!!--!!            )
!!--!!
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!    enddo
!!--!!
!!--!!    ! add the various pieces
!!--!!    ! if-statements are redundant, just put in to be extra safe
!!--!!
!!--!!    do nheft=0,expord
!!--!!          resZZLO = resZZLO + resZZLO_mt(nheft)
!!--!!          resZZNLO = resZZNLO + resZZNLO_mt(nheft)
!!--!!       
!!--!!       if(HiggsExp) then
!!--!!          resHLO = resHLO + resHLO_mt(nheft)
!!--!!          resHNLO = resHNLO + resHNLO_mt(nheft)
!!--!!       endif
!!--!!          resZZHLO = resZZHLO + resZZHLO_mt(nheft)
!!--!!          resZZHNLO = resZZHNLO + resZZHNLO_mt(nheft)
!!--!!    enddo
!!--!!    
!!--!!    first=.false.
!!--!!
!!--!!! define res according to the 'contr' input
!!--!!    if (contr .eq. "bkgd") then                ! |gg->ZZ|^2
!!--!!       if (first) print *, "Doing background |gg->ZZ|^2"
!!--!!       resLO = resZZLO
!!--!!       resNLO = resZZNLO
!!--!!    elseif (contr .eq. "sigl") then            ! |gg->H->ZZ|^2
!!--!!       if (first) print *, "Doing signal |gg->H->ZZ|^2"
!!--!!       resLO = resHLO
!!--!!       resNLO = resHNLO
!!--!!    elseif (contr .eq. "full") then            ! |gg->H->ZZ  +  gg->ZZ|^2
!!--!!       if (first) print *, "Doing full  |gg->H->ZZ + gg->ZZ|^2"
!!--!!       resLO = resZZHLO
!!--!!       resNLO = resZZHNLO
!!--!!    elseif (contr .eq. "intf") then            ! |gg->H->ZZ + gg->ZZ|^2 - |gg->ZZ|^2 - |gg->H->ZZ|^2 = 2*Real( (gg->H->ZZ)*conj(gg->ZZ) )
!!--!!       resLO = resZZHLO - resZZLO - resHLO
!!--!!       resNLO = resZZHNLO - resZZNLO - resHNLO
!!--!!       if (first) print *, "Doing interference 2Real( (gg->H->ZZ)*conj(gg->ZZ) )"
!!--!!    endif
!!--!!
!!--!!    !-- fix pre-factor
!!--!!    resLO = resLO * avegg * colf 
!!--!!    resNLO = resNLO * avegg * colf 
!!--!!
!!--!!!    print *, "resLO",resLO
!!--!!!    print *, "resNLO",resNLO
!!--!!!    stop
!!--!!
!!--!!! poles, as in massless piece
!!--!!! might want to check on this order-by-order in 1/mt exp.
!!--!!    poles(-2) = -2*CA*resLO
!!--!!    poles(-1) = (-two * b0 - two * CA * log(musq/sprod(1,2)))*resLO
!!--!!!    print *, "DP/1l Catani",poles(-2)/resLO
!!--!!!    print *, "SP/1l Catani",poles(-1)/resLO
!!--!!!    stop
!!--!!
!!--!!    return
!!--!!  end subroutine resNLOV_ZZ_heft_swap


end module mod_ampsqNLOV_zz
