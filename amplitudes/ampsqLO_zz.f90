module mod_ampsqLO_zz
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_ampLO_ggzz
  use mod_ampLO_ggh2zz
  use mod_ampLO_ggzz_mass

  implicit none
  private

!  public :: resLO_ZZ_heft_smeft
  public :: resLO_ZZ_heft,resLO_sping_ZZ_heft

contains

! - SMEFT couplings
!  subroutine resLO_ZZ_heft_smeft(p,mass,res)
! full amp. sq. of LO gg->ZZ
! contr=sigl --> |gg -> H -> ZZ|^2
! contr=bkgd --> |gg -> ZZ|^2
! contr=full --> | gg->H->ZZ + gg->ZZ|^2
! contr=intf --> 2*Real((gg->H->ZZ)*conjg(gg->ZZ))
! gg->ZZ through massless and/or massive loops of mass 'mass'
! massive loop may be computed exactly or in 1/mt^2 expansion
!    implicit none
!    real(dp), intent(in) :: p(4,6),mass
!    real(dp), intent(out) :: res
!    real(dp)              :: resZZ,resH,resZZH
!    real(dp)              :: resZZ_mt(0:nheft_max),resH_mt(0:nheft_max),resZZH_mt(0:nheft_max)
!    complex(dp) :: ewampLO_ggzz_all(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!    complex(dp) :: ewampLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
!    complex(dp) :: ampLO_ggh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!    complex(dp) :: ampLO_zzh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
!    complex(dp) :: za(6,6),zb(6,6),res_mt(0:nheft_max),res_mless,res_int(0:nheft_max)
!    complex(dp15) :: zadp(6,6), zbdp(6,6),ewampLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1),ewampLO_ggzz_fullmass(-1:1,-1:1,-1:1,-1:1)
!    real(dp)    :: pdum(4,6),sprod(6,6), pMCFM(6,4),pTZ
!    real(dp15)  :: pMCFMdp(6,4),sproddp(6,6),massdp
!    integer     :: i1,i2,i3,i5,j,nheft,neps
!    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one (N_c^2-1) (10/04/2024)
!    logical, save :: first=.true.

!    res       = zero
!    resH      = zero
!    resZZ     = zero
!    resZZH    = zero
!    resH_mt   = zero
!    resZZ_mt  = zero
!    resZZH_mt = zero

! set up spinors
!    call awayfromzaxis(6,p,pdum)
!    call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)


!    ewampLO_ggzz_mless = czero
!    ewampLO_ggzz_fullmass = czero
!    ewampLO_ggzz_heft  = czero
!    ewampLO_ggzz_all  = czero
!    ampLO_ggh2zz  = czero
!    ampLO_zzh2zz  = czero


! calculate gg -> ZZ amps
!    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then

!       if (mlessloops) then
!          if (first) print *, "including massless loops"
!          ewampLO_ggzz_mless = ewampLO_ggzz(1,2,3,4,5,6,za,zb,sprod)
!       endif
!       if (massloops) then
!          if (first .and. VVExp) print *, "including massive loops in HEFT"
!          if (first .and. .not. VVExp) print *, "including massive loops exactly"


!--- annoyingly, we need to initialize all scalar integrals for the exact mass dependence here
!          call converttoMCFMmom(p,pMCFM)
!          if (kind(p) .ne. 8) then
!             print *, "converting to dp for massive amps", kind(p)
!             pMCFMdp(1:6,1:4) = real(pMCFM(1:6,1:4),kind=dp15)
!             massdp = real(mass,kind=dp15)
!             zadp(1:6,1:6) = cmplx(za(1:6,1:6),kind=dp15)
!             zbdp(1:6,1:6) = cmplx(zb(1:6,1:6),kind=dp15)
!             sproddp(1:6,1:6)=real(sprod(1:6,1:6),kind=dp15)
!          else
!             zadp=za
!             zbdp=zb
!             sproddp=sprod
!             pMCFMdp=pMCFM
!             massdp=mass
!          endif
!          call ZZintegraleval(pMCFMdp,massdp)
!          call getewampsLO_ggzz_heft(zadp,zbdp,sproddp,massdp,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft)


! if pTZ < 0.1 GeV, omit top loos
!          pTZ = sqrt( (p(2,3)+p(2,4))**2 + (p(3,3)+p(3,4))**2 )
!          if ( pTZ .le. 0.1d0) then
!             print *, "warning, massive loop set to zero since pTZ = ",pTZ
!             ewampLO_ggzz_fullmass = czero
!          endif
!       endif

!       if (VVExp) then
!          ewampLO_ggzz_all = ewampLO_ggzz_heft
!          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_all(0,:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
!       else
!          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_fullmass(:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
!       endif
!    endif

! calculate gg->H->ZZ amplitudes. Here we call the smeft amplitude
!    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
!       if (HiggsExp) then
!          if (first) print *, "Using heavy top expansion for Higgs amplitudes"
!          ampLO_ggh2zz =  ewampLO_ggh2zz_heft(1,2,3,4,5,6,za,zb,sprod)
!       else                           ! Higgs with full mass dependence
!          if (first) print *, "Using full mass dependence for Higgs amplitudes"
!          ampLO_ggh2zz(0,:,:,:,:) =  ewampLO_ggh2zz_smeft(1,2,3,4,5,6,za,zb,sprod)
!       endif
!    endif

! add the Higgs and ZZ amplitudes
!    if (contr .eq. "full" .or. contr .eq. "intf") then
!       ampLO_zzh2zz = ewampLO_ggzz_all + ampLO_ggh2zz
!    endif



! square the amps, keeping the order of the 1/mt expansion
!    do i1=-1,1,2
!    do i2=-1,1,2
!    do i3=-1,1,2
!    do i5=-1,1,2

! |H| ^2
!       if (HiggsExp) then
!          resH_mt(0) = resH_mt(0) + abs(ampLO_ggh2zz(0,i1,i2,i3,i5))**2
!          resH_mt(1) = resH_mt(1) + &
!               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(1,i1,i2,i3,i5)))
!          resH_mt(2) = resH_mt(2) + &
!               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(2,i1,i2,i3,i5))) + &
!               abs(ampLO_ggh2zz(1,i1,i2,i3,i5))**2
!          resH_mt(3) = resH_mt(3) + &
!               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(3,i1,i2,i3,i5))) + &
!               2.0_dp*real(ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(2,i1,i2,i3,i5)))
!          resH_mt(4) = resH_mt(4) + &
!               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(4,i1,i2,i3,i5))) + &
!               2.0_dp*real(ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(3,i1,i2,i3,i5))) + &
!               abs(ampLO_ggh2zz(2,i1,i2,i3,i5))**2
!       else
!          resH = resH + abs(ampLO_ggh2zz(0,i1,i2,i3,i5))**2
!       endif

! |ZZ| ^2
!       resZZ_mt(0) = resZZ_mt(0) + abs(ewampLO_ggzz_all(0,i1,i2,i3,i5))**2
!       resZZ_mt(1) = resZZ_mt(1) + &
!             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,i1,i2,i3,i5)))
!       resZZ_mt(2) = resZZ_mt(2) + &
!             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5))) + &
!             abs(ewampLO_ggzz_all(1,i1,i2,i3,i5))**2
!       resZZ_mt(3) = resZZ_mt(3) + &
!            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
!            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5)))
!       resZZ_mt(4) = resZZ_mt(4) + &
!            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(4,i1,i2,i3,i5))) + &
!            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
!            abs(ewampLO_ggzz_all(2,i1,i2,i3,i5))**2

! | ZZ + HZZ |^2
!       resZZH_mt(0) = resZZH_mt(0) + abs(ampLO_zzh2zz(0,i1,i2,i3,i5))**2
!       resZZH_mt(1) = resZZH_mt(1) + &
!             2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,i1,i2,i3,i5)))
!       resZZH_mt(2) = resZZH_mt(2) + &
!             2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5))) + &
!             abs(ampLO_zzh2zz(1,i1,i2,i3,i5))**2
!       resZZH_mt(3) = resZZH_mt(3) + &
!            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
!            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5)))
!       resZZH_mt(4) = resZZH_mt(4) + &
!            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(4,i1,i2,i3,i5))) + &
!            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
!            abs(ampLO_zzh2zz(2,i1,i2,i3,i5))**2

!    enddo
!    enddo
!    enddo
!    enddo

!    do nheft=0,expord
!       resZZ = resZZ + resZZ_mt(nheft)
!       if (HiggsExp) then
!          resH = resH + resH_mt(nheft)
!       endif
!       resZZH = resZZH + resZZH_mt(nheft)
!    enddo

! define res according to the 'contr' input
!    if (contr .eq. "bkgd") then                ! |gg->ZZ|^2
!       if (first) print *, "Doing background |gg->ZZ|^2"
!          res = resZZ
!    elseif (contr .eq. "sigl") then            ! |gg->H->ZZ|^2
!       if (first) print *, "Doing signal |gg->H->ZZ|^2"
!       res = resH
!    elseif (contr .eq. "full") then            ! |gg->H->ZZ  +  gg->ZZ|^2
!       if (first) print *, "Doing full  |gg->H->ZZ + gg->ZZ|^2"
!       res = resZZH
!    elseif (contr .eq. "intf") then            ! |gg->H->ZZ + gg->ZZ|^2 - |gg->ZZ|^2 - |gg->H->ZZ|^2 = 2*Real( (gg->H->ZZ)*conj(gg->ZZ) )
!       res = resZZH - resZZ - resH
!       if (first) print *, "Doing interference 2Real( (gg->H->ZZ)*conj(gg->ZZ) )"
!    endif

!    first=.false.

    !-- fix pre-factor avegg is an average on final states, numerical stuff (10/04/2024)
!    res = res * colf * avegg

!    return
!  end subroutine resLO_ZZ_heft_smeft


  subroutine resLO_ZZ_heft(p,mass,res)
! full amp. sq. of LO gg->ZZ
! contr=sigl --> |gg -> H -> ZZ|^2
! contr=bkgd --> |gg -> ZZ|^2
! contr=full --> | gg->H->ZZ + gg->ZZ|^2
! contr=intf --> 2*Real((gg->H->ZZ)*conjg(gg->ZZ))
! gg->ZZ through massless and/or massive loops of mass 'mass'
! massive loop may be computed exactly or in 1/mt^2 expansion
    implicit none
    real(dp), intent(in) :: p(4,6),mass
    real(dp), intent(out) :: res
    real(dp)              :: resZZ,resH,resZZH
    real(dp)              :: resZZ_mt(0:nheft_max),resH_mt(0:nheft_max),resZZH_mt(0:nheft_max)
    complex(dp) :: ewampLO_ggzz_all(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_ggh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_zzh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: za(6,6),zb(6,6),res_mt(0:nheft_max),res_mless,res_int(0:nheft_max)
    complex(dp15) :: zadp(6,6), zbdp(6,6),ewampLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1),ewampLO_ggzz_fullmass(-1:1,-1:1,-1:1,-1:1)
    real(dp)    :: pdum(4,6),sprod(6,6), pMCFM(6,4),pTZ
    real(dp15)  :: pMCFMdp(6,4),sproddp(6,6),massdp
    integer     :: i1,i2,i3,i5,j,nheft,neps
    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one (N_c^2-1) (10/04/2024)
    logical, save :: first=.true.

    res       = zero
    resH      = zero
    resZZ     = zero
    resZZH    = zero
    resH_mt   = zero
    resZZ_mt  = zero
    resZZH_mt = zero

! set up spinors
    call awayfromzaxis(6,p,pdum)
    call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)


    ewampLO_ggzz_mless = czero
    ewampLO_ggzz_fullmass = czero
    ewampLO_ggzz_heft  = czero
    ewampLO_ggzz_all  = czero
    ampLO_ggh2zz  = czero
    ampLO_zzh2zz  = czero


! calculate gg -> ZZ amps
    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then

       if (mlessloops) then
          if (first) print *, "including massless loops"
          ewampLO_ggzz_mless = ewampLO_ggzz(1,2,3,4,5,6,za,zb,sprod)
       endif
       if (massloops) then
          if (first .and. VVExp) print *, "including massive loops in HEFT"
          if (first .and. .not. VVExp) print *, "including massive loops exactly"


!--- annoyingly, we need to initialize all scalar integrals for the exact mass dependence here
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
          call getewampsLO_ggzz_heft(zadp,zbdp,sproddp,massdp,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft)


! if pTZ < 0.1 GeV, omit top loos
          pTZ = sqrt( (p(2,3)+p(2,4))**2 + (p(3,3)+p(3,4))**2 )
          if ( pTZ .le. 0.1d0) then
             print *, "warning, massive loop set to zero since pTZ = ",pTZ
             ewampLO_ggzz_fullmass = czero
          endif
       endif

       if (VVExp) then
          ewampLO_ggzz_all = ewampLO_ggzz_heft
          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_all(0,:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
       else
          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_fullmass(:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
       endif
    endif

! calculate gg->H->ZZ amplitudes
    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if (HiggsExp) then
          if (first) print *, "Using heavy top expansion for Higgs amplitudes"
          ampLO_ggh2zz =  ewampLO_ggh2zz_heft(1,2,3,4,5,6,za,zb,sprod)
       else                           ! Higgs with full mass dependence
          if (first) print *, "Using full mass dependence for Higgs amplitudes"
          ampLO_ggh2zz(0,:,:,:,:) =  ewampLO_ggh2zz(1,2,3,4,5,6,za,zb,sprod)
       endif
    endif

! add the Higgs and ZZ amplitudes
    if (contr .eq. "full" .or. contr .eq. "intf") then
       ampLO_zzh2zz = ewampLO_ggzz_all + ampLO_ggh2zz
    endif



! square the amps, keeping the order of the 1/mt expansion
    do i1=-1,1,2
    do i2=-1,1,2
    do i3=-1,1,2
    do i5=-1,1,2

! |H| ^2
       if (HiggsExp) then
          resH_mt(0) = resH_mt(0) + abs(ampLO_ggh2zz(0,i1,i2,i3,i5))**2
          resH_mt(1) = resH_mt(1) + &
               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(1,i1,i2,i3,i5)))
          resH_mt(2) = resH_mt(2) + &
               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(2,i1,i2,i3,i5))) + &
               abs(ampLO_ggh2zz(1,i1,i2,i3,i5))**2
          resH_mt(3) = resH_mt(3) + &
               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(3,i1,i2,i3,i5))) + &
               2.0_dp*real(ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(2,i1,i2,i3,i5)))
          resH_mt(4) = resH_mt(4) + &
               2.0_dp*real(ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(4,i1,i2,i3,i5))) + &
               2.0_dp*real(ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(3,i1,i2,i3,i5))) + &
               abs(ampLO_ggh2zz(2,i1,i2,i3,i5))**2
       else
          resH = resH + abs(ampLO_ggh2zz(0,i1,i2,i3,i5))**2
       endif

! |ZZ| ^2
       resZZ_mt(0) = resZZ_mt(0) + abs(ewampLO_ggzz_all(0,i1,i2,i3,i5))**2
       resZZ_mt(1) = resZZ_mt(1) + &
             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,i1,i2,i3,i5)))
       resZZ_mt(2) = resZZ_mt(2) + &
             2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5))) + &
             abs(ewampLO_ggzz_all(1,i1,i2,i3,i5))**2
       resZZ_mt(3) = resZZ_mt(3) + &
            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,i1,i2,i3,i5)))
       resZZ_mt(4) = resZZ_mt(4) + &
            2.0_dp*real(ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(4,i1,i2,i3,i5))) + &
            2.0_dp*real(ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,i1,i2,i3,i5))) + &
            abs(ewampLO_ggzz_all(2,i1,i2,i3,i5))**2

! | ZZ + HZZ |^2
       resZZH_mt(0) = resZZH_mt(0) + abs(ampLO_zzh2zz(0,i1,i2,i3,i5))**2
       resZZH_mt(1) = resZZH_mt(1) + &
             2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,i1,i2,i3,i5)))
       resZZH_mt(2) = resZZH_mt(2) + &
             2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5))) + &
             abs(ampLO_zzh2zz(1,i1,i2,i3,i5))**2
       resZZH_mt(3) = resZZH_mt(3) + &
            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,i1,i2,i3,i5)))
       resZZH_mt(4) = resZZH_mt(4) + &
            2.0_dp*real(ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(4,i1,i2,i3,i5))) + &
            2.0_dp*real(ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,i1,i2,i3,i5))) + &
            abs(ampLO_zzh2zz(2,i1,i2,i3,i5))**2

    enddo
    enddo
    enddo
    enddo

    do nheft=0,expord
       resZZ = resZZ + resZZ_mt(nheft)
       if (HiggsExp) then
          resH = resH + resH_mt(nheft)
       endif
       resZZH = resZZH + resZZH_mt(nheft)
    enddo

! define res according to the 'contr' input
    if (contr .eq. "bkgd") then                ! |gg->ZZ|^2
       if (first) print *, "Doing background |gg->ZZ|^2"
          res = resZZ
    elseif (contr .eq. "sigl") then            ! |gg->H->ZZ|^2
       if (first) print *, "Doing signal |gg->H->ZZ|^2"
       res = resH
    elseif (contr .eq. "full") then            ! |gg->H->ZZ  +  gg->ZZ|^2
       if (first) print *, "Doing full  |gg->H->ZZ + gg->ZZ|^2"
       res = resZZH
    elseif (contr .eq. "intf") then            ! |gg->H->ZZ + gg->ZZ|^2 - |gg->ZZ|^2 - |gg->H->ZZ|^2 = 2*Real( (gg->H->ZZ)*conj(gg->ZZ) )
       res = resZZH - resZZ - resH
       if (first) print *, "Doing interference 2Real( (gg->H->ZZ)*conj(gg->ZZ) )"
    endif

    first=.false.

    !-- fix pre-factor avegg is an average on final states, numerical stuff (10/04/2024)
    res = res * colf * avegg

    return
  end subroutine resLO_ZZ_heft


  subroutine resLO_sping_ZZ_heft(p,mass,res)
! spin correlations
    use mod_ampLO_ggzz
    use mod_ampLO_ggh2zz
    implicit none
    real(dp), intent(in) :: p(4,6),mass
    complex(dp), intent(out) :: res(-1:1,-1:1)
    real(dp)              :: resZZ(-1:1,-1:1),resH(-1:1,-1:1),resZZH(-1:1,-1:1)
    real(dp)              :: resZZ_mt(0:nheft_max,-1:1,-1:1),resH_mt(0:nheft_max,-1:1,-1:1),resZZH_mt(0:nheft_max,-1:1,-1:1)
    complex(dp) :: ewampLO_ggzz_all(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampLO_ggzz_mless(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_ggh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_zzh2zz(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: za(6,6),zb(6,6),res_mt(0:nheft_max,-1:1,-1:1),res_mless,res_int(0:nheft_max)
    complex(dp15) :: zadp(6,6), zbdp(6,6),ewampLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1),ewampLO_ggzz_fullmass(-1:1,-1:1,-1:1,-1:1)
    real(dp)    :: pdum(4,6),sprod(6,6),pMCFM(6,4),pTZ
    real(dp15)  :: pMCFMdp(6,4),sproddp(6,6),massdp
    integer     :: i1,i2,i3,i5,j,nheft,neps,j1
    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one


    res       = zero
    resH      = zero
    resZZ     = zero
    resZZH    = zero
    resH_mt   = zero
    resZZ_mt  = zero
    resZZH_mt = zero

    ewampLO_ggzz_mless    = czero
    ewampLO_ggzz_fullmass = czero
    ewampLO_ggzz_heft     = czero
    ewampLO_ggzz_all      = czero
    ampLO_ggh2zz          = czero
    ampLO_zzh2zz          = czero

! set up spinors
    call awayfromzaxis(6,p,pdum)
    call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)


! calculate gg -> ZZ amps
    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then

       if (mlessloops) then
          ewampLO_ggzz_mless = ewampLO_ggzz(1,2,3,4,5,6,za,zb,sprod)
       endif

       if (massloops) then
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
          call getewampsLO_ggzz_heft(zadp,zbdp,sproddp,massdp,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft)

! if pTZ < 0.1 GeV, omit top loops
          pTZ = sqrt( (p(2,3)+p(2,4))**2 + (p(3,3)+p(3,4))**2 )
          if ( pTZ .le. 0.1d0) then
             print *, "warning, massive loop set to zero since pTZ = ",pTZ
             ewampLO_ggzz_fullmass = czero
          endif
       endif

       if (VVExp) then
          ewampLO_ggzz_all = ewampLO_ggzz_heft
          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_all(0,:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
       else
          ewampLO_ggzz_all(0,:,:,:,:) = ewampLO_ggzz_fullmass(:,:,:,:)+ewampLO_ggzz_mless(:,:,:,:)
       endif
    endif

! calculate gg->H->ZZ amplitudes
    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if (HiggsExp) then      ! Higgs expanded
          ampLO_ggh2zz =  ewampLO_ggh2zz_heft(1,2,3,4,5,6,za,zb,sprod) ! Higgs expanded
       else                           ! Higgs with full mass dependence
          ampLO_ggh2zz(0,:,:,:,:) =  ewampLO_ggh2zz(1,2,3,4,5,6,za,zb,sprod)
       endif
    endif

! add the Higgs and ZZ amplitudes
    if (contr .eq. "full" .or. contr .eq. "intf") then
       ampLO_zzh2zz = ewampLO_ggzz_all + ampLO_ggh2zz
    endif

! square the amps, keeping the order of the 1/mt expansion
!M the factors of 1/mt already exist in these amplitudes
    do i1=-1,1,2
    do j1=-1,1,2
    do i2=-1,1,2
    do i3=-1,1,2
    do i5=-1,1,2

! |H| ^2
       if (HiggsExp) then
          resH_mt(0,i1,j1) = resH_mt(0,i1,j1) + ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(0,j1,i2,i3,i5))
          resH_mt(1,i1,j1) = resH_mt(1,i1,j1) + &
               ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(1,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(0,j1,i2,i3,i5))
          resH_mt(2,i1,j1) = resH_mt(2,i1,j1) + &
               ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(2,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(1,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(2,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(0,j1,i2,i3,i5))
          resH_mt(3,i1,j1) = resH_mt(3,i1,j1) + &
               ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(3,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(2,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(2,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(1,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(3,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(0,j1,i2,i3,i5))
          resH_mt(4,i1,j1) = resH_mt(4,i1,j1) + &
               ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(4,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(1,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(3,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(2,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(2,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(3,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(1,j1,i2,i3,i5)) + &
               ampLO_ggh2zz(4,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(0,j1,i2,i3,i5))
       else
          resH(i1,j1) = resH(i1,j1) + ampLO_ggh2zz(0,i1,i2,i3,i5)*conjg(ampLO_ggh2zz(0,j1,i2,i3,i5))
       endif

! |ZZ| ^2
       resZZ_mt(0,i1,j1) = resZZ_mt(0,i1,j1) + ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(0,j1,i2,i3,i5))
       resZZ_mt(1,i1,j1) = resZZ_mt(1,i1,j1) + &
             ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,j1,i2,i3,i5))+&
             ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(0,j1,i2,i3,i5))
       resZZ_mt(2,i1,j1) = resZZ_mt(2,i1,j1) + &
             ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,j1,i2,i3,i5)) + &
             ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,j1,i2,i3,i5)) + &
             ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(0,j1,i2,i3,i5))
       resZZ_mt(3,i1,j1) = resZZ_mt(3,i1,j1) + &
            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,j1,i2,i3,i5)) + &
            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,j1,i2,i3,i5)) + &
            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,j1,i2,i3,i5)) + &
            ewampLO_ggzz_all(3,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(0,j1,i2,i3,i5))
       resZZ_mt(4,i1,j1) = resZZ_mt(4,i1,j1) + &
            ewampLO_ggzz_all(0,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(4,j1,i2,i3,i5)) + &
            ewampLO_ggzz_all(1,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(3,j1,i2,i3,i5)) + &
            ewampLO_ggzz_all(2,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(2,j1,i2,i3,i5)) + &
            ewampLO_ggzz_all(3,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(1,j1,i2,i3,i5)) + &
            ewampLO_ggzz_all(4,i1,i2,i3,i5)*conjg(ewampLO_ggzz_all(0,j1,i2,i3,i5))

! | ZZ + HZZ |^2
       resZZH_mt(0,i1,j1) = resZZH_mt(0,i1,j1) + ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(0,j1,i2,i3,i5))
       resZZH_mt(1,i1,j1) = resZZH_mt(1,i1,j1) + &
             ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,j1,i2,i3,i5))+&
             ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(0,j1,i2,i3,i5))
       resZZH_mt(2,i1,j1) = resZZH_mt(2,i1,j1) + &
             ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,j1,i2,i3,i5)) + &
             ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,j1,i2,i3,i5)) + &
             ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(0,j1,i2,i3,i5))
       resZZH_mt(3,i1,j1) = resZZH_mt(3,i1,j1) + &
            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,j1,i2,i3,i5)) + &
            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,j1,i2,i3,i5)) + &
            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,j1,i2,i3,i5)) + &
            ampLO_zzh2zz(3,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(0,j1,i2,i3,i5))
       resZZH_mt(4,i1,j1) = resZZH_mt(4,i1,j1) + &
            ampLO_zzh2zz(0,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(4,j1,i2,i3,i5)) + &
            ampLO_zzh2zz(1,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(3,j1,i2,i3,i5)) + &
            ampLO_zzh2zz(2,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(2,j1,i2,i3,i5)) + &
            ampLO_zzh2zz(3,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(1,j1,i2,i3,i5)) + &
            ampLO_zzh2zz(4,i1,i2,i3,i5)*conjg(ampLO_zzh2zz(0,j1,i2,i3,i5))

    enddo
    enddo
    enddo
    enddo
    enddo

    do nheft=0,expord
       resZZ(-1:1,-1:1) = resZZ(-1:1,-1:1) + resZZ_mt(nheft,-1:1,-1:1)
       if (HiggsExp) then
          resH(-1:1,-1:1) = resH(-1:1,-1:1) + resH_mt(nheft,-1:1,-1:1)
       endif
       resZZH(-1:1,-1:1) = resZZH(-1:1,-1:1) + resZZH_mt(nheft,-1:1,-1:1)
    enddo
! define res according to the 'contr' input
    if (contr .eq. "bkgd") then                ! |gg->ZZ|^2
       res(-1:1,-1:1) = resZZ(-1:1,-1:1)
    elseif (contr .eq. "sigl") then            ! |gg->H->ZZ|^2
       res(-1:1,-1:1) = resH(-1:1,-1:1)
    elseif (contr .eq. "full") then            ! |gg->H->ZZ  +  gg->ZZ|^2
       res(-1:1,-1:1) = resZZH(-1:1,-1:1)
    elseif (contr .eq. "intf") then            ! |gg->H->ZZ + gg->ZZ|^2 - |gg->ZZ|^2 - |gg->H->ZZ|^2 = 2*Real( (gg->H->ZZ)*conj(gg->ZZ) )
       res(-1:1,-1:1) = resZZH(-1:1,-1:1) - resZZ(-1:1,-1:1) - resH(-1:1,-1:1)
    endif

    !-- fix pre-factor
    res = res * colf * avegg

    return
  end subroutine resLO_sping_ZZ_heft

end module mod_ampsqLO_zz
