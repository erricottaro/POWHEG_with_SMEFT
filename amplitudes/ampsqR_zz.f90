module mod_ampsqR_zz
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_ampR_gggh2zz
!  use mod_ampR_gggZZ
!  use mod_ampR_gggZZ_mass
!  use mod_amplgggzz_sr

  implicit none 

  private 

  public :: resR_ZZ_heft

contains 

  ! ordering: g g -> g Z[l lb] Z [l lb]
  subroutine resR_ZZ_heft(p,mass,resR)
! full amp. sq. of real emission gg->ZZ+g
! contr=sigl --> |gg -> H(-> ZZ)+g|^2
! contr=bkgd --> |gg -> ZZ+g|^2
! contr=full --> | gg->H(->ZZ)+g + gg->ZZ+g|^2
! contr=intf --> 2*Real((gg->H(->ZZ)+g)*conjg(gg->ZZ+g))
! gg->ZZ+g through massless and/or massive loops of mass 'mass'
! massive loop in 1/mt^2 expansion only
    implicit none
    real(dp), intent(in) :: p(4,7),mass
    real(dp), intent(out) :: resR(-2:0)   ! previously was resR(-2:0), hope it does not break anything
    real(dp) :: resZZ_fsq_mt(0:nheft_max,-2:0),resZZ_dsq_mt(0:nheft_max,-2:0),&
         resZZH_fsq_mt(0:nheft_max,-2:0),resZZH_dsq_mt(0:nheft_max,-2:0)
    real(dp) :: resH_mt(0:nheft_max,-2:0),resZZ_mt(0:nheft_max,-2:0),resZZH_mt(0:nheft_max,-2:0)
    real(dp) :: resZZ(-2:0),resZZH(-2:0),resH(-2:0)
    complex(dp) :: za(7,7),zb(7,7),resR_mt(0:nheft_max)
    complex(dp15) :: amp_gggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1), amp_gggzz_axvec(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewampSRV_mless(-1:1,-1:1,-1:1,-1:1,-1:1), ewampSRA_mless(-1:1,-1:1,-1:1,-1:1,-1:1),ewampSRV_mass(-1:1,-1:1,-1:1,-1:1,-1:1), ewampSRA_mass(-1:1,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: amp123(-1:1,-1:1,-1:1,-1:1,-1:1,-2:0),amp132(-1:1,-1:1,-1:1,-1:1,-1:1,-2:0),&
         amp_mless_fabc(-1:1,-1:1,-1:1,-1:1,-1:1,-2:0),amp_mless_dabc(-1:1,-1:1,-1:1,-1:1,-1:1,-2:0)
    complex(dp) :: amp_gggh2zz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1),amp_gggh2zz(-1:1,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: amp_gggzz_all_dabc(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1,-2:0),&
         amp_gggzz_all_fabc(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1,-2:0),amp_zzh2zz_heft_fabc(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1,-2:0),&
         amp_zzh2zz_heft_dabc(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1,-2:0)
    real(dp)    :: pdum(4,7),sprod(7,7)
    real(dp15)  :: sproddp(7,7),massdp
    complex(dp15) :: zadp(7,7), zbdp(7,7)
    integer     :: i1,i2,i3,i4,i6,j,nheft
    ! we have the result in terms of f(a,b,c)
    real(dp), parameter :: colf = 24.0_dp          ! f^{abc} f^{abc} = Ca (Ca^2 -1)
    real(dp), parameter :: cold = 40.0_dp/3.0_dp   ! d^{abc} d^{abc} = 2 Cf (

    logical, save :: first=.true.

    resZZ_fsq_mt  = zero
    resZZ_dsq_mt  = zero
    resZZH_fsq_mt = zero
    resZZH_dsq_mt = zero
    resH_mt       = zero
    resZZ_mt      = zero
    resZZH_mt     = zero
    resR          = zero
    resZZ         = zero
    resH          = zero
    resZZH        = zero

    amp_gggzz_heft = czero
    amp_gggzz_axvec = czero
    amp123 = czero
    amp132 = czero
    amp_gggh2zz_heft = czero
    amp_gggh2zz = czero
    amp_gggzz_all_dabc = czero
    amp_gggzz_all_fabc = czero
    amp_zzh2zz_heft_fabc = czero
    amp_zzh2zz_heft_dabc = czero

! set up spinors
    call awayfromzaxis(7,p,pdum)
    call spinorur(7,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6),pdum(:,7)/),za,zb,sprod)


! calculate gg -> ZZ amps (eliminate these cases since we are not able to compute background amplitude)
    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then
       print *, "ERROR: SMEFT parameters can be used only to compute signal"
       stop
!       if (massloops) then
!          if (first) print *, "including massive loops in HEFT"
!          if (kind(p) .ne. 8) then
!             massdp = real(mass,kind=dp15)
!             zadp(1:7,1:7) = cmplx(za(1:7,1:7),kind=dp15)
!             zbdp(1:7,1:7) = cmplx(zb(1:7,1:7),kind=dp15)
!             sproddp(1:7,1:7)=real(sprod(1:7,1:7),kind=dp15)
!          else
!             zadp=za
!             zbdp=zb
!             sproddp=sprod
!             massdp=mass
!          endif       
!          call getewampR_gggzz_heft(zadp,zbdp,sproddp,massdp,amp_gggzz_heft,amp_gggzz_axvec)
! add to finite part of ZZ ampl
 !         amp_gggzz_all_fabc(0:nheft_max,:,:,:,:,:,0) = amp_gggzz_all_fabc(0:nheft_max,:,:,:,:,:,0) +  amp_gggzz_heft(0:nheft_max,:,:,:,:,:)
 !         amp_gggzz_all_dabc(0:nheft_max,:,:,:,:,:,0) = amp_gggzz_all_dabc(0:nheft_max,:,:,:,:,:,0) + amp_gggzz_axvec(0:nheft_max,:,:,:,:,:)
! include to make phases of massive and massless loops agree
 !         amp_gggzz_all_fabc = -amp_gggzz_all_fabc
 !         amp_gggzz_all_dabc = -amp_gggzz_all_dabc


! single resonant -- include color factor of 2*sqrt2
 !         call  ewampR_gggzz_srm(1,2,3,za,zb,sprod,expmass,ewampSRV_mass,ewampSRA_mass)
 !         amp_gggzz_all_fabc(0,:,:,:,:,:,0) = amp_gggzz_all_fabc(0,:,:,:,:,:,0) +  ewampSRA_mass(:,:,:,:,:)/sqrt2
 !         amp_gggzz_all_dabc(0,:,:,:,:,:,0) = amp_gggzz_all_dabc(0,:,:,:,:,:,0) +  ewampSRV_mass(:,:,:,:,:)/sqrt2
 !      endif
       
 !      if (mlessloops) then 
 !         if (first) print *, "including massless loops"
 !         call ewampR_gggzz_dr(pdum,1,2,3,amp123,amp132)
! in terms of f^{abc} and d^{abc}, with these written in terms of traces with normalization Tr(t^a t^b) = delta^{ab}
 !         amp_mless_fabc=(amp123-amp132)*ci/two
 !         amp_mless_dabc=(amp123+amp132)/two
! f^{abc}^2 (with 1 norm) = 2 * f^{abc}^2 (with 1/2 norm) (and similarly for dabc^2) -- this returns us to the 1/2 norm
! add to 0th term in 1/mt expanstion of ZZ ampl
 !         amp_mless_fabc = amp_mless_fabc * sqrt(two)
 !         amp_mless_dabc = amp_mless_dabc * sqrt(two)

 !         amp_gggzz_all_fabc(0,:,:,:,:,:,-2:0) = amp_gggzz_all_fabc(0,:,:,:,:,:,-2:0)  + amp_mless_fabc(:,:,:,:,:,-2:0)
 !         amp_gggzz_all_dabc(0,:,:,:,:,:,-2:0) = amp_gggzz_all_dabc(0,:,:,:,:,:,-2:0)  + amp_mless_dabc(:,:,:,:,:,-2:0)
          
! single resonant -- include color factor of 2*sqrt2
 !         call  ewampR_gggzz_sr(1,2,3,za,zb,sprod,ewampSRV_mless,ewampSRA_mless)
 !         amp_gggzz_all_fabc(0,:,:,:,:,:,0) = amp_gggzz_all_fabc(0,:,:,:,:,:,0) +  ewampSRA_mless(:,:,:,:,:)/sqrt2
 !         amp_gggzz_all_dabc(0,:,:,:,:,:,0) = amp_gggzz_all_dabc(0,:,:,:,:,:,0) +  ewampSRV_mless(:,:,:,:,:)/sqrt2
 !      endif
    endif

! calculate gg->H->ZZ amplitudes
    if(contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if(HiggsExp) then
          if (first) print *, "Using heavy top expansion for Higgs amplitudes"
          amp_gggh2zz_heft =  ewampR_gggh2zz_heft(1,2,3,4,5,6,7,za,zb,sprod)
          ! correct color factor 
          amp_gggh2zz_heft = sqrt(two)*amp_gggh2zz_heft 

          ! the sum of both Higgs and ZZ amps -- only for fabc part
          amp_zzh2zz_heft_fabc = amp_gggzz_all_fabc
          amp_zzh2zz_heft_dabc = amp_gggzz_all_dabc
          amp_zzh2zz_heft_fabc(:,:,:,:,:,:,0) = amp_zzh2zz_heft_fabc(:,:,:,:,:,:,0) + amp_gggh2zz_heft(:,:,:,:,:,:)

       else                           ! Higgs with full mass dependence
          if (first) print *, "Using full mass dependence for Higgs amplitudes"
          amp_gggh2zz =  ewampR_gggh2zz(1,2,7,3,4,5,6,za,zb,sprod)

          ! correct color factor 
          amp_gggh2zz = sqrt(two)*amp_gggh2zz 

          ! the sum of both Higgs and ZZ amps -- only for fabc part
          amp_zzh2zz_heft_fabc = amp_gggzz_all_fabc
          amp_zzh2zz_heft_dabc = amp_gggzz_all_dabc
          amp_zzh2zz_heft_fabc(0,:,:,:,:,:,0) = amp_zzh2zz_heft_fabc(0,:,:,:,:,:,0) + amp_gggh2zz(:,:,:,:,:)
       endif
    endif     


! square the amps, keeping the order of the 1/mt expansion
! M the factors of 1/mt already exist in these amplitudes
    do i1=-1,1,2
    do i2=-1,1,2
    do i3=-1,1,2
    do i4=-1,1,2
    do i6=-1,1,2

! | H |^2 ! Unchanged but now includes ep^i index for full result
       if (HiggsExp) then
          resH_mt(0,0) = resH_mt(0,0) + abs(amp_gggh2zz_heft(0,i1,i2,i3,i4,i6))**2
          resH_mt(1,0) = resH_mt(1,0) + &
               2.0_dp*real(amp_gggh2zz_heft(0,i1,i2,i3,i4,i6)*conjg(amp_gggh2zz_heft(1,i1,i2,i3,i4,i6)))
          resH_mt(2,0) = resH_mt(2,0) + &
               2.0_dp*real(amp_gggh2zz_heft(0,i1,i2,i3,i4,i6)*conjg(amp_gggh2zz_heft(2,i1,i2,i3,i4,i6))) + &
               abs(amp_gggh2zz_heft(1,i1,i2,i3,i4,i6))**2
          resH_mt(3,0) = resH_mt(3,0) + &
               2.0_dp*real(amp_gggh2zz_heft(0,i1,i2,i3,i4,i6)*conjg(amp_gggh2zz_heft(3,i1,i2,i3,i4,i6))) + &
               2.0_dp*real(amp_gggh2zz_heft(1,i1,i2,i3,i4,i6)*conjg(amp_gggh2zz_heft(2,i1,i2,i3,i4,i6)))
          resH_mt(4,0) = resH_mt(4,0) + &
               2.0_dp*real(amp_gggh2zz_heft(0,i1,i2,i3,i4,i6)*conjg(amp_gggh2zz_heft(4,i1,i2,i3,i4,i6))) + &
               2.0_dp*real(amp_gggh2zz_heft(1,i1,i2,i3,i4,i6)*conjg(amp_gggh2zz_heft(3,i1,i2,i3,i4,i6))) + &
               abs(amp_gggh2zz_heft(2,i1,i2,i3,i4,i6))**2
       else
          resH_mt(0,0) = resH_mt(0,0) + abs(amp_gggh2zz(i1,i2,i3,i4,i6))**2
       endif
     
! |ZZ| ^2
       resZZ_fsq_mt(0,-2:0) = resZZ_fsq_mt(0,-2:0) + abs(amp_gggzz_all_fabc(0,i1,i2,i3,i4,i6,-2:0))**2
       resZZ_fsq_mt(1,-2:0) = resZZ_fsq_mt(1,-2:0) + &
             2.0_dp*real(amp_gggzz_all_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_fabc(1,i1,i2,i3,i4,i6,-2:0)))
       resZZ_fsq_mt(2,-2:0) = resZZ_fsq_mt(2,-2:0) + &
             2.0_dp*real(amp_gggzz_all_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_fabc(2,i1,i2,i3,i4,i6,-2:0))) + &
             abs(amp_gggzz_all_fabc(1,i1,i2,i3,i4,i6,-2:0))**2
       resZZ_fsq_mt(3,-2:0) = resZZ_fsq_mt(3,-2:0) + &
            2.0_dp*real(amp_gggzz_all_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_fabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_gggzz_all_fabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_fabc(2,i1,i2,i3,i4,i6,-2:0)))
       resZZ_fsq_mt(4,-2:0) = resZZ_fsq_mt(4,-2:0) + &
            2.0_dp*real(amp_gggzz_all_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_fabc(4,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_gggzz_all_fabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_fabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            abs(amp_gggzz_all_fabc(2,i1,i2,i3,i4,i6,-2:0))**2

       resZZ_dsq_mt(0,-2:0) = resZZ_dsq_mt(0,-2:0) + abs(amp_gggzz_all_dabc(0,i1,i2,i3,i4,i6,-2:0))**2
       resZZ_dsq_mt(1,-2:0) = resZZ_dsq_mt(1,-2:0) + &
             2.0_dp*real(amp_gggzz_all_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_dabc(1,i1,i2,i3,i4,i6,-2:0)))
       resZZ_dsq_mt(2,-2:0) = resZZ_dsq_mt(2,-2:0) + &
             2.0_dp*real(amp_gggzz_all_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_dabc(2,i1,i2,i3,i4,i6,-2:0))) + &
             abs(amp_gggzz_all_dabc(1,i1,i2,i3,i4,i6,-2:0))**2
       resZZ_dsq_mt(3,-2:0) = resZZ_dsq_mt(3,-2:0) + &
            2.0_dp*real(amp_gggzz_all_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_dabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_gggzz_all_dabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_dabc(2,i1,i2,i3,i4,i6,-2:0)))
       resZZ_dsq_mt(4,-2:0) = resZZ_dsq_mt(4,-2:0) + &
            2.0_dp*real(amp_gggzz_all_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_dabc(4,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_gggzz_all_dabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_gggzz_all_dabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            abs(amp_gggzz_all_dabc(2,i1,i2,i3,i4,i6,-2:0))**2


! | ZZ + HZZ | ^2
       resZZH_fsq_mt(0,-2:0) = resZZH_fsq_mt(0,-2:0) + abs(amp_zzh2zz_heft_fabc(0,i1,i2,i3,i4,i6,-2:0))**2
       resZZH_fsq_mt(1,-2:0) = resZZH_fsq_mt(1,-2:0) + &
             2.0_dp*real(amp_zzh2zz_heft_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_fabc(1,i1,i2,i3,i4,i6,-2:0)))
       resZZH_fsq_mt(2,-2:0) = resZZH_fsq_mt(2,-2:0) + &
             2.0_dp*real(amp_zzh2zz_heft_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_fabc(2,i1,i2,i3,i4,i6,-2:0))) + &
             abs(amp_zzh2zz_heft_fabc(1,i1,i2,i3,i4,i6,-2:0))**2
       resZZH_fsq_mt(3,-2:0) = resZZH_fsq_mt(3,-2:0) + &
            2.0_dp*real(amp_zzh2zz_heft_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_fabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_zzh2zz_heft_fabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_fabc(2,i1,i2,i3,i4,i6,-2:0)))
       resZZH_fsq_mt(4,-2:0) = resZZH_fsq_mt(4,-2:0) + &
            2.0_dp*real(amp_zzh2zz_heft_fabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_fabc(4,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_zzh2zz_heft_fabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_fabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            abs(amp_zzh2zz_heft_fabc(2,i1,i2,i3,i4,i6,-2:0))**2

       resZZH_dsq_mt(0,-2:0) = resZZH_dsq_mt(0,-2:0) + abs(amp_zzh2zz_heft_dabc(0,i1,i2,i3,i4,i6,-2:0))**2
       resZZH_dsq_mt(1,-2:0) = resZZH_dsq_mt(1,-2:0) + &
             2.0_dp*real(amp_zzh2zz_heft_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_dabc(1,i1,i2,i3,i4,i6,-2:0)))
       resZZH_dsq_mt(2,-2:0) = resZZH_dsq_mt(2,-2:0) + &
             2.0_dp*real(amp_zzh2zz_heft_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_dabc(2,i1,i2,i3,i4,i6,-2:0))) + &
             abs(amp_zzh2zz_heft_dabc(1,i1,i2,i3,i4,i6,-2:0))**2
       resZZH_dsq_mt(3,-2:0) = resZZH_dsq_mt(3,-2:0) + &
            2.0_dp*real(amp_zzh2zz_heft_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_dabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_zzh2zz_heft_dabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_dabc(2,i1,i2,i3,i4,i6,-2:0)))
       resZZH_dsq_mt(4,-2:0) = resZZH_dsq_mt(4,-2:0) + &
            2.0_dp*real(amp_zzh2zz_heft_dabc(0,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_dabc(4,i1,i2,i3,i4,i6,-2:0))) + &
            2.0_dp*real(amp_zzh2zz_heft_dabc(1,i1,i2,i3,i4,i6,-2:0)*conjg(amp_zzh2zz_heft_dabc(3,i1,i2,i3,i4,i6,-2:0))) + &
            abs(amp_zzh2zz_heft_dabc(2,i1,i2,i3,i4,i6,-2:0))**2

    enddo
    enddo
    enddo
    enddo
    enddo

! multiply by the color factors
    resH_mt   = colf * resH_mt
    resZZ_mt  = colf * resZZ_fsq_mt  + cold * resZZ_dsq_mt
    resZZH_mt = colf * resZZH_fsq_mt + cold * resZZH_dsq_mt

! combine all pieces
!    resZZ(-2:0) = resZZ(-2:0) + resZZ_mt(0,-2:0) + resZZ_mt(1,-2:0) + resZZ_mt(2,-2:0) + resZZ_mt(3,-2:0) + resZZ_mt(4,-2:0)
!    if (HiggsExp) then      
!       resH(0) = resH(0) + resH_mt(0,0) + resH_mt(1,0) + resH_mt(2,0) + resH_mt(3,0) + resH_mt(4,0)
!    else
!       resH(0) = resH_mt(0,0)   
!    endif
!    resZZH(-2:0) = resZZH(-2:0) + resZZH_mt(0,-2:0) + resZZH_mt(1,-2:0) + resZZH_mt(2,-2:0) + resZZH_mt(3,-2:0) + resZZH_mt(4,-2:0)

    do nheft=0,expord
       resZZ(-2:0) = resZZ(-2:0) + resZZ_mt(nheft,-2:0)
       resZZH(-2:0) = resZZH(-2:0) + resZZH_mt(nheft,-2:0)
       if (HiggsExp) then
          resH(0) = resH(0) + resH_mt(nheft,0)
       endif
    enddo
    if (.not. HiggsExp) resH(0) = resH_mt(0,0)

! defines res according to the 'contr' input
    if (contr .eq. "bkgd") then
       if (first) print*, "Doing background |gg->ZZ+g|^2"
       resR = resZZ
    elseif (contr .eq. "sigl") then
       if (first) print *, "Doing signal |gg->H+g->ZZ+g|^2"
       resR = resH
    elseif (contr .eq. "full") then
       if (first) print *,"Doing full | gg->H+g->ZZ+g + gg->ZZ+g |^2"
       resR = resZZH
    elseif (contr .eq. "intf") then
       if (first) print *, "Doing interference 2Real( (gg->H+g->ZZ+g)*conj(gg->ZZ+g) )"
       resR = resZZH - resZZ - resH
    endif

    first=.false.

    !-- fix pre-factor -- no color factor, this has already been included 
    resR = avegg * resR
!    print*, "Real amplitude squared", resR
!    stop

    return
  end subroutine resR_ZZ_heft

end module mod_ampsqR_zz
