module mod_ampsqNLOV_ww
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_ampLO_ggww
  use mod_ampNLO_ggww
  use mod_ampLO_ggh2ww
  use mod_ampNLO_ggh2ww

  implicit none
  private
  
  public :: resNLOV_WW
  
contains

  ! full amp. sq. of NLO gg->WW
  ! ordering: g g -> W[l lb] W[l lb]
  ! coefficient of (as/twopi)**3
  subroutine resNLOV_WW(p,E1LL,E1LR,E2LL,E2LR,resLO,resNLO,poles)
    real(dp), intent(in) :: p(4,6)
    complex(dp), intent(in) :: E1LL(9),E1LR(9),E2LL(9),E2LR(9)
    real(dp), intent(out) :: resLO,resNLO,poles(-2:-1)
    real(dp) :: resHLO,resHNLO,resWWLO,resWWNLO,resWWHLO,resWWHNLO
    complex(dp) :: ewamplo_ggww_all(-1:1,-1:1,-1:1,-1:1),ewampnlo_ggww_all(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewamplo_ggww_mless(-1:1,-1:1,-1:1,-1:1),ewampnlo_ggww_mless(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ewamplo_ggww_mass(-1:1,-1:1,-1:1,-1:1),ewampnlo_ggww_mass(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_wwh2ww(-1:1,-1:1,-1:1,-1:1),ampNLO_wwh2ww(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_ggh2ww(-1:1,-1:1,-1:1,-1:1),ampNLO_ggh2ww(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampLO_ggh2ww_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1),ampNLO_ggh2ww_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: pdum(4,6),sprod(6,6),pTV
    integer     :: i1,i2,i3,i5
    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one
    logical, save :: first=.true.

    resLO     = zero
    resNLO    = zero
    resHLO    = zero
    resHNLO   = zero
    resWWLO   = zero
    resWWNLO  = zero
    resWWHLO  = zero
    resWWHNLO = zero
    poles = zero

    ewampLO_ggww_all = czero
    ewampNLO_ggww_all = czero

    ewampLO_ggww_all = czero
    ewampLO_ggww_mless = czero
    ewampLO_ggww_mass = czero
    ewampNLO_ggww_mless = czero
    ewampNLO_ggww_mass = czero
    ewampNLO_ggww_all = czero
    ampLO_wwh2ww = czero
    ampNLO_wwh2ww = czero
    ampLO_ggh2ww = czero
    ampNLO_ggh2ww = czero
    ampLO_ggh2ww_heft = czero
    ampNLO_ggh2ww_heft = czero

    
    ! set up spinors
    call awayfromzaxis(6,p,pdum)
    call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)
    ! calculate gg -> WW amps
    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then

       if (mlessloops) then
          if (first) print *, "including massless loops"
          ewampLO_ggww_mless = ewampLO_ggww(1,2,3,4,5,6,za,zb,sprod)
          call spinorur(6,pdum,za,zb,sprod)
          call getewampsNLO_ggww(za,zb,sprod,E2LL,E2LR,ewampNLO_ggww_mless) !-- g g -> WW notation
          ! phase shift to have agreement with LO results (checked against MCFM)
          ewampNLO_ggww_mless = ci * ewampNLO_ggww_mless

          if (massloops) then   !    
             ewampLO_ggww_mass = ewampLO_ggWWmass(pdum,za,zb,sprod,mt)   ! call one-loop amp, with mt and mb=0

             ! reweight for each helicity
             do i1=-1,1,2
                do i2=-1,1,2
                   do i3=-1,1,2
                      do i5=-1,1,2
                         if (abs(ewampLO_ggww_mless(i1,i2,i3,i5)) .gt. zero) then
                            ewampNLO_ggww_mass(i1,i2,i3,i5) =  ewampNLO_ggww_mless(i1,i2,i3,i5) * ewampLO_ggww_mass(i1,i2,i3,i5)/ewampLO_ggww_mless(i1,i2,i3,i5)
                         endif
                         
                      enddo
                   enddo
                enddo
             enddo
          endif

          ! get "normal" spinors again
          call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)
       endif
      
       if (massloops .and. first) then
             print *, "Massive loops exact at LO"
             print *, "Massive two-loop amps included via reweighting"
       endif

       ewamplo_ggww_all = ewamplo_ggww_mless  + ewamplo_ggww_mass
       ewampnlo_ggww_all = ewampnlo_ggww_mless + ewampnlo_ggww_mass

    endif

    ! calculate gg->H->WW amplitudes
    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
       call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)
       if (HiggsExp) then     
          if (first) print *, "Using heavy top expansion for Higgs amplitudes, only leading term retained"
          ampLO_ggh2ww_heft = ewamplo_ggh2ww_heft(1,2,3,4,5,6,za,zb,sprod)
          call getewamplsNLO_ggh2ww_heft(1,2,3,4,5,6,za,zb,sprod,ampNLO_ggh2ww_heft)
          ampLO_ggh2ww(:,:,:,:) = ampLO_ggh2ww_heft(0,:,:,:,:)
          ampNLO_ggh2ww(:,:,:,:) = ampNLO_ggh2ww_heft(0,:,:,:,:)
       else
          if (first) print *, "Using full mass dependence for Higgs amplitudes"
          ampLO_ggh2ww = ewamplo_ggh2ww(1,2,3,4,5,6,za,zb,sprod)
          call getewamplsNLO_ggh2ww(1,2,3,4,5,6,za,zb,sprod,ampNLO_ggh2ww)
       endif
    endif

    ! add the Higgs and WW amplitudes
    if (contr .eq. "full" .or. contr .eq. "intf") then
          ampLO_wwh2ww = ewampLO_ggww_all + ampLO_ggh2ww
          ampNLO_wwh2ww = ewampNLO_ggww_all + ampNLO_ggh2ww
    endif

    do i1=-1,1,2
    do i2=-1,1,2
    do i3=-1,1,2
    do i5=-1,1,2

       ! |H|^2
       resHLO = resHLO + abs(ampLO_ggh2ww(i1,i2,i3,i5))**2
       resHNLO = resHNLO + 2.0_dp*real(ampLO_ggh2ww(i1,i2,i3,i5)*conjg(ampNLO_ggh2ww(i1,i2,i3,i5)),kind=dp)

       ! |WW|^2
       resWWLO = resWWLO + abs(ewamplo_ggww_all(i1,i2,i3,i5))**2
       resWWNLO = resWWNLO + 2.0_dp*real(ewamplo_ggww_all(i1,i2,i3,i5)*conjg(ewampnlo_ggww_all(i1,i2,i3,i5)),kind=dp)

       ! |W+HW|^2
       resWWHLO = resWWHLO + abs(ampLO_wwh2ww(i1,i2,i3,i5))**2
       resWWHNLO = resWWHNLO + 2.0_dp*real(amplo_wwh2ww(i1,i2,i3,i5)*conjg(ampnlo_wwh2ww(i1,i2,i3,i5)),kind=dp)

    enddo
    enddo
    enddo
    enddo

    first=.false.
    ! define res according to the 'contr' input
    if (contr .eq. "bkgd") then 
       if (first) print *, "Doing background |gg->WW|^2"
       resLO = resWWLO
       resNLO = resWWNLO
    elseif (contr .eq. "sigl") then
       if (first) print *, "Doing signal |gg->H->WW|^2"
       resLO = resHLO
       resNLO = resHNLO
    elseif (contr .eq. "full") then 
       if (first) print *, "Doing full |gg->H->WW + gg->WW|^2"
       resLO = resWWHLO
       resNLO = resWWHNLO
    elseif (contr .eq. "intf") then
       resLO = resWWHLO - resWWLO - resHLO
       resNLO = resWWHNLO - resWWNLO - resHNLO
       if (first) print *, "Doing interference 2Real((gg->H->WW)*conj(gg->WW))"
    endif

    !-- fix pre-factor
    resLO = resLO * avegg * colf 
    resNLO = resNLO * avegg * colf 

 
!    print *, "p1=(/",p(2,1), p(3,1), p(4,1), p(1,1),"/)"
!    print *, "p2=(/",p(2,2), p(3,2), p(4,2), p(1,2),"/)"
!    print *, "p3=(/",p(2,3), p(3,3), p(4,3), p(1,3),"/)"
!    print *, "p4=(/",p(2,4), p(3,4), p(4,4), p(1,4),"/)"
!    print *, "p5=(/",p(2,5), p(3,5), p(4,5), p(1,5),"/)"
!    print *, "p6=(/",p(2,6), p(3,6), p(4,6), p(1,6),"/)"
!    print *, "resLO",resLO
!    print *, "resNLO",resNLO
!    stop
    ! poles, as in massless piece
    ! might want to check on this order-by-order in 1/mt exp.
    poles(-2) = -2*CA*resLO
    poles(-1) = (-two * b0 - two * CA * log(musq/sprod(1,2)))*resLO

    return

  end subroutine resNLOV_WW

end module mod_ampsqNLOV_ww
