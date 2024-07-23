module mod_ampsqLO_ww
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_ampLO_ggww
  use mod_ampLO_ggh2ww

  implicit none
  private

  public :: resLO_WW,resLO_sping_WW

contains

  !-- g(1) g(2) -> [W -> l(3) lb(4)] [W -> l(5) lb(6)]
  !-- return coefficient of (as/twopi)**2
  subroutine resLO_WW(p,mass,res)
    real(dp), intent(in) :: p(4,6),mass
    real(dp), intent(out) :: res
    complex(dp) :: ampww(-1:1,-1:1,-1:1,-1:1),ampww_mass(-1:1,-1:1,-1:1,-1:1),ampww_mless(-1:1,-1:1,-1:1,-1:1),amp_h2ww(-1:1,-1:1,-1:1,-1:1),amp_wwh2ww(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: amp_h2ww_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: pdum(4,6),sprod(6,6)
    real(dp) :: res_WW,res_H,res_HWW
    integer :: i1,i2
    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one
    logical, save :: first=.true.

    res_WW      = zero
    res_HWW     = zero
    res_H       = zero
    res         = zero
    ampww       = czero
    ampww_mless = czero
    ampww_mass  = czero
    amp_h2ww    = czero
    amp_h2ww_heft = czero

    call awayfromzaxis(6,p,pdum)
    call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)
! gg -> WW amplitudes
    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if (mlessloops) then
          if (first) print *, "including massless loops"
          ampww_mless = ewampLO_ggww(1,2,3,4,5,6,za,zb,sprod)
       endif
       if (massloops) then
          if (first .and. VVExp) then
             print *, "massive loops in HEFT not implemented!"
             stop
          endif
          if (.not. VVExp) then
             if (first)  print *, "including massive loops exactly"
             ampww_mass=ewampLO_ggWWmass(pdum,za,zb,sprod,mass)
          endif
       endif
       ampww = ampww_mless + ampww_mass
    endif

! gg -> H -> WW amplitudes
    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if (HiggsExp) then
          if (first) print *, "H->WW through massive loop only implemented for the leading term"
          amp_h2ww_heft = ewamplo_ggh2ww_heft(1,2,3,4,5,6,za,zb,sprod)
          amp_h2ww(:,:,:,:) = amp_h2ww_heft(0,:,:,:,:)
       else
          if (first) print *, "Using full mass dependence for Higgs amplitudes"
          amp_h2ww = ewamplo_ggh2ww(1,2,3,4,5,6,za,zb,sprod)
       endif
    endif

! add the Higgs and WW amplitudes
    if (contr .eq. "full" .or. contr .eq. "intf") then
       amp_wwh2ww = ampww + amp_h2ww
    endif

    do i1=-1,1,2
       do i2=-1,1,2
          if (contr .eq. "bkgd" .or. contr .eq. "intf") then
             res_WW = res_WW + abs(ampww(i1,i2,-1,-1))**2
          endif
          if (contr .eq. "sigl" .or. contr .eq. "intf") then
             res_H = res_H + abs(amp_h2ww(i1,i2,-1,-1))**2
          endif
          if (contr .eq. "full" .or. contr .eq. "intf") then
             res_HWW = res_HWW + abs(amp_wwh2ww(i1,i2,-1,-1))**2
          endif
       enddo
    enddo

! define res according to the 'contr' input
    if (contr .eq. "bkgd") then                ! |gg->WW|^2
       if (first) print *, "Doing background |gg->WW|^2"
       res = res_WW
    elseif (contr .eq. "sigl") then            ! |gg->H->WW|^2
       if (first) print *, "Doing signal |gg->H->WW|^2"
       res = res_H
    elseif (contr .eq. "full") then            ! |gg->H->WW  +  gg->WW|^2
       if (first) print *, "Doing full  |gg->H->WW + gg->WW|^2"
       res = res_HWW
    elseif (contr .eq. "intf") then            ! |gg->H->WW + gg->WW|^2 - |gg->WW|^2 - |gg->H->WW|^2 = 2*Real( (gg->H->WW)*conj(gg->WW) )
       res = res_HWW - res_WW - res_H
       if (first) print *, "Doing interference 2Real( (gg->H->WW)*conj(gg->WW) )"
    endif

    !-- fix pre-factor
    res = res * colf

    res = res*avegg
    first=.false.
    return

  end subroutine resLO_WW



  !-- g(1) g(2) -> [W -> l(3) lb(4)] [W -> l(5) lb(6)]
  !-- return coefficient of (as/twopi)**2
  subroutine resLO_sping_WW(p,mass,res)
    real(dp), intent(in) :: p(4,6),mass
    complex(dp), intent(out) :: res(-1:1,-1:1)
    complex(dp) :: ampww(-1:1,-1:1,-1:1,-1:1),ampww_mass(-1:1,-1:1,-1:1,-1:1),ampww_mless(-1:1,-1:1,-1:1,-1:1),amp_h2ww(-1:1,-1:1,-1:1,-1:1),amp_wwh2ww(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: amp_h2ww_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: pdum(4,6),sprod(6,6)
    complex(dp) :: res_WW(-1:1,-1:1),res_H(-1:1,-1:1),res_HWW(-1:1,-1:1)
    integer :: i1,i2,j1
    real(dp), parameter :: colf = 8.0_dp ! Ca**2-one
    logical, save :: first=.true.

    res_WW      = czero
    res_HWW     = czero
    res_H       = czero
    res         = zero
    ampww       = czero
    ampww_mless = czero
    ampww_mass  = czero
    amp_h2ww    = czero
    amp_h2ww_heft = czero

    call awayfromzaxis(6,p,pdum)
    call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)
! gg -> WW amplitudes
    if (contr .eq. "bkgd" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if (mlessloops) then
          if (first) print *, "including massless loops"
          ampww_mless = ewampLO_ggww(1,2,3,4,5,6,za,zb,sprod)
       endif
       if (massloops) then
          if (first .and. VVExp) then
             print *, "massive loops in HEFT not implemented!"
             stop
          endif
          if (first .and. .not. VVExp) then
             print *, "including massive loops exactly"
             ampww_mass=ewampLO_ggWWmass(pdum,za,zb,sprod,mass)
          endif
       endif
       ampww = ampww_mless + ampww_mass
    endif

! gg -> H -> WW amplitudes
    if (contr .eq. "sigl" .or. contr .eq. "full" .or. contr .eq. "intf") then
       if (HiggsExp) then
          if (first) print *, "H->WW through massive loop only implemented for the leading term"
          amp_h2ww_heft = ewamplo_ggh2ww_heft(1,2,3,4,5,6,za,zb,sprod)
          amp_h2ww(:,:,:,:) = amp_h2ww_heft(0,:,:,:,:)
       else
          if (first) print *, "Using full mass dependence for Higgs amplitudes"
          amp_h2ww = ewamplo_ggh2ww(1,2,3,4,5,6,za,zb,sprod)
       endif
    endif

! add the Higgs and WW amplitudes
    if (contr .eq. "full" .or. contr .eq. "intf") then
       amp_wwh2ww = ampww + amp_h2ww
    endif

    do i1=-1,1,2
       do j1=-1,1,2
          do i2=-1,1,2
             if (contr .eq. "bkgd" .or. contr .eq. "intf") then
                res_WW(i1,j1) = res_WW(i1,j1) + ampww(i1,i2,-1,-1)*conjg(ampww(j1,i2,-1,-1))
             endif
             if (contr .eq. "sigl" .or. contr .eq. "intf") then
                res_H(i1,j1) = res_H(i1,j1) + amp_h2ww(i1,i2,-1,-1)*conjg(amp_h2ww(j1,i2,-1,-1))
             endif
             if (contr .eq. "full" .or. contr .eq. "intf") then
                res_HWW(i1,j1) = res_HWW(i1,j1) + amp_wwh2ww(i1,i2,-1,-1)*conjg(amp_wwh2ww(j1,i2,-1,-1))
             endif
          enddo
       enddo
    enddo

! define res according to the 'contr' input
    if (contr .eq. "bkgd") then                ! |gg->WW|^2
       if (first) print *, "Doing background |gg->WW|^2"
       res = res_WW
    elseif (contr .eq. "sigl") then            ! |gg->H->WW|^2
       if (first) print *, "Doing signal |gg->H->WW|^2"
       res = res_H
    elseif (contr .eq. "full") then            ! |gg->H->WW  +  gg->WW|^2
       if (first) print *, "Doing full  |gg->H->WW + gg->WW|^2"
       res = res_HWW
    elseif (contr .eq. "intf") then            ! |gg->H->WW + gg->WW|^2 - |gg->WW|^2 - |gg->H->WW|^2 = 2*Real( (gg->H->WW)*conj(gg->WW) )
       res = res_HWW - res_WW - res_H
       if (first) print *, "Doing interference 2Real( (gg->H->WW)*conj(gg->WW) )"
    endif

    !-- fix pre-factor
    res = res * colf

    res = res*avegg
    first=.false.
    return

  end subroutine resLO_sping_WW

end module mod_ampsqLO_ww


