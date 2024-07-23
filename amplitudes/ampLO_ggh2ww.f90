module mod_ampLO_ggh2ww
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_func_for_h2vv
  implicit none
  private

  public :: ewamplo_ggh2ww
  public :: ewamplo_ggh2ww_heft
  public :: ewamptilde_ggh2ww !-- amplitude without the ff massive form factor

contains

  !-- 0 -> g(1) g(2) [Wp -> l(3) lb(4)] [Wm -> l(5) lb(6)]
  !-- label: h1,h2,h3,h5,coefficient of as/twopi * delta(a,b)
  !-- sign to match MCFM
  function ewampLO_ggh2ww(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: ewampLO_ggh2ww(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    real(dp) :: s12
    complex(dp) :: ff,ff_t,ff_b,spinamp(-1:1,-1:1,-1:1,-1:1)

    ewampLO_ggh2ww = czero

    s12 = sprod(j1,j2)

    !-- finite mt-mb
    ff_t = getff(s12,expmass**2)
    if (hwithb) then
       ff_b = getff(s12,mbsq)
    else
       ff_b = czero
    endif
    ff = ff_t + ff_b

    spinamp = ewamptilde_ggh2ww(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    ewampLO_ggh2ww(-1,-1,-1,-1) = ff * spinamp(-1,-1,-1,-1)
    ewampLO_ggh2ww(+1,+1,-1,-1) = ff * spinamp(+1,+1,-1,-1)

    return

  end function ewampLO_ggh2ww

  !-- same in the 1/mt expansion, first index is the order
  function ewampLO_ggh2ww_heft(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: ewampLO_ggh2ww_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    real(dp) :: s12,mfsq,ff(0:nheft_max)
    complex(dp) :: spinamp(-1:1,-1:1,-1:1,-1:1)

    ewampLO_ggh2ww_heft = czero

    !-- adapted from MCFM, src/ZZ/getggHZZamps.f
    s12 = sprod(j1,j2)

    !-- finite mt
    ff = getff_heft(s12,expmass)

    spinamp = ewamptilde_ggh2ww(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    ewampLO_ggh2ww_heft(:,-1,-1,-1,-1) = ff(:) * spinamp(-1,-1,-1,-1)
    ewampLO_ggh2ww_heft(:,+1,+1,-1,-1) = ff(:) * spinamp(+1,+1,-1,-1)

    return

  end function ewampLO_ggh2ww_heft

  !-- spinor structure and coupling, without the ff form factor
  !-- 0 -> g(1) g(2) [Wp -> l(3) lb(4)] [Wm -> l(5) lb(6)]
  !-- label: h1,h2,h3,h5,coefficient of as/twopi * delta(a,b)
  !-- sign to match MCFM
  function ewamptilde_ggh2ww(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: ewamptilde_ggh2ww(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    real(dp) :: prefactor_coupl
    real(dp) :: s12,mfsq
    complex(dp) :: prod(-1:1,-1:1),decay(-1:1,-1:1)
    complex(dp) :: prop34,prop56,prefac


    prefactor_coupl = gwsq**2/six !-- see docs/higgsamp.pdf
    ewamptilde_ggh2ww = czero

    !-- adapted from MCFM, src/ZZ/getggHZZamps.f
    s12 = sprod(j1,j2)

    prod(-1,-1) = za(j1,j2)**2
    prod(+1,+1) = zb(j1,j2)**2

    !-- H -> WW, ignore interference for the time being
    prop34 = one/(sprod(j3,j4)-mwsq + ci * mw*gaw)
    prop56 = one/(sprod(j5,j6)-mwsq + ci * mw*gaw)

    prefac = ci * prefactor_coupl * prop34 * prop56 / (s12 - mhsq + ci * mh * gah)

    decay(-1,-1) = za(j3,j5)*zb(j6,j4)

    ewamptilde_ggh2ww(-1,-1,-1,-1) = prod(-1,-1) * decay(-1,-1) * prefac
    ewamptilde_ggh2ww(+1,+1,-1,-1) = prod(+1,+1) * decay(-1,-1) * prefac

    return

  end function ewamptilde_ggh2ww

end module mod_ampLO_ggh2ww
