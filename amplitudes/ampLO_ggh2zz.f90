module mod_ampLO_ggh2zz
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_func_for_h2vv
  implicit none
  private

  public :: ewamplo_ggh2zz_smeft !-- amplitude with SMEFT couplings
  public :: ewamplo_ggh2zz
  public :: ewamplo_ggh2zz_heft
  public :: ewamptilde_ggh2zz !-- amplitude without the ff massive form factor

contains
  ! -- with SMEFT couplings
  function ewampLO_ggh2zz_smeft(j1,j2,j3,j4,j5,j6,za,zb,sprod)
	complex(dp) :: ewampLO_ggh2zz_smeft(-1:1,-1:1,-1:1,-1:1) !amplitude is an array of 4 arrays of dimension 3. 0th index is not used, -1 and 1 are helicities (date: 30/03/24)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    real(dp) :: s12
    complex(dp) :: ff,ff_t,ff_b,spinamp(-1:1,-1:1,-1:1,-1:1)
    integer :: h3,h5

    ewampLO_ggh2zz_smeft = czero ! means complex zero (date: 23/03/24)

    s12 = sprod(j1,j2) ! invariant mass of the system, sprod is scalar product (02/04/2024)

    !-- finite mt-mb
    ff_t = getff(s12,expmass**2) 
	if(hwithb) then
		ff_b = getff(s12, mbsq)
	else
		ff_b = czero
	endif
	ff = three/two*cggh + ct*ff_t + cb*ff_b
    
    spinamp = ewamptilde_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    do h3 = -1,1,2
    do h5 = -1,1,2
       ewampLO_ggh2zz_smeft(-1,-1,h3,h5) = ff * spinamp(-1,-1,h3,h5)
       ewampLO_ggh2zz_smeft(+1,+1,h3,h5) = ff * spinamp(+1,+1,h3,h5)
    enddo
    enddo

    return
	
  end function ewamplo_ggh2zz_smeft  

  !-- 0 -> g(1) g(2) [l(3) lb(4)] [l(5) lb(6)]
  !-- label: h1,h2,h3,h5,coefficient of as/twopi * delta(a,b)
  !-- sign to match MCFM
  !-- it is simpler to implement the smeft couplings directly in the original function
  function ewampLO_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: ewampLO_ggh2zz(-1:1,-1:1,-1:1,-1:1) !amplitude is an array of 4 arrays of dimension 3. 0th index is not used, -1 and 1 are helicities (date: 30/03/24)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    real(dp) :: s12
    complex(dp) :: ff,ff_t,ff_b,spinamp(-1:1,-1:1,-1:1,-1:1)
    integer :: h3,h5

    ewampLO_ggh2zz = czero ! means complex zero (date: 23/03/24)

    s12 = sprod(j1,j2) ! invariant mass of the system, sprod is scalar product (02/04/2024)

    !-- finite mt-mb
    ff_t = getff(s12,expmass**2)
    if (hwithb) then
       ff_b = getff(s12,mbsq)
    else
       ff_b = czero
    endif
    ff = three/two*cggh + ct*ff_t + cb*ff_b

    spinamp = ewamptilde_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    do h3 = -1,1,2
    do h5 = -1,1,2
       ewampLO_ggh2zz(-1,-1,h3,h5) = ff * spinamp(-1,-1,h3,h5)
       ewampLO_ggh2zz(+1,+1,h3,h5) = ff * spinamp(+1,+1,h3,h5)
    enddo
    enddo

    return

  end function ewampLO_ggh2zz

  !-- same in the 1/mt expansion, first index is the order
  function ewampLO_ggh2zz_heft(j1,j2,j3,j4,j5,j6,za,zb,sprod)
    complex(dp) :: ewampLO_ggh2zz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    real(dp) :: s12,ff(0:nheft_max)
    complex(dp) :: spinamp(-1:1,-1:1,-1:1,-1:1)
    integer :: h3,h5

    ewampLO_ggh2zz_heft = czero

    s12 = sprod(j1,j2)

    !-- finite mt
    ff = getff_heft(s12,expmass**2)


    spinamp = ewamptilde_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    do h3 = -1,1,2
    do h5 = -1,1,2
       ewampLO_ggh2zz_heft(:,-1,-1,h3,h5) = ff(:) * spinamp(-1,-1,h3,h5)
       ewampLO_ggh2zz_heft(:,+1,+1,h3,h5) = ff(:) * spinamp(+1,+1,h3,h5)
    enddo
    enddo

    return

  end function ewampLO_ggh2zz_heft

  !-- spinor and ew structures, no mass form factor
  !-- 0 -> g(1) g(2) [l(3) lb(4)] [l(5) lb(6)]
  !-- label: h1,h2,h3,h5,coefficient of as/twopi * delta(a,b)
  !-- sign to match MCFM
  function ewamptilde_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod) ! j1 and j2 are gg 4-momenta and j3-6 are 4l 4-momenta, I guess (02/04/2024)
    complex(dp) :: ewamptilde_ggh2zz(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6), zb(6,6) !dont know what these matrices are (02/04/2024)
    real(dp), intent(in) :: sprod(6,6)
    real(dp) :: prefactor_coupl
    complex(dp) :: prod(-1:1,-1:1),decay(-1:1,-1:1)
    complex(dp) :: prop34,prop56,prefac
    real(dp) :: s12
    integer :: h3,h5

    prefactor_coupl = gwsq**2/12.0_dp/cosW2**2 !-- see docs/higgsamp.pdf
    ewamptilde_ggh2zz = czero

    !-- result taken from MCFM, src/ZZ/getggHZZamps.f
    s12 = sprod(j1,j2)

    prod(-1,-1) = za(j1,j2)**2 ! what the hell are za and zb (03/04/24)
    prod(+1,+1) = zb(j1,j2)**2

    !-- H -> ZZ, ignore interference for the time being
    prop34 = one/(sprod(j3,j4)-mzsq + ci * mz*gaz) ! gaz and gah are widths (02/04/24)
    prop56 = one/(sprod(j5,j6)-mzsq + ci * mz*gaz)

    prefac = ci * prefactor_coupl * prop34 * prop56 / (s12 - mhsq + ci * mh * gah)

    decay(-1,-1) = za(j3,j5)*zb(j6,j4) * Llep * Llep ! Llep and Rlep are related to left and right handed leptons (03/04/24)
    decay(+1,-1) = za(j4,j5)*zb(j6,j3) * Rlep * Llep
    decay(-1,+1) = za(j3,j6)*zb(j5,j4) * Llep * Rlep
    decay(+1,+1) = za(j4,j6)*zb(j5,j3) * Rlep * Rlep

    do h3=-1,1,2
    do h5=-1,1,2
       ewamptilde_ggh2zz(-1,-1,h3,h5) = prod(-1,-1) * decay(h3,h5) * prefac
       ewamptilde_ggh2zz(+1,+1,h3,h5) = prod(+1,+1) * decay(h3,h5) * prefac
    enddo
    enddo

    return

  end function ewamptilde_ggh2zz


end module mod_ampLO_ggh2zz
