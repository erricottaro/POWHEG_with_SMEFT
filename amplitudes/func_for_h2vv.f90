module mod_func_for_h2vv
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: getff,getff_heft

contains

  !-- top quark triangle amplitude, -> 1 in the heavy top approximation
  !-- this is 3/8/s Mbar[s] in higgsamp.pdf
  function getff(s12,mfsq)
    complex(dp) :: getff
    real(dp) :: s12,mfsq

    getff = three/s12 * mfsq * (two - sI3(s12,mfsq)*(one-four*mfsq/s12))

    return

  end function getff

  !-- the same, but as an expansion in 1/mt
  function getff_heft(s12,mfsq)
    integer :: nterm
    real(dp) :: getff_heft(0:nheft_max)
    real(dp) :: s12,mfsq,x

    x = s12/four/mfsq

    getff_heft(0) = one
    getff_heft(1) = 7/30.0_dp * x
    getff_heft(2) = 2/21.0_dp * x**2
    getff_heft(3) = 26/525.0_dp * x**3
    getff_heft(4) = 512/17325.0_dp * x**4
    getff_heft(5) = 1216/63063.0_dp * x**5

    return

  end function getff_heft

  !-- scalar triangle multiplied by s
  function sI3(s12,mfsq) ! -2*f(tau) 
    complex(dp) :: sI3
    real(dp) :: s12,mfsq,invtau

    invtau = four*mfsq/s12 !-- s12 is the invariant mass? (30/03/24)

    if (invtau.ge.one) then !-- below threshold
       sI3 = -two * asin(one/sqrt(invtau))**2
    else !-- above threshold
       sI3 = half * ( log((one+sqrt(one-invtau))/(one-sqrt(one-invtau))) - ci*pi)**2
    endif

    return

  end function sI3

end module mod_func_for_h2vv
