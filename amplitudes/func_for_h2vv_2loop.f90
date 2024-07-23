module mod_func_for_h2vv_2loop
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: getff_2l,getff_2l_heft

contains

  !-- two-loop massive coefficient c1_mh -> 11/2
  function getff_2l(s12,mfsq)
    complex(dp) :: getff_2l
    complex(dp15) :: x_dp15,HPL1,HPL2,HPL3,HPL4
    real(dp) :: s12,mfsq,tau
    complex(dp) :: t(9),x,onemx

    tau = s12/four/mfsq

    if (tau.lt.1) then
       x = cmplx(one-two*tau,two*sqrt(tau*(one-tau)),kind=dp)
    else
       x = cmplx((sqrt(one-one/tau)-one)/(sqrt(one-one/tau)+one),zero,kind=dp)
    endif

    !-- get basis functions, using x -> x + I delta
    x_dp15 = x
    t(1) = HPL1(0,x_dp15)
    t(2) = HPL1(1,x_dp15) 
    t(3) = HPL2(0,-1,x_dp15)
    t(4) = HPL2(0,1,x_dp15)
    t(5) = HPL3(0,0,1,x_dp15)
    t(6) = HPL3(0,0,-1,x_dp15)
    t(7) = HPL4(0,0,0,-1,x_dp15)
    t(8) = HPL4(0,0,0,1,x_dp15)
    t(9) = HPL4(1,0,-1,0,x_dp15)

    onemx = one-x
	! - z3 is Riemann zeta function evaluated in 3 (08/04/2024)
    getff_2l = (-94*x)/onemx**2 + (2*pisq**2*x*(1 + x - 3*x**2 - 3*x**3))/(5.0_dp*onemx**5) + &
         (2*x*(31 + 34*x + 31*x**2)*zeta3)/onemx**4 - (24*x*(1 + x)*t(1))/onemx**3 + &
         (2*x*(11 + 11*x - 43*x**2 - 43*x**3)*zeta3*t(1))/onemx**5 + &
         (3*x*(3 + 22*x + 3*x**2)*t(1)**2)/(2.0_dp*onemx**4) + &
         (pisq*x*(1 + x - 17*x**2 - 17*x**3)*t(1)**2)/(6.0_dp*onemx**5) + &
         (x*(6 + 59*x + 58*x**2 + 33*x**3)*t(1)**3)/(3.0_dp*onemx**5) + &
         (x*(5 + 5*x - 13*x**2 - 13*x**3)*t(1)**4)/(24.0_dp*onemx**5) - &
         (x*(47 + 66*x + 47*x**2)*t(1)**2*t(2))/onemx**4 + &
         (16*x*(1 + x + x**2 + x**3)*t(1)**2*t(3))/onemx**5 + &
         (2*x*(51 + 74*x + 51*x**2)*t(1)*t(4))/onemx**4 - &
         (2*x*(5 + 5*x + 23*x**2 + 23*x**3)*t(1)**2*t(4))/onemx**5 - &
         (2*x*(55 + 82*x + 55*x**2)*t(5))/onemx**4 + &
         (4*x*(1 + x)*(23 + 41*x**2)*t(1)*(t(5) - t(6)))/onemx**5 + &
         (36*x*(5 + 5*x + 11*x**2 + 11*x**3)*t(7))/onemx**5 - &
         (36*x*(5 + 5*x + 7*x**2 + 7*x**3)*t(8))/onemx**5 + &
         (x*(1 + x)**2*((-4*pisq*t(1))/3.0_dp + 108*zeta3*t(2) + 6*pisq*t(1)*t(2) - &
         6*t(1)**3*t(2) - 32*t(1)*t(3) - 6*pisq*t(4) + 64*t(6) + 72*t(9)))/onemx**4 

    return
    
  end function getff_2l

  !-- the same, but as an expansion in 1/mt
  function getff_2l_heft(s12,mfsq)
    integer :: nterm
    real(dp) :: getff_2l_heft(0:nheft_max)
    real(dp) :: s12,mfsq,tau,Lsm,Lum

    tau = s12/four/mfsq

    getff_2l_heft(0) = 11.0_dp/2.0_dp
    getff_2l_heft(1) = (1237._dp*tau)/540.0_dp
    getff_2l_heft(2) = (17863._dp*tau**2)/14175.0_dp
    getff_2l_heft(3) = (157483._dp*tau**3)/198450.0_dp
    getff_2l_heft(4) = (1273388._dp*tau**4)/2338875.0_dp
    getff_2l_heft(5) = (97424769412._dp*tau**5)/245827456875.0_dp

    return
    
  end function getff_2l_heft

end module mod_func_for_h2vv_2loop
