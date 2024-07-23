module mod_consts_dp
  use mod_types; use common_def
  implicit none
  private

  ! --- numerical constants
  real(dp), public, parameter :: pi =&
       & 3.141592653589793238462643383279502884197_dp
  real(dp), public, parameter :: twopi =&
       & 6.283185307179586476925286766559005768394_dp
  real(dp), public, parameter :: zeta2 =&
       & 1.644934066848226436472415166646025189219_dp
  real(dp), public, parameter :: zeta3 =&
       & 1.202056903159594285399738161511449990765_dp
  real(dp), parameter, public :: pisq =&
       & 9.869604401089358618834490999876151135314_dp
  real(dp), parameter, public :: eulergamma =&
       & 0.577215664901532860606512090082402431042_dp
  real(dp), public, parameter :: half  = 0.5_dp, two = 2.0_dp
  real(dp), public, parameter :: zero  = 0.0_dp, one = 1.0_dp
  real(dp), public, parameter :: mone = -1._dp
  real(dp), public, parameter :: three = 3._dp
  real(dp), public, parameter :: four  = 4._dp
  real(dp), public, parameter :: five  = 5._dp
  real(dp), public, parameter :: six  =  6._dp
  real(dp), public, parameter :: sqrt2 = &
       &1.4142135623730950488016887242096980785696718753769_dp
  real(dp), public, parameter :: msqrt2 = &
       &-1.4142135623730950488016887242096980785696718753769_dp
  real(dp), public, parameter :: sqrt3 = &
       &1.7320508075688772935274463415058723669428_dp
  complex(dp), public, parameter :: ci = cmplx(zero,1.0_dp,kind=dp)
  complex(dp), public, parameter :: czero = cmplx(0.0_dp,0.0_dp,kind=dp)
  complex(dp), public, parameter :: cone = cmplx(1.0_dp,0.0_dp,kind=dp)

  ! --- color factors
  real(dp), public, parameter :: CF = four/three
  real(dp), public, parameter :: CA = three
  real(dp), public, parameter :: Tr = half
  real(dp), public, parameter :: xn = 3.0_dp ! number of colors (03/04/2024)
  real(dp), public, parameter :: xnsq = 9.0_dp
  real(dp), public, parameter :: V = xnsq-one !-- for MCFM routines

  ! --- averaging factors
  real(dp), public, parameter :: aveqq = one/four/xnsq
  real(dp), public, parameter :: avegg = one/four/(xnsq-one)**2
  real(dp), public, parameter :: aveqg = one/four/(xnsq-one)/xn

  ! --- couplings and parameters
  real(dp), public :: mt = 173.2_dp
  real(dp), public :: mtsq != mt**2
  real(dp), public :: mb = 4.5_dp
  real(dp), public :: mbsq != mb**2
  real(dp), public :: mh = 125.0_dp
  real(dp), public :: mhsq != mh**2
  real(dp), public :: mz = 91.1876_dp
  real(dp), public :: mzsq != mz**2
  real(dp), public :: mw = 80.398_dp
  real(dp), public :: mwsq != mw**2
  !real(dp), public, parameter :: GaW = (two*xn+three)*Gf*MW**3/(six*sqrt2*Pi) !-- LO W width
  real(dp), public :: GaW = 2.1054_dp !-- MCFM W width
  real(dp), public :: GaZ = 2.4952_dp !-- MCFM Z width
  !real(dp), public, parameter :: GaH = 0.0041650_dp !-- MCFM H width
  real(dp), public :: GaH = 0.0041650_dp !-- MCFM H width
  real(dp), public, parameter :: Gf = 1.16639d-5
  real(dp), public :: gwsq != four * sqrt2 * MW**2 * Gf
  real(dp), public :: vev != one/sqrt(Gf*sqrt2)
  real(dp), public :: sinW2 != 0.2226459_dp
  real(dp), public :: cosW2 != one-sinW2
  real(dp), public :: sw != sqrt(sinW2)
  real(dp), public :: cw != sqrt(cosW2)
  !--- SMEFT couplings
  real(dp), public :: cggh
  real(dp), public :: cb
  real(dp), public :: ct
  !--
  real(dp), public, parameter :: Qup =  2.0_dp/3.0_dp
  real(dp), public, parameter :: Qdn = -1.0_dp/3.0_dp
  real(dp), public :: Vup != 1.0_dp/2.0_dp - 4.0_dp/3.0_dp*sinW2
  real(dp), public, parameter :: Aup = 1.0_dp/2.0_dp
  real(dp), public :: Vdn != -1.0_dp/2.0_dp + 2.0_dp/3.0_dp*sinW2
  real(dp), public, parameter :: Adn = -1.0_dp/2.0_dp
  real(dp), public :: Lup != Vup + Aup
  real(dp), public :: Rup != Vup - Aup
  real(dp), public :: Ldn != Vdn + Adn
  real(dp), public :: Rdn != Vdn - Adn
  !--
  real(dp), public, parameter :: Qel = -one
  real(dp), public, parameter :: Vnu = half
  real(dp), public, parameter :: Anu = half
  real(dp), public :: Vel != -half + 2.0_dp*sinW2
  real(dp), public, parameter :: Ael = -half
  real(dp), public :: Lel != Vel + Ael
  real(dp), public :: Rel != Vel - Ael
  real(dp), public, parameter :: Lnu = Vnu + Anu
  real(dp), public, parameter :: Rnu = Vnu - Anu
  real(dp), public :: Llep
  real(dp), public :: Rlep

  !--

  !-- quark name definition, according to LHAPDF
  integer, public, parameter :: gl = 0
  integer, public, parameter :: dn = 1
  integer, public, parameter :: up = 2
  integer, public, parameter :: st = 3
  integer, public, parameter :: ch = 4
  integer, public, parameter :: bt = 5
  integer, public, parameter :: tp = 6

  ! --- active flavors and QCD beta function
  real(dp), public :: nf != four
  real(dp), public :: nflav != nf   ! this is MCFM nflavors
  real(dp), public :: b0!=(xn*11.0_dp - 2.0_dp*nflav)/6.0_dp     ! QCD beta-function
  real(dp), public :: b1!=(17.0_dp*Ca**2 - 5.0_dp*Ca*nf - 3.0_dp*Cf*nf)/6.0_dp

  ! -- color structures
  real(dp), public, parameter :: CaCf = CA*CF
  real(dp), public, parameter :: CfCf = CF**2
  real(dp), public, parameter :: CaCa = CA**2
  real(dp), public, parameter :: CaTr = Ca/two
  real(dp), public, parameter :: CfTr = Cf/two
  real(dp), public :: CfNf != CF*nf
  real(dp), public :: CaNf != CA*nf
  real(dp), public :: NfNf != nf*nf

  ! --- buffers
  real(dp), public, parameter :: tolb = 1d-8
  real(dp), public, parameter :: tolb_lpsep = 1d-8
  real(dp), public, parameter :: fbuff = 1d-15
  real(dp), public, parameter :: lpsepbuff = 1d-40


  ! --- the RR/LL part
  integer, public, parameter :: LL = 1
  integer, public, parameter :: RR = 2

  !--- max number of terms in the 1/mt expansion
  integer, public, parameter :: nheft_max = 5
  integer, public, parameter :: expord = nheft_max !-- order of 1/mt expansion to keep
  real(dp), public :: expmass

  !--- process switches
  logical, public :: hwithb !-- finite b-quark effects in the Higgs amplitude
  character, public :: proc*(5)
  character, public :: contr*(4)
  logical, public   :: HiggsExp
  logical, public   :: VVExp
  logical, public   :: mlessloops
  logical, public   :: massloops


end module mod_consts_dp
