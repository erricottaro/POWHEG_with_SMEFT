module common_def
  use mod_types
  implicit none
  private

  real(dp), public :: s

  real(dp), public :: mu
  real(dp), public :: musq

  real(dp), public :: buff

! --- the number of generations
  integer, public :: nup = 2
  integer, public :: ndn = 3
  integer, public :: ng = 2

end  module common_def
