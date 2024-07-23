module mod_ampNLO_ggh2zz
  use mod_types; use common_def; use mod_consts_dp
  use mod_ampLO_ggh2zz
  use mod_auxfunctions
  use mod_func_for_h2vv
  use mod_func_for_h2vv_2loop
  implicit none
  private
  
  logical :: polecheck = .false.

  public :: getewamplsLONLO_ggh2zz,getewamplsLONLO_ggh2zz_heft

contains

  !-- g(1) g(2) -> [Z->l(3) lb(4)] [Z->l(5) lb(6)]
  !-- returns both the LO and NLO amplitudes, with delta(a,b) asontwopi[asontwopi**2] factored out
  !-- finite part according to qT --- CHECK THE I BETA AND SIMILAR PHASES
  subroutine getewamplsLONLO_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod,amp1l,amp2l)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: amp1l(-1:1,-1:1,-1:1,-1:1),amp2l(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ff1l,ff2l,ff_t,ff_b,spinamp(-1:1,-1:1,-1:1,-1:1)
    integer :: h3,h5
    real(dp) :: s12

    s12 = sprod(j1,2)

    if (polecheck) then
       print *, "*********************************************************************"
       print *, "* WARNING: substituting double pole values for two-loop amplitudes! *"
       print *, "*********************************************************************"
    endif


    !-- 1-loop, finite mt-mb
    ff_t = getff(s12,expmass**2)
    if (hwithb) then
       ff_b = getff(s12,mbsq)
    else
       ff_b = czero
    endif
    !-- adding SMEFT coefficients (23/06/2024)
    ff1l =three/two*cggh + ct*ff_t + cb*ff_b

    !-- 2-loop, finite mt-mb (expmass is the mass of the top (09/04/2024))
    ff_t = getff_2l(s12,expmass**2)
    if (hwithb) then
       ff_b = getff_2l(s12,mbsq)
    else
       ff_b = czero
    endif
    ff2l = ct*ff_t + cb*ff_b

    !-- this way, ff2l -> 11/2. Now add the ''right'' pisq
    ff2l = ff2l + Ca * pisq/two * ff1l

    if (polecheck) then
       ff2l = -CA * ff1l
    endif

    !-- now get the amplitudes
    spinamp = ewamptilde_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    do h3 = -1,1,2
    do h5 = -1,1,2
       amp1l(-1,-1,h3,h5) = ff1l * spinamp(-1,-1,h3,h5)
       amp1l(+1,+1,h3,h5) = ff1l * spinamp(+1,+1,h3,h5)

       amp2l(-1,-1,h3,h5) = ff2l * spinamp(-1,-1,h3,h5)
       amp2l(+1,+1,h3,h5) = ff2l * spinamp(+1,+1,h3,h5)
    enddo
    enddo

    return
    
  end subroutine getewamplsLONLO_ggh2zz

  subroutine getewamplsLONLO_ggh2zz_heft(j1,j2,j3,j4,j5,j6,za,zb,sprod,amp1l,amp2l)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: amp1l(0:nheft_max,-1:1,-1:1,-1:1,-1:1),amp2l(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: spinamp(-1:1,-1:1,-1:1,-1:1)
    real(dp) :: s12,ff1l(0:nheft_max),ff2l(0:nheft_max)
    integer :: h3,h5

    s12 = sprod(j1,2)

    !-- 1-loop, finite mt-mb
    ff1l = getff_heft(s12,expmass**2)

    !-- 2-loop, finite mt-mb
    ff2l = getff_2l_heft(s12,expmass**2)
    
    !-- this way, ff2l -> 11/2. Now add the ''right'' pisq
    ff2l = ff2l + Ca * pisq/two * ff1l

    !-- now get the amplitudes
    spinamp = ewamptilde_ggh2zz(j1,j2,j3,j4,j5,j6,za,zb,sprod)

    do h3 = -1,1,2
    do h5 = -1,1,2
       amp1l(:,-1,-1,h3,h5) = ff1l(:) * spinamp(-1,-1,h3,h5)
       amp1l(:,+1,+1,h3,h5) = ff1l(:) * spinamp(+1,+1,h3,h5)

       amp2l(:,-1,-1,h3,h5) = ff2l(:) * spinamp(-1,-1,h3,h5)
       amp2l(:,+1,+1,h3,h5) = ff2l(:) * spinamp(+1,+1,h3,h5)
    enddo
    enddo

    return
    
  end subroutine getewamplsLONLO_ggh2zz_heft

  !-- g(1) g(2) -> H -> [Z -> l(3) lb(4)] [Z -> l(5) lb(6)]
  !-- return coefficient of (as/twopi)**3 (for the 2-loop), finite part according to QT
  !-- PHASES OF THE IMAGINARY PART MUST BE CHECKED
  ! subroutine resNLOV_h2zz(p,res)
  !   real(dp), intent(in) :: p(4,6)
  !   real(dp), intent(out) :: res
  !   complex(dp) :: ampzztr(-1:1,-1:1,-1:1,-1:1),ampzz1l(-1:1,-1:1,-1:1,-1:1)
  !   complex(dp) :: za(6,6),zb(6,6),fg1l
  !   real(dp) :: pdum(4,6),sprod(6,6)
  !   integer :: i1,i2,i3,i5
  !   real(dp), parameter :: colf = 8.0_dp ! Ca**2-one

  !   res = zero

  !   call awayfromzaxis(6,p,pdum)
  !   call spinorur(6,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6)/),za,zb,sprod)

  !   call getewamplsLONLO_ggh2zz(1,2,3,4,5,6,za,zb,sprod,ampzztr,ampzz1l)

  !   do i1=-1,1,2
  !   do i2=-1,1,2
  !   do i3=-1,1,2
  !   do i5=-1,1,2
  !      res = res + ampzztr(i1,i2,i3,i5)*conjg(ampzz1l(i1,i2,i3,i5)) &
  !                + ampzz1l(i1,i2,i3,i5)*conjg(ampzztr(i1,i2,i3,i5))
  !   enddo
  !   enddo
  !   enddo
  !   enddo

  !   !-- fix pre-factor
  !   res = res * colf * avegg

  !   return

  ! end subroutine resNLOV_h2zz

end module mod_ampNLO_ggh2zz
