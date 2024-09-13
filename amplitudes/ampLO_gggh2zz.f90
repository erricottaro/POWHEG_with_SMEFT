module mod_ampR_gggh2zz
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_func_for_h2vv_real
  implicit none
  private

  public :: resR_H2ZZ
  public :: ewampR_gggh2zz,ewampR_gggh2zz_heft

contains

  !-- g(1) g(2) -> H -> [Z -> l(3) lb(4)] [Z -> l(5) lb(6)]
  !-- return coefficient of (as/twopi)**2
  subroutine resR_H2ZZ(p,res)
    real(dp), intent(in) :: p(4,7)
    real(dp), intent(out) :: res
    complex(dp) :: ampzz(-1:1,-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: za(7,7),zb(7,7)
    real(dp) :: pdum(4,7),sprod(7,7)
    integer :: i2,i3,i4,i6
    real(dp), parameter :: colf = four * CA**2 * CF

    res = zero

    call awayfromzaxis(7,p,pdum)
    call spinorur(7,(/-pdum(:,1),-pdum(:,2),pdum(:,3),pdum(:,4),pdum(:,5),pdum(:,6),pdum(:,7)/),za,zb,sprod)

    ampzz = ewampR_gggh2zz(1,2,3,4,5,6,7,za,zb,sprod)

    do i2=-1,1,2
    do i3=-1,1,2
    do i4=-1,1,2
    do i6=-1,1,2
       res = res + abs(ampzz(+1,i2,i3,i4,i6))**2
    enddo
    enddo
    enddo
    enddo

    res = res * two !-- other helicities

    !-- fix pre-factor
    res = res * colf 

    res = res*avegg

    return

  end subroutine resR_H2ZZ

  !-- 0 -> g(1) g(2) g(3) H -> [l(4) lb(5)] [l(6) lb(7)]
  !-- label: h1,h2,h3,h4,h6 coefficient of as/twopi * gs * (sqrt2 f(a,b,c))
  !-- same notation of Raoul's interference paper, see higgsamp_real.pdf
  function ewampR_gggh2zz(j1,j2,j3,j4,j5,j6,j7,za,zb,sprod)
    complex(dp) :: ewampR_gggh2zz(-1:1,-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
    complex(dp), intent(in) :: za(7,7),zb(7,7)
    real(dp), intent(in) :: sprod(7,7)
    complex(dp) :: ampprod(-1:1,-1:1,-1:1),decay(-1:1,-1:1)
    complex(dp) :: propH,prop45,prop67,prefac
    integer :: h1,h2,h3,h4,h6
    real(dp) :: qsq






    

    ewampR_gggh2zz = czero

    qsq = sprod(j1,j2)+sprod(j1,j3)+sprod(j2,j3)

    propH = one/(qsq - mhsq + ci * mh * gah)
    prop45 = one/(sprod(j4,j5)-mzsq + ci * mz*gaz)
    prop67 = one/(sprod(j6,j7)-mzsq + ci * mz*gaz)

!   here introduce Z rescaling cz
    prefac = -ci * gwsq**2/12.0_dp/cosW2**2*propH*prop45*prop67*cz

    call get_ampprod_h2vv_real(j1,j2,j3,j4,j5,j6,j7,za,zb,sprod,ampprod)

    decay(-1,-1) = za(j4,j6)*zb(j7,j5) * Lel * Lel 
    decay(+1,-1) = za(j5,j6)*zb(j7,j4) * Rel * Lel 
    decay(-1,+1) = za(j4,j7)*zb(j6,j5) * Lel * Rel 
    decay(+1,+1) = za(j5,j7)*zb(j6,j4) * Rel * Rel 

    do h1=-1,1,2
    do h2=-1,1,2
    do h3=-1,1,2
    do h4=-1,1,2
    do h6=-1,1,2
       ewampR_gggh2zz(h1,h2,h3,h4,h6) = prefac * ampprod(h1,h2,h3) * decay(h4,h6)
    enddo
    enddo
    enddo
    enddo
    enddo

    return

  end function ewampR_gggh2zz

  !-- same but as an expansion in 1/fourmtsq. First index is the expansion order
  function ewampR_gggh2zz_heft(j1,j2,j3,j4,j5,j6,j7,za,zb,sprod)
    complex(dp) :: ewampR_gggh2zz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1,-1:1)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
    complex(dp), intent(in) :: za(7,7),zb(7,7)
    real(dp), intent(in) :: sprod(7,7)
    complex(dp) :: ampprod(0:nheft_max,-1:1,-1:1,-1:1),decay(-1:1,-1:1)
    complex(dp) :: propH,prop45,prop67,prefac
    integer :: h1,h2,h3,h4,h6
    real(dp) :: qsq

    ewampR_gggh2zz_heft = czero

    qsq = sprod(j1,j2)+sprod(j1,j3)+sprod(j2,j3)

    propH = one/(qsq - mhsq + ci * mh * gah)
    prop45 = one/(sprod(j4,j5)-mzsq + ci * mz*gaz)
    prop67 = one/(sprod(j6,j7)-mzsq + ci * mz*gaz)

    prefac = -ci * gwsq**2/12.0_dp/cosW2**2*propH*prop45*prop67

    call get_ampprod_h2vv_real_heft(j1,j2,j3,j4,j5,j6,j7,za,zb,sprod,ampprod)

    decay(-1,-1) = za(j4,j6)*zb(j7,j5) * Lel * Lel 
    decay(+1,-1) = za(j5,j6)*zb(j7,j4) * Rel * Lel 
    decay(-1,+1) = za(j4,j7)*zb(j6,j5) * Lel * Rel 
    decay(+1,+1) = za(j5,j7)*zb(j6,j4) * Rel * Rel 

    do h1=-1,1,2
    do h2=-1,1,2
    do h3=-1,1,2
    do h4=-1,1,2
    do h6=-1,1,2
       ewampR_gggh2zz_heft(:,h1,h2,h3,h4,h6) = prefac * ampprod(:,h1,h2,h3) * decay(h4,h6)
    enddo
    enddo
    enddo
    enddo
    enddo

    return

  end function ewampR_gggh2zz_heft
  
end module mod_ampR_gggh2zz
