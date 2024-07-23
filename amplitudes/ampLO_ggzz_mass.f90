module mod_ampLO_ggzz_mass
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions

  implicit none
  private

  public :: getewampsLO_ggzz_mass_eps
  public :: getewampsLO_ggzz_heft

contains


  subroutine getewampsLO_ggzz_heft(za,zb,sprod,mass,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft)
! wrapper routine to return only the finite part of the 1loop ampl.
    complex(dp15), intent(in) :: za(6,6),zb(6,6)
    real(dp15), intent(in) :: sprod(6,6),mass
    complex(dp15), intent(out) :: ewampLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)
    complex(dp15), intent(out) :: ewampLO_ggzz_fullmass(-1:1,-1:1,-1:1,-1:1)
    complex(dp15)              :: ewampLO_ggzz_heft_eps(0:nheft_max,-1:1,-1:1,-1:1,-1:1,0:2)


    call getewampsLO_ggzz_mass_eps(za,zb,sprod,mass,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft_eps)
    ewampLO_ggzz_heft(0:nheft_max,-1:1,-1:1,-1:1,-1:1)=ewampLO_ggzz_heft_eps(0:nheft_max,-1:1,-1:1,-1:1,-1:1,0)

    return

  end subroutine getewampsLO_ggzz_heft


  !-- 0 -> g(1) g(2) [Z->l(3) lb(4)] [Z->l(5) lb(6)]
  !-- returns the LO amplitudes, with delta(a,b) asontwopi[asontwopi**2] factored out
  !-- first index is the expansion index (0=mt->infty limit)
  !-- next indices are h1,h2,h3,h5
  !-- final index is expansion in eps
  subroutine getewampsLO_ggzz_mass_eps(za,zb,sprod,mass,ewampLO_ggzz_fullmass,ewampLO_ggzz_heft_eps)
    complex(dp15), intent(in) :: za(6,6),zb(6,6)
    real(dp15), intent(in) :: sprod(6,6),mass
    complex(dp15), intent(out) :: ewampLO_ggzz_heft_eps(0:nheft_max,-1:1,-1:1,-1:1,-1:1,0:2)
    complex(dp15), intent(out) :: ewampLO_ggzz_fullmass(-1:1,-1:1,-1:1,-1:1)
    integer :: i1,i2,h1,h2,h3,h5,nheft,neps,hh1,hh2,hh3,hh5
    complex(dp15) :: A_lo_ax(20,0:5,-2:2),A_lo_vec(20,0:5,-2:2),T(20,-1:1,-1:1,-1:1,-1:1)
    complex(dp15) :: prop34,prop56
    complex(dp15) :: Mloop_uptype(-1:1,-1:1,-1:1,-1:1)
    real(dp15) :: twosinwcosw,cl1(-1:1),cl2(-1:1),cvec(1:2),cax(1:2),q1,q2,Qu,Qd
    complex(dp15) :: ewampLO_ggzz_ax(0:nheft_max,-1:1,-1:1,-1:1,-1:1,0:2),ewampLO_ggzz_vec(0:nheft_max,-1:1,-1:1,-1:1,-1:1,0:2)
    complex(dp15) :: spinorweight(-1:1,-1:1,-1:1,-1:1),ewampLO_ggzz_summed(-1:1,-1:1,-1:1,-1:1)
    complex(dp15) :: MCFMtop(-1:1,-1:1,-1:1,-1:1),MCFMup(-1:1,-1:1,-1:1,-1:1),MCFMdn(-1:1,-1:1,-1:1,-1:1)
    complex(dp15) :: AmtLLa(2,2,2,2),AmtLRa(2,2,2,2),Amt_vec,Amt_ax,AmtLL(-1:1,-1:1,-1:1,-1:1),AmtLR(-1:1,-1:1,-1:1,-1:1)

    ewampLO_ggzz_heft_eps = czero
    ewampLO_ggzz_ax = czero
    ewampLO_ggzz_vec = czero

! --- first get the exact mass-dependence in the ampl from the MCFM routines:
    call ggZZmassamp_new(za,zb,sprod,mass,AmtLLa,AmtLRa)


! -- ampl checks against MCFM
!    include "spinorweight.f90"
!    if (abs(expmass-500d0) .le. 1d-15) then
!       include "Checks/massexpLO/MCFM/ggZZ_ampcheck_pt1_mt500.out"
!    elseif (abs(expmass-1000d0) .le. 1d-15) then
!       include "Checks/massexpLO/MCFM/ggZZ_ampcheck_pt1_mt1000.out"
!    elseif (abs(expmass-2000d0) .le. 1d-15) then
!       include "Checks/massexpLO/MCFM/ggZZ_ampcheck_pt1_mt2000.out"
!    elseif (abs(expmass-5000d0) .le. 1d-15) then
!       include "Checks/massexpLO/MCFM/ggZZ_ampcheck_pt3_mt5000.out"
!    endif
! -- end checks



    !-- populate the form factors seperately for vector and axial
    call ffaxLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_lo_ax)
    call ffvecLO_ggzz_heft(1,2,3,4,5,6,mass,mu,za,zb,sprod,A_lo_vec)
    !-- populate the helicity array
    call heli_zz_heft(1,2,3,4,5,6,za,zb,sprod,T)


    !-- build the amplitudes
    do nheft=0,nheft_max
    do h1=-1,+1,2
    do h2=-1,+1,2
    do h3=-1,+1,2
    do h5=-1,+1,2
    do i1=1,20
    do neps=0,0
       ewampLO_ggzz_ax(nheft,h1,h2,h3,h5,neps) = ewampLO_ggzz_ax(nheft,h1,h2,h3,h5,neps) + &
            A_lo_ax(i1,nheft,neps)*T(i1,h1,h2,h3,h5)

       ewampLO_ggzz_vec(nheft,h1,h2,h3,h5,neps) = ewampLO_ggzz_vec(nheft,h1,h2,h3,h5,neps) + &
            A_lo_vec(i1,nheft,neps)*T(i1,h1,h2,h3,h5)
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo


! now include couplings, as done in ampLO_ggzz


 !--- propagator factors

    prop34=sprod(3,4)/(sprod(3,4)-mzsq + ci*mz*gaz)
    prop56=sprod(5,6)/(sprod(5,6)-mzsq + ci*mz*gaz)
    !-- couplings, see src/Need/couplz.f and coupling.f
    !-- xw = sin(tw)^2, sin2w = 2*sin(tw)*cos(tw)
    twosinwcosw=two*sw*cw

!    cl1(-1)=Llep/(twosinwcosw)
!    cl1(+1)=Rlep/(twosinwcosw)
!    cl2(-1)=Llep/(twosinwcosw)
!    cl2(+1)=Rlep/(twosinwcosw)
!
!    cvec(up)=Vup/(twosinwcosw)!*sqrt(gwsq*sinW2)
!    cvec(dn)=Vdn/(twosinwcosw)!*sqrt(gwsq*sinW2)
!    cax(up)=Aup/(twosinwcosw)!*sqrt(gwsq*sinW2)
!    cax(dn)=Adn/(twosinwcosw)!*sqrt(gwsq*sinW2)

    cl1(-1)=Llep/(2.0_dp*cw)
    cl1(+1)=Rlep/(2.0_dp*cw)
    cl2(-1)=Llep/(2.0_dp*cw)
    cl2(+1)=Rlep/(2.0_dp*cw)

    cvec(up)=Vup/(2.0_dp*cw)!*sqrt(gwsq*sinW2)
    cvec(dn)=Vdn/(2.0_dp*cw)!*sqrt(gwsq*sinW2)
    cax(up)=Aup/(2.0_dp*cw)!*sqrt(gwsq*sinW2)
    cax(dn)=Adn/(2.0_dp*cw)!*sqrt(gwsq*sinW2)

    q1=-one
    q2=-one
    Qu=Qup
    Qd=Qdn

    !--- dress vector and axial amplitudes with appropriate couplings
    do h1=-1,1,2
    do h2=-1,1,2
    do h3=-1,1,2
    do h5=-1,1,2
!-- exact mass-dependence in the loop
!-- indices translate from MCFM helicity conventions to ours
!      Amt_vec=2d0*(AmtLL(h1,h2,h3,h5)+AmtLR(h1,h2,h3,h5))
!      Amt_ax =2d0*(AmtLL(h1,h2,h3,h5)-AmtLR(h1,h2,h3,h5))
       hh1=(h1+3)/2
       hh2=(h2+3)/2
       hh3=(h3+3)/2
       hh5=(h5+3)/2
      Amt_vec=2d0*(AmtLLa(hh1,hh2,hh3,hh5)+AmtLRa(hh1,hh2,hh3,hh5))
      Amt_ax =2d0*(AmtLLa(hh1,hh2,hh3,hh5)-AmtLRa(hh1,hh2,hh3,hh5))

      ewampLO_ggzz_fullmass(h1,h2,h3,h5)=ci*( &
           Amt_vec*(Qu*q1*sinW2+cvec(up)*cl1(h3)*prop34)*(Qu*q2*sinW2+cvec(up)*cl2(h5)*prop56) &
           +Amt_ax*(cax(up)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56))

    do neps=0,2

       !--- internal loops of top quarks
!       print *, "2*ax", 2d0*abs(ewampLO_ggzz_ax(1,h1,h2,h3,h5)),2d0*abs(ewampLO_ggzz_ax(2,h1,h2,h3,h5))
!       print *, "2*ax", 2d0*abs(ewampLO_ggzz_ax(1,h1,h2,h3,h5))
!       print *, "2*ax, 2*vec", 2d0*abs(ewampLO_ggzz_ax(1,h1,h2,h3,h5)+ewampLO_ggzz_ax(2,h1,h2,h3,h5)+ewampLO_ggzz_ax(3,h1,h2,h3,h5)+ewampLO_ggzz_ax(4,h1,h2,h3,h5)),2d0*abs(ewampLO_ggzz_vec(2,h1,h2,h3,h5)+ewampLO_ggzz_vec(3,h1,h2,h3,h5)+ewampLO_ggzz_vec(4,h1,h2,h3,h5))
!       print *, "2*vec",2d0*abs(ewampLO_ggzz_vec(2,h1,h2,h3,h5))
!       print *, "2*vec",2d0*abs(ewampLO_ggzz_vec(2,h1,h2,h3,h5)+ewampLO_ggzz_vec(3,h1,h2,h3,h5)+ewampLO_ggzz_vec(4,h1,h2,h3,h5))

       ewampLO_ggzz_heft_eps(:,h1,h2,h3,h5,neps)=ewampLO_ggzz_vec(:,h1,h2,h3,h5,neps)*( &
            (Qu*q1*sinW2+cvec(up)*cl1(h3)*prop34)*(Qu*q2*sinW2+cvec(up)*cl2(h5)*prop56)) &
            +ewampLO_ggzz_ax(:,h1,h2,h3,h5,neps)*(cax(up)*cl1(h3)*prop34)*(cax(up)*cl2(h5)*prop56)

    enddo
    enddo
    enddo
    enddo
    enddo

! -- change normalizations from Tr[TaTb]=1/2*delta_[ab] (used by Matthew)
! -- to Tr[TaTb]=delta_[ab] (used in code)
! -- should include factor of sqrt(2) for each Ta
!    ewampLO_ggzz_eps=2d0*ewampLO_ggzz_eps


! normalization to agree with MCFM
    ewampLO_ggzz_heft_eps=4d0*ewampLO_ggzz_heft_eps
!    print *, "************ multiplying our result by 4 *****************"


! -- now add in phase to agree with MCFM
    ewampLO_ggzz_heft_eps=ci*ewampLO_ggzz_heft_eps



!! --R.R. -- check the phases using spinorweights
!    ewampLO_ggzz_summed=czero
!    do nheft=0,nheft_max
!       ewampLO_ggzz_summed(:,:,:,:) = ewampLO_ggzz_summed(:,:,:,:) + ewampLO_ggzz_heft_eps(nheft,:,:,:,:,0)
!    enddo
!    print *, "top amplitudes divided by spinorweights"
!    do h1=-1,1,2
!    do h2=-1,1,2
!    do h3=-1,1,2
!    do h5=-1,1,2
!       print *, h1,h2,h3,h5
!       print *, "our results/(2 sin^2 theta_W)/spinorweight",ewampLO_ggzz_fullmass(h1,h2,h3,h5)/2d0/sinW2**2/spinorweight(h1,h2,h3,h5)*2d0
!       print *, "our results/(2 sin^2 theta_W)/spinorweight",ewampLO_ggzz_summed(h1,h2,h3,h5)/2d0/sinW2**2/spinorweight(h1,h2,h3,h5
!       print *, "MCFM/spinorweight",MCFMtop(h1,h2,h3,h5)
!       print *, "ratio",ewampLO_ggzz_summed(h1,h2,h3,h5)/2d0/sinW2**2/spinorweight(h1,h2,h3,h5)/MCFMtop(h1,h2,h3,h5)

!    enddo
!    enddo
!    enddo
!    enddo


! --Overall  normalization
!    ewampLO_ggzz_eps = ewampLO_ggZZ_eps * two * (gwsq * sinW2)**2
    ewampLO_ggzz_heft_eps  = ewampLO_ggZZ_heft_eps  * gwsq**2
    ewampLO_ggzz_fullmass = ewampLO_ggZZ_fullmass * gwsq**2
!-- given the normalization of the MCFM massive ampl, we need an extra factor of 2
    ewampLO_ggzz_fullmass = ewampLO_ggZZ_fullmass * 2.0_dp


  end subroutine getewampsLO_ggzz_mass_eps

end module mod_ampLO_ggzz_mass


