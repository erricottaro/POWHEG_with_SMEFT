! -- massless and massive form factors for two triangles linked by an offshell gluon
! -- massive FF is expanded about mass->infty
! -- called in mod_DoubleTri.f90
! -- R.R. Nov. 2015

    subroutine F2mless(k1sq,k2sq,k3sq,F2)
      use mod_types; use mod_consts_dp
!-- Form factor for massless triangle
!-- Taken from Eq. 2.14 of Hagiwara, Kuruma, Yamada, Nucl. Phys. B358 (1991),80-96
!-- excl. factor of 1/pi^2
       implicit none
       real(dp),intent(in)     :: k1sq,k2sq,k3sq
       complex(dp),intent(out) :: F2
       real(dp)                :: k1Dk2

        k1Dk2=(k3sq-k1sq-k2sq)/2.0_dp
        F2=1.0_dp/2.0_dp/k1Dk2 * (1.0_dp/2.0_dp + k3sq/4.0_dp/k1Dk2*log(k2sq/k3sq))

    end subroutine F2mless

    subroutine F2mass(k1sq,k2sq,k3sq,mass,F2m)
      use mod_types; use mod_consts_dp
!-- Form factor for massive triangle, expanded about mass->infty
       real(dp),intent(in)     :: k1sq,k2sq,k3sq,mass
!       complex(dp),intent(out) :: F2
       complex(dp),intent(out)  :: F2m(0:5)
       integer                 :: nheft

       F2m(0)=czero

       F2m(1) = -1.0_dp/6.0_dp

       F2m(2) = (-2.0_dp*k2sq - 2.0_dp*k3sq)/30.0_dp + (k2sq + 2.0_dp*k3sq)/90.0_dp + (3.0_dp*k2sq + 3.0_dp*k3sq)/90.0_dp

       F2m(3) = (-3.0_dp*k2sq**2 - 3.0_dp*k2sq*k3sq - 3.0_dp*k3sq**2)/315.0_dp +&
            (6.0_dp*k2sq**2 + 6.0_dp*k2sq*k3sq + 6.0_dp*k3sq**2)/1260.0_dp + (3.0_dp*k2sq**2 + 6.0_dp*k2sq*k3sq + 9.0_dp*k3sq**2)/2520.0_dp

       F2m(4) = (-4.0_dp*k2sq**3 - 4.0_dp*k2sq**2*k3sq - 4.0_dp*k2sq*k3sq**2 - 4.0_dp*k3sq**3)/2520.0_dp +&
             (2.0_dp*k2sq**3 + 4.0_dp*k2sq**2*k3sq + 6.0_dp*k2sq*k3sq**2 + 8.0_dp*k3sq**3)/12600.0_dp +&
            -  (10.0_dp*k2sq**3 + 10.0_dp*k2sq**2*k3sq + 10.0_dp*k2sq*k3sq**2 + 10.0_dp*k3sq**3)/12600.0_dp

       F2m(5) = (-10.0_dp*k2sq**4 - 10.0_dp*k2sq**3*k3sq - 10.0_dp*k2sq**2*k3sq**2 - 10.0_dp*k2sq*k3sq**3 - 10.0_dp*k3sq**4)/34650.0_dp +&
            -  (15.0_dp*k2sq**4 + 15.0_dp*k2sq**3*k3sq + 15.0_dp*k2sq**2*k3sq**2 + 15.0_dp*k2sq*k3sq**3 + 15.0_dp*k3sq**4)/103950.0_dp +&
            (5.0_dp*k2sq**4 + 10.0_dp*k2sq**3*k3sq + 15.0_dp*k2sq**2*k3sq**2 + 20.0_dp*k2sq*k3sq**3 + 25.0_dp*k3sq**4)/207900.0_dp

!     F2=czero
!     do nheft=1,5
!        F2 = F2 + F2m(nheft)*mass**(-2*nheft)
!     enddo

! include the mass factors
     do nheft=0,5
        F2m(nheft)=F2m(nheft)*mass**(-2*nheft)
     enddo

    end subroutine F2mass
