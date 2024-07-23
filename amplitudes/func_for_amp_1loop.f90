module mod_func_for_amp_1loop
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: Lnrat,L0,L1,Lsm1_2me,Lsm1_2mh,I3m,dilog2

contains

  function klog(x)
    real(dp), intent(in)  :: x
    complex(dp) :: klog

    if (x.gt.zero) then
       klog = log(x)
    else
       klog = log(abs(x)) - ci*pi
    endif

  end function klog

  function Lnrat(x,y)
    real(dp), intent(in) :: x,y
    complex(dp) :: Lnrat

    Lnrat = klog(x) - klog(y)

  end function Lnrat

  function L0(r1,r2)
    implicit none
    real(dp), intent(in) :: r1,r2
    complex(dp) :: L0

    L0 = Lnrat(r1,r2)/(one-r1/r2)

  end function L0

  function L1(r1,r2)
    implicit none
    real(dp), intent(in) :: r1,r2
    complex(dp) :: L1

    L1 = (L0(r1,r2)+one)/(one-r1/r2)

  end function L1

  function dilog2(xIn)
    implicit none
    real(dp), intent(in):: xIn
    real(dp):: x,z,z2
    complex(dp):: dilog2, Li2tmp,Const
    real(dp):: Fact
    real(dp), parameter:: Pi26=1.64493406684822643647241516665_dp
    real(dp), parameter:: Pi23=3.28986813369645287294483033329_dp
    real(dp), parameter:: B1= -0.25_dp
    real(dp), parameter:: B2=  2.7777777777777778d-02
    real(dp), parameter:: B4= -2.7777777777777778d-04
    real(dp), parameter:: B6=  4.7241118669690098d-06
    real(dp), parameter:: B8= -9.1857730746619635d-08
    real(dp), parameter:: B10= 1.8978869988970999d-09
    real(dp), parameter:: B12=-4.0647616451442255d-11
    real(dp), parameter:: B14= 8.9216910204564525d-13
    real(dp), parameter:: B16=-1.9939295860721075d-14
    real(dp), parameter:: B18= 4.5189800296199181d-16

    Const=cmplx(0.0_dp,0.0_dp)
    Fact =1.0_dp

    x = xIn
    if ( x.gt.1.0_dp.and.x.lt.2.0_dp  ) then
       Fact = -1.0_dp
       Const = Pi23-0.5_dp*log(x)**2-pi*log(x)*ci
       x = 1.0_dp/x
       write(6,*) 'error in dilog2'
       stop
    elseif(x.gt.2.0_dp) then
       Fact = -1.0_dp
       Const = Pi23-0.5_dp*log(x)**2-pi*log(x)*ci
       x = 1.0_dp/x
    elseif ( x.lt.(-1.0_dp) ) then
       Fact = -1.0_dp
       Const =-Pi26-0.5_dp*log(-x)**2
       x = 1.0_dp/x
    elseif ( x.eq.1.0_dp) then
       dilog2 = Pi26
       return
    endif

    if ( x.gt.0.5_dp ) then
       Fact = -1.0_dp*Fact
       Const = Const + Pi26 - log(x)*log(1.0_dp-x)
       x = 1.0_dp-x
    endif

    z = -log(1.0_dp-x)
    Li2tmp = z
    z2 = z*z
    Li2tmp = Li2tmp + B1  * z2
    z = z*z2
    Li2tmp = Li2tmp + B2  * z
    z = z*z2
    Li2tmp = Li2tmp + B4  * z
    z = z*z2
    Li2tmp = Li2tmp + B6  * z
    z = z*z2
    Li2tmp = Li2tmp + B8  * z
    z = z*z2
    Li2tmp = Li2tmp + B10 * z
    z = z*z2
    Li2tmp = Li2tmp + B12 * z
    z = z*z2
    Li2tmp = Li2tmp + B14 * z
    z = z*z2
    Li2tmp = Li2tmp + B16 * z
    z = z*z2
    Li2tmp = Li2tmp + B18 * z

    dilog2 = Fact*Li2tmp + Const

    return

  end function dilog2

  function kdilog(r1,r2)  ! this is really just dilog(1-r1/r2) but with
    implicit none         ! particular imaginary part
    real(dp), intent(in) :: r1, r2
    real(dp) :: x
    complex(dp) :: kdilog

    x = r1/r2

    if (x.gt.zero) then
       kdilog = dilog2(one-x)
    else
       kdilog = pisq/6.0_dp - dilog2(x) &
            - log(one-x)*Lnrat(r1,r2)
    endif

  end function kdilog

  function Lsm1_2me(s,t,m1sq,m3sq)
    implicit none
    complex(dp) :: Lsm1_2me
    !---- formula taken from G.~Duplancic and B~Nizic [arXiv:hep-ph/0006249 v2]
    !---- Eq 71
    !---- Lsm1_2me notation follows from
    !----  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
    !----  %``Dimensionally regulated pentagon integrals,''
    !----  Nucl.\ Phys.\ B {\bf 412}, 751 (1994)
    !----  [arXiv:hep-ph/9306240].
    !----  %%CITATION = HEP-PH 9306240;%%
    !----  Eqs. (I.13)
    !---- analytic continuation has been checked by calculating numerically.
    integer j
    real(dp) :: s,t,m1sq,m3sq,arg(4),omarg(4),f2me,htheta
    complex(dp) :: Li2(4),wlog(4)

    !--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
    htheta(s)=half+half*sign(one,s)

    f2me=(s+t-m1sq-m3sq)/(s*t-m1sq*m3sq)

    arg(1)=f2me*s
    arg(2)=f2me*t
    arg(3)=f2me*m1sq
    arg(4)=f2me*m3sq

    do j=1,4
       omarg(j)=one-arg(j)
       wlog(j)=log(abs(arg(j))) &
            +(ci*pi)*cmplx(htheta(-arg(j))*sign(one,f2me))
       if (omarg(j) .gt. one) then
          Li2(j)=cmplx(zeta2-dilog2(arg(j))) &
               -wlog(j)*cmplx(log(omarg(j)))
       else
          Li2(j)=cmplx(dilog2(omarg(j)))
       endif
    enddo
    Lsm1_2me=Li2(1)+Li2(2)-Li2(3)-Li2(4)

    return
  end function Lsm1_2me


  !-- this one has problem with the imaginary part
  function Lsm1_2me_OLD(s,t,m12,m22)
    implicit none
    complex(dp) :: Lsm1_2me_OLD
    real(dp) :: s,t,m12,m22

    Lsm1_2me_OLD = -kdilog(-m12,-s)-kdilog(-m12,-t)-kdilog(-m22,-s)-kdilog(-m22,-t) &
         + kdilog(-m12*m22,-s*t) - half*Lnrat(-s,-t)**2

  end function Lsm1_2me_OLD

  function Lsm1_2mh(s,t,m12,m22)
    implicit none
    complex(dp) :: Lsm1_2mh
    real(dp) :: s,t,m12,m22

    Lsm1_2mh = -kdilog(-m12,-t)-kdilog(-m22,-t)-half*Lnrat(-s,-t)**2  &
         + half*Lnrat(-s,-m12)*Lnrat(-s,-m22)  &
         + (half*(s-m12-m22)+m12*m22/t)*I3m(s,m12,m22)

  end function Lsm1_2mh

  function I3m(s1,s2,s3)
    !     This is the function I3m, a massless triangle with all three external
    !     lines offshell defined in BDK
    !     %\cite{Bern:1997sc}
    !     \bibitem{Bern:1997sc}
    !     Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
    !     %``One-loop amplitudes for e+ e- to four partons,''
    !     Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
    !     [arXiv:hep-ph/9708239].
    !     %%CITATION = HEP-PH 9708239;%%
    !     defined in their equation II.9
    !     \int da_1 da_2 da_3 /(-a_1*a_2*s1-a_2*a_3*s2-a_3*a_1*s3)
    implicit none
    real(dp), intent(in)  :: s1,s2,s3
    real(dp) :: smax,smid,smin,del3,rtdel3
    complex(dp) :: I3m
    real(dp) ::  flag

    smax=max(s1,s2,s3)
    smin=min(s1,s2,s3)
    smid=s1+s2+s3-smax-smin
    del3=s1**2+s2**2+s3**2-two*(s1*s2+s2*s3+s3*s1)

    if (del3 .gt. zero) then
       rtdel3=sqrt(del3)
       if (smax .lt. zero) then
          !---case all negative
          flag=zero
          i3m=i3m1b(smax,smid,smin,rtdel3,flag)
       elseif (smin .gt. zero) then
          !---case all positive
          flag=zero
          i3m=-i3m1b(-smin,-smid,-smax,rtdel3,flag)
       elseif ((smid .lt. zero) .and. (smin .lt. 0)) then
          !---case two negative and one positive
          flag=one
          i3m=i3m1b(smin,smid,smax,rtdel3,flag)
       elseif ((smax .gt. zero).and.(smid .gt. zero)) then
          !---case two positive and one negative
          flag=-one
          i3m=-i3m1b(-smax,-smid,-smin,rtdel3,flag)
       endif
    elseif (del3 .lt. zero) then
       rtdel3=sqrt(-del3)
       if (smax .lt. zero) then
          !---case all negative
          i3m=cmplx(i3m1a(+s1,+s2,+s3,rtdel3),kind=dp)
       elseif (smin .gt. zero) then
          !---case all positive
          i3m=-cmplx(i3m1a(-s1,-s2,-s3,rtdel3),kind=dp)
       endif
    endif

    return
  end  function I3m


  function I3m1a(s1,s2,s3,rtmdel)
    implicit none
    !     symmetric form of Lu and Perez
    !     %\cite{Lu:1992ny}
    !     \bibitem{Lu:1992ny}
    !     H.~J.~Lu and C.~A.~Perez,
    !     %``Massless one loop scalar three point integral and associated Clausen,
    !     %Glaisher and L functions,''
    !     SLAC-PUB-5809
    real(dp) :: I3m1a
    real(dp) :: s1,s2,s3,rtmdel
    real(dp) :: d1,d2,d3,arg1,arg2,arg3

    d1=s1-s2-s3
    d2=s2-s3-s1
    d3=s3-s1-s2

    arg1=two*atan(rtmdel/d1)
    arg2=two*atan(rtmdel/d2)
    arg3=two*atan(rtmdel/d3)
    i3m1a=two/rtmdel*(clausen(arg1)+clausen(arg2)+clausen(arg3))

  end function I3m1a


  function I3m1b(s1,s2,s3,rtdel,flag)
    implicit none
    complex(dp) :: I3m1b
    !     form of Ussyukina and Davydychev
    ! %\cite{Usyukina:1994iw}
    ! \bibitem{Usyukina:1994iw}
    !   N.~I.~Usyukina and A.~I.~Davydychev,
    !   %``New results for two loop off-shell three point diagrams,''
    !  Phys.\ Lett.\ B {\bf 332}, 159 (1994)
    !  [arXiv:hep-ph/9402223].
    !  %%CITATION = HEP-PH 9402223;%%
    real(dp), intent(in) ::  s1,s2,s3,rtdel,flag
    real(dp) ::  x,y,rho,argx,argy,argdlx,argdly,d3,temp,xlog,ylog,rat
    d3=s3-s1-s2
    x=s1/s3
    y=s2/s3
    rat=half*(d3+rtdel)/s3
    if (abs(rat) .lt. 1d-3) rat=two*s1*s2/(s3*(d3-rtdel))
    rho=one/rat
    argx=rho*x
    argy=rho*y
    argdlx=-argx
    argdly=-argy

    if ((argdlx .gt. 1.0_dp) .or. (argdly .gt. 1.0_dp)) then
       write(6,*) 'problems with call of I3m1b'
       write(6,*) 'argdlx',argdlx
       write(6,*) 'argdly',argdly
       stop
    endif

    xlog=log(abs(argx))
    ylog=log(abs(argy))
    temp=xlog*ylog+pisq/three+(ylog-xlog)*log((one+argy)/(one+argx)) &
         +two*(dilog2(argdlx)+dilog2(argdly))
    I3m1b=cmplx(temp-abs(flag)*pisq,kind=dp) &
         +ci*pi*cmplx(flag*(xlog+ylog),kind=dp)
    I3m1b=-I3m1b/rtdel

  end function I3m1b

  function clausen(theta)
    complex(dp15) :: xspenz
    real(dp), intent(in) :: theta
    real(dp) :: theta1, z
    real(dp) :: claus
    real(dp) :: clausen

    ! Statement function
    claus(z) = imag(xspenz(cmplx(exp(ci*z),kind=dp15)))

    if (theta.gt.two*pi/three) then
       theta1 = pi - theta
       clausen = claus(theta1)-half*claus(two*theta1)
    else
       clausen = claus(theta)
    endif

    return

  end function clausen

end module mod_func_for_amp_1loop
