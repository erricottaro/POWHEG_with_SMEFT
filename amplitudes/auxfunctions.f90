
module mod_auxfunctions
  use mod_types; use mod_consts_dp; use common_def
  !  use mod_datacuts

  implicit none
  private

  public :: aa, bb, sprodu, spinoru, spinorur, scr, sc3, sc, &
       boost, boostz, NF1, get_NaN, delta2, give2vect, delta3, give3vect, &
       give3to4vect, give2to4vect, give1to4vect, cyclicswap, cyclicswap2

  public :: rotate, rotate_xy, rotate_yz,  rotate_xy_c, rotate_yz_c

  public :: awayfromzaxis, getPolVectors, eps, converttoMCFMmom

contains

  subroutine awayfromzaxis(n,pin,pout)
    integer, intent(in) :: n
    real(dp), intent(in) :: pin(4,n)
    real(dp), intent(out) :: pout(4,n)

    pout(1,:) = pin(1,:)
    pout(2,:) = pin(4,:)
    pout(3,:) = pin(2,:)
    pout(4,:) = pin(3,:)

    return

  end subroutine awayfromzaxis

  subroutine backtozaxis(epsin,epsout)
    complex(dp), intent(in) :: epsin(0:3)
    complex(dp), intent(out) :: epsout(0:3)
    complex(dp) :: tmp(0:3)

    tmp=epsin

    epsout(0) = tmp(0)
    epsout(3) = tmp(1)
    epsout(1) = tmp(2)
    epsout(2) = tmp(3)

    return

  end subroutine backtozaxis


  subroutine getPolVectors(p, eps1, eps2)
    real(dp), intent(in) :: p(4,6)
    complex(dp), intent(out) :: eps1(0:3, -1:1), eps2(0:3, -1:1)
    complex(dp) :: tmp(0:3)
    real(dp) :: pdum(4,6)

    call awayfromzaxis(6,p,pdum)
    ! get back to traditional axes
    call backtozaxis(eps(-1,pdum(:,1),pdum(:,2)),eps1(:,-1))
    call backtozaxis(eps(1,pdum(:,1),pdum(:,2)),eps1(:,1))
    call backtozaxis(eps(-1,pdum(:,2),pdum(:,1)),eps2(:,-1))
    call backtozaxis(eps(1,pdum(:,2),pdum(:,1)),eps2(:,1))

    return
  end subroutine getPolVectors


  subroutine give1to4vect(k1,v)
    implicit none
    integer i
    complex(dp), intent(in)  ::  k1(4)
    complex(dp) ::  v1(4),v2(4),v3(4),v4(4)
    complex(dp) ::  vaux(4),vaux1(4)
    complex(dp) ::  v2aux(4)
    complex(dp) ::  v(4,4)
    complex(dp) ::  d2,d1,ko1,ko2,ko3
    complex(dp) ::  a11,a12,a22,dnorm
    complex(dp) ::  a1v,a2v,a3v,a33



    d1 = sc(k1,k1)


    v1=k1/d1

    v2aux = (/1.3_dp, 1.7_dp, 2.3_dp, 3.4_dp/)

    !          v2aux = (/2.1_dp, 1.8_dp, 0.45_dp, 0.56_dp/)

    !        v2aux = (/1.3_dp, 1.7_dp,2.4_dp,3.5_dp/)

    a1v = sc(v2aux,v1)

    v2aux =v2aux -a1v*k1


    dnorm =  sc(v2aux,v2aux)
    dnorm=sqrt(dnorm)


    v2  = v2aux/dnorm


    a12 =  sc(v1,v2)
    a11 =  sc(v1,v1)
    a22 =  sc(v2,v2)

    d2=a12**2 - a11*a22
    ko1= a22/d2
    ko2= a11/d2
    ko3=-a12/d2

    vaux = (/1.5_dp, 13.2_dp, 14.1_dp, 9.1_dp/)


    a1v = sc(v1,vaux)
    a2v = sc(v2,vaux)


    vaux=vaux+ko1*a1v*v1+ko2*a2v*v2 + ko3*(a1v*v2+a2v*v1)


    dnorm =  sc(vaux,vaux)
    dnorm =  sqrt(dnorm)



    v3 = vaux/dnorm

    vaux1 = (/2.4_dp, 1.3_dp, 3.45_dp, 0.56_dp/)


    a1v = sc(v1,vaux1)
    a2v = sc(v2,vaux1)


    vaux =vaux1 + ko1*a1v*v1 + ko2*a2v*v2  + ko3*(a1v*v2  + a2v*v1 )


    a3v = sc(vaux,v3)
    a33 = sc(v3,v3)

    vaux = vaux  - a3v/a33*v3


    dnorm = sc(vaux,vaux)
    dnorm=sqrt(dnorm)


    v4 = vaux/dnorm


    v(1,1:4) = v1(1:4)
    v(2,1:4) = v2(1:4)
    v(3,1:4) = v3(1:4)
    v(4,1:4) = v4(1:4)


    return
  end subroutine give1to4vect



  subroutine give2to4vect(k1,k2,v)
    implicit none
    integer i
    complex(dp), intent(in) ::  k1(4),k2(4)
    complex(dp), intent(out) :: v(4,4)
    complex(dp) ::  v1(4),v2(4),v3(4),v4(4)
    complex(dp) ::  vaux(4), d2,ko1,ko2,ko3, a11,a12,a22,dnorm, a1v,a2v,a3v,a33, vaux1(4)


    call delta2(k1,k2,d2)

    call give2vect(k2,k1,v1)
    call give2vect(k1,k2,v2)

    v1 = v1/d2
    v2 = v2/d2

    a11 = sc(k1,k1)
    a12 = sc(k1,k2)
    a22 = sc(k2,k2)

    d2=a12**2-a11*a22
    ko1=a22/d2
    ko2=a11/d2
    ko3=-a12/d2

    !        vaux = (/1.3_dp, 1.7_dp,2.4_dp,3.5_dp/)

    vaux = (/2.3_dp, 1.7_dp,1.4_dp,0.5_dp/)

    a1v = sc(k1,vaux)
    a2v = sc(k2,vaux)


    vaux = vaux + ko1*a1v*k1  +ko2*a2v*k2  + ko3*(a1v*k2 +a2v*k1 )


    dnorm = sc(vaux,vaux)
    dnorm = sqrt(dnorm)

    v3 =vaux/dnorm


    !       obtaining the forth vector

    !        vaux = (/2.1_dp, 1.2_dp, 3.4_dp, 0.5_dp/)

    vaux = (/100.1_dp, 1.2_dp, 100.4_dp, 0.5_dp/)


    a1v = sc(k1,vaux)
    a2v = sc(k2,vaux)


    vaux =vaux +ko1*a1v*k1 +ko2*a2v*k2 + ko3*(a1v*k2 + a2v*k1 )


    a3v = sc(vaux,v3)
    a33 = sc(v3,v3)

    vaux = vaux  - a3v/a33*v3

    dnorm = sc(vaux,vaux)
    dnorm= sqrt(dnorm)

    v4  = vaux/dnorm

    v(1,:) = v1
    v(2,:) = v2
    v(3,:) = v3
    v(4,:) = v4


    return
  end subroutine give2to4vect




  subroutine give3to4vect(k1,k2,k3,v)
    implicit none
    complex(dp), intent(in)  ::  k1(4),k2(4),k3(4)
    complex(dp), intent(out)  ::  v(4,4)
    complex(dp) ::  v1(4),v2(4),v3(4),v4(4)
    complex(dp) ::  vaux(4),vauxa(4),v3a(4)
    complex(dp) ::  d3,d2,ko1,ko2,ko3
    complex(dp) ::  a11,a12,a22,a13,a23,dnorm,a2v,a1v
    complex(dp) ::  a3v
    integer i

    call delta3(k1,k2,k3,d3)
    call give3vect(k2,k3,k1,v1)
    call give3vect(k1,k3,k2,v2)
    call give3vect(k1,k2,k3,v3)


    v1 = v1/d3
    v2 = v2/d3
    v3 = v3/d3

    !       I need to construct a vector in a plane transverse to k1,k1,k3
    !       To do that, I first create the metric tensor of the plane, transverse to (k1,k2):
    !       this is given by three coefficients:

    a12 = sc(k1,k2)
    a11 = sc(k1,k1)
    a22 = sc(k2,k2)

    a13 = sc(k1,k3)
    a23 = sc(k2,k3)

    d2=a12**2 -a11*a22
    ko1=a22/d2
    ko2=a11/d2
    ko3=-a12/d2

    !        after that, I take a k3 vector an convolute it with that metric

    vaux =k3 +ko1*a13*k1 +ko2*a23*k2 + ko3*(a13*k2  +a23*k1)

    dnorm =  sc(vaux,vaux)
    dnorm=sqrt(dnorm)

    v3a=vaux/dnorm


    vaux(1) =0.5_dp
    vaux(2)= 2.7_dp
    vaux(3)= 3.2_dp
    vaux(4)= 4.3_dp

    a1v =  sc(k1,vaux)
    a2v =  sc(k2,vaux)


    vauxa=vaux+ko1*a1v*k1+ko2*a2v*k2 + ko3*(a1v*k2+a2v*k1)

    a3v = sc(vauxa,v3a)

    vaux = vauxa - a3v*v3a

    dnorm = sqrt(sc(vaux,vaux))

    v4=vaux/dnorm


    v(1,1:4) = v1(1:4)
    v(2,1:4) = v2(1:4)
    v(3,1:4) = v3(1:4)
    v(4,1:4) = v4(1:4)


    return
  end subroutine give3to4vect



  subroutine delta2(k1,k2,d2)
    implicit none
    complex(dp), intent(in) ::  k1(4), k2(4)
    complex(dp), intent(out) :: d2
    complex(dp) :: tk1(4), tk2(4)
    integer i

    tk1 = -k1
    tk2 = -k2

    tk1(1) = k1(1)
    tk2(1) = k2(1)


    d2=-tk1(3)*tk2(1)*k1(1)*k2(3) -tk1(3)*tk2(4)*k1(4)*k2(3) &
         -tk1(1)*tk2(2)*k1(2)*k2(1) -tk1(3)*tk2(2)*k1(2)*k2(3) &
         -tk1(2)*tk2(4)*k1(4)*k2(2) -tk1(2)*tk2(1)*k1(1)*k2(2) &
         -tk1(2)*tk2(3)*k1(3)*k2(2) -tk1(4)*tk2(1)*k1(1)*k2(4) &
         -tk1(4)*tk2(2)*k1(2)*k2(4) &
         +tk1(4)*tk2(1)*k1(4)*k2(1) &
         -tk1(1)*tk2(3)*k1(3)*k2(1) &
         -tk1(1)*tk2(4)*k1(4)*k2(1) &
         -tk1(4)*tk2(3)*k1(3)*k2(4) &
         +tk1(2)*tk2(1)*k1(2)*k2(1) &
         +tk1(3)*tk2(4)*k1(3)*k2(4) &
         +tk1(3)*tk2(1)*k1(3)*k2(1) &
         +tk1(1)*tk2(3)*k1(1)*k2(3)   &
         +tk1(4)*tk2(3)*k1(4)*k2(3)  &
         +tk1(2)*tk2(4)*k1(2)*k2(4)  &
         +tk1(3)*tk2(2)*k1(3)*k2(2)  &
         +tk1(4)*tk2(2)*k1(4)*k2(2)  &
         +tk1(1)*tk2(2)*k1(1)*k2(2) &
         +tk1(1)*tk2(4)*k1(1)*k2(4) &
         +tk1(2)*tk2(3)*k1(2)*k2(3)


    return
  end subroutine delta2



  subroutine give2vect(k2,k1,v)
    implicit none
    complex(dp), intent(in)  ::  k1(4), k2(4)
    complex(dp)  ::   tk2(4)
    complex(dp), intent(out) ::   v(4)
    integer :: i


    tk2 = -k2
    tk2(1) = k2(1)



    v(1)=-tk2(3)*k1(3)*k2(1)-tk2(4)*k1(4)*k2(1) &
         -tk2(2)*k1(2)*k2(1)+tk2(2)*k1(1)*k2(2) &
         +tk2(3)*k1(1)*k2(3)+tk2(4)*k1(1)*k2(4)

    v(2)=-tk2(1)*k1(1)*k2(2)-tk2(4)*k1(4)*k2(2)  &
         -tk2(3)*k1(3)*k2(2)+tk2(4)*k1(2)*k2(4) &
         +tk2(3)*k1(2)*k2(3)+tk2(1)*k1(2)*k2(1)

    v(3)=-tk2(2)*k1(2)*k2(3)-tk2(1)*k1(1)*k2(3) &
         -tk2(4)*k1(4)*k2(3)+tk2(1)*k1(3)*k2(1) &
         +tk2(2)*k1(3)*k2(2)+tk2(4)*k1(3)*k2(4)

    v(4)= tk2(1)*k1(4)*k2(1)-tk2(1)*k1(1)*k2(4) &
         +tk2(3)*k1(4)*k2(3)-tk2(3)*k1(3)*k2(4)  &
         +tk2(2)*k1(4)*k2(2)-tk2(2)*k1(2)*k2(4)

    return
  end subroutine give2vect



  subroutine delta3(k1,k2,k3,d3)
    implicit none
    complex(dp), intent(in)  ::  k1(4), k2(4), k3(4)
    complex(dp), intent(out) :: d3
    complex(dp) ::  tk1(4), tk2(4), tk3(4)
    complex(dp) ::  s3,s2,s1
    complex(dp) ::  k11,k12,k13,k14
    complex(dp) ::  k21,k22,k23,k24
    complex(dp) ::  k31,k32,k33,k34
    complex(dp) ::  tk11,tk12,tk13,tk14
    complex(dp) ::  tk21,tk22,tk23,tk24
    complex(dp) ::  tk31,tk32,tk33,tk34
    integer i


    tk1 = -k1
    tk2 = -k2
    tk3 = -k3

    tk1(1) = k1(1)
    tk2(1) = k2(1)
    tk3(1) = k3(1)

    k11=k1(1)
    k12=k1(2)
    k13=k1(3)
    k14=k1(4)

    k21=k2(1)
    k22=k2(2)
    k23=k2(3)
    k24=k2(4)


    k31=k3(1)
    k32=k3(2)
    k33=k3(3)
    k34=k3(4)


    tk11=tk1(1)
    tk12=tk1(2)
    tk13=tk1(3)
    tk14=tk1(4)

    tk21=tk2(1)
    tk22=tk2(2)
    tk23=tk2(3)
    tk24=tk2(4)


    tk31=tk3(1)
    tk32=tk3(2)
    tk33=tk3(3)
    tk34=tk3(4)

    s3 = tk12*tk23*tk31*k11*k22*k33+tk12*tk21*tk33*k13*k22*k31-tk13*tk21*tk34*k11*k23*k34 &
         +tk13*tk21*tk32*k13*k21*k32+tk11*tk22*tk34*k11*k22*k34+tk14*tk23*tk32*k13*k22*k34 &
         +tk13*tk24*tk31*k14*k21*k33+tk12*tk24*tk31*k11*k22*k34+tk11*tk24*tk32*k11*k24*k32 &
         +tk14*tk22*tk33*k12*k23*k34-tk12*tk24*tk33*k14*k22*k33-tk13*tk22*tk31*k11*k22*k33 &
         -tk13*tk22*tk34*k14*k22*k33-tk14*tk22*tk31*k12*k24*k31-tk12*tk24*tk33*k13*k24*k32 &
         +tk14*tk22*tk31*k12*k21*k34-tk12*tk23*tk34*k14*k23*k32-tk13*tk22*tk34*k13*k24*k32

    s2 = s3-tk11*tk24*tk33*k11*k23*k34+tk13*tk22*tk34*k14*k23*k32+tk11*tk23*tk32*k12*k21*k33 &
         +tk13*tk24*tk32*k12*k23*k34+tk13*tk22*tk31*k13*k22*k31-tk11*tk24*tk33*k13*k24*k31 &
         -tk13*tk24*tk32*k13*k22*k34-tk11*tk22*tk33*k12*k21*k33+tk11*tk23*tk34*k13*k24*k31 &
         +tk14*tk21*tk32*k12*k24*k31+tk12*tk24*tk33*k12*k24*k33+tk11*tk24*tk32*k14*k22*k31  &
         +tk12*tk21*tk33*k11*k23*k32-tk12*tk24*tk31*k12*k21*k34+tk12*tk24*tk31*k12*k24*k31  &
         -tk12*tk21*tk34*k14*k21*k32-tk12*tk23*tk31*k12*k21*k33+tk14*tk21*tk32*k14*k21*k32


    s3 = s2-tk12*tk21*tk34*k12*k24*k31+tk13*tk21*tk32*k11*k22*k33-tk14*tk22*tk31*k14*k21*k32  &
         +tk13*tk21*tk34*k11*k24*k33+tk14*tk23*tk32*k14*k23*k32+tk12*tk24*tk33*k13*k22*k34        &
         -tk11*tk24*tk32*k12*k24*k31+tk12*tk24*tk33*k14*k23*k32+tk14*tk22*tk33*k13*k24*k32        &
         -tk11*tk24*tk32*k11*k22*k34+tk12*tk23*tk31*k13*k21*k32+tk13*tk24*tk32*k14*k22*k33       &
         +tk12*tk24*tk31*k14*k21*k32-tk11*tk24*tk32*k14*k21*k32+tk11*tk22*tk34*k14*k21*k32       &
         +tk12*tk21*tk34*k12*k21*k34+tk11*tk24*tk33*k13*k21*k34

    s1 = s3+tk13*tk21*tk32*k12*k23*k31-tk14*tk22*tk33*k13*k22*k34+tk12*tk21*tk33*k12*k21*k33  &
         -tk11*tk24*tk33*k14*k21*k33+tk11*tk24*tk33*k11*k24*k33+tk11*tk22*tk33*k12*k23*k31         &
         +tk13*tk21*tk34*k13*k21*k34+tk14*tk21*tk32*k11*k22*k34+tk14*tk23*tk31*k11*k24*k33        &
         +tk14*tk23*tk31*k13*k21*k34-tk13*tk22*tk31*k12*k23*k31-tk14*tk23*tk32*k12*k23*k34        &
         -tk14*tk23*tk31*k11*k23*k34-tk12*tk24*tk31*k14*k22*k31+tk14*tk21*tk33*k11*k23*k34    &
         -tk14*tk21*tk33*k11*k24*k33-tk14*tk22*tk33*k12*k24*k33+tk13*tk21*tk34*k14*k23*k31 &
         +tk12*tk23*tk31*k12*k23*k31

    s3 = s1+tk14*tk23*tk31*k14*k23*k31-tk14*tk22*tk31*k11*k22*k34 &
         +tk13*tk24*tk31*k13*k24*k31+tk14*tk21*tk33*k13*k24*k31-tk13*tk21*tk32*k12*k21*k33 &
         -tk14*tk23*tk31*k13*k24*k31+tk12*tk23*tk34*k13*k24*k32+tk13*tk24*tk32*k13*k24*k32 &
         +tk14*tk21*tk33*k14*k21*k33-tk14*tk21*tk32*k12*k21*k34-tk11*tk23*tk32*k13*k21*k32 &
         -tk12*tk21*tk34*k11*k22*k34-tk13*tk24*tk31*k11*k24*k33+tk13*tk22*tk34*k12*k24*k33 &
         +tk13*tk22*tk34*k13*k22*k34-tk11*tk22*tk34*k14*k22*k31+tk11*tk22*tk34*k12*k24*k31

    s2 = s3+tk14*tk23*tk32*k12*k24*k33-tk12*tk23*tk31*k11*k23*k32   &
         -tk13*tk21*tk32*k11*k23*k32+tk11*tk22*tk33*k11*k22*k33-tk11*tk22*tk33*k13*k22*k31 &
         +tk12*tk21*tk34*k11*k24*k32+tk11*tk23*tk32*k13*k22*k31-tk14*tk21*tk32*k14*k22*k31 &
         -tk11*tk23*tk32*k11*k22*k33-tk11*tk23*tk34*k14*k23*k31+tk14*tk22*tk31*k14*k22*k31 &
         -tk11*tk22*tk33*k11*k23*k32-tk12*tk21*tk33*k12*k23*k31+tk11*tk23*tk34*k11*k23*k34 &
         +tk11*tk24*tk33*k14*k23*k31-tk11*tk23*tk32*k12*k23*k31  &
         +tk11*tk22*tk33*k13*k21*k32+tk13*tk22*tk31*k12*k21*k33

    s3 = s2+tk12*tk21*tk34*k14*k22*k31-tk13*tk22*tk31*k13*k21*k32+tk13*tk22*tk31*k11*k23*k32  &
         +tk14*tk22*tk31*k11*k24*k32-tk12*tk23*tk31*k13*k22*k31-tk14*tk23*tk31*k14*k21*k33  &
         -tk12*tk24*tk31*k11*k24*k32-tk13*tk24*tk31*k14*k23*k31-tk13*tk24*tk31*k13*k21*k34  &
         +tk13*tk24*tk31*k11*k23*k34+tk11*tk23*tk32*k11*k23*k32  &
         -tk13*tk21*tk32*k13*k22*k31-tk14*tk21*tk32*k11*k24*k32-tk14*tk23*tk32*k13*k24*k32 &
         -tk14*tk23*tk32*k14*k22*k33+tk11*tk24*tk32*k12*k21*k34-tk13*tk24*tk32*k14*k23*k32 &
         -tk13*tk24*tk32*k12*k24*k33

    d3 = s3-tk12*tk21*tk33*k13*k21*k32-tk12*tk21*tk33*k11*k22*k33  &
         -tk14*tk21*tk33*k14*k23*k31-tk14*tk21*tk33*k13*k21*k34  &
         -tk14*tk22*tk33*k14*k23*k32+tk14*tk22*tk33*k14*k22*k33-tk12*tk24*tk33*k12*k23*k34 &
         -tk13*tk21*tk34*k13*k24*k31-tk13*tk21*tk34*k14*k21*k33-tk11*tk22*tk34*k11*k24*k32  &
         -tk11*tk22*tk34*k12*k21*k34-tk13*tk22*tk34*k12*k23*k34+tk11*tk23*tk34*k14*k21*k33  &
         -tk11*tk23*tk34*k11*k24*k33-tk11*tk23*tk34*k13*k21*k34+tk12*tk23*tk34*k14*k22*k33  &
         -tk12*tk23*tk34*k12*k24*k33-tk12*tk23*tk34*k13*k22*k34+tk12*tk23*tk34*k12*k23*k34


    return
  end subroutine delta3


  subroutine give3vect(k2,k3,k1,v)
    implicit none
    complex(dp), intent(in)  ::  k1(4), k2(4), k3(4)
    complex(dp)  ::  tk2(4), tk3(4)
    complex(dp), intent(out)  ::  v(4)
    integer i

    tk2 = -k2
    tk3 = -k3

    tk2(1) = k2(1)
    tk3(1) = k3(1)


    v(1)=tk2(2)*tk3(3)*k1(1)*k2(2)*k3(3) &
         +tk2(2)*tk3(4)*k1(2)*k2(4)*k3(1)    &
         -tk2(4)*tk3(2)*k1(1)*k2(2)*k3(4)  &
         +tk2(2)*tk3(3)*k1(3)*k2(1)*k3(2) &
         -tk2(3)*tk3(2)*k1(1)*k2(2)*k3(3) &
         -tk2(2)*tk3(4)*k1(4)*k2(2)*k3(1) &
         -tk2(3)*tk3(2)*k1(3)*k2(1)*k3(2) &
         +tk2(4)*tk3(2)*k1(4)*k2(2)*k3(1) &
         -tk2(2)*tk3(4)*k1(1)*k2(4)*k3(2) &
         +tk2(3)*tk3(4)*k1(3)*k2(4)*k3(1) &
         -tk2(2)*tk3(3)*k1(3)*k2(2)*k3(1) &
         -tk2(3)*tk3(4)*k1(3)*k2(1)*k3(4) &
         +tk2(4)*tk3(2)*k1(2)*k2(1)*k3(4) &
         +tk2(3)*tk3(4)*k1(1)*k2(3)*k3(4) &
         -tk2(2)*tk3(3)*k1(1)*k2(3)*k3(2) &
         +tk2(3)*tk3(4)*k1(4)*k2(1)*k3(3) &
         -tk2(4)*tk3(3)*k1(4)*k2(1)*k3(3) &
         -tk2(4)*tk3(2)*k1(2)*k2(4)*k3(1) &
         -tk2(2)*tk3(4)*k1(2)*k2(1)*k3(4) &
         -tk2(3)*tk3(4)*k1(1)*k2(4)*k3(3) &
         -tk2(4)*tk3(2)*k1(4)*k2(1)*k3(2) &
         +tk2(3)*tk3(2)*k1(2)*k2(1)*k3(3) &
         -tk2(3)*tk3(4)*k1(4)*k2(3)*k3(1) &
         -tk2(4)*tk3(3)*k1(1)*k2(3)*k3(4) &
         -tk2(2)*tk3(3)*k1(2)*k2(1)*k3(3) &
         +tk2(4)*tk3(3)*k1(4)*k2(3)*k3(1) &
         +tk2(3)*tk3(2)*k1(1)*k2(3)*k3(2) &
         +tk2(2)*tk3(3)*k1(2)*k2(3)*k3(1) &
         -tk2(3)*tk3(2)*k1(2)*k2(3)*k3(1)  &
         +tk2(2)*tk3(4)*k1(1)*k2(2)*k3(4)  &
         +tk2(2)*tk3(4)*k1(4)*k2(1)*k3(2) &
         +tk2(4)*tk3(2)*k1(1)*k2(4)*k3(2) &
         +tk2(3)*tk3(2)*k1(3)*k2(2)*k3(1) &
         -tk2(4)*tk3(3)*k1(3)*k2(4)*k3(1) &
         +tk2(4)*tk3(3)*k1(1)*k2(4)*k3(3) &
         +tk2(4)*tk3(3)*k1(3)*k2(1)*k3(4)

    v(2)=tk2(4)*tk3(3)*k1(2)*k2(4)*k3(3)  &
         -tk2(3)*tk3(4)*k1(2)*k2(4)*k3(3)  &
         +tk2(1)*tk3(3)*k1(2)*k2(1)*k3(3) &
         +tk2(4)*tk3(3)*k1(3)*k2(2)*k3(4) &
         -tk2(3)*tk3(4)*k1(3)*k2(2)*k3(4) &
         -tk2(3)*tk3(1)*k1(2)*k2(1)*k3(3) &
         -tk2(4)*tk3(3)*k1(2)*k2(3)*k3(4) &
         +tk2(4)*tk3(1)*k1(2)*k2(4)*k3(1) &
         +tk2(4)*tk3(1)*k1(4)*k2(1)*k3(2) &
         +tk2(1)*tk3(3)*k1(1)*k2(3)*k3(2) &
         +tk2(1)*tk3(4)*k1(2)*k2(1)*k3(4) &
         -tk2(3)*tk3(1)*k1(1)*k2(3)*k3(2) &
         -tk2(4)*tk3(3)*k1(4)*k2(2)*k3(3) &
         +tk2(3)*tk3(4)*k1(4)*k2(2)*k3(3) &
         -tk2(1)*tk3(4)*k1(4)*k2(1)*k3(2) &
         +tk2(3)*tk3(1)*k1(2)*k2(3)*k3(1) &
         -tk2(1)*tk3(3)*k1(2)*k2(3)*k3(1) &
         -tk2(4)*tk3(1)*k1(4)*k2(2)*k3(1) &
         -tk2(4)*tk3(1)*k1(2)*k2(1)*k3(4) &
         -tk2(1)*tk3(3)*k1(3)*k2(1)*k3(2) &
         +tk2(3)*tk3(1)*k1(3)*k2(1)*k3(2) &
         +tk2(3)*tk3(4)*k1(3)*k2(4)*k3(2) &
         +tk2(3)*tk3(1)*k1(1)*k2(2)*k3(3) &
         -tk2(3)*tk3(1)*k1(3)*k2(2)*k3(1) &
         -tk2(1)*tk3(4)*k1(1)*k2(2)*k3(4) &
         +tk2(1)*tk3(4)*k1(4)*k2(2)*k3(1) &
         -tk2(4)*tk3(3)*k1(3)*k2(4)*k3(2) &
         +tk2(1)*tk3(3)*k1(3)*k2(2)*k3(1) &
         +tk2(1)*tk3(4)*k1(1)*k2(4)*k3(2) &
         -tk2(1)*tk3(3)*k1(1)*k2(2)*k3(3) &
         +tk2(4)*tk3(1)*k1(1)*k2(2)*k3(4) &
         +tk2(3)*tk3(4)*k1(2)*k2(3)*k3(4)  &
         -tk2(4)*tk3(1)*k1(1)*k2(4)*k3(2) &
         -tk2(3)*tk3(4)*k1(4)*k2(3)*k3(2) &
         -tk2(1)*tk3(4)*k1(2)*k2(4)*k3(1) &
         +tk2(4)*tk3(3)*k1(4)*k2(3)*k3(2)

    v(3)=tk2(2)*tk3(1)*k1(1)*k2(3)*k3(2)  &
         +tk2(2)*tk3(4)*k1(2)*k2(4)*k3(3) &
         -tk2(1)*tk3(4)*k1(1)*k2(3)*k3(4) &
         -tk2(4)*tk3(1)*k1(3)*k2(1)*k3(4) &
         -tk2(4)*tk3(2)*k1(2)*k2(4)*k3(3) &
         -tk2(1)*tk3(4)*k1(4)*k2(1)*k3(3) &
         -tk2(2)*tk3(4)*k1(2)*k2(3)*k3(4) &
         -tk2(4)*tk3(1)*k1(4)*k2(3)*k3(1) &
         +tk2(2)*tk3(4)*k1(4)*k2(3)*k3(2) &
         +tk2(2)*tk3(1)*k1(2)*k2(1)*k3(3) &
         -tk2(1)*tk3(2)*k1(2)*k2(1)*k3(3) &
         +tk2(1)*tk3(4)*k1(4)*k2(3)*k3(1)+tk2(4)*tk3(1)*k1(1)*k2(3)*k3(4) &
         +tk2(1)*tk3(4)*k1(1)*k2(4)*k3(3)+tk2(4)*tk3(1)*k1(4)*k2(1)*k3(3) &
         -tk2(1)*tk3(2)*k1(1)*k2(3)*k3(2)+tk2(1)*tk3(2)*k1(3)*k2(1)*k3(2) &
         -tk2(2)*tk3(4)*k1(3)*k2(4)*k3(2)-tk2(2)*tk3(1)*k1(3)*k2(1)*k3(2) &
         -tk2(1)*tk3(4)*k1(3)*k2(4)*k3(1)+tk2(1)*tk3(2)*k1(2)*k2(3)*k3(1) &
         +tk2(4)*tk3(1)*k1(3)*k2(4)*k3(1)-tk2(2)*tk3(4)*k1(4)*k2(2)*k3(3) &
         +tk2(1)*tk3(2)*k1(1)*k2(2)*k3(3)-tk2(2)*tk3(1)*k1(1)*k2(2)*k3(3) &
         +tk2(4)*tk3(2)*k1(2)*k2(3)*k3(4)-tk2(4)*tk3(1)*k1(1)*k2(4)*k3(3) &
         +tk2(2)*tk3(1)*k1(3)*k2(2)*k3(1)-tk2(1)*tk3(2)*k1(3)*k2(2)*k3(1) &
         -tk2(4)*tk3(2)*k1(4)*k2(3)*k3(2)+tk2(4)*tk3(2)*k1(3)*k2(4)*k3(2) &
         -tk2(2)*tk3(1)*k1(2)*k2(3)*k3(1)+tk2(2)*tk3(4)*k1(3)*k2(2)*k3(4) &
         -tk2(4)*tk3(2)*k1(3)*k2(2)*k3(4)+tk2(4)*tk3(2)*k1(4)*k2(2)*k3(3)  &
         +tk2(1)*tk3(4)*k1(3)*k2(1)*k3(4)

    v(4)=  &
         tk2(2)*tk3(1)*k1(1)*k2(4)*k3(2)+tk2(3)*tk3(1)*k1(3)*k2(1)*k3(4)  &
         -tk2(2)*tk3(3)*k1(4)*k2(3)*k3(2)+tk2(3)*tk3(2)*k1(4)*k2(3)*k3(2) &
         -tk2(1)*tk3(2)*k1(2)*k2(1)*k3(4)-tk2(1)*tk3(2)*k1(1)*k2(4)*k3(2) &
         -tk2(1)*tk3(3)*k1(3)*k2(1)*k3(4)+tk2(1)*tk3(2)*k1(2)*k2(4)*k3(1) &
         -tk2(2)*tk3(3)*k1(2)*k2(4)*k3(3)+tk2(2)*tk3(1)*k1(2)*k2(1)*k3(4) &
         +tk2(1)*tk3(3)*k1(3)*k2(4)*k3(1)-tk2(3)*tk3(1)*k1(3)*k2(4)*k3(1) &
         +tk2(3)*tk3(1)*k1(1)*k2(4)*k3(3)-tk2(1)*tk3(3)*k1(1)*k2(4)*k3(3) &
         +tk2(2)*tk3(3)*k1(2)*k2(3)*k3(4)-tk2(2)*tk3(1)*k1(2)*k2(4)*k3(1) &
         -tk2(3)*tk3(2)*k1(2)*k2(3)*k3(4)-tk2(3)*tk3(1)*k1(4)*k2(1)*k3(3) &
         +tk2(3)*tk3(2)*k1(2)*k2(4)*k3(3)+tk2(2)*tk3(3)*k1(3)*k2(4)*k3(2) &
         -tk2(1)*tk3(3)*k1(4)*k2(3)*k3(1)-tk2(3)*tk3(1)*k1(1)*k2(3)*k3(4) &
         +tk2(1)*tk3(2)*k1(4)*k2(1)*k3(2)+tk2(3)*tk3(1)*k1(4)*k2(3)*k3(1) &
         +tk2(1)*tk3(3)*k1(4)*k2(1)*k3(3)-tk2(2)*tk3(1)*k1(1)*k2(2)*k3(4) &
         +tk2(1)*tk3(3)*k1(1)*k2(3)*k3(4)-tk2(3)*tk3(2)*k1(3)*k2(4)*k3(2) &
         -tk2(3)*tk3(2)*k1(4)*k2(2)*k3(3)+tk2(2)*tk3(3)*k1(4)*k2(2)*k3(3) &
         +tk2(3)*tk3(2)*k1(3)*k2(2)*k3(4)-tk2(2)*tk3(3)*k1(3)*k2(2)*k3(4) &
         +tk2(2)*tk3(1)*k1(4)*k2(2)*k3(1)+tk2(1)*tk3(2)*k1(1)*k2(2)*k3(4) &
         -tk2(1)*tk3(2)*k1(4)*k2(2)*k3(1)-tk2(2)*tk3(1)*k1(4)*k2(1)*k3(2)

    return
  end subroutine give3vect



  subroutine rotate_xy(theta,pin,pout)
    implicit none
    real(dp), intent(in) :: theta
    real(dp), intent(in) :: pin(1:4)
    real(dp), intent(out) :: pout(1:4)


    pout = pin

    pout(2) =  cos(theta)*pin(2) + sin(theta)*pin(3)
    pout(3) = -sin(theta)*pin(2) + cos(theta)*pin(3)

  end subroutine rotate_xy


  subroutine rotate_yz(theta,pin,pout)
    implicit none
    real(dp), intent(in) :: theta
    real(dp), intent(in) :: pin(1:4)
    real(dp), intent(out) :: pout(1:4)

    pout = pin

    pout(3) =  cos(theta)*pin(3) + sin(theta)*pin(4)
    pout(4) = -sin(theta)*pin(3) + cos(theta)*pin(4)

  end subroutine rotate_yz




  subroutine rotate_xy_c(theta,pin,pout)
    implicit none
    real(dp), intent(in) :: theta
    complex(dp), intent(in) :: pin(1:4)
    complex(dp), intent(out) :: pout(1:4)


    pout = pin

    pout(2) =  cos(theta)*pin(2) + sin(theta)*pin(3)
    pout(3) = -sin(theta)*pin(2) + cos(theta)*pin(3)

  end subroutine rotate_xy_c


  subroutine rotate_yz_c(theta,pin,pout)
    implicit none
    real(dp), intent(in) :: theta
    complex(dp), intent(in) :: pin(1:4)
    complex(dp), intent(out) :: pout(1:4)

    pout = pin

    pout(3) =  cos(theta)*pin(3) + sin(theta)*pin(4)
    pout(4) = -sin(theta)*pin(3) + cos(theta)*pin(4)

  end subroutine rotate_yz_c


  !------- this subroutine will remove j component of a vector n, by rotating it
  !------- in the i-j plane and will output all momenta in the rotated reference
  !------- frame

  !------- ``k'' labels the particle whose vector is rotated

  subroutine rotate(k,i,j,pin,nin,pepin,nepin,pout,nout,pepout,nepout)
    implicit none
    integer, intent(in) :: i,j,k
    real(dp), intent(in) ::  pin(4,8),nin(3,8),pepin(8),nepin(8)
    real(dp), intent(out) :: pout(4,8),nout(3,8),pepout(8),nepout(8)
    real(dp) :: nnorm, cphi, sphi,n(3)

    pout = pin
    nout = nin
    pepout = pepin
    nepout =  nepin

    n(:) = nin(:,k)

    if (j.ne.5) then

       if (n(j).eq.zero) then
          return
       else
          nnorm = sqrt(n(i)**2 + n(j)**2)

          if (n(i).lt.zero) then
             cphi = -n(i)/nnorm
             sphi = -n(j)/nnorm
          else
             cphi = n(i)/nnorm
             sphi = n(j)/nnorm
          endif

          pout(i+1,:) =  cphi*pin(i+1,:) + sphi*pin(j+1,:)
          pout(j+1,:) =  -sphi*pin(i+1,:) + cphi*pin(j+1,:)

          nout(i,:) =  cphi*nin(i,:) + sphi*nin(j,:)
          nout(j,:) = -sphi*nin(i,:) + cphi*nin(j,:)

       endif

    elseif(j.eq.5) then

       if (nepin(k).eq.zero) then
          return
       else

          nnorm = sqrt(n(i)**2 + nepin(k)**2)

          if (n(i).lt.zero) then
             cphi = -n(i)/nnorm
             sphi = -nepin(k)/nnorm
          else
             cphi = n(i)/nnorm
             sphi = nepin(k)/nnorm
          endif

          pout(i+1,:) =  cphi*pin(i+1,:) + sphi*pepin(:)
          pepout(:) =   -sphi*pin(i+1,:) + cphi*pepin(:)

          nout(i,:) =  cphi*nin(i,:) + sphi*nepin(:)
          nepout(:) = -sphi*nin(i,:) + cphi*nepin(:)

       endif

    endif


  end subroutine rotate


  function NF1(x3,x4,x5)
    implicit none
    real(dp), intent(in) :: x3,x4,x5
    real(dp) :: NF1

    NF1 = 1.0_dp + x4*(1.0_dp-2.0_dp*x3) &
         - 2.0_dp*(1.0_dp-2.0_dp*x5)*sqrt(abs(x4*(1.0_dp-x3)*(1.0_dp-x3*x4)))

  end function NF1

  ! take momenta p_in in frame in which particle one is at rest with mass
  ! "mass" and convert to frame in which particle one has fourvector p1
  subroutine boost(mass,p1,p_in,p_out)
    use mod_types; use mod_consts_dp; use common_def
    implicit none
    real(dp), intent(in)  ::  mass,p1(4),p_in(4)
    real(dp), intent(out) ::  p_out(4)
    ! ------------------------------------
    real(dp) ::  gam,beta(2:4),bdotp
    integer :: j,k

    gam=p1(1)/mass
    bdotp=zero
    do j=2,4
       beta(j)=-p1(j)/p1(1)
       bdotp=bdotp+p_in(j)*beta(j)
    enddo
    p_out(1)=gam*(p_in(1)-bdotp)
    do k=2,4
       p_out(k)=p_in(k)+gam*beta(k)*(gam/(gam+one)*bdotp-p_in(1))
    enddo
  end subroutine boost

  !-- boost along the z-axis from the com frame to the lab frame according to
  !-- p1 = sqrt(s*xa*xb)/2 * (1,0,0,+1) --> sqrt(s)/2 * xa * (1,0,0,+1)
  !-- p2 = sqrt(s*xa*xb)/2 * (1,0,0,-1) --> sqrt(s)/2 * xb * (1,0,0,-1)
  subroutine boostz(xa,xb,npart,pin,pout)
    real(dp), intent(in) :: xa,xb
    integer, intent(in) :: npart
    real(dp), intent(in) :: pin(4,npart)
    real(dp), intent(out) :: pout(4,npart)
    real(dp) :: beta,gamma
    integer :: i

    beta = (xb-xa)/(xb+xa)
    gamma = (xa+xb)/two/sqrt(xa*xb)

    pout = pin

    do i = 1,npart
       pout(1,i) = gamma * (pin(1,i)-beta*pin(4,i))
       pout(4,i) = gamma * (pin(4,i)-beta*pin(1,i))
    enddo

    return

  end subroutine boostz


  !--- scalar products, 2*pi.pj for massless particles
  subroutine sprodu(j,p,sprod)
    implicit none
    integer, intent(in) :: j
    real(dp), intent(in) :: p(4,j)
    real(dp), intent(out) :: sprod(j,j)
    integer :: i1, i2

    sprod = zero

    do i1 = 1, j
       do i2 = i1+1, j
          sprod(i1,i2) = two * scr(p(:,i1),p(:,i2))
          sprod(i2,i1) = sprod(i1,i2)
       enddo
    enddo

    return

  end subroutine sprodu


  !--- spinor products

  subroutine spinorur(j,p,zza,zzb,ss1)
    implicit none
    integer, intent(in) :: j
    real(dp), intent(in) :: p(4,j)
    complex(dp), intent(out) :: zza(j,j),zzb(j,j)
    real(dp), intent(out) :: ss1(j,j)
    complex(dp) :: pdum(4,j),ss1dum(j,j)

    pdum = p

    call spinoru(j,pdum,zza,zzb,ss1dum)
    ss1 = real(ss1dum,kind=dp)

  end subroutine spinorur


  subroutine spinoru(j,p,zza,zzb,ss1)
    implicit none
    integer, intent(in) :: j
    complex(dp), intent(out)  :: zza(j,j),zzb(j,j)
    complex(dp), intent(out) :: ss1(j,j)
    complex(dp), intent(in) :: p(1:4,j)
    integer :: i1,i2



    zza = czero
    zzb = czero
    ss1 = czero


    do i1=1,j
       do i2=i1+1,j

          zza(i1,i2) = aa(p(:,i1),p(:,i2))
          zzb(i1,i2) = bb(p(:,i1),p(:,i2))

          zza(i2,i1) = -zza(i1,i2)
          zzb(i2,i1) = -zzb(i1,i2)

          ss1(i1,i2) = zza(i1,i2)*zzb(i2,i1)
          ss1(i2,i1) = ss1(i1,i2)


       enddo
    enddo


    return
  end subroutine spinoru


  !-- real
  subroutine spinoru1(j,p,zza,zzb,ss1)
    implicit none
    integer, intent(in) :: j
    complex(dp), intent(out)  :: zza(j,j),zzb(j,j)
    real(dp), intent(out) :: ss1(j,j)
    real(dp), intent(in) :: p(1:4,j)
    integer :: i1,i2

    zza = czero
    zzb = czero
    ss1 = zero


    do i1=1,j
       do i2=i1+1,j

          zza(i1,i2) = aar(p(:,i1),p(:,i2))
          zzb(i1,i2) = bbr(p(:,i1),p(:,i2))

          zza(i2,i1) = -zza(i1,i2)
          zzb(i2,i1) = -zzb(i1,i2)

          print *, zza(i1,i2)*zzb(i2,i1)
          ss1(i1,i2) = real(zza(i1,i2)*zzb(i2,i1),kind=dp)
          ss1(i2,i1) = ss1(i1,i2)

       enddo
    enddo


    return
  end subroutine spinoru1

  function aar(p1,p2)
    real(dp), intent(in) :: p1(4), p2(4)
    complex(dp) :: aar, eta(1:2)

    eta = one

    if (p1(1).lt.zero) eta(1) = -ci
    if (p2(1).lt.zero) eta(2) = -ci

    aar = (p1(2)+ci*p1(3))*(p2(1)+p2(4)) - (p1(1)+p1(4))*(p2(2)+ci*p2(3))
    aar = eta(1)*eta(2)*aar/sqrt( abs(p1(1)+p1(4)) )/sqrt( abs(p2(1)+p2(4)) )

  end function aar

  function bbr(p1,p2)
    real(dp), intent(in) :: p1(4), p2(4)
    complex(dp) :: bbr, eta(1:2)

    if (p1(1).lt.zero) eta(1) = ci
    if (p2(1).lt.zero) eta(2) = ci

    bbr =  (p1(1)+p1(4))*(p2(2)-ci*p2(3)) - (p1(2)-ci*p1(3))*(p2(1)+p2(4))
    bbr = eta(1)*eta(2)*bbr/sqrt(abs(p1(1)+p1(4)))/sqrt(abs(p2(1)+p2(4)))

  end function bbr

  ! -- spinor helicity routines

  ! | p > spinor (R, +)

  function spa(p)
    complex(dp), intent(in) :: p(4)
    complex(dp) :: spa(2)

    spa = (/sqrt(p(1)+p(4)),(p(2)+ci*p(3))/sqrt(p(1)+p(4))/)

  end function spa


  ! [ p | spinor conjugate +

  function bspb(p)
    complex(dp), intent(in) :: p(4)
    complex(dp) :: bspb(2)

    bspb = (/ sqrt(p(1)+p(4)), (p(2)-ci*p(3))/sqrt(p(1)+p(4))/)

  end function bspb



  ! | p ] spinor L, ``-''

  function spb(p)
    complex(dp), intent(in) :: p(4)
    complex(dp) :: spb(2)

    spb = (/sqrt(p(1)-p(4)),-(p(2)+ci*p(3))/sqrt(p(1)-p(4)) /)

  end function spb


  ! < p | spinor conjugate `-'

  function bspa(p)
    complex(dp), intent(in) :: p(4)
    complex(dp) :: bspa(2)

    bspa = (/ sqrt(p(1)-p(4)), -(p(2)-ci*p(3))/sqrt(p(1)-p(4)) /)

  end function bspa


  !--- spinor products


  ! <1 2 >

  function aa(p1,p2)
    complex(dp), intent(in) :: p1(4), p2(4)
    complex(dp) :: aa



    aa = (p1(2)+ci*p1(3))*(p2(1)+p2(4)) - (p1(1)+p1(4))*(p2(2)+ci*p2(3))
    aa = aa/sqrt( p1(1)+p1(4) )/sqrt(p2(1)+p2(4))



  end function aa


  ! [1 2]

  function bb(p1,p2)
    complex(dp), intent(in) :: p1(4), p2(4)
    complex(dp) :: bb

    bb =  (p1(1)+p1(4))*(p2(2)-ci*p2(3)) - (p1(2)-ci*p1(3))*(p2(1)+p2(4))
    bb = bb/sqrt(p1(1)+p1(4))/sqrt(p2(1)+p2(4))


  end function bb


  ! Polarization vector for gluons

  function eps(h,p,q)
    integer, intent(in) :: h
    real(dp), intent(in) :: p(4),q(4)
    real(dp) :: hl, pp, pm, qp, norm
    complex(dp) :: ppTh, qpTh, eps(4)

    if ((h.ne.1).and.(h.ne.-1)) then
       stop "hel must be 1 or -1"
    endif
    hl = real(h, kind=dp)
    pp = p(1)+p(4)
    pm = p(1)-p(4)
    ppTh = cmplx(p(2),-hl*p(3),kind=dp) ! pTb if h=1, pT if h=-1
    qp = q(1)+q(4)
    qpTh = cmplx(q(2),hl*q(3),kind=dp)  ! qT if h=1, qTb if h=-1

    eps(1) = pp*qp+ppTh*qpTh
    eps(2) = qp*ppTh+pp*qpTh
    eps(3) = -hl*ci*(-qp*ppTh+pp*qpTh)
    eps(4) = pp*qp-ppTh*qpTh

    norm = one/(sqrt2*(qp*conjg(ppTh) - pp*qpTh))

    eps = eps*norm

!    print *, eps

  end function eps

  !--- Spinors

  ! -- ubar spinor, massless
  function ubar0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: ubar0(4)
    real(dp)    :: p0,px,py,pz
    real(dp) ::  theta, phi, rrr
    complex(dp) :: cp0
    logical :: flag_nan

    flag_nan = .false.

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)


    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in ubar -> p0=0'
       !       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = py/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = asin(rrr)
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif


    call get_NaN(theta,flag_nan)
    if (flag_nan) then
       write(6,*) 'ubar-th', p0,px,py,pz
       stop
    endif


    call get_NaN(phi,flag_nan)
    if (flag_nan) then
       write(6,*) 'ubar-phi', p0,px,py,pz
       stop
    endif


    if (i.eq.1) then
       ubar0(1)=czero
       ubar0(2)=czero
       ubar0(3)=sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       ubar0(4)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
    elseif(i.eq.-1) then
       ubar0(1)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       ubar0(2)=-sqrt2*sqrt(cp0)*cmplx(abs(cos(theta/two)),0.0_dp,kind=dp)
       ubar0(3)=czero
       ubar0(4)=czero
    else
       stop 'ubar0: i out of range'
    endif

  end function ubar0



  ! -- u0  spinor, massless
  function u0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: u0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    logical :: flag_nan

    flag_nan = .false.

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in v0 -> p0=0'
       !       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (py/p0.gt.zero)  phi = pi/two
       if (py/p0.lt.zero)  phi = three*pi/two
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (px/p0.gt.zero)  phi = zero
       if (px/p0.lt.zero)  phi = pi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif

    call get_NaN(theta,flag_nan)
    if (flag_nan) then
       write(6,*) 'th-v0', p0,px,py,pz
       write(6,*) 'ratio', pz/p0, acos(pz/p0)
       stop
    endif


    call get_NaN(phi,flag_nan)
    if (flag_nan) then
       write(6,*) 'ph-v0', p0,px,py,pz, px/p0/sin(theta)
       stop
    endif

    if (i.eq.-1) then
       u0(1)=czero
       u0(2)=czero
       u0(3)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
       u0(4)=-sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
    elseif (i.eq.1) then
       u0(1)= sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       u0(2)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       u0(3)=czero
       u0(4)=czero
    else
       stop 'u0: i out of range'
    endif

  end function u0


  ! -- v0  spinor, massless
  function v0(p,i)
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: v0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    logical :: flag_nan

    flag_nan = .false.

    p0=real(p(1),kind=dp)
    px=real(p(2),kind=dp)
    py=real(p(3),kind=dp)
    pz=real(p(4),kind=dp)

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in v0 -> p0=0'
       !       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if (pz/p0.gt.0.0_dp) theta = zero
       if (pz/p0.lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (py/p0.gt.zero)  phi = pi/two
       if (py/p0.lt.zero)  phi = three*pi/two
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (px/p0.gt.zero)  phi = zero
       if (px/p0.lt.zero)  phi = pi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif

    call get_NaN(theta,flag_nan)
    if (flag_nan) then
       write(6,*) 'th-v0', p0,px,py,pz
       write(6,*) 'ratio', pz/p0, acos(pz/p0)
       stop
    endif

    call get_NaN(phi,flag_nan)
    if (flag_nan) then
       write(6,*) 'ph-v0', p0,px,py,pz, px/p0/sin(theta)
       stop
    endif

    if (i.eq.1) then
       v0(1)=czero
       v0(2)=czero
       v0(3)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
       v0(4)=-sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
    elseif (i.eq.-1) then
       v0(1)= sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       v0(2)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       v0(3)=czero
       v0(4)=czero
    else
       stop 'v0: i out of range'
    endif

  end function v0



  !---- THESE ARE POLARIZATION/OTHER ROUTINES

  ! -- massless vector polarization subroutine
  function pol_mless(p,i,outgoing)
    real(dp), intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------
    integer :: pol
    real(dp) :: p0,px,py,pz
    real(dp) :: pv,ct,st,cphi,sphi
    complex(dp) :: pol_mless(4)

    !^^^IFmp
    !    p0=(p(1)+conjg(p(1)))/two
    !    px=(p(2)+conjg(p(2)))/two
    !    py=(p(3)+conjg(p(3)))/two
    !    pz=(p(4)+conjg(p(4)))/two
    !^^^ELSE
    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)
    !^^^END


    pv=sqrt(abs(p0**2))

    if (pv.eq.zero) then
       write(6,*) 'error in pol_mless - pv=0'
       !       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       ct = 1.0_dp
       st = 0.0_dp
       cphi = 1.0_dp
       sphi = 0.0_dp
    elseif(px.eq.zero.and.py.ne.zero) then
       ct=pz/pv
       st=sqrt(abs(1.0_dp-ct**2))
       sphi = 1.0_dp
       cphi = 0.0_dp
    elseif(py.eq.zero.and.px.ne.zero) then
       ct=pz/pv
       st=sqrt(abs(1.0_dp-ct**2))
       sphi = 0.0_dp
       cphi = 1.0_dp
    else
       ct=pz/pv
       st=sqrt(abs(1.0_dp-ct**2))
       cphi= px/pv/st
       sphi= py/pv/st
    endif

    !    if (st < tolb) then
    !       cphi=1.0_dp
    !       sphi=0.0_dp
    !    else
    !       cphi= px/pv/st
    !       sphi= py/pv/st
    !    endif

    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0_dp) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing
    if (present(outgoing)) then
       if (outgoing) pol = -pol
    endif

    pol_mless(1)=czero
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2

  end function pol_mless


  function pol_mless2(p,i,out)
    integer, intent(in)       :: i
    real(dp), intent(in) :: p(4)
    character(len=*), intent(in):: out
    complex(dp)             :: pol_mless2(4)
    ! -------------------------------------

    if (out == 'out') then
       pol_mless2 = pol_mless(p,i,outgoing=.true.)
    else
       pol_mless2 = pol_mless(p,i,outgoing=.false.)
    endif
  end function pol_mless2


  !--------massive vector boson polarization routine

  function pol_mass(p,m,i,outgoing)
    integer, intent(in)       :: i
    complex(dp), intent(in) :: p(4)
    real(dp),  intent(in)   :: m
    logical, intent(in),optional :: outgoing
    complex(dp)             :: pol_mass(4)
    ! -------------------------------------
    real(dp) :: p0,px,py,pz, pv
    real(dp) :: ct,st,cphi,sphi
    integer :: pol

    !^^^IFmp
    !    p0=(p(1)+conjg(p(1)))/two
    !    px=(p(2)+conjg(p(2)))/two
    !    py=(p(3)+conjg(p(3)))/two
    !    pz=(p(4)+conjg(p(4)))/two
    !^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    !^^^END


    ! i=0 is longitudinal polarization

    ! -- distinguish between positive and negative energies
    if ( p0 > zero) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing
    if (present(outgoing)) then
       if (outgoing) pol = -pol
    endif

    pv= sqrt(abs(p0**2 - m**2))

    if (pv.gt.1d-8) then  ! for ``moving W'

       ct= pz/pv
       st= sqrt(abs(one-ct**2))

       if (st < tolb) then
          cphi=one
          sphi=zero
       else
          cphi= px/pv/st
          sphi= py/pv/st
       endif


       if(pol == -1.or.pol == 1) then
          pol_mass(1)=czero
          pol_mass(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
          pol_mass(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
          pol_mass(4)=-st/sqrt2
       elseif (pol == 0) then
          pol_mass(1)= pv/m
          pol_mass(2)= p0/m/pv*px
          pol_mass(3)= p0/m/pv*py
          pol_mass(4)= p0/m/pv*pz
       else
          stop 'pol_mass: pol out of range'
       endif

    else

       pol_mass = czero

       if (pol == -1.or.pol == 1) then

          pol_mass(2)=one/sqrt2
          pol_mass(3)=+ ci*pol/sqrt2

       elseif(pol == 0) then

          pol_mass(4) = (1.0_dp,zero)

       endif

    endif

  end function pol_mass

  function pol_mass2(p,m,i,out)
    integer, intent(in)       :: i
    complex(dp), intent(in) :: p(4)
    real(dp),  intent(in)   :: m
    character(len=*), intent(in):: out
    complex(dp)             :: pol_mass2(4)
    ! -------------------------------------

    if (out == 'out') then
       pol_mass2 = pol_mass(p,m,i,outgoing=.true.)
    else
       pol_mass2 = pol_mass(p,m,i,outgoing=.false.)
    endif
  end function pol_mass2



  !--- scalar products

  function sc(p1,p2)
    complex(dp), intent(in) :: p1(4), p2(4)
    complex(dp)             :: sc
    sc = p1(1)*p2(1)
    sc = sc - sum(p1(2:4)*p2(2:4))
  end function sc


  function scr(p1,p2)   !scalar product of real vectors
    real(dp), intent(in) :: p1(4), p2(4)
    real(dp) :: scr
    scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
  end function scr


  function sca(p1,p2)   !scalar product of real vectors, returned as
    real(dp), intent(in) :: p1(4), p2(4)  ! a complex number
    real(dp) :: foo
    complex(dp) :: sca
    foo = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
    sca = cmplx(foo,kind=dp)
  end function sca


  function scc(p1,eg)
    real(dp), intent(in) :: p1(:)
    complex(dp), intent(in) :: eg(:)
    complex(dp)             :: scc
    integer :: sizemin

    sizemin=min(size(p1),size(eg))

    scc = p1(1)*eg(1)
    scc = scc - sum(p1(2:sizemin)*eg(2:sizemin))

  end function scc


  function sc3(n1,n2)
    real(dp), intent(in) :: n1(3),n2(3)
    real(dp) :: sc3

    sc3 = sum(n1*n2)

  end function sc3


  function spb2(sp,v)
    complex(dp), intent(in) :: sp(:),v(:)
    complex(dp) :: spb2(size(sp))
    complex(dp) :: x0(4,4),xx(4,4),xy(4,4)
    complex(dp) :: xz(4,4),x5(4,4)
    complex(dp) :: y1,y2,y3,y4,bp,bm,cp,cm
    integer :: i,i1,i2,i3,Dv,Ds,imax

    Ds = size(sp)

    if (Ds == 4) then
       Dv = 4
    elseif (Ds == 8) then
       Dv = 6
    elseif (Ds == 16) then
       Dv = 8
    else
       stop 'spb2:Dv not allowed'
    endif

    imax = Ds/4

    do i=1,imax
       i1= 1+4*(i-1)
       i2=i1+3

       y1=sp(i1)
       y2=sp(i1+1)
       y3=sp(i1+2)
       y4=sp(i1+3)

       x0(1,i)=y3
       x0(2,i)=y4
       x0(3,i)=y1
       x0(4,i)=y2

       xx(1,i) = y4
       xx(2,i) = y3
       xx(3,i) = -y2
       xx(4,i) = -y1

       xy(1,i)=ci*y4
       xy(2,i)=-ci*y3
       xy(3,i)=-ci*y2
       xy(4,i)=ci*y1

       xz(1,i)=y3
       xz(2,i)=-y4
       xz(3,i)=-y1
       xz(4,i)=y2

       x5(1,i)=y1
       x5(2,i)=y2
       x5(3,i)=-y3
       x5(4,i)=-y4
    enddo

    if (Dv.eq.4) then

       do i=1,4
          spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)
       enddo

    elseif (Dv.eq.6) then
       bp = (v(5)+ci*v(6))
       bm=(v(5)-ci*v(6))

       do i=1,4

          spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
               -v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)

          i1 = i+4
          spb2(i1)= v(1)*x0(i,2)-v(2)*xx(i,2) &
               -v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)
       enddo

    elseif (Dv.eq.8) then
       bp=(v(5)+ci*v(6))
       bm=(v(5)-ci*v(6))
       cp=(v(7)+ci*v(8))
       cm=(v(7)-ci*v(8))

       do i=1,4

          spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
               -v(3)*xy(i,1)-v(4)*xz(i,1) &
               +bm*x5(i,2)-cm*x5(i,3)

          i1 = i+4

          spb2(i1) = v(1)*x0(i,2)-v(2)*xx(i,2) &
               -v(3)*xy(i,2)-v(4)*xz(i,2) &
               -bp*x5(i,1)+cm*x5(i,4)

          i2 = i1+4

          spb2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3) &
               -v(3)*xy(i,3)-v(4)*xz(i,3) &
               +bm*x5(i,4)+cp*x5(i,1)

          i3=i2+4

          spb2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4) &
               -v(3)*xy(i,4)-v(4)*xz(i,4) &
               -bp*x5(i,3)-cp*x5(i,2)

       enddo
    else
       stop 'spb2: Dv out of bound'
    endif

  end function spb2


  function spi2(v,sp)
    complex(dp), intent(in) :: sp(:),v(:)
    complex(dp) :: spi2(size(sp))
    complex(dp) :: x0(4,4),xx(4,4),xy(4,4)
    complex(dp) :: xz(4,4),x5(4,4)
    complex(dp) ::  y1,y2,y3,y4,bp,bm,cp,cm
    integer :: i,i1,i2,i3,imax,Dv,Ds

    Ds = size(sp)

    if (Ds == 4) then
       Dv = 4
    elseif (Ds == 8) then
       Dv = 6
    elseif (Ds == 16) then
       Dv = 8
    else
       stop 'spi2:Dv not allowed'
    endif


    imax = Ds/4

    do i=1,imax
       i1= 1+4*(i-1)
       i2=i1+3

       y1=sp(i1)
       y2=sp(i1+1)
       y3=sp(i1+2)
       y4=sp(i1+3)

       x0(1,i)=y3
       x0(2,i)=y4
       x0(3,i)=y1
       x0(4,i)=y2


       xx(1,i) = -y4
       xx(2,i) = -y3
       xx(3,i) = y2
       xx(4,i) = y1


       xy(1,i)=ci*y4
       xy(2,i)=-ci*y3
       xy(3,i)=-ci*y2
       xy(4,i)=ci*y1

       xz(1,i)=-y3
       xz(2,i)=y4
       xz(3,i)=y1
       xz(4,i)=-y2

       x5(1,i)=y1
       x5(2,i)=y2
       x5(3,i)=-y3
       x5(4,i)=-y4

    enddo

    if(Dv.eq.4) then

       do i=1,4

          spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
               -v(3)*xy(i,1)-v(4)*xz(i,1)
       enddo

    elseif (Dv.eq.6) then
       bp = (v(5)+ci*v(6))
       bm=(v(5)-ci*v(6))


       do i=1,4

          spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
               -v(3)*xy(i,1)-v(4)*xz(i,1) &
               -bp*x5(i,2)

          i1=i+4

          spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2) &
               -v(3)*xy(i,2)-v(4)*xz(i,2) &
               +bm*x5(i,1)

       enddo

    elseif (Dv.eq.8) then

       bp = (v(5)+ci*v(6))
       bm=(v(5)-ci*v(6))
       cp=(v(7)+ci*v(8))
       cm=(v(7)-ci*v(8))

       do i=1,4

          spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)&
               -v(3)*xy(i,1)-v(4)*xz(i,1)&
               -bp*x5(i,2)+ cp*x5(i,3)

          i1=i+4

          spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)&
               -v(3)*xy(i,2)-v(4)*xz(i,2)&
               +bm*x5(i,1)-cp*x5(i,4)

          i2=i1+4

          spi2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)&
               -v(3)*xy(i,3)-v(4)*xz(i,3)&
               -bp*x5(i,4)-cm*x5(i,1)

          i3=i2+4

          spi2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)&
               -v(3)*xy(i,4)-v(4)*xz(i,4)&
               +bm*x5(i,3)+cm*x5(i,2)


       enddo

    else
       stop 'spi2: Dv out of bounds'
    end if

  end function spi2


  function  psp1(sp1,sp2)
    complex(dp), intent(in) :: sp1(:)
    complex(dp), intent(in) :: sp2(:)
    complex(dp) :: psp1

    psp1 = sum(sp1(1:)*sp2(1:))

  end function psp1








  subroutine get_NaN(value,flag_nan)
    implicit none
    real(dp), intent(in)  :: value
    logical, intent(out) :: flag_nan

    flag_nan = .false.

    if (.not.value.le.0.0_dp .and. .not.value.gt.0.0_dp ) then
       flag_nan =.true.
    endif

  end subroutine get_NaN

  subroutine cyclicswap(amp,swaptype)
    ! -- perform a cyclic swap on h1,h2,h3.
    ! -- and on leptons h4,h5
    implicit none
    complex(dp) ::  amp(2,2,2,2,2),swapamp(2,2,2,2,2)
    integer :: swaptype

    swapamp=amp
    if (swaptype .eq. 1) then
       ! -- swap 1,2,3 --> 2,3,1
       swapamp(1,1,2,:,:)=amp(2,1,1,:,:)
       swapamp(1,2,1,:,:)=amp(1,1,2,:,:)
       swapamp(2,1,1,:,:)=amp(1,2,1,:,:)
       swapamp(1,2,2,:,:)=amp(2,1,2,:,:)
       swapamp(2,2,1,:,:)=amp(1,2,2,:,:)
       swapamp(2,1,2,:,:)=amp(2,2,1,:,:)
    elseif (swaptype .eq. 2) then
       ! -- swap 1,2,3 --> 3,1,2
       swapamp(1,1,2,:,:)=amp(1,2,1,:,:)
       swapamp(1,2,1,:,:)=amp(2,1,1,:,:)
       swapamp(2,1,1,:,:)=amp(1,1,2,:,:)
       swapamp(1,2,2,:,:)=amp(2,2,1,:,:)
       swapamp(2,2,1,:,:)=amp(2,1,2,:,:)
       swapamp(2,1,2,:,:)=amp(1,2,2,:,:)
       ! -- lepton swap 4,5 <--> 5,4
    elseif (swaptype .eq. 3)  then
       swapamp(:,:,:,1,2)=amp(:,:,:,2,1)
       swapamp(:,:,:,2,1)=amp(:,:,:,1,2)

       ! ---- old stuff, keep for now
       !      elseif (swaptype .eq. 3) then
       !c -- swap 1,2,3 <--> 2,1,3
       !         swapamp(1,1,2)=amp(1,1,2)
       !         swapamp(1,2,1)=amp(2,1,1)
       !         swapamp(2,1,1)=amp(1,2,1)
       !         swapamp(1,2,2)=amp(2,1,2)
       !         swapamp(2,2,1)=amp(2,2,1)
       !         swapamp(2,1,2)=amp(1,2,2)
       !      elseif (swaptype .eq. 4) then
       !c -- first swap 1<-->2 and then swap 1,2,3 <--> 3,1,2
       !c -- overall : 1,2,3 <--> 2,1,3 <--> 3,2,1
       !         swapamp(1,1,2)=amp(2,1,1)
       !         swapamp(1,2,1)=amp(1,2,1)
       !         swapamp(2,1,1)=amp(1,1,2)
       !         swapamp(1,2,2)=amp(2,2,1)
       !         swapamp(2,2,1)=amp(1,2,2)
       !         swapamp(2,1,2)=amp(2,1,2)
       !      elseif (swaptype .eq. 5) then
       !c -- first swap 1<-->2 and then swap 1,2,3 <--> 2,3,1
       !c -- overall : 1,2,3 <--> 2,1,3 <--> 1,3,2
       !         swapamp(1,1,2)=amp(1,2,1)
       !         swapamp(1,2,1)=amp(1,1,2)
       !         swapamp(2,1,1)=amp(2,1,1)
       !         swapamp(1,2,2)=amp(1,2,2)
       !         swapamp(2,2,1)=amp(2,1,2)
       !         swapamp(2,1,2)=amp(2,2,1)
    endif

    amp=swapamp


    return
  end subroutine cyclicswap



  subroutine cyclicswap2(amp,swaptype)
    ! -- perform a cyclic swap on h1,h2,h3.
    implicit none
    complex(dp) :: amp(2,2,2),swapamp(2,2,2)
    integer :: swaptype
    swapamp=amp
    if (swaptype .eq. 1) then
       ! -- swap 1,2,3 --> 2,3,1
       swapamp(1,1,2)=amp(2,1,1)
       swapamp(1,2,1)=amp(1,1,2)
       swapamp(2,1,1)=amp(1,2,1)
       swapamp(1,2,2)=amp(2,1,2)
       swapamp(2,2,1)=amp(1,2,2)
       swapamp(2,1,2)=amp(2,2,1)
    elseif (swaptype .eq. 2) then
       ! -- swap 1,2,3 --> 3,1,2
       swapamp(1,1,2)=amp(1,2,1)
       swapamp(1,2,1)=amp(2,1,1)
       swapamp(2,1,1)=amp(1,1,2)
       swapamp(1,2,2)=amp(2,2,1)
       swapamp(2,2,1)=amp(2,1,2)
       swapamp(2,1,2)=amp(1,2,2)
    endif

    amp=swapamp


    return
  end subroutine cyclicswap2


 subroutine converttoMCFMmom(pin,pout)
    implicit none
    real(dp), intent(in)  :: pin(4,6)
    real(dp), intent(out) :: pout(6,4)
    integer               :: npart

    do npart=1,6
       pout(npart,4)=pin(1,npart)           ! E
       pout(npart,1)=pin(2,npart)           ! px
       pout(npart,2)=pin(3,npart)           ! py
       pout(npart,3)=pin(4,npart)           ! pz
    enddo
    pout(1,:)=-pout(1,:)
    pout(2,:)=-pout(2,:)
  end subroutine converttoMCFMmom

end module mod_auxfunctions
