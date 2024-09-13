module mod_func_for_h2vv_real
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: get_ampprod_h2vv_real,get_ampprod_h2vv_real_heft

contains

  !-- square bracket in higgsamp_real.pdf
  !-- index: h1,h2,h3
  !-- 0 -> g(1) g(2) g(3) [H], without couplings
  subroutine get_ampprod_h2vv_real(j1,j2,j3,j4,j5,j6,j7,za,zb,sprod,ampprod)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
    complex(dp), intent(in) :: za(7,7),zb(7,7)
    real(dp), intent(in) :: sprod(7,7)
    complex(dp), intent(out) :: ampprod(-1:1,-1:1,-1:1)
    complex(dp) :: mya2a,mya2b,mya2c,mya4
    complex(dp) :: mya2a_t,mya2b_t,mya2c_t,mya4_t
    complex(dp) :: mya2a_b,mya2b_b,mya2c_b,mya4_b
    real(dp) :: hmass2

    call geta2a4(sprod(j1,j2),sprod(j1,j3),sprod(j2,j3),expmass**2,mya2a_t,mya2b_t,mya2c_t,mya4_t)
	! print *, expmass
    if (hwithb) then
       ! print *, "USING MASSIVE BOTTOM", mbsq 
       call geta2a4(sprod(j1,j2),sprod(j1,j3),sprod(j2,j3),mbsq,mya2a_b,mya2b_b,mya2c_b,mya4_b)
    else
       mya2a_b = czero
       mya2b_b = czero
       mya2c_b = czero
       mya4_b = czero
    endif
    ! -- here bottom and top contributions are summed (here we insert cggh, ct and cb)
    mya2a = cggh + ct*mya2a_t + cb*mya2a_b
    mya2b = cggh + ct*mya2b_t + cb*mya2b_b
    mya2c = cggh + ct*mya2c_t + cb*mya2c_b
    mya4 = cggh + ct*mya4_t + cb*mya4_b
    
    hmass2 = sprod(j1,j2)+sprod(j1,j3)+sprod(j2,j3)

    ampprod(+1,+1,+1) = hmass2**2/za(j1,j2)/za(j2,j3)/za(j1,j3)*mya4
    ampprod(+1,+1,-1) = zb(j2,j1)**3/zb(j3,j1)/zb(j3,j2)*mya2a
    ampprod(+1,-1,+1) = zb(j3,j1)**3/zb(j2,j1)/zb(j3,j2)*mya2c
    ampprod(+1,-1,-1) = za(j2,j3)**3/za(j1,j2)/za(j1,j3)*mya2b
    ampprod(-1,-1,-1) = hmass2**2/zb(j2,j1)/zb(j3,j2)/zb(j3,j1)*mya4
    ampprod(-1,-1,+1) = za(j1,j2)**3/za(j2,j3)/za(j1,j3)*mya2a
    ampprod(-1,+1,-1) = za(j1,j3)**3/za(j1,j2)/za(j2,j3)*mya2c
    ampprod(-1,+1,+1) = zb(j3,j2)**3/zb(j2,j1)/zb(j3,j1)*mya2b

    return

  end subroutine get_ampprod_h2vv_real

  !-- same in the mt->infty limit, first index is the expansion in 1/fourmtsq
  subroutine get_ampprod_h2vv_real_heft(j1,j2,j3,j4,j5,j6,j7,za,zb,sprod,ampprod)
    integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
    complex(dp), intent(in) :: za(7,7),zb(7,7)
    real(dp), intent(in) :: sprod(7,7)
    complex(dp), intent(out) :: ampprod(0:nheft_max,-1:1,-1:1,-1:1)
    real(dp) :: mya2a(0:nheft_max),mya2b(0:nheft_max),mya2c(0:nheft_max),mya4(0:nheft_max)
    real(dp) :: hmass2

    call geta2a4_heft(sprod(j1,j2),sprod(j1,j3),sprod(j2,j3),expmass**2,mya2a,mya2b,mya2c,mya4)

    hmass2 = sprod(j1,j2)+sprod(j1,j3)+sprod(j2,j3)
    
    ampprod(:,+1,+1,+1) = hmass2**2/za(j1,j2)/za(j2,j3)/za(j1,j3)*mya4(:)
    ampprod(:,+1,+1,-1) = zb(j2,j1)**3/zb(j3,j1)/zb(j3,j2)*mya2a(:)
    ampprod(:,+1,-1,+1) = zb(j3,j1)**3/zb(j2,j1)/zb(j3,j2)*mya2c(:)
    ampprod(:,+1,-1,-1) = za(j2,j3)**3/za(j1,j2)/za(j1,j3)*mya2b(:) 
    ampprod(:,-1,-1,-1) = hmass2**2/zb(j2,j1)/zb(j3,j2)/zb(j3,j1)*mya4(:)
    ampprod(:,-1,-1,+1) = za(j1,j2)**3/za(j2,j3)/za(j1,j3)*mya2a(:)
    ampprod(:,-1,+1,-1) = za(j1,j3)**3/za(j1,j2)/za(j2,j3)*mya2c(:)
    ampprod(:,-1,+1,+1) = zb(j3,j2)**3/zb(j2,j1)/zb(j3,j1)*mya2b(:)


    return

  end subroutine get_ampprod_h2vv_real_heft
    

  !-- result: a functions normalized to 1 in the mt->infty limit
  subroutine geta2a4(ss,tt,uu,mfsq,mya2a,mya2b,mya2c,mya4)
    real(dp), intent(in) :: ss,tt,uu,mfsq
    complex(dp), intent(out) :: mya2a,mya2b,mya2c,mya4
    complex(dp) :: a2a,a2b,a2c,a4
    real(dp) :: hmass2

    hmass2=ss+tt+uu
    a2a = ehsva2(ss,tt,uu,mfsq) !-- 123
    a2b = ehsva2(uu,ss,tt,mfsq) !-- 231
    a2c = ehsva2(tt,uu,ss,mfsq) !-- 312
   
    a4 = ehsva4(ss,tt,uu,mfsq)

    mya2a = -three*a2a*hmass2**2/ss**2
    mya2b = -three*a2b*hmass2**2/uu**2
    mya2c = -three*a2c*hmass2**2/tt**2

    mya4 = -three*a4

    return

  end subroutine geta2a4

  !-- same in the 1/mt expansion
  subroutine geta2a4_heft(ss,tt,uu,mfsq,mya2a,mya2b,mya2c,mya4)
    real(dp), intent(in) :: ss,tt,uu,mfsq
    real(dp), intent(out) :: mya2a(0:nheft_max),mya2b(0:nheft_max),mya2c(0:nheft_max),mya4(0:nheft_max)
    real(dp) :: hmass2,sh,th,uh,fourmfsq

    hmass2=ss+tt+uu
    sh = ss/hmass2
    th = tt/hmass2
    uh = uu/hmass2
    fourmfsq = 4.0_dp * mfsq

    mya2a = mya2_heft(sh,th,uh,hmass2,fourmfsq) !-- 123
    mya2b = mya2_heft(uh,sh,th,hmass2,fourmfsq) !-- 231
    mya2c = mya2_heft(th,uh,sh,hmass2,fourmfsq) !-- 312
    
    mya4 = mya4_heft(sh,th,uh,hmass2,fourmfsq)

    return

  end subroutine geta2a4_heft
    
  !-- expansion in 1/fourmfsq
  function mya2_heft(sh,th,uh,hmass2,fourmfsq)
    real(dp) :: sh,th,uh,hmass2,fourmfsq
    real(dp) :: mya2_heft(0:nheft_max)

    mya2_heft(0) = one
    mya2_heft(1) = (7*hmass2)/(30.0_dp*fourmfsq)
    mya2_heft(2) = (hmass2**2*(10*sh**2 + 10*th**2 + 23*th*uh + 10*uh**2 + 20*sh*(th + uh)))/(105.0_dp*fourmfsq**2)
    mya2_heft(3) = (2*hmass2**3*(39*sh**3 + 117*sh**2*(th + uh) + &
         9*sh*(13*th**2 + 30*th*uh + 13*uh**2) + &
         13*(3*th**3 + 10*th**2*uh + 10*th*uh**2 + 3*uh**3)))/(1575.0_dp*fourmfsq**3)
    mya2_heft(4) = (8*hmass2**4*(64*sh**4 + 64*th**4 + 275*th**3*uh + 410*th**2*uh**2 + &
         275*th*uh**3 + 64*uh**4 + 256*sh**3*(th + uh) + &
         sh**2*(384*th**2 + 873*th*uh + 384*uh**2) + &
         16*sh*(16*th**3 + 55*th**2*uh + 55*th*uh**2 + 16*uh**3)))/ &
         (17325.0_dp*fourmfsq**4)
    mya2_heft(5) = (16*hmass2**5*(380*sh**5 + 1900*sh**4*(th + uh) + &
         152*sh**3*(25*th**2 + 56*th*uh + 25*uh**2) + &
         20*sh**2*(190*th**3 + 651*th**2*uh + 651*th*uh**2 + 190*uh**3) + &
         sh*(1900*th**4 + 8540*th**3*uh + 13041*th**2*uh**2 + 8540*th*uh**3 + &
         1900*uh**4) + 4*(95*th**5 + 497*th**4*uh + 973*th**3*uh**2 + &
         973*th**2*uh**3 + 497*th*uh**4 + 95*uh**5)))/(315315.0_dp*fourmfsq**5)

    return

  end function mya2_heft

  !-- expansion in 1/fourmfsq
  function mya4_heft(sh,th,uh,hmass2,fourmfsq)
    real(dp) :: sh,th,uh,hmass2,fourmfsq
    real(dp) :: mya4_heft(0:nheft_max)

    mya4_heft(0) = one
    mya4_heft(1) = (hmass2*(7*sh**3 + 21*sh**2*(th + uh) + 7*(th + uh)**3 + &
         3*sh*(7*th**2 + 12*th*uh + 7*uh**2)))/(30.0_dp*fourmfsq)
    mya4_heft(2) = (hmass2**2*(10*sh**3 + 30*sh**2*(th + uh) + 10*(th + uh)**3 + &
         sh*(30*th**2 + 57*th*uh + 30*uh**2)))/(105.0_dp*fourmfsq**2)
    mya4_heft(3) = (2*hmass2**3*(39*sh**5 + 195*sh**4*(th + uh) + 390*sh**3*(th + uh)**2 + &
         39*(th + uh)**5 + 10*sh**2* &
         (39*th**3 + 121*th**2*uh + 121*th*uh**2 + 39*uh**3) + &
         5*sh*(39*th**4 + 156*th**3*uh + 242*th**2*uh**2 + 156*th*uh**3 + 39*uh**4)) &
         )/(1575.0_dp*fourmfsq**3)
    mya4_heft(4) = (8*hmass2**4*(64*sh**6 + 384*sh**5*(th + uh) + 64*(th + uh)**6 + &
         15*sh**4*(64*th**2 + 129*th*uh + 64*uh**2) + &
         5*sh**3*(256*th**3 + 803*th**2*uh + 803*th*uh**2 + 256*uh**3) + &
         5*sh**2*(192*th**4 + 803*th**3*uh + 1230*th**2*uh**2 + 803*th*uh**3 + &
         192*uh**4) + sh*(384*th**5 + 1935*th**4*uh + 4015*th**3*uh**2 + &
         4015*th**2*uh**3 + 1935*th*uh**4 + 384*uh**5)))/(17325.0_dp*fourmfsq**4) 
    mya4_heft(5) = (16*hmass2**5*(380*sh**7 + 2660*sh**6*(th + uh) + 380*(th + uh)**7 + &
         84*sh**5*(95*th**2 + 192*th*uh + 95*uh**2) + &
         28*sh**4*(475*th**3 + 1491*th**2*uh + 1491*th*uh**2 + 475*uh**3) + &
         7*sh**3*(1900*th**4 + 8020*th**3*uh + 12261*th**2*uh**2 + 8020*th*uh**3 + &
         1900*uh**4) + 21*sh**2* &
         (380*th**5 + 1988*th**4*uh + 4087*th**3*uh**2 + 4087*th**2*uh**3 + &
         1988*th*uh**4 + 380*uh**5) + &
         28*sh*(95*th**6 + 576*th**5*uh + 1491*th**4*uh**2 + 2005*th**3*uh**3 + &
         1491*th**2*uh**4 + 576*th*uh**5 + 95*uh**6)))/(315315.0_dp*fourmfsq**5)

    return

  end function mya4_heft
    
  !-- taken from MCFM, src/Httbar,ehsv.f
  !-- note that MCFM uses internal cli2
  !-- I use chaplin, the result should be real so there should be no phase problem

  function ehsva4(ss,tt,uu,mfsq)
    !-- ehsv:EqnA.8
    complex(dp) :: ehsva4
    real(dp) :: ss,tt,uu,mfsq

    ehsva4=ehsvb4(ss,tt,uu,mfsq)+ehsvb4(uu,ss,tt,mfsq)+ehsvb4(tt,uu,ss,mfsq)

    return 

  end function ehsva4

  function ehsva2(ss,tt,uu,mfsq)
    !-- ehsv:EqnA.9
    complex(dp) :: ehsva2
    real(dp) :: ss,tt,uu,mfsq

    ehsva2=ehsvb2(ss,tt,uu,mfsq)+ehsvb2(ss,uu,tt,mfsq)
    
    return 
  
  end function ehsva2

  function ehsvb4(ss,tt,uu,mfsq)
    !-- ehsv:EqnA.10
    complex(dp) :: ehsvb4
    real(dp) :: hmass2,ss,tt,uu,mfsq

    hmass2=ss+tt+uu
    !--- The Fermilab preprint has w2(s), but it makes no difference due
    !--- to symmetrization in ehsva4 above      
    ehsvb4=mfsq/hmass2*(-2.0_dp/3.0_dp +(mfsq/hmass2-0.25_dp)*(w2(tt,mfsq)-w2(hmass2,mfsq)+w3(ss,tt,uu,hmass2,mfsq)))

    return 
    
  end function ehsvb4

  function ehsvb2(ss,tt,uu,mfsq)
    !-- ehsv:EqnA.11
    complex(dp) :: ehsvb2
    real(dp) :: hmass2,ss,tt,uu,mfsq
    
    hmass2=ss+tt+uu
    
    ehsvb2=mfsq/hmass2**2*(ss*(uu-ss)/(ss+uu) &
         +2.0_dp*uu*tt*(uu+2.0_dp*ss)/(ss+uu)**2*(w1(tt,mfsq)-w1(hmass2,mfsq)) &
         +(mfsq-0.25_dp*ss) &
         *(0.5_dp*w2(ss,mfsq)+0.5_dp*w2(hmass2,mfsq)-w2(tt,mfsq)+w3(ss,tt,uu,hmass2,mfsq)) &
         +ss**2*(2.0_dp*mfsq/(ss+uu)**2-0.5_dp/(ss+uu))*(w2(tt,mfsq)-w2(hmass2,mfsq)) &
         +0.5_dp*uu*tt/ss*(w2(hmass2,mfsq)-2.0_dp*w2(tt,mfsq)) &
         +0.125_dp*(ss-12.0_dp*mfsq-4.0_dp*uu*tt/ss)*w3(tt,ss,uu,hmass2,mfsq))

    return 
   
  end function ehsvb2

  !-- this one would be for the qg amplitude
  function ehsva5(ss,tt,uu,mfsq)
    !-- ehsv:EqnA.14
    complex(dp) :: ehsva5
    real(dp) :: hmass2,ss,tt,uu,mfsq
    
    hmass2=ss+tt+uu
    
    ehsva5=mfsq/hmass2*(4.0_dp+4.0_dp*ss/(uu+tt)*(w1(ss,mfsq)-w1(hmass2,mfsq)) &
         +(1.0_dp-4.0_dp*mfsq/(uu+tt))*(w2(ss,mfsq)-w2(hmass2,mfsq)))
    
    return 
    
  end function ehsva5
  
  function w1(ss,mfsq)
    !-- ehsv:EqnA.19
    complex(dp) :: w1
    real(dp) :: ss,mfsq,rat,temp
    
    rat=4.0_dp*mfsq/ss
    temp=sqrt(abs(1.0_dp/rat))
    
    if (rat .lt. zero) then
       w1=2.0_dp*sqrt(1.0_dp-rat)*asinh(temp)
    elseif (rat .gt. 1.0_dp) then
       w1=2.0_dp*sqrt(rat-1.0_dp)*asin(temp)
    else 
       temp=2.0_dp*acosh(temp)
       w1=sqrt(1.0_dp-rat)*cmplx(temp,-pi,kind=dp)
    endif
    
    return
    
  end function w1
  
  function w2(ss,mfsq)
    !-- ehsv:EqnA.20
    complex(dp) :: w2
    real(dp) :: ss,mfsq,rat,tempr,tempi
    
    rat=ss/(4.0_dp*mfsq)
    tempr=sqrt(abs(rat))
    
    if (rat .lt. zero) then
       tempr=asinh(tempr)
       w2=4.0_dp*tempr**2
    elseif (rat .gt. 1.0_dp) then
       tempr=acosh(tempr)
       tempi=-4.0_dp*tempr*pi
       tempr=+4.0_dp*tempr**2-pi**2
       w2=cmplx(tempr,tempi,kind=dp)
    else 
       tempr=asin(tempr)
       w2=-4.0_dp*tempr**2
    endif
    
    return
    
  end function w2
  
  function w3(ss,tt,uu,varg,mfsq)
    !-- ehsv:EqnA.17
    complex(dp) :: w3
    real(dp) :: ss,tt,uu,varg,mfsq
    
    w3=i3(ss,tt,uu,varg,mfsq)-i3(ss,tt,uu,ss,mfsq)-i3(ss,tt,uu,uu,mfsq)
    
    return
    
  end function w3
  
  function i3(ss,tt,uu,varg,mfsq)
    !-- ehsv:EqnA.21
    complex(dp) :: i3,cli2,zth,zph
    real(dp) :: ss,tt,uu,varg,mfsq,rat,al,be,ga,r,theta,phi
    real(dp) :: arg1,arg2,arg3,arg4,arg
    real(dp15) :: HPL2real
    complex(dp15) :: HPL2
    
    rat=4.0_dp*mfsq/varg
    
    if (rat .lt. zero) then
       be=0.5_dp*(1.0_dp+sqrt(1.0_dp+4.0_dp*tt*mfsq/(uu*ss)))
       ga=0.5_dp*(1.0_dp+sqrt(1.0_dp-rat))
       arg1=ga/(ga+be-1.0_dp)
       arg2=(ga-1.0_dp)/(ga+be-1.0_dp)
       arg3=(be-ga)/be
       arg4=(be-ga)/(be-1.0_dp)
       i3=2.0_dp/(2.0_dp*be-1.0_dp) &
            !*(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4) &
            *(-HPL2real(0,1,real(arg1,kind=dp15),0d0)+HPL2real(0,1,real(arg2,kind=dp15),0d0)+HPL2real(0,1,real(arg3,kind=dp15),0d0)-HPL2real(0,1,real(arg4,kind=dp15),0d0) &
            +0.5_dp*(log(be)**2-log(be-1.0_dp)**2) &
            +log(ga)*log((ga+be-1.0_dp)/be) &
            +log(ga-1.0_dp)*log((be-1.0_dp)/(ga+be-1.0_dp)))
    elseif (rat .gt. 1.0_dp) then
       be=0.5_dp*(1.0_dp+sqrt(1.0_dp+4.0_dp*tt*mfsq/(uu*ss)))
       al=sqrt(rat-1.0_dp)           
       r=sqrt((al**2+1.0_dp)/(al**2+(2.0_dp*be-1.0_dp)**2))
       arg=r*(al**2+2.0_dp*be-1.0_dp)/(1.0_dp+al**2)
       if (arg .ge. 1.0_dp) then
          phi=zero
       else
          phi=acos(arg)
       endif
       arg=r*(al**2-2.0_dp*be+1.0_dp)/(1.0_dp+al**2)
       if (arg .ge. 1.0_dp) then
          theta=zero
       else
          theta=acos(arg)
       endif
       zth=r*cmplx(cos(theta),sin(theta),kind=dp)
       zph=r*cmplx(cos(phi),sin(phi),kind=dp)
       i3=2.0_dp/(2.0_dp*be-1.0_dp) &
            !.     *(2.0_dp*dble(cli2(zth))-2.0_dp*dble(cli2(zph))
            *(2.0_dp*real(HPL2(0,1,cmplx(zth,kind=dp15)),kind=dp)-2.0_dp*real(HPL2(0,1,cmplx(zph,kind=dp15)),kind=dp) &
            +(phi-theta)*(phi+theta-pi))
    else
       be=0.5_dp*(1.0_dp+sqrt(1.0_dp+4.0_dp*tt*mfsq/(uu*ss)))
       ga=0.5_dp*(1.0_dp+sqrt(1.0_dp-rat))
       arg1=ga/(ga+be-1.0_dp)
       arg2=(ga-1.0_dp)/(ga+be-1.0_dp)
       arg3=ga/(ga-be)
       arg4=(ga-1.0_dp)/(ga-be)
       
       i3=2.0_dp/(2.0_dp*be-1.0_dp) &
            !*(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4) &
            *(-HPL2real(0,1,real(arg1,kind=dp15),0d0)+HPL2real(0,1,real(arg2,kind=dp15),0d0)+HPL2real(0,1,real(arg3,kind=dp15),0d0)-HPL2real(0,1,real(arg4,kind=dp15),0d0) &
            +log(ga/(1.0_dp-ga))*log((ga+be-1.0_dp)/(be-ga)) &
            -ci*pi*log((ga+be-1.0_dp)/(be-ga)))
    endif
    
    return
    
  end function i3
  
  function acosh(y)
    real(dp) :: acosh,y
    
    acosh=log(y+sqrt(y**2-1.0_dp))
    
    return
    
  end function acosh
  
  function asinh(y)
    real(dp) :: asinh,y
    
    asinh=log(y+sqrt(y**2+1.0_dp))
    
    return
    
  end function asinh

end module mod_func_for_h2vv_real
