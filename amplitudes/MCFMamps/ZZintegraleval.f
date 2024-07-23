      subroutine ZZintegraleval(k,mass)
      use mod_types; use mod_consts_dp; use common_def
C---- Calculate integrals for the gg->ZZ process
c---- with a loop of quarks with mass mass
C-----Author: Keith Ellis, September 2013
      implicit none
      include 'ZZdlabels.f'
      include 'ZZclabels.f'
      include 'blabels.f'
      include 'ggZZintegrals.f'
      real(dp)  :: k(6,4)
      integer nu
      real(dp)  :: p1(4),p2(4),p3(4),p4(4),mass,mass2,
     & s12,s13,s23,s14,s24,s34,p3sq,p4sq
      complex(dp)  :: qlI4,qlI3,qlI2,qlI4_6explicit
C      complex(dp)  :: qlI4_6



c--- initialize QCDLoop on first pass
!      if (qlfirst) then
!        call qlinit
!        qlfirst=.false.
!      endif

      p1(:)=k(1,:)
      p2(:)=k(2,:)
      p3(:)=k(3,:)+k(4,:)
      p4(:)=k(5,:)+k(6,:)

      p3sq=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4sq=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      s12=(p1(4)+p2(4))**2
      s13=(p1(4)+p3(4))**2
      s23=(p2(4)+p3(4))**2
      s14=(p1(4)+p4(4))**2
      s24=(p2(4)+p4(4))**2
      s34=(p3(4)+p4(4))**2
      do nu=1,3
      s12=s12-(p1(nu)+p2(nu))**2
      s13=s13-(p1(nu)+p3(nu))**2
      s23=s23-(p2(nu)+p3(nu))**2
      s14=s14-(p1(nu)+p4(nu))**2
      s24=s24-(p2(nu)+p4(nu))**2
      s34=s34-(p3(nu)+p4(nu))**2
      enddo 
      mass2=mass**2

c--- 4-d integrals
      intD0(d1_2_34)=qlI4(zero,zero,p3sq,p4sq,s12,s23,mass2,mass2,mass2,
     &     mass2,musq,0)
      intD0(d2_1_34)=qlI4(zero,zero,p3sq,p4sq,s12,s13,mass2,mass2,mass2,
     &     mass2,musq,0)
      intD0(d1_34_2)=qlI4(zero,p3sq,zero,p4sq,s13,s23,mass2,mass2,mass2,
     &     mass2,musq,0)

      intC0(c1_34)=qlI3(zero,p3sq,s13,mass2,mass2,mass2,musq,0)
      intC0(c2_56)=qlI3(zero,p4sq,s24,mass2,mass2,mass2,musq,0)
      intC0(c2_34)=qlI3(zero,p3sq,s23,mass2,mass2,mass2,musq,0)
      intC0(c1_56)=qlI3(zero,p4sq,s14,mass2,mass2,mass2,musq,0)
      intC0(c1_2)=qlI3(zero,zero,s12,mass2,mass2,mass2,musq,0)
      intC0(c12_34)=qlI3(p3sq,p4sq,s12,mass2,mass2,mass2,musq,0)
      
      intB0(b12)=qlI2(s12,mass2,mass2,musq,0)
      intB0(b34)=qlI2(p3sq,mass2,mass2,musq,0)
      intB0(b56)=qlI2(p4sq,mass2,mass2,musq,0)
      intB0(b134)=qlI2(s13,mass2,mass2,musq,0)
      intB0(b234)=qlI2(s23,mass2,mass2,musq,0)
      intB0(rat)=1d0

c--- 6-d integrals
c      intD0(d6_1_2_34)=qlI4_6(p3sq,p4sq,s12,s23,mass2,musq)
c      intD0(d6_2_1_34)=qlI4_6(p3sq,p4sq,s12,s13,mass2,musq)
c       write(6,*) 'Old intD0(d6_1_2_34)',intD0(d6_1_2_34)
c       write(6,*) 'Old intD0(d6_2_1_34)',intD0(d6_2_1_34)

c--- 6-d integrals
      intD0(d6_1_2_34)=qlI4_6explicit(p3sq,p4sq,s12,s23,mass2,musq,
     & intD0(d1_2_34),intC0(c2_34),intC0(c12_34),intC0(c1_56),
     & intC0(c1_2))
      intD0(d6_2_1_34)=qlI4_6explicit(p3sq,p4sq,s12,s13,mass2,musq,
     & intD0(d2_1_34),intC0(c1_34),intC0(c12_34),intC0(c2_56),
     & intC0(c1_2))
c       write(6,*) 'New intD0(d6_1_2_34)',intD0(d6_1_2_34)
c       write(6,*) 'New intD0(d6_2_1_34)',intD0(d6_2_1_34)

c--- For checking by isolating integrals one by one
c      intB0(:)=czero

c      intC0(c1_34)=czero
c      intC0(c2_56)=czero
c      intC0(c1_2)=czero
c      intC0(c1_56)=czero
c      intC0(c2_34)=czero
c      intC0(c12_34)=czero

c      intD0(d1_2_34)=cone
c      intD0(d2_1_34)=cone
c      intD0(d1_34_2)=czero

c      if (docheck) then
c        write(6,*) 'Integraleval:resetting integrals for debugging'
c        write(6,*) 'intC0(c1_34)',intC0(c1_34)
c        write(6,*) 'intC0(c2_56)',intC0(c2_56)
c        write(6,*) 'intC0(c1_2)',intC0(c1_2)
c        write(6,*) 'intC0(c1_56)',intC0(c1_56)
c        write(6,*) 'intC0(c2_34)',intC0(c2_34)
c        write(6,*) 'intC0(c12_34)',intC0(c12_34)
c        write(6,*) 'intD0(d1_2_34)',intD0(d1_2_34)
c        write(6,*) 'intD0(d2_1_34)',intD0(d2_1_34)
c        write(6,*) 'intD0(d1_34_2)',intD0(d1_34_2)
c      endif
      
      return
      end





      function qlI4_6explicit(s34,s56,s12,s234,masssq,
     & musq,D0,C0_1,C0_2,C0_3,C0_4)
      use mod_types; use mod_consts_dp
c--- Evaluates the 6-dimensional scalar "hard" box integral
c--- formula taken from sixdim.prc
c--- Differs from qlI4_6 in that values of scalar integrals are passed in
c--- as arguments
      implicit none
      complex(dp)  :: qlI4_6explicit
      real(dp)  :: s34,s56,s12,s234,s134,musq,masssq,Y
      complex(dp)  :: D0,C0_1,C0_2,C0_3,C0_4,xD0,xC0(4)

c      D0=qlI4(0d0,0d0,s34,s56,s12,s234,masssq,masssq,masssq,masssq,musq,0)
c      C0_1=qlI3(0d0,s34,s234,masssq,masssq,masssq,musq,0)
c      C0_2=qlI3(s12,s34,s56,masssq,masssq,masssq,musq,0)
c      C0_3=qlI3(0d0,s234,s56,masssq,masssq,masssq,musq,0)
c      C0_4=qlI3(0d0,0d0,s12,masssq,masssq,masssq,musq,0)
      
      s134=s56+s34-s12-s234
      Y=((s134-s34)*(s234-s34)-s12*s34)
            
      xD0=-s12*s234-4d0*masssq/s234*Y
      xC0(1)=s234-s34
      xC0(2)=s12-s34-s56+2d0*s34*s56/s234
      xC0(3)=s234-s56
      xC0(4)=s12

      qlI4_6explicit=
     & (xD0*D0+xC0(1)*C0_1+xC0(2)*C0_2+xC0(3)*C0_3+xC0(4)*C0_4)
     & *s234/(2d0*Y)

      return
      end
  

      function qlI4_6(s34,s56,s12,s234,masssq,musq)
      use mod_types; use mod_consts_dp
c--- Evaluates the 6-dimensional scalar "hard" box integral
c--- formula taken from sixdim.prc
      implicit none
      complex(dp)  :: qlI4_6
      real(dp)  :: s34,s56,s12,s234,s134,musq,masssq,Y
      complex(dp)  :: D0,C0(4),xD0,xC0(4),qlI3,qlI4

      D0=qlI4(0d0,0d0,s34,s56,s12,s234,masssq,masssq,masssq,masssq,musq,
     & 0)
      C0(1)=qlI3(0d0,s34,s234,masssq,masssq,masssq,musq,0)
      C0(2)=qlI3(s12,s34,s56,masssq,masssq,masssq,musq,0)
      C0(3)=qlI3(0d0,s234,s56,masssq,masssq,masssq,musq,0)
      C0(4)=qlI3(0d0,0d0,s12,masssq,masssq,masssq,musq,0)
      
      s134=s56+s34-s12-s234
      Y=((s134-s34)*(s234-s34)-s12*s34)
C      Gram=s12/4d0*Y
            
      xD0=-s12*s234-4d0*masssq/s234*Y
      xC0(1)=s234-s34
      xC0(2)=s12-s34-s56+2d0*s34*s56/s234
      xC0(3)=s234-s56
      xC0(4)=s12

c--- DEBUG: isolate only box
c      C0(1)=dcmplx(0d0,0d0)
c      C0(3)=dcmplx(0d0,0d0)
c      C0(4)=dcmplx(0d0,0d0)
c      D0=dcmplx(1d0,0d0)
c--- DEBUG: isolate only box

      qlI4_6=s234/2d0/Y
     & *(xD0*D0+xC0(1)*C0(1)+xC0(2)*C0(2)+xC0(3)*C0(3)+xC0(4)*C0(4))

c      write(6,*) 'D0=',D0
c      write(6,*) 'xD0=',xD0
c      write(6,*) xD0*D0
c      write(6,*) xC0(1)*C0(1)
c      write(6,*) xC0(2)*C0(2)
c      write(6,*) xC0(3)*C0(3)
c      write(6,*) xC0(4)*C0(4)

      return
      end
  


