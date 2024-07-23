      subroutine Acalc(k1,k2,k3,k4,k5,k6,sprod,mass,A)
      use mod_types; use mod_consts_dp
      implicit none 
      include 'ZZclabels.f'
      include 'ZZdlabels.f'
      include 'ggZZintegrals.f'
      complex(dp)  :: A(6)
      real(dp)  :: mass,masssq,s12,s134,s234,s134bar,s234bar,p3sq,Y
      complex(dp)  :: za(6,6),zb(6,6)
      real(dp)  :: sprod(6,6)
      integer k1,k2,k3,k4,k5,k6
      masssq=mass**2
      s12=sprod(k1,k2)
      s134bar=sprod(k1,k3)+sprod(k1,k4)
      s234bar=sprod(k2,k3)+sprod(k2,k4)
      p3sq=sprod(k3,k4)
      s234=s234bar+p3sq
      s134=s134bar+p3sq
      Y=s134bar*s234bar-s12*p3sq

      A(1)=+0.5d0*masssq/s12*(
     & +2d0*intC0(c1_34)*s134bar
     & +2d0*intC0(c2_56)*(-s234bar-s12)
     & +2d0*intC0(c2_34)*s234bar
     & +2d0*intC0(c1_56)*(-s12-s134bar)
     & -2d0*Y*intD0(d1_34_2)
     & +s12*(s12-4*masssq)
     & *(intD0(d1_2_34)+intD0(d2_1_34)+intD0(d1_34_2)))

      A(2)=2d0*masssq*(intD0(d6_1_2_34)+intD0(d6_2_1_34)+intC0(c12_34)
     & +masssq*(intD0(d1_2_34)+intD0(d2_1_34)-intD0(d1_34_2)))


c      A(2)=masssq/Y*(
c     & +intC0(c1_34)*s134*s134bar
c     & +intC0(c2_56)*s134*(-s234bar-s12)
c     & +intC0(c2_34)*s234*s234bar
c     & +intC0(c1_56)*s234*(-s134bar-s12)
c     & +intC0(c1_2)*s12*(s234+s134)
c     & +intC0(c12_34)*(2d0*Y-((s134bar+s234bar)**2-4d0*s12*p3sq))
c     & -intD0(d2_1_34)*s12*s134**2-intD0(d1_2_34)*s12*s234**2
c     & -2d0*Y*masssq*(intD0(d1_2_34)+intD0(d1_34_2)+intD0(d2_1_34)))

      A(3)=0.5d0*masssq*s12*(intD0(d1_2_34)-intD0(d2_1_34)
     &     -intD0(d1_34_2))
      A(4)=0.5d0*masssq*s12*(intD0(d2_1_34)-intD0(d1_2_34)
     &     -intD0(d1_34_2))

      A(5)=0.5d0*masssq*s12/(s134*s234)*(
     & +2d0*s134*intD0(d6_1_2_34)+2d0*s234*intD0(d6_2_1_34)
     & -s134*s234*intD0(d1_34_2)
     & +4d0*masssq*(s134*intD0(d1_2_34)+s234*intD0(d2_1_34))
     & +2d0*(s234+s134)*intC0(c12_34))

c      A(5)=0.5d0*masssq*s12/Y*(
c     & +intC0(c1_34)*2d0*s134bar
c     & -intC0(c1_56)*2d0*(s12+s134bar)
c     & +intC0(c1_2)*2d0*s12
c     & -intD0(d1_2_34)*s12*s234-intD0(d2_1_34)*s12*s134-Y*intD0(d1_34_2))


      A(6)=czero
c      A(6)=masssq*s12/Y*(
c     & -intC0(c1_34)*s134bar
c     & +intC0(c1_56)*(s12+s134bar)
c     & +intC0(c2_34)*s234bar
c     & -intC0(c2_56)*(s12+s234bar))
c      write(6,*) 'A(6)',A(6)

      return
      end
