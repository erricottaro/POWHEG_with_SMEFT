      subroutine init_couplings
      use  mod_amplitudes_parameters, only: set_amplitudes_parameters
      implicit none
      include 'pwhg_par.h'
      include 'pwhg_flg.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_physpar.h'
      include 'cvecbos.h'
c Avoid multiple calls to this subroutine. The parameter file is opened
c but never closed ...
      logical called
      data called/.false./
      save called
      integer ewscheme
      double precision  Two, Four, Rt2, Pi
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter( Rt2   = 1.414213562d0 )
      parameter( Pi = 3.14159265358979323846d0 )
      real * 8 powheginput
      external powheginput
      real * 8 mass_low,mass_high
      integer nup,ndown,ngen
      ! - SMEFT parameters
      !double precision  ph_cggh, ph_cb, ph_ct, ph_cz
      real * 8 pwhg_alphas
      external pwhg_alphas
      integer j
      logical mt_expansion

      if(called) then
         return
      else
         called=.true.
         call ffini
      endif

      if(powheginput("#verytinypars").eq.1) then
         par_isrtinycsi = 1d-12
         par_isrtinyy = 1d-12
         par_fsrtinycsi = 1d-12
         par_fsrtinyy = 1d-12
      endif

      physpar_ml(1)=0.511d-3
      physpar_ml(2)=0.1057d0
      physpar_ml(3)=1.777d0
      physpar_mq(1)=0.33d0     ! up
      physpar_mq(2)=0.33d0     ! down
      physpar_mq(3)=0.50d0     ! strange
      physpar_mq(4)=1.50d0     ! charm
      physpar_mq(5)=4.80d0     ! bottom

      call powheginputstring("proc",flg_proc)
      if (flg_proc.eq."ZZ") then
         idvecbos1 = 23
         idvecbos2 = 23
         iddecay1=powheginput("vdecaymodeV1")
         iddecay2=powheginput("vdecaymodeV2")
         iddecay12=-iddecay1
         iddecay22=-iddecay2
      else if (flg_proc.eq."WW") then
         idvecbos1 = 24
         idvecbos2 = -24
c        The process should be W+W-
         iddecay1=min(powheginput("vdecaymodeV1"),powheginput("vdecaymodeV2"))
         iddecay2=max(powheginput("vdecaymodeV1"),powheginput("vdecaymodeV2"))
         if(iddecay1* iddecay2 > 0) then
                print *, 'the decay modes must have opposite sign for W+W-, check the powheg.input'
                stop
         endif
         iddecay12=-sign(1,iddecay1)*(abs(iddecay1)+1)
         iddecay22=-sign(1,iddecay2)*(abs(iddecay2)+1)
      else
         write (*,*) ,"process ",flg_proc," not allowed"
         call pwhg_exit(-1)
      endif

      call powheginputstring("contr",flg_contr)
      if ((flg_contr.ne."full").and.(flg_contr.ne."only_h").and.
     $     (flg_contr.ne."no_h").and.(flg_contr.ne."interf_h")) then
         write (*,*) ,"contribution ",flg_contr," not allowed"
         call pwhg_exit(-1)
      endif


      write(*,*) "======================================"
      write(*,*) " gg4l running for process: ", flg_proc
      write(*,*) "  and contribution: ", flg_contr
      write(*,*) "======================================"

 
 

c the only parameters relevant for this process are set
c via powheginput.

      ph_Zmass=powheginput('#zmass')
      if(ph_Zmass<0d0) ph_Zmass = 91.1876d0
      ph_tmass=powheginput('#tmass')
      if(ph_tmass<0d0) ph_tmass=172.5d0
      ph_Wmass=powheginput('#wmass')
      if(ph_Wmass<0d0) ph_Wmass = 80.398d0
      ph_bmass=powheginput('#bmass')
      if(ph_bmass<0d0) ph_bmass = 4.75d0
      ph_Hmass=powheginput('#hmass')
      if(ph_Hmass<0d0) ph_Hmass = 125d0

      ph_Zwidth=powheginput('#zwidth')
      if(ph_Zwidth<0d0) ph_Zwidth=2.4952d0
      ph_Wwidth=powheginput('#wwidth')
      if(ph_Wwidth<0d0) ph_Wwidth=2.04807d0
      ph_Hwidth=powheginput('#hwidth')
      if(ph_Hwidth<0d0) ph_Hwidth=0.41650d-2 

      ph_WmWw = ph_Wmass * ph_Wwidth
      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_Wmass2 = ph_Wmass**2
      ph_Zmass2 = ph_Zmass**2


      ph_alphaem = powheginput("#alphaem")
      if (ph_alphaem.lt.0d0) ph_alphaem = 1d0/132.3384d0
      ph_sthw2 = powheginput("#sthw2")
      if (ph_sthw2.lt.0d0) ph_sthw2 = 0.2226459
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_unit_e = sqrt(4*pi*ph_alphaem)
      ph_gfermi = pi*ph_alphaem/(sqrt(2d0)*ph_sthw2*ph_Wmass2)

c set up masses and widths for resonance damping factors
      physpar_pdgmasses(22) = 0
      physpar_pdgmasses(23) = ph_Zmass
      physpar_pdgmasses(24) = ph_Wmass
      physpar_pdgmasses(25) = ph_Hmass
      physpar_pdgmasses(5)  = ph_bmass
      physpar_pdgmasses(6)  = ph_tmass
      physpar_pdgwidths(23) = ph_Zwidth
      physpar_pdgwidths(24) = ph_Wwidth
      physpar_pdgwidths(25) = ph_Hwidth
      physpar_pdgwidths(6) = ph_twidth
      physpar_pdgwidths(22) = 0

      do j=1,physpar_npdg
         physpar_pdgmasses(-j) = physpar_pdgmasses(j)
         physpar_pdgwidths(-j) = physpar_pdgwidths(j)
      enddo

      physpar_phspmasses = physpar_pdgmasses
      physpar_phspwidths = physpar_pdgwidths

c     Initialize couplings and parameters in amplitudes library
      nup=powheginput("#nup")
      if (nup.lt.0)  nup=2
      ndown=powheginput("#ndown")
      if (ndown.lt.0) ndown=3
      ngen=powheginput("#ngen")
      if (ngen.lt.0) ngen=2
      massiveloops=powheginput('#massiveloops') /= 0
      write(*,*), " check that b-mass is non-zero and massive loops are T",ph_bmass, massiveloops
c     Initialize SMEFT couplings
	  ph_cggh=powheginput("#cggh")
	  if(ph_cggh<0d0) ph_cggh=1.0d0
	  ph_cb=powheginput("#cb")
	  if(ph_cb<0d0) ph_cb=1.0d0 
	  ph_ct=powheginput("#ct")
	  if(ph_ct<0d0) ph_ct=1.0d0
	  ph_cz=powheginput("#cz")
	  if(ph_cz<0d0) ph_cz=1.0d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Set here the number of light flavours
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      st_nlight=5
      if (nup+ndown.ne.st_nlight) then
         write(*,*) "WARNING: mismatch between st_nlight and nup+ndown",
     $        st_nlight, nup+ndown
      endif

c     Consistency checks:
      if((ph_bmass>0d0).and.(.not.massiveloops)) then
         write(*,*) "ERROR: bottom mass must be set to zero"
         write(*,*) "when only massless loops are considered"
         call exit(-1)
      endif
	! - adding 
      call set_amplitudes_parameters(ph_Zmass,ph_Zwidth,ph_Wmass
     $     ,ph_Wwidth,ph_sthw2, ph_unit_e, ph_Hmass, ph_Hwidth, ph_tmass, ph_bmass,
     $     st_nlight, nup, ndown, ngen, iddecay1, iddecay2, flg_proc
     $     ,flg_contr, ph_cggh, ph_cb, ph_ct, ph_cz)

*********************************************************
* Print out of all the relevant couplings
*   so that we have them in a log file
*********************************************************
      write(*,*) '**************************************************'
      write(*,*) '* init_couplings.f writeout'
      write(*,*) '**************************************************'
      write(*,*) 'alpha=', ph_alphaem
      write(*,*) 'gfermi=', ph_gfermi
      write(*,*) 'bmass=', ph_bmass
      write(*,*) 'zmass=', ph_Zmass
      write(*,*) 'zwidth=', ph_Zwidth
      write(*,*) 'wmass=', ph_Wmass
      write(*,*) 'wwidth=', ph_Wwidth
      write(*,*) 'tmass=', ph_tmass
      write(*,*) 'twidth=', ph_twidth
      write(*,*) 'hmass=',  ph_hmass
      write(*,*) 'hwidth=', ph_hwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) 'n up = ',nup
      write(*,*) 'n down  = ',ndown
      write(*,*) 'alphas(MZ) = ', pwhg_alphas(ph_Zmass**2
     $        ,st_lambda5MSB,st_nlight)
      write(*,*) 'cggh = ', ph_cggh
      write(*,*) 'cb = ', ph_cb
      write(*,*) 'ct = ', ph_ct
      write(*,*) 'cz = ', ph_cz
      write(*,*) '**************************************************'
      end

