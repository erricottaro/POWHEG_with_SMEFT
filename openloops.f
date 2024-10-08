!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!                 ===POWHEG+OpenLoops Interface===
!
! For independent implementation you can download OpenLoops seperately from:
!                http://openloops.hepforge.org
!
! and follow the instruction for installation. Afterwards, add the
! following lines at the approbirate place in the Makefile of the
! POWHEG-BOX:
!
! ///////////////////////////////////////////////////////////////
! OLPATH=<PATH_TO_OPENLOOPS>
! ifdef OLPATH
! FF+= -I$(OLPATH)/lib_src/openloops/mod
! LIBS+= -Wl,-rpath=$(OLPATH)/lib -L$(OLPATH)/lib  -lopenloops
! USER+=openloops.o
! endif
! ///////////////////////////////////////////////////////////////
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      module openloops_powheg
        use openloops, only: set_parameter
        implicit none

        type reg_processes
          integer :: id
          integer :: nleg
          integer, allocatable :: proc(:)
          integer :: amptype
          integer :: ncolb
          integer, allocatable :: flowbasis(:,:,:)
          character(len=50) :: approx
        end type
        type(reg_processes), allocatable, save :: processes(:)

        contains

! Register & Loopup processes
        function register(flav,amptype,intid,approx)
          use openloops, only: register_process, start,
     &                         tree_colbasis_dim, tree_colourflow
          use ol_parameters_decl_dp, only: approximation
          implicit none
          integer, intent(in) :: flav(:)
          integer, intent(in) :: amptype
          integer, intent(out), optional :: intid
          integer :: register, i, nreg, flavreg(size(flav))
          integer :: ncolb, colelemsz, nhel_out
          character(len=*), intent(in),  optional :: approx
          character(len=50) :: approxbak

          type(reg_processes), allocatable :: processes_bak(:)


          if (allocated(processes)) then
            do i = 1,size(processes)
              if ( processes(i)%nleg == size(flav) .and. processes(i)%amptype == amptype) then
                if (all(processes(i)%proc == flav)) then
                  if (present(approx)) then
                    if (trim(processes(i)%approx) == trim(approx)) then
                      register = processes(i)%id
                      if (present(intid)) intid = i
                      return
                    else
                      goto 11
                    end if
                  end if
                  register = processes(i)%id
                  if (present(intid)) intid = i
                  return
                end if
              end if
            end do
          end if

11       continue

          ! register unknown process
          do i=1, size(flav)
            if (flav(i) == 0 ) then
                flavreg(i) = 21
            else
                flavreg(i) = flav(i)
            end if
          end do
          if (present(approx)) then
            approxbak = approximation
            call set_parameter("approx", trim(approx))
          end if
          register =  register_process(flavreg,amptype)
          if (present(approx)) then
            call set_parameter("approx", approxbak)
          end if

          if ( register < 0 ) then
            print*, "[POWHEG-BOX+OpenLoops] Process not found!"
            stop    ! process not found -> fatal stop
          else
            if (allocated(processes)) then
              nreg = size(processes)
              allocate(processes_bak(nreg))
              processes_bak = processes
              deallocate(processes)
              allocate(processes(nreg+1))
              processes(1:nreg) = processes_bak
              deallocate(processes_bak)
            else
              nreg=0
              allocate(processes(1))
            end if
            processes(nreg+1)%id = register
            processes(nreg+1)%nleg = size(flav)
            allocate(processes(nreg+1)%proc(processes(nreg+1)%nleg))
            processes(nreg+1)%proc = flav
            processes(nreg+1)%amptype = amptype
            if (amptype .eq. 12) then
!              call loop_colbasis_dim(processes(nreg+1)%id,
!     &          processes(nreg+1)%ncolb,colelemsz,nhel_out) ! obtain
!     dimension of colour flow basis
!              allocate(processes(nreg+1)%flowbasis(2,
!     &               processes(nreg+1)%nleg,processes(nreg+1)%ncolb))
!              call loop_colourflow(processes(nreg+1)%id,
!     &           processes(nreg+1)%flowbasis)! obtain colour flow basis
            else
              call tree_colbasis_dim(processes(nreg+1)%id,
     &                  processes(nreg+1)%ncolb,colelemsz,nhel_out) ! obtain dimension of colour flow basis
              allocate(processes(nreg+1)%flowbasis(2,
     &                  processes(nreg+1)%nleg,processes(nreg+1)%ncolb))
              call tree_colourflow(processes(nreg+1)%id,
     &                  processes(nreg+1)%flowbasis)! obtain colour flow basis
            end if
            if (present(approx)) then
              processes(nreg+1)%approx = trim(approx)
            else
              processes(nreg+1)%approx = ""
            end if
          end if
        end function register

! Born
        subroutine openloops_born(p,bflav,born,bornjk,bmunu)
          use openloops, only: evaluate_tree, evaluate_ccmatrix, evaluate_scpowheg
          implicit none
          include 'pwhg_st.h'
          double precision,intent(in) :: p(:,:)
          integer,intent(in) :: bflav(:)
          double precision,intent(out) :: born
          double precision,intent(out),optional :: bornjk(:,:)
          double precision,intent(out),optional :: bmunu(:,:,:)
          double precision :: bornew
          double precision :: bornmunu(0:3,0:3)
          integer j, id
          id = register(bflav,1)
          call set_parameter("alpha_qcd", st_alpha)   !alpha_qcd

          if (.not. present(bornjk) ) then ! only born
            call evaluate_tree(id, p, born)
          else ! born && colour-correlated born
            bornjk = 0
!            call evaluate_ccmatrix(id, p, born, bornjk, bornew)
            call evaluate_ccmatrix(id, p, born, bornjk(1:size(bflav),1:size(bflav)), bornew)
            bornjk = -1*bornjk
          end if

          if (present(bmunu)) then  ! spin-correlated born
            bmunu = 0
            do j=1,size(bflav)
              if (bflav(j) == 0) then
                bornmunu = 0
                call evaluate_scpowheg(id, p, j, born, bornmunu)
                bmunu(:,:,j) = bornmunu
              else
                bmunu(:,:,j) = 0
              end if
            end do
          end if
        end subroutine openloops_born

! Loop^2
        subroutine openloops_loop2(p,lflav,res,loopjk,lmunu,approx)
          use openloops, only: evaluate_loop2, evaluate_ccmatrix2, evaluate_scpowheg2
          implicit none
          include 'pwhg_st.h'
          double precision,intent(in) :: p(:,:)
          integer,intent(in) :: lflav(:)
          double precision,intent(out) :: res
          double precision,intent(out),optional :: loopjk(:,:)
          double precision,intent(out),optional :: lmunu(:,:,:)
          character(len=*), intent(in),  optional :: approx
          double precision :: acc
          double precision :: bornew
          double precision :: loopmunu(0:3,0:3)
          integer j, id

          id = register(lflav,12,approx=trim(approx))
          call set_parameter("alpha_qcd", st_alpha)   !alpha_qcd


          if (.not. present(loopjk) ) then ! only born
            call evaluate_loop2(id, p, res, acc)
          else ! loop^2 && colour-correlated loop^2
            loopjk = 0
            call evaluate_ccmatrix2(id, p, res, loopjk(1:size(lflav),1:size(lflav)), bornew)
            loopjk = -1*loopjk
          end if

          if (present(lmunu)) then  ! spin-correlated loop^2
            lmunu = 0
            if (res /= 0) then
              do j=1,size(lflav)
                if (lflav(j) == 0) then
                  loopmunu = 0
                  call evaluate_scpowheg2(id, p, j, res, loopmunu)
                  lmunu(:,:,j) = loopmunu
                else
                  lmunu(:,:,j) = 0
                end if
              end do
            end if
          end if
        end subroutine openloops_loop2

! Real
        subroutine openloops_real(p,rflav,res)
          use openloops, only: evaluate_tree
          implicit none
          include 'pwhg_st.h'
          include 'pwhg_math.h'
          double precision, intent(in) :: p(:,:)
          integer, intent(in) :: rflav(:)
          double precision, intent(out) :: res
          integer id
          id = register(rflav,1)
          call set_parameter("alpha_qcd", st_alpha)   !alpha_qcd
          call evaluate_tree(id, p, res)
          res=res/(st_alpha/(2d0*pi))
        end subroutine openloops_real

! Loop^2 real this is the one I will probably need to modify (03/09/24)
        subroutine openloops_loop2real(p,lflav,res,acc,approx)
          use openloops, only: evaluate_loop2
          implicit none
          include 'pwhg_st.h'
          include 'pwhg_math.h'
          double precision,intent(in) :: p(:,:)
          integer,intent(in) :: lflav(:)
          double precision,intent(out) :: res
          double precision,intent(out),optional :: acc
          character(len=*), intent(in),  optional :: approx
          integer :: id
          double precision :: olacc
          id = register(lflav,12,approx=trim(approx))
          call set_parameter("alpha_qcd", st_alpha)   !alpha_qcd
          ! print *, st_alpha
          call evaluate_loop2(id, p, res, olacc)
          res=res/(st_alpha/(2d0*pi))
          if(present(acc)) then
            acc=olacc
          end if
        end subroutine openloops_loop2real

! Virtual
        subroutine openloops_virtual(p,vflav,virtual, born)
          use openloops, only: evaluate_loop
          implicit none
          include 'pwhg_st.h'
          include 'pwhg_math.h'
          double precision, intent(in) ::  p(:,:)
          integer, intent(in) :: vflav(:)
          double precision, intent(out) ::  virtual
          double precision, intent(out), optional ::  born
          double precision :: m2L0, m2L1(0:2), acc
          integer :: id
          id = register(vflav,11)
          call set_parameter("alpha_qcd", st_alpha)   ! alpha_qcd
          call set_parameter("mu", sqrt(st_muren2))   ! renormalization scale
          call evaluate_loop(id, p, m2L0, m2L1, acc)
          virtual = m2L1(0)/(st_alpha/(2d0*pi))
          if (present(born)) born = m2L0
        end subroutine openloops_virtual

! Tree
        subroutine openloops_tree(p,flav,res)
          use openloops, only: evaluate_tree
          implicit none
          include 'pwhg_st.h'
          double precision,intent(in) :: p(:,:)
          integer,intent(in) :: flav(:)
          double precision,intent(out) :: res
          integer :: id
          id = register(flav,1)
          call set_parameter("alpha_qcd", st_alpha)   !alpha_qcd
          call evaluate_tree(id, p, res)
        end subroutine openloops_tree

! Color flow
      subroutine openloops_colour(p,flav,color,ifl)
        use openloops, only: evaluate_tree_colvect2, tree_colbasis_dim, tree_colourflow
        implicit none
        include 'pwhg_st.h'
        include 'pwhg_math.h'
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: flav(:)
        integer, intent(out) ::  color(:,:)
        integer, intent(out) :: ifl
        integer :: id, intid, randomflow, i
        integer :: ncolb, colelemsz, nhel_out
        double precision, allocatable :: m2arr(:)
        double precision random, toss, cumm2arr
        external random
        id = register(flav,1,intid)
        allocate(m2arr(processes(intid)%ncolb))
        call set_parameter("alpha_qcd", st_alpha)   ! alpha_qcd
        call evaluate_tree_colvect2(id, p, m2arr) ! obtain squared amplitudes for individual colour flows
        m2arr  = m2arr/sum(m2arr)   ! normalize
        toss = random() ! get random number
        cumm2arr=0
        do i=1,size(m2arr)
          cumm2arr = cumm2arr + m2arr(i)
          if (toss < cumm2arr) then
            randomflow = i
            exit
          end if
        end do
        color(:,1:2) = processes(intid)%flowbasis(:,1:2,randomflow) ! return colour flow basis element with largest probability
        color(2,3:) = processes(intid)%flowbasis(1,3:,randomflow) ! switch final state ordering
        color(1,3:) = processes(intid)%flowbasis(2,3:,randomflow)
        call map_flow_basis(color,size(flav))
        ifl =  randomflow
        deallocate(m2arr)
      end subroutine openloops_colour

      subroutine map_flow_basis(color,n)
      implicit none
      integer, intent(inout) ::  color(2,n)
      integer, intent(in) :: n
      integer :: i, j
      color = color
      do i=1,2
        do j=1,n
          if (color(i,j)/=0) color(i,j)=color(i,j)+500
        end do
      end do
      end subroutine map_flow_basis

      subroutine openloops_borncolour(p,bflav,color)
        implicit none
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: bflav(:)
        integer, intent(out) ::  color(:,:)
        integer ifl
        call openloops_colour(p,bflav,color,ifl)
      end subroutine openloops_borncolour

      subroutine openloops_realcolour(p,rflav,color)
        implicit none
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: rflav(:)
        integer, intent(out) ::  color(:,:)
        integer ifl
        call openloops_colour(p,rflav,color,ifl)
      end subroutine openloops_realcolour

      subroutine openloops_borncolour_ifl(p,bflav,color,ifl)
        implicit none
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: bflav(:)
        integer, intent(out) ::  color(:,:), ifl
        call openloops_colour(p,bflav,color,ifl)
      end subroutine openloops_borncolour_ifl

      subroutine openloops_realcolour_ifl(p,rflav,color,ifl)
        implicit none
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: rflav(:)
        integer, intent(out) ::  color(:,:), ifl
        call openloops_colour(p,rflav,color,ifl)
      end subroutine openloops_realcolour_ifl


! Initialize OpenLoops parameters
      subroutine openloops_init()
        use openloops, only: start
        use openloops_blha, only: olp_printparameter
        implicit none
        real * 8 powheginput
        external powheginput
        include 'PhysPars.h'
        include 'pwhg_st.h'

! Physical parameters
		! write(*,*) "Initializing openloops parameters"
		! write(*,*) ph_cz, ph_ct, ph_cb
		! write(*,*) ph_cz*ph_Zmass, ph_ct*ph_tmass, ph_cb*ph_bmass 
        call set_parameter("alpha_qcd", st_alpha)
        call set_parameter("ew_scheme", 2)
        call set_parameter("alpha_qed", ph_alphaem)
        call set_parameter("z_mass",  ph_Zmass)
        call set_parameter("z_width", ph_Zwidth)
        call set_parameter("w_mass",  ph_Wmass)
        call set_parameter("w_width",  ph_Wwidth)
        call set_parameter("h_mass",  ph_Hmass)
        call set_parameter("h_width", ph_Hwidth)
        call set_parameter("t_mass",  ph_tmass)
        call set_parameter("t_width", ph_twidth)
        call set_parameter("b_mass",  ph_bmass)
!        if (ph_bwidth   /= 0) call set_parameter("b_width",  ph_bwidth)
!        if (ph_cmass    /= 0) call set_parameter("c_mass",  ph_cmass)
!        if (ph_cwidth   /= 0) call set_parameter("c_width",  ph_cwidth)
!        if (ph_emass    /= 0) call set_parameter("e_mass", ph_emass)
!        if (ph_ewidth   /= 0) call set_parameter("e_width", ph_ewidth)
!        if (ph_mumass   /= 0) call set_parameter("mu_mass", ph_mumass)
!        if (ph_muwidth  /= 0) call set_parameter("mu_width", ph_muwidth)
!        if (ph_taumass  /= 0) call set_parameter("tau_mass", ph_taumass)
!        if (ph_tauwidth /= 0) call set_parameter("tau_width", ph_tauwidth)

! OpenLoops parameters
        ! verbosity level
      if (powheginput('#ol_verbose')>0) then
        call set_parameter("verbose", powheginput('#ol_verbose'))
      end if

! switch off complex mass scheme -> weak mixing angle is real. Should be switched on for processes with intermediate resonances (default)
!      call set_parameter("use_cms",0)

! Initialize parameters
        call start

! Print OL parameters to file
!        call openloops_printparameter()

      end subroutine openloops_init

! Initialize OpenLoops parameters
      subroutine openloops_printparameter()
        use openloops_blha, only: olp_printparameter
        implicit none
        call olp_printparameter("parameters.ol")
      end subroutine openloops_printparameter

      end module openloops_powheg

! Fortran77 Interfaces
      subroutine openloops_born(p,bflav,born,bornjk,bmunu)
        use openloops_powheg, only: ol_born => openloops_born
        implicit none
        include 'nlegborn.h'
        double precision,intent(in) :: p(0:3,nlegbornexternal)
        integer,intent(in) :: bflav(nlegbornexternal)
        double precision,intent(out) :: born
        double precision,intent(out),optional :: bornjk(nlegbornexternal,nlegbornexternal)
        double precision,intent(out),optional :: bmunu(0:3,0:3,nlegbornexternal)
        call ol_born(p,bflav,born,bornjk,bmunu)
      end subroutine openloops_born

      subroutine openloops_real(p,rflav,amp2real)
        use openloops_powheg, only: ol_real => openloops_real
        implicit none
        include 'nlegborn.h'
        double precision, intent(in) ::p(0:3,nlegbornexternal+1)
        integer, intent(in) :: rflav(nlegbornexternal+1)
        double precision, intent(out) :: amp2real
        call ol_real(p,rflav,amp2real)
      end subroutine openloops_real

      subroutine openloops_virtual(p,vflav,virtual)
        use openloops_powheg, only: ol_virtual => openloops_virtual
        implicit none
        include 'nlegborn.h'
        double precision, intent(in) ::  p(0:3,nlegbornexternal)
        integer, intent(in) :: vflav(nlegbornexternal)
        double precision, intent(out) ::  virtual
        call ol_virtual(p,vflav,virtual)
      end subroutine openloops_virtual

      subroutine openloops_loop2(p,lflav,res,loopjk,lmunu)
        use openloops_powheg, only: ol_loop2 => openloops_loop2
        implicit none
        include 'nlegborn.h'
        double precision,intent(in) :: p(0:3,nlegbornexternal)
        integer,intent(in) :: lflav(nlegbornexternal)
        double precision,intent(out) :: res
        double precision,intent(out),optional :: loopjk(nlegbornexternal,nlegbornexternal)
        double precision,intent(out),optional :: lmunu(0:3,0:3,nlegbornexternal)
        call ol_loop2(p,lflav,res,loopjk=loopjk,lmunu=lmunu)
      end subroutine openloops_loop2

      subroutine openloops_loop2real(p,lflav,res,acc)
        use openloops_powheg, only: ol_loop2real => openloops_loop2real
        implicit none
        include 'nlegborn.h'
        double precision,intent(in) :: p(0:3,nlegbornexternal+1)
        integer,intent(in) :: lflav(nlegbornexternal+1)
        double precision,intent(out) :: res
        double precision,intent(out),optional :: acc
        call ol_loop2real(p,lflav,res,acc=acc)
      end subroutine openloops_loop2real

      subroutine openloops_borncolour_ifl(p,bflav,color,ifl)
        use openloops_powheg, only: ol_colour => openloops_colour
        implicit none
        include 'nlegborn.h'
        double precision, intent(in) ::  p(0:3,nlegbornexternal)
        integer, intent(in) :: bflav(nlegbornexternal)
        integer, intent(out) ::  color(2,nlegbornexternal), ifl
        call ol_colour(p,bflav,color,ifl)
      end subroutine openloops_borncolour_ifl

      subroutine openloops_realcolour_ifl(p,rflav,color,ifl)
        use openloops_powheg, only: ol_colour => openloops_colour
        implicit none
        include 'nlegborn.h'
        double precision, intent(in) ::  p(0:3,nlegbornexternal+1)
        integer, intent(in) :: rflav(nlegbornexternal+1)
        integer, intent(out) ::  color(2,nlegbornexternal+1), ifl
        call ol_colour(p,rflav,color,ifl)
      end subroutine openloops_realcolour_ifl

      subroutine openloops_borncolour(p,bflav,color)
        use openloops_powheg, only: ol_colour => openloops_colour
        implicit none
        include 'nlegborn.h'
        double precision, intent(in) ::  p(0:3,nlegbornexternal)
        integer, intent(in) :: bflav(nlegbornexternal)
        integer, intent(out) ::  color(2,nlegbornexternal)
        integer ifl
        call ol_colour(p,bflav,color,ifl)
      end subroutine openloops_borncolour

      subroutine openloops_realcolour(p,rflav,color)
        use openloops_powheg, only: ol_colour => openloops_colour
        implicit none
        include 'nlegborn.h'
        double precision, intent(in) ::  p(0:3,nlegbornexternal+1)
        integer, intent(in) :: rflav(nlegbornexternal+1)
        integer, intent(out) ::  color(2,nlegbornexternal+1)
        integer ifl
        call ol_colour(p,rflav,color,ifl)
      end subroutine openloops_realcolour

      subroutine openloops_init()
        use openloops_powheg, only: ol_init => openloops_init
        implicit none
        call ol_init()
      end subroutine openloops_init

      subroutine openloops_printparameter()
        use openloops_powheg, only: ol_printparameter => openloops_printparameter
        implicit none
        call ol_printparameter()
      end subroutine openloops_printparameter


