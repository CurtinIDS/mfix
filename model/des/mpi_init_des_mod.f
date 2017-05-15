!----------------------------------------------------------------------!
!  Module: MPI_INIT_DES                                                !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for setting up DES MPI communications.   !
!                                                                      !
!----------------------------------------------------------------------!
      module mpi_init_des

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use parallel_mpi
      use mpi_utility
      use discretelement
      use desgrid
      use compar
      use physprop
      use sendrecv
      use des_bc
      use desmpi_wrapper
      use sendrecvnode
      use mfix_pic
      use des_thermo
      use run, only: ENERGY_EQ,ANY_SPECIES_EQ
      use param, only: DIMENSION_N_s
      use des_rxns

      use desmpi


      contains


!----------------------------------------------------------------------!
!  Module: DESMPI_INIT                                                 !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Allocates and initializes variables used by MPI send/recv  !
!  calls. Sets flags related to periodic boundaries and processor      !
!  interfaces.                                                         !
!                                                                      !
!----------------------------------------------------------------------!
      subroutine desmpi_init

      use particle_filter, only: DES_INTERP_GARG

      use desmpi, only: iGhostPacketSize
      use desmpi, only: iParticlePacketSize
      use desmpi, only: iPairPacketSize
! Particle orientation
      use discretelement, only: PARTICLE_ORIENTATION
      use discretelement, only: DES_USR_VAR_SIZE

      use error_manager

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lfaces
      integer :: lmaxlen1,lmaxlen2,lmaxarea,lmaxghostpar,ii

      DOUBLE PRECISION, PARAMETER :: ONEMBo8 = 131072.0

!-----------------------------------------------

! Calculate the size of ghost particle packets:
      iGhostPacketSize = 15 + DES_USR_VAR_SIZE
      IF(ENERGY_EQ) &
         iGhostPacketSize = iGhostPacketSize + 1

! Calculate the size of particle packets.
      iParticlePacketSize = 30 + DES_USR_VAR_SIZE
      IF(ENERGY_EQ) &
         iParticlePacketSize = iParticlePacketSize + 1 + DIMENSION_N_s
      IF(DO_OLD) &
         iParticlePacketSize = iParticlePacketSize + 15
      IF(DO_OLD .AND. ENERGY_EQ) &
         iParticlePacketSize = iParticlePacketSize + 1
      IF(MPPIC) &
         iParticlePacketSize = iParticlePacketSize + 1
      IF(PARTICLE_ORIENTATION) &
         iParticlePacketSize = iParticlePacketSize + 3
      IF(DES_EXPLICITLY_COUPLED) THEN
         iParticlePacketSize = iParticlePacketSize + 3
         IF(ENERGY_EQ)iParticlePacketSize = iParticlePacketSize + 1
         IF(ANY_SPECIES_EQ)iParticlePacketSize = iParticlePacketSize + 1
      ENDIF

! Calculate the size of neighbor data
      iPairPacketSize = 11

! Calculate the initial size of send and recv buffer based on max_pip,
! total cells max number of boundary cells, and ghost par packet size.
      lfaces = dimn*2
      lmaxlen1 = dg_iend2-dg_istart2+1
      lmaxlen2 = dg_jend2-dg_jstart2+1
      if (do_K) then
         lmaxlen1 = max(lmaxlen1,dg_kend2 -dg_kstart2+1)
         lmaxlen2 = max(lmaxlen2,dg_kend2 -dg_kstart2+1)
      else
         lmaxlen1 = max(lmaxlen1,lmaxlen2)
         lmaxlen2 = 1
      end if


! Note: 10 is added for buffer and is required for send and recv indices
      lmaxarea = lmaxlen1*lmaxlen2 + 10

! Random value. This gets resized by DESMPI_CHECK_SENDRECVBUF
      lmaxghostpar = 100
      imaxbuf = lmaxghostpar*lmaxarea*iGhostPacketSize

      WRITE(ERR_MSG, 1000) iMAXBUF/ONEMBo8, ONEMBo8/iGhostPacketSize,  &
         ONEMBo8/iParticlePacketSize, ONEMBo8/iPairPacketSize
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT(/'DES MPI send/recv buffer: ',F7.1,' MB',/' o  ',F6.0,1X, &
         'Ghost Particles/MB',/' o  ',F6.0,1X,'Particles/MB',/' o  ',  &
         F6.0,1X,'Neighbor Pairs/MB')


      allocate (dsendbuf(2));
      allocate (drecvbuf(2));
      do ii=1, size(dsendbuf)
         allocate (dsendbuf(ii)%facebuf(imaxbuf));
         allocate (drecvbuf(ii)%facebuf(imaxbuf));
      end do

      allocate (isendindices(lmaxarea,lfaces)); isendindices=0
      allocate (irecvindices(lmaxarea,lfaces)); irecvindices=0

      allocate (isendreq(lfaces)); isendreq=0
      allocate (irecvreq(lfaces)); irecvreq=0
      allocate (isendcnt(lfaces)); isendcnt=0

      allocate (dcycl_offset(lfaces,dimn)); dcycl_offset=0.0
      allocate (ineighproc(lfaces)); ineighproc=0
      allocate (iexchflag(lfaces)); iexchflag=.FALSE.

! allocate variables related to scattergather
      allocate(iscattercnts(0:numpes-1)); iscattercnts=0
      allocate(igathercnts(0:numpes-1));  igathercnts=0
      allocate(idispls(0:numpes-1)); idispls=0

! set the communication flags
      CALL DESMPI_SETCOMM

!      call des_dbgmpi(1)

      end subroutine desmpi_init


!------------------------------------------------------------------------
! Subroutine       : desmpi_setcomm
! Purpose          : sets the flags required for interprocessor communication
!------------------------------------------------------------------------
      subroutine desmpi_setcomm()
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lijkproc,liproc,ljproc,lkproc
      integer li,lj,lk,lis,lie,ljs,lje,lks,lke,lcount,lface
      integer listart1,liend1,ljstart1,ljend1,lkstart1,lkend1
      integer listart2,liend2,ljstart2,ljend2,lkstart2,lkend2
!-----------------------------------------------

! set flags for interprocessor boundaries and set the corresponding to proc
      lijkproc = mype
      liproc = iofproc(lijkproc)
      ljproc = jofproc(lijkproc)
      lkproc = kofproc(lijkproc)
      iexchflag(:) = .false.
      ineighproc(:) = 0
      if(liproc.gt.0) then
         iexchflag(2) = .true.
         ineighproc(2)=procijk(liproc-1,ljproc,lkproc)
      end if
      if(liproc.lt.nodesi-1) then
         iexchflag(1) = .true.
         ineighproc(1)=procijk(liproc+1,ljproc,lkproc)
      end if
      if(ljproc.gt.0) then
         iexchflag(4)= .true.
         ineighproc(4)=procijk(liproc,ljproc-1,lkproc)
      end if
      if(ljproc.lt.nodesj-1) then
         iexchflag(3) = .true.
         ineighproc(3)=procijk(liproc,ljproc+1,lkproc)
      end if
      if(lkproc.gt.0) then
         iexchflag(6)=.true.
         ineighproc(6)=procijk(liproc,ljproc,lkproc-1)
      end if
      if(lkproc.lt.nodesk-1) then
         iexchflag(5) =.true.
         ineighproc(5)=procijk(liproc,ljproc,lkproc+1)
      end if

!set flags for cyclic boundary conditions and corresponding to proc
      dcycl_offset(:,:) = 0
      if (des_periodic_walls_x) then
         if(liproc.eq.0) then
            iexchflag(2)=.true.
            ineighproc(2)= procijk(nodesi-1,ljproc,lkproc)
            dcycl_offset(2,1)= xlength
         end if
         if(liproc.eq.nodesi-1) then
            iexchflag(1)=.true.
            ineighproc(1)= procijk(0,ljproc,lkproc)
            dcycl_offset(1,1)=-xlength
         end if
      end if
      if (des_periodic_walls_y) then
         if(ljproc.eq.0) then
            iexchflag(4)=.true.
            ineighproc(4)= procijk(liproc,nodesj-1,lkproc)
            dcycl_offset(4,2)= ylength
         end if
         if(ljproc.eq.nodesj-1) then
            iexchflag(3)=.true.
            ineighproc(3)= procijk(liproc,0,lkproc)
            dcycl_offset(3,2)=-ylength
         end if
      end if
      if (des_periodic_walls_z) then
         if(lkproc.eq.0) then
            iexchflag(6)=.true.
            ineighproc(6)= procijk(liproc,ljproc,nodesk-1)
            dcycl_offset(6,3)=zlength
         end if
         if(lkproc.eq.nodesk-1) then
            iexchflag(5)=.true.
            ineighproc(5)= procijk(liproc,ljproc,0)
            dcycl_offset(5,3)=-zlength
         end if
      end if

      listart1=dg_istart1; liend1=dg_iend1
      listart2=dg_istart2; liend2=dg_iend2

      ljstart1=dg_jstart1; ljend1=dg_jend1
      ljstart2=dg_jstart2; ljend2=dg_jend2

      lkstart1=dg_kstart1; lkend1=dg_kend1
      lkstart2=dg_kstart2; lkend2=dg_kend2

! Extend the domain indices to account for mass inlets and outlets. Do
! not extend the domain for periodic walls becuase 1) they should not
! include inflows or outflows and 2) they are already expanded.
      IF(.NOT.DES_PERIODIC_WALLS_X) THEN
         if(listart1.eq.dg_imin1) listart1 = dg_imin1-1
         if(liend1.eq.dg_imax1) liend1 = dg_imax1+1
      ENDIF

      IF(.NOT.DES_PERIODIC_WALLS_Y) THEN
         if(ljstart1.eq.dg_jmin1) ljstart1 = dg_jmin1-1
         if(ljend1.eq.dg_jmax1) ljend1 = dg_jmax1+1
      ENDIF

      IF(DO_K .AND. .NOT.DES_PERIODIC_WALLS_Z) THEN
         if(lkstart1.eq.dg_kmin1) lkstart1 = dg_kmin1-1
         if(lkend1.eq.dg_kmax1) lkend1 = dg_kmax1+1
      ENDIF

! set the ghost cell indices for e-w, n-s and t-b
! for east and west faces
      lks = lkstart1
      lke = lkend1
      ljs = ljstart1
      lje = ljend1

!east face
      lface = 1
      li = liend1
      lcount = 1
      do lk = lks,lke
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li+1,lj,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

!west face
      lface = 2
      li = listart1
      lcount = 1
      do lk = lks,lke
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li-1,lj,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

! for north and south faces
      lks = lkstart1
      lke = lkend1
      lis = listart2
      lie = liend2

!north face
      lface = 3
      lj = ljend1
      lcount = 1
      do lk = lks,lke
         do li = lis,lie
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj+1,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1


!south face
      lface = 4
      lj = ljstart1
      lcount = 1
      do lk = lks,lke
         do li = lis,lie
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj-1,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

! for top and bottom
      if (no_k) return
      lis = listart2
      lie = liend2
      ljs = ljstart2
      lje = ljend2

!top face
      lface = 5
      lk = lkend1
      lcount = 1
      do li = lis,lie
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj,lk+1)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1


!bottom face
      lface = 6
      lk = lkstart1
      lcount = 1
      do li = lis,lie
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj,lk-1)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

      return
      end subroutine desmpi_setcomm


!------------------------------------------------------------------------
! Subroutine       : des_scatter_particle
! Purpose          : Scatters the particles read from input file or
!                    generated based on specified volume fraction
! Comments         : In the main thread (pe_io) particles will be seperated
!                    based on the location and it will be packed in drootbuf
!                    Each proc receives its particles along with particle count
!------------------------------------------------------------------------
      subroutine des_scatter_particle

      use mpi_comm_des, only: desmpi_scatterv
      use des_allocate, only: particle_grow

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lcurpar,lproc,lbuf,lpacketsize
      integer lproc_parcnt(0:numpes-1),lpar_proc(particles)
!-----------------------------------------------

      integer :: rdimn

      rdimn = merge(2,3, NO_K)

! set the packet size for transfer
      lpacketsize = 2*rdimn + 2

! build the send buffer in PE_IO proc
! first pass to get the count of particles
      lpar_proc(:) =-1
      lproc_parcnt(:) = 0
      if(myPE.eq.pe_io) then
         if (no_k) then
            do lcurpar = 1,particles
               do lproc= 0,numpes-1
                  if (   dpar_pos(lcurpar,1).ge.xe(istart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,1).lt.xe(iend1_all(lproc))     &
                   .and. dpar_pos(lcurpar,2).ge.yn(jstart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,2).lt.yn(jend1_all(lproc))) then
                     lpar_proc(lcurpar) = lproc
                     lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                     exit
                  endif
               enddo
               if (lpar_proc(lcurpar).eq.-1) then
                  WRITE(*,500) lcurpar
                  call des_mpi_stop
               endif
            enddo
         else
            do lcurpar = 1,particles
               do lproc= 0,numpes-1
                  if (   dpar_pos(lcurpar,1).ge.xe(istart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,1).lt.xe(iend1_all(lproc))     &
                   .and. dpar_pos(lcurpar,2).ge.yn(jstart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,2).lt.yn(jend1_all(lproc))     &
                   .and. dpar_pos(lcurpar,3).ge.zt(kstart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,3).lt.zt(kend1_all(lproc))) then
                     lpar_proc(lcurpar) = lproc
                     lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                     exit
                  end if
               end do
               if (lpar_proc(lcurpar).eq.-1) then
                  WRITE(*,501) lcurpar
                  call des_mpi_stop
               endif
            enddo
         endif  ! if (no_k)
      endif ! if (my_pe.eq.pe_io)
      call bcast(lproc_parcnt(0:numpes-1),pe_io)

! second pass: set and allocate scatter related variables
      pip = lproc_parcnt(mype)
      call PARTICLE_GROW(pip)
      max_pip = max(pip,max_pip)
      iscr_recvcnt = pip*lpacketsize
      allocate (dprocbuf(iscr_recvcnt))
      if (mype.eq.pe_io) then
         allocate (drootbuf(particles*lpacketsize))
      else
         allocate (drootbuf(10))
      endif

! in the IO processor build the drootbuffer and idispls required
! for mpi communication
      if(mype.eq.pe_io) then
         idispls(0) = 0
         iscattercnts(0) = lproc_parcnt(0)*lpacketsize
         do lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)
            iscattercnts(lproc) = lproc_parcnt(lproc)*lpacketsize
         end do
         lproc_parcnt(:) = 0
         do lcurpar = 1,particles
            lproc = lpar_proc(lcurpar)
            lbuf = idispls(lproc)+lproc_parcnt(lproc)*lpacketsize+1
            drootbuf(lbuf) = dpar_rad(lcurpar); lbuf = lbuf + 1
            drootbuf(lbuf) = dpar_den(lcurpar); lbuf = lbuf + 1
            drootbuf(lbuf:lbuf+rdimn-1) = dpar_pos(lcurpar,1:rdimn); lbuf = lbuf + rdimn
            drootbuf(lbuf:lbuf+rdimn-1) = dpar_vel(lcurpar,1:rdimn); lbuf = lbuf + rdimn
            lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
         enddo
      endif
      call desmpi_scatterv(ptype=2)

! unpack the particles in each processor and set the pip
      do lcurpar = 1,pip
         lbuf = (lcurpar-1)*lpacketsize+1
         des_radius(lcurpar) = dprocbuf(lbuf); lbuf = lbuf+1
         ro_sol(lcurpar) = dprocbuf(lbuf); lbuf = lbuf+1
         des_pos_new(lcurpar,1:rdimn) = dprocbuf(lbuf:lbuf+rdimn-1); lbuf = lbuf+rdimn
         des_vel_new(lcurpar,1:rdimn) = dprocbuf(lbuf:lbuf+rdimn-1); lbuf = lbuf+rdimn
         call set_normal(lcurpar)
      enddo
      deallocate (dprocbuf,drootbuf)

 500  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 501  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')

      RETURN
      end subroutine des_scatter_particle



!------------------------------------------------------------------------
! Subroutine       : DES_RESTART_GHOST
! Purpose          : restart file contains neighbour information in terms
!                    global id. This routine converts the global id into
!                    local particle number
!                    steps
!                    1. Exchange the ghost particles (does not involve any neighbour info)
!                    2. loop through particles neighbour and contact list
!                       and convert global numbers to local numbers
! Parameters       : none
!------------------------------------------------------------------------

      subroutine DES_RESTART_GHOST

      use mpi_comm_des, only: desmpi_sendrecv_init
      use mpi_comm_des, only: desmpi_sendrecv_wait

      use mpi_funs_des, only: desmpi_check_sendrecvbuf
      use mpi_pack_des, only: desmpi_pack_ghostpar
      use mpi_unpack_des, only: desmpi_unpack_ghostpar

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer linter,lface
      integer ii
!-----------------------------------------------
! set do_nsearch true so that the ghost cell will be updated
      do_nsearch = .true.
      call desgrid_pic(plocate=.true.)
      call desmpi_check_sendrecvbuf(check_global=.true.)

!call ghost particle exchange in E-W, N-S, T-B order

      do ii=1, size(dsendbuf)
         dsendbuf(ii)%facebuf(1) = 0
         drecvbuf(ii)%facebuf(1) = 0
      end do

      ighost_updated(:) = .false.
      ispot = 1
      do linter = 1,dimn
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface))cycle
            call desmpi_pack_ghostpar(lface)
            call desmpi_sendrecv_init(lface)
         enddo
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface)) cycle
            call desmpi_sendrecv_wait(lface)
            call desmpi_unpack_ghostpar(lface)
         enddo
! update pic required as particles in ghost cell can move between ghost cells
         do lface = linter*2-1,linter*2
            if(dsendbuf(1+mod(lface,2))%facebuf(1).gt.0.or.drecvbuf(1+mod(lface,2))%facebuf(1).gt.0) then
               call desgrid_pic(plocate=.false.)
               exit
            endif
         enddo
      enddo
      call des_mpi_barrier

      end subroutine DES_RESTART_GHOST

      end module mpi_init_des
