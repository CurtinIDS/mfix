!----------------------------------------------------------------------!
!  Module: MPI_UNPACK_DES                                              !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for unpacking real and ghost particles   !
!     from the MPI recv buffers.                                       !
!----------------------------------------------------------------------!
      MODULE MPI_UNPACK_DES

      PRIVATE
      PUBLIC :: DESMPI_UNPACK_PARCROSS, DESMPI_UNPACK_GHOSTPAR

      interface unpack_dbuf
         module procedure unpack_db0 ! real scalars
         module procedure unpack_db1 ! real arrays
         module procedure unpack_i0  ! integer scalars
         module procedure unpack_i1  ! integer arrays
         module procedure unpack_l0  ! logical scalars
      end interface unpack_dbuf

      CONTAINS

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_UNPACK_GHOSTPAR                                  !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Unpacks ghost particle from the recv buffer.                !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_UNPACK_GHOSTPAR(pface)

! Global Variables:
!---------------------------------------------------------------------//
! Size of Particle data packet
      use desmpi, only: iGhostPacketSize
! Index of last particle added to this process.
      use desmpi, only: iSPOT
! Flag indicating that the ghost particle was updated
      use discretelement, only: iGHOST_UPDATED
! The MPI receive buffer
      use desmpi, only: dRECVBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Dimensions of DES grid
      use desgrid, only: DG_IJKSIZE2
! DES grid cell containing each particle: current/previous
      use discretelement, only: DG_PIJK, DG_PIJKPRV
! The global ID for each particle
      use discretelement, only: iGLOBAL_ID
! Particle positions: current/previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! Particle tangential velocities: current/previous
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Particle rotational velocities: current/previous
      use discretelement, only: OMEGA_NEW, OMEGA_OLD
! Particle tempertures
      use des_thermo, only: DES_T_s
! Particle radius, volume
      use discretelement, only: DES_RADIUS, PVOL
! Map to fluid grid cells and solids phase (I,J,K,IJK,M)
      use discretelement, only: PIJK
! Flag to send/recv old (previous) values
      use discretelement, only: DO_OLD
! Number of particles on the process (max particle array size)
      use discretelement, only: PIP
! Number of ghost particles on the current process
      use discretelement, only: iGHOST_CNT
! User-defined variables for each particle.
      use discretelement, only: DES_USR_VAR, DES_USR_VAR_SIZE

      use des_allocate

! Global Constants:
!---------------------------------------------------------------------//
      use constant, only: PI
! Dimension of particle spatial arrays.
      use discretelement, only: DIMN
      use discretelement, only: max_pip

      use functions, only: is_nonexistent
      use functions, only: is_normal,  set_normal
      use functions, only: is_exiting, set_exiting, set_exiting_ghost
      use functions, only: set_ghost

      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lcurpar,lparid,lprvijk,lparijk,lparcnt,ltot_ind
      integer :: lbuf,llocpar,lnewcnt,lpicloc
      logical,dimension(:),allocatable :: lfound
      integer,dimension(:),allocatable :: lnewspot,lnewpic
      logical :: tmp
!......................................................................!

! unpack the particles:
! if it already exists update the position
! if not and do_nsearch is true then add to the particle array

      lparcnt = drecvbuf(1+mod(pface,2))%facebuf(1)
      lnewcnt = lparcnt
      allocate (lfound(lparcnt),lnewspot(lparcnt),lnewpic(dg_ijksize2))
      lfound(:) = .false.
      lnewspot(:) =0
      lnewpic = 0

      do lcurpar = 1,lparcnt
         lbuf = (lcurpar-1)*iGhostPacketSize+ibufoffset

! 1) Global ID
         call unpack_dbuf(lbuf,lparid,pface)
! 2) DES Grid IJK
         call unpack_dbuf(lbuf,lparijk,pface)
! 3) DES Grid IJK - Previous
         call unpack_dbuf(lbuf,lprvijk,pface)

! Determine if this particle already exists on this process as a
! ghost particle. If so, (lfound), then the current infomration is
! updated on the current process. Otherwise (.NOT.lfound) a new
! ghost particle is created on this process.
         lfound(lcurpar) = locate_par(lparid,lprvijk,llocpar)
         if (lparijk .ne. lprvijk .and. .not.lfound(lcurpar)) then
            lfound(lcurpar) = locate_par(lparid,lparijk,llocpar)
         endif

         if(lfound(lcurpar)) then
! Store the local variables
            dg_pijk(llocpar) = lparijk
            dg_pijkprv(llocpar) = lprvijk

! 4) Radious
            call unpack_dbuf(lbuf,des_radius(llocpar),pface)
! 5) Phase index
            call unpack_dbuf(lbuf,pijk(llocpar,5),pface)
! 6) Position
            call unpack_dbuf(lbuf,des_pos_new(llocpar,1:dimn),pface)
! 7) Translational Velocity
            call unpack_dbuf(lbuf,des_vel_new(llocpar,1:dimn),pface)
! 8) Rotational Velocity
            call unpack_dbuf(lbuf,omega_new(llocpar,1:3),pface)
! 9) Exiting particle flag
            call unpack_dbuf(lbuf,tmp,pface)
            if (tmp) call set_exiting_ghost(llocpar)
! 10) Temperature
            IF(ENERGY_EQ) &
               call unpack_dbuf(lbuf,des_t_s(llocpar),pface)
! 11) User Variables
            IF(DES_USR_VAR_SIZE > 0) &
               call unpack_dbuf(lbuf,des_usr_var(:,llocpar),pface)

! Calculate the volume of the ghost particle.
            PVOL(llocpar) = (4.0D0/3.0D0)*PI*DES_RADIUS(llocpar)**3
! Flag that the ghost particle was updated.
            ighost_updated(llocpar) = .true.
            lnewcnt = lnewcnt-1

! Copy the current value to the previous value if needed.
            IF (DO_OLD) THEN
               des_pos_old(llocpar,:)= des_pos_new(llocpar,:)
               des_vel_old(llocpar,:)= des_vel_new(llocpar,:)
               omega_old(llocpar,:)= omega_new(llocpar,:)
            ENDIF

         else
            lnewpic(lparijk) = lnewpic(lparijk) + 1
         endif
      enddo

! iAdd new ghost particles
      if(lnewcnt > 0) then
         call PARTICLE_GROW(pip+lnewcnt)
         ighost_cnt = ighost_cnt + lnewcnt
         pip = pip + lnewcnt
         max_pip = max(pip,max_pip)
         do lcurpar = 1,lparcnt
            if(lfound(lcurpar)) cycle
            lbuf = (lcurpar-1)*iGhostPacketSize+ibufoffset

!  1) Global particle ID
            call unpack_dbuf(lbuf,lparid,pface)
!  2) DES grid IJK
            call unpack_dbuf(lbuf,lparijk,pface)
!  3) DES grid IJK - Previous
            call unpack_dbuf(lbuf,lprvijk,pface)
! Locate the first open space in the particle array.
            do while(.not.is_nonexistent(ispot))
               ispot = ispot + 1
            enddo
! Set the flags for the ghost particle and store the local variables.
            call set_ghost(ispot)
            iglobal_id(ispot)  = lparid
            dg_pijk(ispot) = lparijk
            dg_pijkprv(ispot) = lprvijk
!  4) Radius
            call unpack_dbuf(lbuf,des_radius(ispot),pface)
!  5) Phase index
            call unpack_dbuf(lbuf,pijk(ispot,5),pface)
!  6) Position
            call unpack_dbuf(lbuf,des_pos_new(ispot,1:dimn),pface)
!  7) Translational velocity
            call unpack_dbuf(lbuf,des_vel_new(ispot,1:dimn),pface)
!  8) Rotational velocity
            call unpack_dbuf(lbuf,omega_new(ispot,1:dimn),pface)
!  9) Exiting particle flag
            call unpack_dbuf(lbuf,tmp,pface)
            if (tmp) call set_exiting_ghost(ispot)
! 10) Temperature.
            IF(ENERGY_EQ) &
               call unpack_dbuf(lbuf,des_t_s(ispot),pface)
! 11) User varaible
            IF(DES_USR_VAR_SIZE > 0)&
               call unpack_dbuf(lbuf,des_usr_var(:,ispot),pface)

            ighost_updated(ispot) = .true.
            lnewspot(lcurpar) = ispot

            PVOL(ispot) = (4.0D0/3.0D0)*PI*DES_RADIUS(ispot)**3

            IF (DO_OLD) THEN
               des_pos_old(ispot,1:dimn) = des_pos_new(ispot,1:dimn)
               des_vel_old(ispot,1:dimn) = des_vel_new(ispot,1:dimn)
               omega_old(ispot,1:3) = omega_new(ispot,1:3)
            ENDIF
         enddo
      endif

!deallocate temporary variablies
      deallocate (lfound,lnewspot,lnewpic)

      end subroutine desmpi_unpack_ghostpar

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_UNPACK_PARCROSS                                  !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Unpacks real particle from the recv buffer.                 !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_UNPACK_PARCROSS(pface)

! Global Variables:
!---------------------------------------------------------------------//
! Size of ghost particle data packet
      use desmpi, only: iParticlePacketSize
! Index of last particle added to this process.
      use desmpi, only: iSPOT
! The MPI receive buffer
      use desmpi, only: dRECVBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Runtime flag for solving species equations
      use run, only: ANY_SPECIES_EQ
! Runtime flag for MPPIC solids
      use mfix_pic, only: MPPIC
! DES grid cell containing each particle: current/previous
      use discretelement, only: DG_PIJK, DG_PIJKPRV
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! The statistical weight of each particle.
      use mfix_pic, only: DES_STAT_WT
! The global ID for each particle
      use discretelement, only: iGLOBAL_ID
! Particle positions: current/previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! Particle tangential velocities: current/previous
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Particle rotational velocities: current/previous
      use discretelement, only: OMEGA_NEW, OMEGA_OLD
!Particle orientation
      use discretelement, only: PARTICLE_ORIENTATION,ORIENTATION
! Particle radius, volume, density, mass
      use discretelement, only: DES_RADIUS, PVOL, RO_SOL, PMASS
! Previous value for particle acceleration (tangential/rotational)
      use discretelement, only: DES_ACC_OLD, ROT_ACC_OLD
! Particle species composition
      use des_rxns, only: DES_X_s
! Particle tempertures.
      use des_thermo, only: DES_T_s
! Force arrays acting on the particle
      use discretelement, only: FC, TOW
! One of the moment of inertia
      use discretelement, only: OMOI
! Map to fluid grid cells and solids phase (I,J,K,IJK,M)
      use discretelement, only: PIJK
! Flag to send/recv old (previous) values
      use discretelement, only: DO_OLD
! Number of particles on the process (max particle array size)
      use discretelement, only: PIP
! Number of ghost particles on the current process
      use discretelement, only: iGHOST_CNT
! Flag indicating the the fluid-particle drag is explicitly coupled.
      use discretelement, only: DES_EXPLICITLY_COUPLED
! Explicit fluid-particle drag force
      use discretelement, only: DRAG_FC
! Explict convection and HOR
      use des_thermo, only: CONV_Qs, RXNS_Qs
! User-defined variables for each particle.
      use discretelement, only: DES_USR_VAR, DES_USR_VAR_SIZE
! Neighbor collision history information
      use discretelement, only: PFT_NEIGHBOR
! Dimension of particle spatial arrays.
      use discretelement, only: DIMN
! The ID of the current process
      use compar, only: myPE

! Module Procedures:
!---------------------------------------------------------------------//
      use des_allocate
      use desmpi_wrapper, only: DES_MPI_STOP
      use discretelement, only: max_pip
      use functions, only: IS_NORMAL, IS_NONEXISTENT
      use functions, only: SET_ENTERING, SET_EXITING, SET_NORMAL

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lcurpar,lparcnt,llocpar,lparid,lparijk,lprvijk
      integer :: lneigh,lcontactindx,lcontactid,lcontact,&
                 lneighid,lneighijk
      logical :: lfound
      integer :: lbuf,lcount
      logical :: lneighfound
      integer :: cc,kk,num_neighborlists_sent,nn

      logical :: tmp
!......................................................................!

! loop through particles and locate them and make changes
      lparcnt = drecvbuf(1+mod(pface,2))%facebuf(1)

! if mppic make sure enough space available
      call PARTICLE_GROW(pip+lparcnt)
      max_pip = max(pip+lparcnt,max_pip)

      do lcurpar =1,lparcnt

         lfound = .false.
         lbuf = (lcurpar-1)*iParticlePacketSize + ibufoffset
! 1) Global ID
         call unpack_dbuf(lbuf,lparid,pface)
! 2) DES Grid IJK
         call unpack_dbuf(lbuf,lparijk,pface)
! 3) DES grid IJK - previous
         call unpack_dbuf(lbuf,lprvijk,pface)

! PIC particles are always 'new' to the receiving process. Find the
! first available array position and store the global ID. Increment
! the PIP counter to include the new particle.
         IF(MPPIC) THEN
            DO WHILE(.NOT.IS_NONEXISTENT(ISPOT))
               ISPOT = ISPOT + 1
            ENDDO
            lLOCPAR = iSPOT
            iGLOBAL_ID(lLOCPAR) = lPARID
            PIP = PIP + 1

! A DEM particle should already exist on the current processor as a
! ghost particle. Match the sent particle to the local ghost particle
! by matching the global IDs. Decrement the iGHOST_CNT counter to
! account for the switch from ghost to real particle.
         ELSE
            lFOUND  = LOCATE_PAR(lPARID,lPRVIJK,lLOCPAR)
            IF (.NOT. lFOUND) THEN
               lFOUND = exten_locate_par(lPARID, lPARIJK, lLOCPAR)
               IF(.NOT.lFOUND) THEN
                  WRITE(*,1000) iNEIGHPROC(PFACE), MYPE, lPARID
                  CALL DES_MPI_STOP
               ENDIF
            ENDIF
            iGHOST_CNT = iGHOST_CNT - 1
         ENDIF

 1000 FORMAT(2/1X,72('*'),/1x,'From: DESMPI_UNPACK_PARCROSS: ',/       &
         ' Error 1000: Unable to match particles crossing processor ', &
         'boundaries.',/3x,'Source Proc: ',I9,' ---> Destination ',    &
         'Proc: ', I9,/3x,'Global Particle ID: ',I12,/1x,72('*'))

! convert the local particle from ghost to existing and update its position
         call set_normal(llocpar)
         dg_pijk(llocpar) = lparijk
         dg_pijkprv(llocpar) = lprvijk
! 4) Radius
         call unpack_dbuf(lbuf,des_radius(llocpar),pface)
! 5-9) Fluid cell I,J,K,IJK, and solids phase index
         call unpack_dbuf(lbuf,pijk(llocpar,:),pface)
! 10) Entering particle flag.
         call unpack_dbuf(lbuf,tmp,pface)
         if (tmp) CALL SET_ENTERING(llocpar)
! 11) Exiting particle flag.
         call unpack_dbuf(lbuf,tmp,pface)
         if (tmp) CALL SET_EXITING(llocpar)
! 12) Density
         call unpack_dbuf(lbuf,ro_sol(llocpar),pface)
! 13) Volume
         call unpack_dbuf(lbuf,pvol(llocpar),pface)
! 14) Mass
         call unpack_dbuf(lbuf,pmass(llocpar),pface)
! 15) 1/Moment of Inertia
         call unpack_dbuf(lbuf,omoi(llocpar),pface)
! 16) Position with cyclic shift
         call unpack_dbuf(lbuf,des_pos_new(llocpar,:),pface)
! 17) Translational velocity
         call unpack_dbuf(lbuf,des_vel_new(llocpar,:),pface)
! 18) Rotational velocity
         call unpack_dbuf(lbuf,omega_new(llocpar,:),pface)
! 19) Accumulated translational forces
         call unpack_dbuf(lbuf,fc(llocpar,:),pface)
! 20) Accumulated torque forces
         call unpack_dbuf(lbuf,tow(llocpar,:),pface)
         IF(ENERGY_EQ) THEN
! 21) Temperature
            call unpack_dbuf(lbuf,des_t_s(llocpar),pface)
! 22) Species composition
            call unpack_dbuf(lbuf,des_x_s(llocpar,:),pface)
          ENDIF
! 23) User defined variable
         IF(DES_USR_VAR_SIZE > 0) &
            call unpack_dbuf(lbuf,des_usr_var(:,llocpar),pface)
! 24) Particle orientation
         IF(PARTICLE_ORIENTATION) &
            call unpack_dbuf(lbuf,orientation(:,llocpar),pface)

! -- Higher order integration variables
         IF (DO_OLD) THEN
! 25) Position (previous)
            call unpack_dbuf(lbuf,des_pos_old(llocpar,:),pface)
! 26) Translational velocity (previous)
            call unpack_dbuf(lbuf,des_vel_old(llocpar,:),pface)
! 27) Rotational velocity (previous)
            call unpack_dbuf(lbuf,omega_old(llocpar,:),pface)
! 28) Translational acceleration (previous)
            call unpack_dbuf(lbuf,des_acc_old(llocpar,:),pface)
! 29) Rotational acceleration (previous)
            call unpack_dbuf(lbuf,rot_acc_old(llocpar,:),pface)
         ENDIF

         IF(DES_EXPLICITLY_COUPLED) THEN
! 30) Explicit drag force
            call unpack_dbuf(lbuf,drag_fc(llocpar,:),pface)
! 31) Explicit convective heat transfer
            IF(ENERGY_EQ) call unpack_dbuf(lbuf,conv_qs(llocpar),pface)
! 32) Explicit heat of reaction
            IF(ANY_SPECIES_EQ) call unpack_dbuf(lbuf,rxns_qs(llocpar),pface)
         ENDIF

! 33) Statistical weight
         IF(MPPIC) call unpack_dbuf(lbuf,des_stat_wt(llocpar),pface)

      end do

! 34) Number of pair datasets
      lbuf = lparcnt*iParticlePacketSize + ibufoffset
      call unpack_dbuf(lbuf,num_neighborlists_sent,pface)

      do nn = 1, num_neighborlists_sent
! 35) Global ID of packed particle.
         call unpack_dbuf(lbuf,lparid,pface)
! 36) DES grid IJK of cell receiving the particle.
         call unpack_dbuf(lbuf,lparijk,pface)

! Locate the particle on the current process.
         if (.not. locate_par(lparid,lparijk,llocpar)) then
            if (.not. exten_locate_par(lparid,lparijk,llocpar)) then
               print *,"at buffer location",lbuf," pface = ",pface
               print *,"COULD NOT FIND PARTICLE ",lparid," IN IJK ",lparijk
               call des_mpi_stop
            endif
         endif
! 37) Global ID of neighbor particle.
         call unpack_dbuf(lbuf,lneighid,pface)
! 38) DES grid IJK of cell containing the neighbor particle.
         call unpack_dbuf(lbuf,lneighijk,pface)

! Locate the neighbor particle on the current process.
         if (.not. locate_par(lneighid,lneighijk,lneigh)) then
            if (.not. exten_locate_par(lneighid,lparijk,lneigh)) then
               print *,"  "
               print *,"  "
               print *," fail on  ", myPE
               print *,"at buffer location",lbuf," pface = ",pface
               print *,"COULD NOT FIND NEIGHBOR ",lneighid," IN IJK ",lneighijk
               call des_mpi_stop
            endif
         endif

! If the neighbor particle is a 'real' particle on this processor, then
! the pair data may already exist. Check before adding it.
! Create a new neighbor pair if it was not matched to an exiting pair.
          cc = add_pair(llocpar,lneigh)
! 39) Tangential collision history.
         call unpack_dbuf(lbuf,pft_neighbor(:,cc),pface)
      enddo

      END SUBROUTINE desmpi_unpack_parcross

!----------------------------------------------------------------------!
! Function: LOCATE_PAR                                                 !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Return the local index of the particle matching the passed  !
!    global ID. The function returns TRUE if the particle is matched,  !
!    otherwise it returns FALSE.                                       !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION LOCATE_PAR(pGLOBALID, pIJK, pLOCALNO)

      use discretelement, only: iGLOBAL_ID
      use desgrid, only: DG_IJKStart2, DG_IJKEnd2
      use derived_types, only: dg_pic

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Global ID of the particle
      INTEGER, INTENT(IN) :: pGlobalID
! IJK of DES grid cell containing the particle
      INTEGER, INTENT(IN) :: pIJK
! Local ID for the particle.
      INTEGER, INTENT(OUT) :: pLocalNO

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: lpicloc, lcurpar

! Initialize the result.
      locate_par = .false.

! Verify that the passied IJK value is within a valid range.
      if(pIJK < dg_ijkstart2 .or. pIJK > dg_ijkend2)  RETURN

! Loop the the particles in DES grid cell pIJK. Return to the calling
! routine if the passed global ID matches the global ID of one of
! the local particles.
      DO lpicloc = 1,dg_pic(pijk)%isize
         lcurpar = dg_pic(pijk)%p(lpicloc)
         IF(iGLOBAL_ID(lcurpar) == pGlobalID) THEN
            plocalno = lcurpar
            locate_par = .true.
            RETURN
         ENDIF
      ENDDO

      RETURN
      end function locate_par

!----------------------------------------------------------------------!
! Function: EXTEN_LOCATE_PAR                                           !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Return the local index of the particle matching the passed  !
!    global ID. The function returns TRUE if the particle is matched,  !
!    otherwise it returns FALSE.                                       !
!------------------------------------------------------------------------
      LOGICAL FUNCTION EXTEN_LOCATE_PAR(pGlobalID, pIJK, pLocalNO)

      use derived_types, only: dg_pic
      use discretelement, only: iGLOBAL_ID
      use desgrid, only: DG_IJKStart2, DG_IJKEnd2
      use desgrid, only: dg_Iof_LO, DG_Jof_LO, DG_Kof_LO
      use geometry, only: NO_K

      use desgrid, only: dg_funijk

      implicit none

! Dummy variables:
!---------------------------------------------------------------------//
! The global ID of the particle to be matched locally
      INTEGER, INTENT(IN) :: pGlobalId
! The DES grid cell index expected to contain the particle.
      INTEGER, INTENT(IN) :: pIJK
! The local ID of the matching particle.
      INTEGER, INTENT(OUT) :: pLocalNo

! Local variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: lijk, li, lj, lk, lic, ljc, lkc, lkoffset
      INTEGER :: lpicloc,lcurpar

      exten_locate_par = .false.

      lic = dg_iof_lo(pijk)
      ljc = dg_jof_lo(pijk)
      lkc = dg_kof_lo(pijk)
      lkoffset = merge(0, 1, NO_K)
      DO  lk = lkc-lkoffset,lkc+lkoffset
      DO  lj = ljc-1,ljc+1
      DO  li = lic-1,lic+1
         lijk = dg_funijk(li,lj,lk)
         IF (lijk .lt. dg_ijkstart2 .or. lijk .gt. dg_ijkend2) CYCLE
         DO lpicloc = 1, dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)
            IF (iglobal_id(lcurpar) .eq. pglobalid) THEN
               plocalno = lcurpar
               exten_locate_par = .true.
               RETURN
            END IF
         END DO
      END DO
      END DO
      END DO

      RETURN
      END FUNCTION EXTEN_LOCATE_PAR

!----------------------------------------------------------------------!
! Unpack subroutine for single real variables                          !
!----------------------------------------------------------------------!
      subroutine unpack_db0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf)
      lbuf = lbuf + 1

      return
      end subroutine unpack_db0

!----------------------------------------------------------------------!
! Unpack subroutine for real arrays                                    !
!----------------------------------------------------------------------!
      subroutine unpack_db1(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1)
      lbuf = lbuf + lsize

      return
      end subroutine unpack_db1

!----------------------------------------------------------------------!
! Unpack subroutine for single integer variables                       !
!----------------------------------------------------------------------!
      subroutine unpack_i0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf)
      lbuf = lbuf + 1

      return
      end subroutine unpack_i0

!----------------------------------------------------------------------!
! Unpack subroutine for integer arrays                                 !
!----------------------------------------------------------------------!
      subroutine unpack_i1(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1)
      lbuf = lbuf + lsize

      return
      end subroutine unpack_i1

!----------------------------------------------------------------------!
! Unpack subroutine for logical variables                              !
!----------------------------------------------------------------------!
      subroutine unpack_l0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      logical, intent(inout) :: idata

      idata = merge(.true.,.false.,0.5<drecvbuf(1+mod(pface,2))%facebuf(lbuf))
      lbuf = lbuf + 1

      return
      end subroutine unpack_l0

      end module mpi_unpack_des
