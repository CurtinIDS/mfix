!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM_MI                                           !
!  Author: J.Musser                                   Date: 23-Nov-09  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM_MI

! Modules
!---------------------------------------------------------------------//
      USE bc
      USE compar
      USE constant
      USE des_allocate
      USE des_bc
      USE discretelement
      USE error_manager
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run
      USE toleranc
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: BCV
      INTEGER BCV_I      ! BC loop counter
      INTEGER M, MM           ! Mass phase loop counter
      INTEGER RANGE_TOP, RANGE_BOT ! Dummy values
      INTEGER PHASE_CNT        ! Number of solid phases at bc
      INTEGER PHASE_LIST(DIM_M) ! List of phases used in current bc

! the number of particles injected in a solids time step
      DOUBLE PRECISION NPMpSEC(DIM_M) ! For solid phase m
      DOUBLE PRECISION NPpSEC
      DOUBLE PRECISION NPpDT        ! Total for BC
      DOUBLE PRECISION SCALED_VAL
      DOUBLE PRECISION MAX_DIA ! Max diameter of incoming particles at bc

      DOUBLE PRECISION :: EPs_ERR
      DOUBLE PRECISION :: VOL_FLOW

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag
      LOGICAL :: FATAL

! Temp inlet velocity for solids phase M
      DOUBLE PRECISION VEL_TMP(DIM_M)
      DOUBLE PRECISION EPs_TMP(DIM_M)

! Minimum/maximum solids velocity at inlet.  Also used in the iterative
! steps as the starting and ending velocities
      DOUBLE PRECISION  MAX_VEL

      DOUBLE PRECISION  MINIPV, MAXIPV
      INTEGER :: OCCUPANTS

!......................................................................!

      CALL INIT_ERR_MSG("SET_BC_DEM_MI")

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,2x,'DEM inlet count: ',I4)") DEM_BCMI

! Allocate the MI data structures for NEW and RES2 cases. Allocation for
! RES1 cases is done prior to reading the RES file.
      IF(RUN_TYPE /= 'RESTART_1') CALL ALLOCATE_DEM_MI

! Loop over BCs that flagged for DEM mass inflow.
      DO BCV_I = 1, DEM_BCMI

! Get the user defined BC ID.
         BCV = DEM_BCMI_MAP(BCV_I)

         if(dFlag) write(*,"(2/,'Setting DEM_MI:',I3)") BCV_I

! The number of mass phases at this inlet.  While a system may be
! polydisperse, the inlet could consist of a single mass phase
         PHASE_CNT = 0
! The mass phase indices of incoming particles at this inlet
         PHASE_LIST(:) = -1
! The max diameter of incoming particles at this inlet
         MAX_DIA = ZERO

! Determine if the inlet is mono or polydisperse
         DO M=1, MMAX + DES_MMAX
            IF(SOLIDS_MODEL(M) /= 'DEM') CYCLE
            IF(BC_ROP_s(BCV,M) == UNDEFINED) CYCLE
            IF(COMPARE(BC_ROP_s(BCV,M),ZERO)) CYCLE
            PHASE_CNT = PHASE_CNT + 1
            PHASE_LIST(PHASE_CNT) = M
            MAX_DIA = MAX(MAX_DIA,D_P0(M))
         ENDDO
! Set the polydispersity flag.
         DEM_MI(BCV_I)%POLYDISPERSE = (PHASE_CNT > 1)

! Layout the feed pattern.
         CALL LAYOUT_MI_DEM(BCV, BCV_I, MAX_DIA)

! Initialize
         MAX_VEL = ZERO
         NPMpSEC(:) = ZERO

! Calculate the individual velocities for each solid phase
         DO MM = 1, PHASE_CNT
            M = PHASE_LIST(MM)

! Pull off the BC velocity normal to the flow plane.
            SELECT CASE(BC_PLANE(BCV))
            CASE('N','S'); VEL_TMP(M) = abs(BC_V_s(BCV,M))
            CASE('E','W'); VEL_TMP(M) = abs(BC_U_s(BCV,M))
            CASE('T','B'); VEL_TMP(M) = abs(BC_W_s(BCV,M))
            END SELECT

! Check for min/max inlet velocity
            MAX_VEL = MAX(ABS(VEL_TMP(M)), MAX_VEL)

! Calculate the number of particles of mass phase M are injected per
! second for each solid phase present at the boundary.

! Use the mass flow rate if defined.
            IF(BC_MASSFLOW_s(BCV,M) == ZERO) THEN
               NPMPSEC(M) = ZERO

            ELSEIF(BC_MASSFLOW_s(BCV,M) /= UNDEFINED) THEN
               NPMPSEC(M) = BC_MASSFLOW_s(BCV,M) / &
                  (RO_s0(M)*(PI/6.d0*D_P0(M)**3))

! Otherwise use the volumetric flow rate if defined.
            ELSEIF(BC_VOLFLOW_S(BCV,M) == ZERO) THEN
               NPMPSEC(M) = ZERO

            ELSEIF(BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
               NPMpSEC(M) = BC_VOLFLOW_s(BCV,M) / (PI/6.d0*D_P0(M)**3)

! Calculate volumetric flow rate to convert to particle count. BC_AREA
! was already corrected for cut cells and velocity was recalculated
! to ensure user-specified mass or volumetric flow rates.
            ELSE
               VOL_FLOW = VEL_TMP(M) * BC_AREA(BCV) * BC_EP_S(BCV,M)
! Calculate the number of particles of mass phase M are injected per
! second for each solid phase present at the boundary
               NPMpSEC(M) = VOL_FLOW / (PI/6.d0*D_P0(M)**3)
! Write some debugging information if needed.
            ENDIF
            if(dFlag) write(*,1100) M, VEL_TMP(M), NPMpSEC(M)
         ENDDO

! Total number of particles at BCV injected per second
         NPpSEC = sum(NPMpSEC)

 1100 FORMAT(/2x,'Conversion Info: Phase ',I2,/4x,'Velocity: ',g12.5,/ &
         4X,'NPMpSEC = ',F11.1)

         if(dFlag) write(*,"(/2x,'Max Velocity:',3x,g12.5)") MAX_VEL
         if(dFlag) write(*,"( 2x,'NPpSEC:',3x,F11.1)") NPpSEC

! The number of total particles per solid time step DTSOLID
         NPpDT = NPpSEC * DTSOLID

! Inject one particle every PI_FACTOR solids time steps.
         IF(NPpDT == ZERO)THEN
            PI_COUNT(BCV_I) = 0
            PI_FACTOR(BCV_I) = UNDEFINED_I
         ELSEIF(NPpDT .LT. 1.0)THEN
            PI_COUNT(BCV_I) = 1
            PI_FACTOR(BCV_I) = FLOOR(real(1.d0/NPpDT))
! Inject PI_COUNT particles every soilds time step.
         ELSE
            PI_COUNT(BCV_I) = CEILING(real(NPpDT))
            PI_FACTOR(BCV_I) = 1
         ENDIF

         OCCUPANTS = DEM_MI(BCV_I)%OCCUPANTS

         IF(PI_COUNT(BCV_I) > 0) THEN
! Calculate the minimum inlet velocity. The cutoff is associated with
! square packing of disks on a plane.
            MINIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I)) * dble(  &
               FLOOR( real(OCCUPANTS)/real(PI_COUNT(BCV_I)))))
! Calculate the velocity needed to ensure that half the inlet is free.
! Inlets with velocities greater than this value can be randomly seeded,
! otherwise, particles are seeded in according to the grid.
            MAXIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I)) * dble(  &
             FLOOR(CEILING(real(OCCUPANTS)/2.0)/real(PI_COUNT(BCV_I)))))
         ELSE
            MINIPV = -UNDEFINED
            MAXIPV =  UNDEFINED
         ENDIF

         if(dFlag) write(*,"(/2x,'MaxIPV:',3x,g12.5)") MAXIPV
         if(dFlag) write(*,"( 2x,'MinIPV:',3x,g12.5)") MINIPV

         IF(MAX_VEL < MINIPV) THEN
            WRITE(ERR_MSG,1110) BCV, MAX_VEL, MINIPV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1110 FORMAT('Error 1110: Solids velocity for BC ',I3,' is too low ', &
         'to satisfy DEM',/'inlet packing restrictions. Potential ',  &
         'solutions:',//,' > If the BC velocities (BC_U_s, BC_V_s, ', &
         'BC_W_s) are defined, specify',/3x,'a larger value for the ',&
         'velocity normal to the flow plane.',//' > If MASSFLOW or ', &
         'VOLFLOW are defined, decrease the solids volume',/3x,       &
         'fraction to increase solids velocity.',//2x,'Max user-',    &
         'specified BC velocity:   ',g12.5,/2x,'Minimum required ',   &
         'solids Velocity: ',g12.5)

! Set all BC solids velocities to the largest velocity and recalculate
! BC_EP_s to determine the magnitude of the change.
         EPs_ERR = ZERO
         EPs_TMP = ZERO
         DO MM = 1, PHASE_CNT
            M = PHASE_LIST(MM)

            IF(BC_MASSFLOW_s(BCV,M) == ZERO .OR. &
               BC_VOLFLOW_S(BCV,M) == ZERO) THEN
! If no solids, set the tmp EPs to zero and clear solid velocity
               EPs_TMP(M) = 0.0d0
               SELECT CASE(BC_PLANE(BCV))
               CASE('N','S'); BC_V_s(BCV,M) =  ZERO
               CASE('E','W'); BC_U_s(BCV,M) =  ZERO
               CASE('T','B'); BC_W_s(BCV,M) =  ZERO
               END SELECT

            ELSE
               EPs_TMP(M) = BC_EP_s(BCV,M) * (VEL_TMP(M) / MAX_VEL)
! Over-write the current BC value.
               SELECT CASE(BC_PLANE(BCV))
               CASE('N'); BC_V_s(BCV,M) =  abs(MAX_VEL)
               CASE('S'); BC_V_s(BCV,M) = -abs(MAX_VEL)
               CASE('E'); BC_U_s(BCV,M) =  abs(MAX_VEL)
               CASE('W'); BC_U_s(BCV,M) = -abs(MAX_VEL)
               CASE('T'); BC_W_s(BCV,M) =  abs(MAX_VEL)
               CASE('B'); BC_W_s(BCV,M) = -abs(MAX_VEL)
               END SELECT

            ENDIF
            EPs_ERR = EPs_ERR + (BC_EP_s(BCV,M) - EPs_TMP(M))
         ENDDO

! If the net change in solids volume fraction is greatere than 0.01,
! flag this as an error and exit. >> Let the user fix the input.
         IF(.NOT.COMPARE(EPs_ERR,SMALL_NUMBER)) THEN
            IF(EPs_ERR > 0.01) THEN
               WRITE(ERR_MSG,1200) BCV, MAX_VEL, EPs_ERR
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               FATAL = .TRUE.

! Report the amount of changes imposed on the BC in setting a
! uniform inlet velocity.
            ELSE
               WRITE(ERR_MSG,1205) BCV, MAX_VEL, EPs_ERR
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               FATAL = .FALSE.
            ENDIF

            WRITE(ERR_MSG, 1210)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            DO MM = 1, PHASE_CNT
            M = PHASE_LIST(MM)
               WRITE(ERR_MSG,1211) M, BC_EP_s(BCV,M), EPs_TMP(M), &
                  (BC_EP_s(BCV,M)-EPs_TMP(M))
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDDO
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=FATAL)
         ENDIF


 1200 FORMAT('Error 1200: Unable to impose a uniform solids velocity ',&
         'on BC',I3,'.',/'Setting all solids to the highest velocity ',&
         'results in a large error',/'in the solids volume fraction. ',&
         'Please correct the mfix.dat file.',2/5X,'Max Inlet Velocity',&
         ':',1X,ES11.4,/5x,'Total BC_EP_s Error:',ES11.4)


 1205 FORMAT('Warning 1201: Uniform solids velocity imposed on BC',I3, &
         '.',2/,5X,'Uniform Inlet Velocity:',1X,ES11.4,/5X,'Total ',   &
         'BC_EP_s Error:',4X,ES11.4,/' ')

 1210 FORMAT(/,5X,'|',11('-'),'|',3(14('-'),'|'),/5X,'|',3X,'Phase',3X,&
         '|',4X,'BC_EP_s',3X,'|',2X,'Calculated',2X,'|',3X,'ABS ',   &
         'Error',2X,'|',/5X,'|',11('-'),'|',3(14('-'),'|'))

 1211 FORMAT(5X,'|',4X,I2,5X,'| ',1X,ES11.4,1X,'|',1X,ES11.4,2X,'|',   &
         1X,ES11.4,1X,' |',/5X,'|',11('-'),'|',3(14('-'),'|'))


! For polydisperse inlets, construct the DES_POLY_LAYOUT array
         IF(DEM_MI(BCV_I)%POLYDISPERSE) THEN
            RANGE_BOT = 1
            DO MM=1,PHASE_CNT - 1
               M = PHASE_LIST(MM)
               SCALED_VAL = dble(NUMFRAC_LIMIT)*(NPMpSEC(M)/NPpSEC)
               RANGE_TOP = FLOOR(SCALED_VAL) + (RANGE_BOT-1)
               DEM_BC_POLY_LAYOUT(BCV_I,RANGE_BOT:RANGE_TOP) = M
               RANGE_BOT = RANGE_TOP+1
            ENDDO

            M = PHASE_LIST(PHASE_CNT)
            DEM_BC_POLY_LAYOUT(BCV_I,RANGE_BOT:NUMFRAC_LIMIT) = M
! For monodisperse inlets, store the single mass phase used
         ELSE
            DEM_BC_POLY_LAYOUT(BCV_I,:) = PHASE_LIST(1)
         ENDIF


! Calculate des mass inlet time; time between injection.  If the run
! type is RESTART_1, DES_MI_TIME will be picked up from the restart file
! with an updated value.
         IF(RUN_TYPE /= 'RESTART_1') &
            DEM_MI_TIME(BCV_I) = TIME + dble(PI_FACTOR(BCV_I)) * DTSOLID


         WRITE(ERR_MSG,1000) BCV, NPpDT, int(NPpDT/DTSOLID), &
            PI_FACTOR(BCV_I), PI_COUNT(BCV_I), DEM_MI_TIME(BCV_I)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT(2/,2X,'For mass inlet BC: ', I3,/,&
         4X,'No. particles injected per solids time step = ', ES15.8,/,&
         4X,'No. particles injected per second = ', I10,/,&
         4X,'PI_FACTOR = ', I10,' PI_COUNT = ', I5,/,&
         4X,'start DES_MI_TIME = ', ES15.8,/'  ')

      ENDDO

!
      CALL SET_DEM_BCMI_IJK


      CALL FINL_ERR_MSG


      RETURN
      END SUBROUTINE SET_BC_DEM_MI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM                                              !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_DEM_BCMI_IJK

! Modules
!---------------------------------------------------------------------//
      use discretelement, only: min_radius, max_radius
      use bc, only: BC_PLANE
      use bc, only: BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t
      use des_bc, only: DEM_BCMI, DEM_BCMI_MAP, DEM_BCMI_IJK
      use des_bc, only: DEM_BCMI_IJKSTART, DEM_BCMI_IJKEND
      use desgrid, only: dg_IMIN1, dg_IMAX1
      use desgrid, only: dg_JMIN1, dg_JMAX1
      use desgrid, only: dg_KMIN1, dg_KMAX1
      use desgrid, only: DG_FUNIJK
      use desgrid, only: IofPOS, JofPOS, KofPOS
      use desgrid, only: dg_is_ON_myPE_plus1layers
      use error_manager
      use funits, only: DMP_LOG
      use mpi_utility
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER, ALLOCATABLE :: LOC_DEM_BCMI_IJK(:)
      INTEGER :: BCV, BCV_I
      INTEGER :: LC
      INTEGER :: MAX_CELLS
      INTEGER :: BND1, BND2
      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag
      INTEGER :: I,J,K,IJK
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t
!......................................................................!

      CALL INIT_ERR_MSG("SET_DEM_BCMI_IJK")

      dFlag = (DMP_LOG .AND. setDBG)

      if(dFlag) write(*,"(2/,2x,'From: SET_DEM_BCMI_IJK')")

! Loop over all inflow BCs to get an approximate count of the number
! of fluid cells that are adjacent to them.
      MAX_CELLS = 0
      DO BCV_I=1, DEM_BCMI
         BCV = DEM_BCMI_MAP(BCV_I)

! Set the search area a little bigger than the inlet area.
         if(dFlag) WRITE(*,"(/2x,'Adding cells for BCV: ',I3)") BCV
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S')
            BND1 = min(IofPOS(BC_X_e(BCV))+1,dg_IMAX1) - &
               max(IofPOS(BC_X_w(BCV))-1,dg_IMIN1)
            BND2 = min(KofPOS(BC_Z_t(BCV))+1,dg_KMAX1) - &
               max(KofPOS(BC_Z_b(BCV))-1,dg_KMIN1)

         CASE('E','W')
            BND1 = min(JofPOS(BC_Y_n(BCV))+1,dg_JMAX1) - &
               max(JofPOS(BC_Y_s(BCV))-1,dg_JMIN1)
            BND2 = min(KofPOS(BC_Z_t(BCV))+1,dg_KMAX1) - &
               max(KofPOS(BC_Z_b(BCV))-1,dg_KMIN1)

         CASE('T','B')
            BND1 = min(IofPOS(BC_X_e(BCV))+1,dg_IMAX1) - &
               max(IofPOS(BC_X_w(BCV))-1,dg_IMIN1)
            BND2 = min(JofPOS(BC_Y_n(BCV))+1,dg_JMAX1) - &
               max(JofPOS(BC_Y_s(BCV))-1,dg_JMIN1)
         END SELECT

         MAX_CELLS = MAX_CELLS + (BND1 + 1)*(BND2 + 1)
         if(dFlag) WRITE(*,"(4x,'Plane:   ',A)") BC_PLANE(BCV)
         if(dFlag) WRITE(*,"(4x,'Cells: ', I8)") (BND1+1)*(BND2+1)
      ENDDO

      if(dFlag) write(*,"(2x,'Max Cells: ',I8)") MAX_CELLS

! Allocate an array to hold the IJK values. This array should be
! more than enough to store all the IJKs.
      allocate( LOC_DEM_BCMI_IJK(MAX_CELLS) )


! Loop over the IJKs for each BC and store only the IJKs that you
! own as well as the start/end array positions for each BC.
      LC = 1
      DO BCV_I = 1, DEM_BCMI

         DEM_BCMI_IJKSTART(BCV_I) = LC
         BCV = DEM_BCMI_MAP(BCV_I)

         if(dFlag) write(*,"(/2x,'Searching for fluid cells:',I3)") BCV

         I_w = max(IofPOS(BC_X_w(BCV))-1,dg_IMIN1)
         I_e = min(IofPOS(BC_X_e(BCV))+1,dg_IMAX1)
         J_s = max(JofPOS(BC_Y_s(BCV))-1,dg_JMIN1)
         J_n = min(JofPOS(BC_Y_n(BCV))+1,dg_JMAX1)
         K_b = max(KofPOS(BC_Z_b(BCV))-1,dg_KMIN1)
         K_t = min(KofPOS(BC_Z_t(BCV))+1,dg_KMAX1)

! Depending on the flow plane, the 'common' index needs set to reference
! the fluid cell.
         SELECT CASE (BC_PLANE(BCV))
         CASE('N')
            J_s = JofPOS(BC_Y_s(BCV)+MIN_RADIUS)
            J_n = JofPOS(BC_Y_n(BCV)+MAX_RADIUS)
         CASE('S')
            J_s = JofPOS(BC_Y_s(BCV)-MAX_RADIUS)
            J_n = JofPOS(BC_Y_n(BCV)-MIN_RADIUS)
         CASE('E')
            I_w = IofPOS(BC_X_w(BCV)+MIN_RADIUS)
            I_e = IofPOS(BC_X_e(BCV)+MAX_RADIUS)
         CASE('W')
            I_w = IofPOS(BC_X_w(BCV)-MAX_RADIUS)
            I_e = IofPOS(BC_X_e(BCV)-MIN_RADIUS)
         CASE('T')
            K_b = KofPOS(BC_Z_b(BCV)+MIN_RADIUS)
            K_t = KofPOS(BC_Z_t(BCV)+MAX_RADIUS)
         CASE('B')
            K_b = KofPOS(BC_Z_b(BCV)-MAX_RADIUS)
            K_t = KofPOS(BC_Z_t(BCV)-MIN_RADIUS)
         END SELECT

         if(dFlag) then
            write(*,"(4x,'Search bounds: ')")
            write(*,"(6x,'I_w/I_e:',2(2x,I6))") I_w, I_e
            write(*,"(6x,'J_s/J_n:',2(2x,I6))") J_s, J_n
            write(*,"(6x,'K_b/K_t:',2(2x,I6))") K_b, K_t
         endif

! Store the IJKs.
         DO K = K_b, K_t
         DO J = J_s, J_n
         DO I = I_w, I_e
! Skip cells that this rank does not own or are considered dead.
            IF(.NOT.dg_is_ON_myPE_plus1layers(I,J,K))CYCLE

            IJK = DG_FUNIJK(I,J,K)
            LOC_DEM_BCMI_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO

         DEM_BCMI_IJKEND(BCV_I) = LC-1

         IF(dFLAG) write(*,1111) BCV, BCV_I,                           &
            DEM_BCMI_IJKSTART(BCV_I), DEM_BCMI_IJKEND(BCV_I)

      ENDDO

 1111 FORMAT(/2x,'DEM Mass Inflow:',/4x,'BC:',I4,3x,'MAP:',I4,         &
         /4x,'IJKSTART:',I6,/4x,'IJKEND:  ',I6)


! Allocate the global store arrary array. This changes across MPI ranks.
      IF(LC > 1) THEN
         allocate( DEM_BCMI_IJK(LC-1) )
         DEM_BCMI_IJK(1:LC-1) = LOC_DEM_BCMI_IJK(1:LC-1)
      ELSE
         allocate( DEM_BCMI_IJK(1) )
         DEM_BCMI_IJK(1) = LOC_DEM_BCMI_IJK(1)
      ENDIF

      deallocate(LOC_DEM_BCMI_IJK)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_DEM_BCMI_IJK
