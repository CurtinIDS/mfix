MODULE CUT_CELL_PREPROC
   USE exit, ONLY: mfix_exit
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CUT_CELL_PREPROCESSING                                 C
!  Purpose: Perform the cut-cell preprocessing stage:                  C
!           Identify cut cells, define face areas, and volumes         C
!           Set flags                                                  C
!           Compute Interpolations factors                             C
!           Compute Non-orthogonality Corrections terms                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CUT_CELL_PREPROCESSING

      USE cdist
      USE compar
      USE constant
      USE cutcell
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE quadric
      USE run
      USE sendrecv
      USE toleranc
      USE vtk

      IMPLICIT NONE

      INTEGER :: SAFE_MODE_COUNT
      DOUBLE PRECISION :: CPU_PP_START,CPU_PP_END

      IF(.NOT.CG_HEADER_WAS_PRINTED) CALL PRINT_CG_HEADER

      CALL CPU_TIME (CPU_PP_START)

      CALL OPEN_CUT_CELL_FILES

      CALL ALLOCATE_CUT_CELL_ARRAYS

      CALL DEFINE_QUADRICS

      CALL SET_3D_CUT_CELL_FLAGS

      CALL GATHER_DATA

!======================================================================
! Gather Data  and writes surface(s) defined by all cut cells
!======================================================================

      IF(WRITE_VTK_FILES.AND.(.NOT.BDIST_IO)) THEN
         CALL WRITE_CUT_SURFACE_VTK
      ENDIF

      CALL SET_3D_CUT_U_CELL_FLAGS
      CALL SET_3D_CUT_V_CELL_FLAGS
      IF(DO_K) CALL SET_3D_CUT_W_CELL_FLAGS

      CALL SET_3D_CUT_CELL_TREATMENT_FLAGS

      CALL GET_3D_ALPHA_U_CUT_CELL
      CALL GET_3D_ALPHA_V_CUT_CELL
      IF(DO_K) CALL GET_3D_ALPHA_W_CUT_CELL

      CALL SET_GHOST_CELL_FLAGS

      CALL SET_ODXYZ_U_CUT_CELL
      CALL SET_ODXYZ_V_CUT_CELL
      IF(DO_K) CALL SET_ODXYZ_W_CUT_CELL

      CALL GET_U_MASTER_CELLS
      CALL GET_V_MASTER_CELLS
      IF(DO_K) CALL GET_W_MASTER_CELLS

      CALL SEND_RECEIVE_CUT_CELL_VARIABLES

      CALL GET_DISTANCE_TO_WALL

      CALL PRINT_GRID_STATISTICS

      CALL CG_GET_BC_AREA

      CALL FLOW_TO_VEL(.FALSE.)

      CALL CG_FLOW_TO_VEL

      CALL CONVERT_CG_MI_TO_PS

      CALL CPU_TIME (CPU_PP_END)

      IF(myPE == PE_IO) THEN
         WRITE(*,20)'CARTESIAN GRID PRE-PROCESSING COMPLETED IN ',CPU_PP_END - CPU_PP_START, ' SECONDS.'
         WRITE(*,10)'============================================================================='
      ENDIF

      IF(myPE == PE_IO) THEN

         SAFE_MODE_COUNT = SUM(CG_SAFE_MODE)

         IF(SAFE_MODE_COUNT>0) THEN


            WRITE(*,10)'######################################################################'
            WRITE(*,10)'######################################################################'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                              \/                                  ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##  ===>   WARNING: RUNNING CARTESIAN GRID IN SAFE MODE !  <===     ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##  SAFE MODE ACTIVATED FOR :                                       ##'
            IF(CG_SAFE_MODE(1)==1) WRITE(*,10)'##                            - All scalar quantities               ##'
            IF(CG_SAFE_MODE(3)==1) WRITE(*,10)'##                            - X-Velocity (Gas and Solids)         ##'
            IF(CG_SAFE_MODE(4)==1) WRITE(*,10)'##                            - Y-Velocity (Gas and Solids)         ##'
            IF(CG_SAFE_MODE(5)==1) WRITE(*,10)'##                            - Z-Velocity (Gas and Solids)         ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##                              /\                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'######################################################################'
            WRITE(*,10)'######################################################################'

         ENDIF
      ENDIF

      RETURN

10    FORMAT(1X,A)
20    FORMAT(1X,A,F8.2,A)

      END SUBROUTINE CUT_CELL_PREPROCESSING

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_FLOW_TO_VEL                                         C
!  Purpose: Convert flow to velocity bc's                              C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CG_FLOW_TO_VEL

      USE bc
      USE compar
      USE constant
      USE cutcell
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE physprop
      USE quadric
      USE run
      USE scales
      USE sendrecv
      USE toleranc
      USE vtk

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!     loop/variable indices
      INTEGER :: M, BCV
!     Volumetric flow rate computed from mass flow rate
      DOUBLE PRECISION :: VOLFLOW
!     Solids phase volume fraction
      DOUBLE PRECISION :: EPS
!     Average molecular weight
      DOUBLE PRECISION :: MW
!
!-----------------------------------------------

      DO BCV = 1, DIMENSION_BC

         IF (BC_TYPE_ENUM(BCV)==CG_MI) THEN

            IF(BC_VELMAG_g(BCV)==UNDEFINED) THEN
!
!           If gas mass flow is defined convert it to volumetric flow
!
               IF (BC_MASSFLOW_G(BCV) /= UNDEFINED) THEN
                  IF (RO_G0 /= UNDEFINED) THEN
                     VOLFLOW = BC_MASSFLOW_G(BCV)/RO_G0
                  ELSE
                     IF (BC_P_G(BCV)/=UNDEFINED .AND. BC_T_G(BCV)/=UNDEFINED) &
                        THEN
                        IF (MW_AVG == UNDEFINED) THEN
                           MW = CALC_MW(BC_X_G,DIMENSION_BC,BCV,NMAX(0),MW_G)
                        ELSE
                           MW = MW_AVG
                        ENDIF
                        VOLFLOW = BC_MASSFLOW_G(BCV)/EOSG(MW,(BC_P_G(BCV)-P_REF), &
                                                 BC_T_G(BCV))
                     ELSE
                        IF (BC_TYPE_ENUM(BCV) == CG_MO) THEN
                           IF (BC_MASSFLOW_G(BCV) == ZERO) THEN
                              VOLFLOW = ZERO
                           ENDIF
                        ELSE
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1020) BCV
                           call mfix_exit(myPE)
                        ENDIF
                     ENDIF
                  ENDIF
!
!             If volumetric flow is also specified compare both
!
                  IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
                     IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_G(BCV))) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, VOLFLOW, BC_VOLFLOW_G(BCV)
                        call mfix_exit(myPE)
                     ENDIF
                  ELSE
                     BC_VOLFLOW_G(BCV) = VOLFLOW
                  ENDIF
               ENDIF
!
!           If gas volumetric flow is defined convert it to velocity
!
               IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
                  IF (BC_EP_G(BCV) /= UNDEFINED) THEN
                     BC_VELMAG_g(BCV) = BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV))
                  ELSE
                     RETURN                      !Error caught in Check_data_07
                  ENDIF
               ENDIF

            ENDIF


!
!  Do flow conversions for solids phases
!
            DO M = 1, MMAX

               IF(BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
!
!             If solids mass flow is defined convert it to volumetric flow
!
                  IF (BC_MASSFLOW_S(BCV,M) /= UNDEFINED) THEN
                     IF (RO_S0(M) /= UNDEFINED) THEN
                        VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S0(M)
                     ELSE
                        RETURN                   !  This error will be caught in a previous routine
                     ENDIF
!
!               If volumetric flow is also specified compare both
!
                     IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
                        IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_S(BCV,M))) THEN
                           IF(DMP_LOG)WRITE(UNIT_LOG,1200)BCV,VOLFLOW,M,BC_VOLFLOW_S(BCV,M)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        BC_VOLFLOW_S(BCV,M) = VOLFLOW
                     ENDIF
                  ENDIF

                  IF (BC_ROP_S(BCV,M)==UNDEFINED .AND. MMAX==1) BC_ROP_S(BCV,M)&
                        = (ONE - BC_EP_G(BCV))*RO_S0(M)
                  IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
                     IF (BC_ROP_S(BCV,M) /= UNDEFINED) THEN
                        EPS = BC_ROP_S(BCV,M)/RO_S0(M)
                        IF (EPS /= ZERO) THEN
                           BC_VELMAG_s(BCV,M) = BC_VOLFLOW_S(BCV,M)/(BC_AREA(BCV)*EPS)
                        ELSE
                           IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN
                              BC_VELMAG_s(BCV,M) = ZERO
                           ELSE
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1250) BCV, M
                              call mfix_exit(myPE)
                           ENDIF
                        ENDIF
                     ELSE
                        IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN
                           BC_VELMAG_s(BCV,M) = ZERO
                        ELSE
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1260) BCV, M
                           call mfix_exit(myPE)
                        ENDIF
                     ENDIF
                  ENDIF

               ENDIF
            END DO
         ENDIF
      END DO

100         FORMAT(1X,A,I8)
110         FORMAT(1X,A,A)
120         FORMAT(1X,A,F14.8,/)
130         FORMAT(1X,A,I8,F14.8,/)

 1000 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_g) = ',G14.7,/1X,70('*')/)


 1020 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,&
         '  BC_P_g, BC_T_g, and BC_X_g',/' should be specified',/1X,70('*')/)


 1200 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_s',I1,') = ',G14.7,/1X,70('*')/)

 1250 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Non-zero vol. or mass flow specified with BC_ROP_s',&
         I1,' = 0.',/1X,70('*')/)
 1260 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' BC_ROP_s',I1,' not specified',/1X,70('*')/)
      RETURN

      END SUBROUTINE CG_FLOW_TO_VEL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONVERT_CG_MI_TO_PS                                    C
!  Purpose: Convert CG_MI BCs to Point sources                         C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 06-Jan-14  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CONVERT_CG_MI_TO_PS

      USE bc
      USE compar
      USE constant
      USE cutcell
      USE eos, only: EOSG
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE physprop
      USE ps
      USE quadric
      USE run
      USE scales
      USE sendrecv
      USE toleranc
      USE vtk

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!     loop/variable indices
      INTEGER :: IJK, M, NN, BCV
!
      INTEGER :: iproc
      INTEGER :: NPS,PSV

!-----------------------------------------------
!

!      print*,'Entering test',MyPE
#ifdef MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

! Each procesor waits for its turn to find cells where to add a point source and updates the list of point sources

      do iproc = 0,NumPEs-1
         if (MyPE==iproc) Then

! First, find how many point sources are already defined. This could be regular PS from mfix.dat or new ones
! coming from the convertion of CG_MI to PS
               NPS = 0

               PS_LP: do PSV = 1, DIMENSION_PS
                  if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
                  NPS = PSV
               enddo PS_LP

!               print *,'Last PS=',NPS

! Next loop through all cells, and when a cut-cell with CG_MI is found, add a point source in this cell

            DO IJK = ijkstart3, ijkend3
               BCV = BC_ID(IJK)
               IF(BCV>0) THEN
                  IF(CG_MI_CONVERTED_TO_PS(BCV).AND.INTERIOR_CELL_AT(IJK).AND.VOL(IJK)>ZERO) THEN

                     NPS = NPS + 1

!                     print*,MyPE,NPS

                     PS_DEFINED(NPS) = .TRUE.

                     POINT_SOURCE = .TRUE.

                     PS_I_w(NPS) = I_OF(IJK)
                     PS_I_e(NPS) = I_OF(IJK)
                     PS_J_s(NPS) = J_OF(IJK)
                     PS_J_n(NPS) = J_OF(IJK)
                     PS_K_b(NPS) = K_OF(IJK)
                     PS_K_t(NPS) = K_OF(IJK)

                     PS_VOLUME(NPS) = VOL(IJK)

                     PS_MASSFLOW_g(NPS) = BC_MASSFLOW_g(BCV) * VOL(IJK) / BC_VOL(BCV)

                     PS_T_g(NPS)    = BC_T_g(BCV)

                     IF(BC_U_g(BCV)==UNDEFINED) THEN
                        PS_U_g(NPS)    = Normal_S(IJK,1)
                     ELSE
                        PS_U_g(NPS)    = BC_U_g(BCV)
                     ENDIF

                     IF(BC_V_g(BCV)==UNDEFINED) THEN
                        PS_V_g(NPS)    = Normal_S(IJK,2)
                     ELSE
                        PS_V_g(NPS)    = BC_V_g(BCV)
                     ENDIF

                     IF(BC_W_g(BCV)==UNDEFINED) THEN
                        PS_W_g(NPS)    = Normal_S(IJK,3)
                     ELSE
                        PS_W_g(NPS)    = BC_W_g(BCV)
                     ENDIF

                     DO NN=1,NMAX(0)
                        PS_X_g(NPS,NN)    = BC_X_g(BCV,NN)
                     ENDDO

                     DO M=1, MMAX
                        PS_MASSFLOW_s(NPS,M) = BC_MASSFLOW_s(BCV,M) * VOL(IJK) / BC_VOL(BCV)

                        PS_T_s(NPS,1)  = BC_T_s(BCV,M)

                        IF(BC_U_s(BCV,M)==UNDEFINED) THEN
                           PS_U_s(NPS,M)    = Normal_S(IJK,1)
                        ELSE
                           PS_U_s(NPS,M)    = BC_U_s(BCV,M)
                        ENDIF

                        IF(BC_V_s(BCV,M)==UNDEFINED) THEN
                           PS_V_s(NPS,M)    = Normal_S(IJK,2)
                        ELSE
                           PS_V_s(NPS,M)    = BC_V_s(BCV,M)
                        ENDIF

                        IF(BC_W_s(BCV,M)==UNDEFINED) THEN
                           PS_W_s(NPS,M)    = Normal_S(IJK,3)
                        ELSE
                           PS_W_s(NPS,M)    = BC_W_s(BCV,M)
                        ENDIF


                        DO NN=1,NMAX(M)
                           PS_X_s(NPS,M,NN)    = BC_X_s(BCV,M,NN)
                        ENDDO

                     ENDDO

!                     print*,'PS created:',NPS,PS_MASSFLOW_g(NPS),PS_VOLUME(NPS),BC_VOL(BCV)
                  ENDIF
               ENDIF

            ENDDO  ! IJK Loop

         endif  ! Work done by each processor in same order as rank

#ifdef MPI
         CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
         call bcast(POINT_SOURCE,iproc)
         call bcast(PS_DEFINED,iproc)
         call bcast(PS_I_w,iproc)
         call bcast(PS_I_e,iproc)
         call bcast(PS_J_s,iproc)
         call bcast(PS_J_n,iproc)
         call bcast(PS_K_b,iproc)
         call bcast(PS_K_t,iproc)
         call bcast(PS_MASSFLOW_g,iproc)
         call bcast(PS_U_g,iproc)
         call bcast(PS_V_g,iproc)
         call bcast(PS_W_g,iproc)
         call bcast(PS_X_g,iproc)
         call bcast(PS_T_g,iproc)
         call bcast(PS_MASSFLOW_s,iproc)
         call bcast(PS_U_s,iproc)
         call bcast(PS_V_s,iproc)
         call bcast(PS_W_s,iproc)
         call bcast(PS_X_s,iproc)
         call bcast(PS_T_s,iproc)
         call bcast(PS_VOLUME,iproc)

      enddo

#ifdef MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
!      print*,'Leaving test',MyPE
!      call mfix_exit(myPE)

      RETURN

      END SUBROUTINE CONVERT_CG_MI_TO_PS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONVERT_CG_MI_TO_PS_PE                                 C
!  Purpose: Convert CG_MI BCs to Point sources                         C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 06-Jan-14  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CONVERT_CG_MI_TO_PS_PE

      USE bc
      USE compar
      USE constant
      USE cutcell
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE physprop
      USE ps
      USE quadric
      USE run
      USE scales
      USE sendrecv
      USE toleranc
      USE vtk

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!     loop/variable indices
      INTEGER :: IJK, BCV
!
      INTEGER :: NPS,PSV
!
!-----------------------------------------------
!

! Find the last Point source that is defined. New point sources
! will be added after that.

!     print*,'setting bc_type to CG_NSW and exiting'
!      DO BCV = 1, DIMENSION_BC
!         IF (BC_TYPE_ENUM(BCV) == 'CG_MI') THEN
!            BC_TYPE_ENUM(BCV) = 'CG_NSW'
!            print*,'Converted CG_MI to CG_FSW for BC#',BCV
!         ENDIF
!      ENDDO
!      RETURN

      NPS = 0

      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         NPS = PSV
      enddo PS_LP

      print *,'Last PS=',NPS
!      read(*,*)

! Loop though each cell. When a CG_MI is found convert it to a single point source
! and change the BC_TYPE to Free-slip

      DO IJK = ijkstart3, ijkend3
         BCV = BC_ID(IJK)
         IF(BCV>0) THEN
            IF(CG_MI_CONVERTED_TO_PS(BCV).AND.INTERIOR_CELL_AT(IJK).AND.VOL(IJK)>ZERO) THEN

               NPS = NPS + 1

               PS_DEFINED(NPS) = .TRUE.

               POINT_SOURCE = .TRUE.

               PS_I_w(NPS) = I_OF(IJK)
               PS_I_e(NPS) = I_OF(IJK)
               PS_J_s(NPS) = J_OF(IJK)
               PS_J_n(NPS) = J_OF(IJK)
               PS_K_b(NPS) = K_OF(IJK)
               PS_K_t(NPS) = K_OF(IJK)

               PS_MASSFLOW_g(NPS) = BC_MASSFLOW_g(BCV) * VOL(IJK) / BC_VOL(BCV)

               PS_VOLUME(NPS) = VOL(IJK)

               PS_T_g(NPS)    = BC_T_g(BCV)

               IF(BC_U_g(NPS)==UNDEFINED) THEN
                  PS_U_g(NPS)    = Normal_S(IJK,1)
               ELSE
                  PS_U_g(NPS)    = BC_U_g(NPS)
               ENDIF

               IF(BC_V_g(NPS)==UNDEFINED) THEN
                  PS_V_g(NPS)    = Normal_S(IJK,2)
               ELSE
                  PS_V_g(NPS)    = BC_V_g(NPS)
               ENDIF

               IF(BC_W_g(NPS)==UNDEFINED) THEN
                  PS_W_g(NPS)    = Normal_S(IJK,3)
               ELSE
                  PS_W_g(NPS)    = BC_W_g(NPS)
               ENDIF

! This is a temporary setting for the solids phase and will need to be generalalized
               PS_MASSFLOW_s(NPS,1) = 0.0

               PS_T_s(NPS,1)  = 298.0

               PS_U_s(NPS,1)    = Normal_S(IJK,1)
               PS_V_s(NPS,1)    = Normal_S(IJK,2)
               PS_W_s(NPS,1)    = Normal_S(IJK,3)

               PS_U_s(NPS,1)    = ZERO
               PS_V_s(NPS,1)    = ZERO
               PS_W_s(NPS,1)    = ZERO

!               IF(Normal_S(IJK,2)/=ONE) print*,'Not vertical'
!               IF(Normal_S(IJK,2)==ONE) print*,'    vertical'

!               IF(Normal_S(IJK,2)/=ONE) PS_MASSFLOW_g(NPS) = ZERO

!               IF(.NOT.CUT_CELL_AT(IJK)) THEN
!                  print*,'turn off PS :not a scalar cut cell'
!                  PS_MASSFLOW_g(NPS) = ZERO
!               ENDIF
!               IF(.NOT.CUT_U_CELL_AT(IJK)) THEN
!                  print*,'turn off PS :not a u cut cell'
!                  PS_MASSFLOW_g(NPS) = ZERO
!               ENDIF
!               IF(.NOT.CUT_V_CELL_AT(IJK)) THEN
!                  print*,'turn off PS :not a v cut cell'
!                  PS_MASSFLOW_g(NPS) = ZERO
!               ENDIF
!               IF(.NOT.CUT_W_CELL_AT(IJK)) THEN
!                  print*,'turn off PS :not a w cut cell'
!                  PS_MASSFLOW_g(NPS) = ZERO
!               ENDIF


!                  PS_MASSFLOW_g(NPS) = ZERO

               print*,'PS created:',NPS,PS_MASSFLOW_g(NPS),PS_VOLUME(NPS),PS_I_w(NPS),PS_J_n(NPS),PS_K_b(NPS), &
                    INTERIOR_CELL_AT(IJK),PS_U_g(NPS),PS_V_g(NPS),PS_W_g(NPS), &
                    CUT_CELL_AT(IJK),Normal_S(IJK,1),Normal_S(IJK,2),Normal_S(IJK,3)
!               ENDIF

!               PS_DEFINED(NPS) = .FALSE.
            ENDIF
         ENDIF
      ENDDO

!      DO BCV = 1, DIMENSION_BC
!         IF (BC_TYPE_ENUM(BCV) == 'CG_MI') THEN
!            BC_TYPE_ENUM(BCV) = 'CG_NSW'
!            print*,'Converted CG_MI to CG_FSW for BC#',BCV
!         ENDIF
!      ENDDO

100         FORMAT(1X,A,I8)
110         FORMAT(1X,A,A)
120         FORMAT(1X,A,F14.8,/)
130         FORMAT(1X,A,I8,F14.8,/)


 1000 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_g) = ',G14.7,/1X,70('*')/)


 1020 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,&
         '  BC_P_g, BC_T_g, and BC_X_g',/' should be specified',/1X,70('*')/)


 1200 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_s',I1,') = ',G14.7,/1X,70('*')/)

 1250 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Non-zero vol. or mass flow specified with BC_ROP_s',&
         I1,' = 0.',/1X,70('*')/)
 1260 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' BC_ROP_s',I1,' not specified',/1X,70('*')/)
      RETURN

      END SUBROUTINE CONVERT_CG_MI_TO_PS_PE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PRINT_CG_HEADER                                        C
!  Purpose: Display Cartesian-Grid Header on screen                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE PRINT_CG_HEADER

      USE param
      USE param1
      USE parallel
      USE compar
      USE cutcell

      IMPLICIT NONE




      IF(myPE == PE_IO) THEN

         WRITE(*,10)'============================================================================='
         WRITE(*,10)'  ____          ___________    _   ____      ____    __________   __________  '
         WRITE(*,10)' |    \        /           |  (_)  \   \    /   /   |          | |          | '
         WRITE(*,10)' |     \      /      ______|  ___   \   \  /   /    |    ______| |    ______| '
         WRITE(*,10)' |      \    /      |______  |   |   \   \/   /     |   |        |   |        '
         WRITE(*,10)' |       \  /              | |   |    \      /  === |   |        |   |  ____  '
         WRITE(*,10)' |   |\   \/   /|    ______| |   |    /      \      |   |        |   | |_   | '
         WRITE(*,10)' |   | \      / |   |        |   |   /   /\   \     |   |______  |   |___|  | '
         WRITE(*,10)' |   |  \    /  |   |        |   |  /   /  \   \    |          | |          | '
         WRITE(*,10)' |___|   \__/   |___|        |___| /___/    \___\   |__________| |__________| '
         WRITE(*,10)'                                                                              '
         WRITE(*,10)'============================================================================='
         WRITE(*,10)'MFIX WITH CARTESIAN GRID IMPLEMENTATION.'


         IF(RE_INDEXING) THEN
            WRITE(*,10)'RE-INDEXING IS TURNED ON.'
!            IF(ADJUST_PROC_DOMAIN_SIZE) THEN
!               WRITE(*,10)'EACH PROCESSOR DOMAIN SIZE WILL BE ADJUSTED FOR BETTER LOAD BALANCING.'
!            ELSE
!               WRITE(*,10)'WARNING: PROCESSOR DOMAIN SIZE WILL BE UNIFORMLY DISTRIBUTED.'
!               WRITE(*,10)'THIS COULD RESULT IN VERY POOR LOAD BALANCING.'
!            ENDIF
         ELSE
            WRITE(*,10)'RE-INDEXING IS TURNED OFF.'
         ENDIF

      ENDIF

      CG_HEADER_WAS_PRINTED = .TRUE.

10    FORMAT(1X,A)

      RETURN

      END SUBROUTINE PRINT_CG_HEADER

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EVAL_F                                                 C
!  Purpose: Evaluate the function f(x,y,z) defining the boundary       C
!                                                                      C
!  Author: Jeff Dietiker                               Date: 21-FEB-08 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!  Modified: ##                                        Date: ##-###-## C
!  Purpose: ##                                                         C
!                                                                      C
!  Literature/Document References: ##                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
        SUBROUTINE EVAL_F(METHOD,x1,x2,x3,Q,f,CLIP_FLAG)

      USE compar
      USE cutcell
      USE parallel
      USE param1, only: undefined
      USE quadric
      USE quadric
      USE sendrecv

      IMPLICIT NONE

      DOUBLE PRECISION x1,x2,x3
      DOUBLE PRECISION f
      INTEGER :: Q,Q_ID,BCID
      LOGICAL :: CLIP_FLAG
      CHARACTER (LEN = 7) :: METHOD
      CHARACTER(LEN=9) :: GR

      INTEGER :: GROUP,GS,P

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F_G

      ALLOCATE(F_G(N_GROUP,0:DIM_QUADRIC))

      SELECT CASE(METHOD)

         CASE('QUADRIC')

            F_G = - UNDEFINED

            DO GROUP = 1, N_GROUP
               GS = GROUP_SIZE(GROUP)
               GR = TRIM(GROUP_RELATION(GROUP))

               DO P = 1 , GS
                  Q_ID = GROUP_Q(GROUP,P)
                  CALL GET_F_QUADRIC(x1,x2,x3,Q_ID,F_G(GROUP,P),CLIP_FLAG)
               ENDDO
               IF(GR == 'AND') THEN
                  F_G(GROUP,0) = MAXVAL(F_G(GROUP,1:GS))
               ELSEIF(GR == 'OR') THEN
                  F_G(GROUP,0) = MINVAL(F_G(GROUP,1:GS))
               ELSEIF(GR == 'PIECEWISE') THEN
                  CALL REASSSIGN_QUADRIC(x1,x2,x3,GROUP,Q_ID)
!                  CLIP_FLAG=.FALSE.
                  CALL GET_F_QUADRIC(x1,x2,x3,Q_ID,F_G(GROUP,0),CLIP_FLAG)
!                  CLIP_FLAG=.TRUE.
               ENDIF

            ENDDO

            f = F_G(1,0)

            DO GROUP = 2, N_GROUP

               GR = TRIM(RELATION_WITH_PREVIOUS(GROUP))

               IF(GR =='AND') THEN
                  f = DMAX1(f,F_G(GROUP,0))
               ELSEIF(GR =='OR') THEN
                  f = DMIN1(f,F_G(GROUP,0))
               ENDIF

            ENDDO


         CASE('POLYGON')

            CALL EVAL_POLY_FCT(x1,x2,x3,Q,f,CLIP_FLAG,BCID)

         CASE('USR_DEF')

            CALL EVAL_USR_FCT(x1,x2,x3,Q,f,CLIP_FLAG)

!         CASE('STL')

!            CALL EVAL_STL_FCT(x1,x2,x3,Q,f,CLIP_FLAG,BCID)

         CASE DEFAULT

            WRITE(*,*)'ERROR IN SUBROUTINE EVAL_F.'
            WRITE(*,*)'UNKNOWN METHOD:',METHOD
            WRITE(*,*)'ACCEPTABLE METHODS:'
            WRITE(*,*)'QUADRIC'
            WRITE(*,*)'POLYGON'
            WRITE(*,*)'USR_DEF'
!            WRITE(*,*)'STL'
            CALL MFIX_EXIT(myPE)

      END SELECT

      DEALLOCATE(F_G)

      RETURN
      END SUBROUTINE EVAL_F

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: intersect_line                                         C
!  Purpose: Finds the intersection between the quadric surface ,       C
!           and the line (xa,ya,za) and (xb,yb,zb).                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE INTERSECT_LINE(METHOD,xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric

      IMPLICIT NONE
      DOUBLE PRECISION:: x1,y1,z1,x2,y2,z2,x3,y3,z3
      DOUBLE PRECISION:: xa,ya,za,xb,yb,zb,xc,yc,zc
      INTEGER :: Q_ID,niter
      DOUBLE PRECISION :: x_intersection
      DOUBLE PRECISION :: f1,f2,f3,fa,fb
      DOUBLE PRECISION :: t1,t2,t3
      LOGICAL :: CLIP_FLAG,CLIP_FLAG1,CLIP_FLAG2,CLIP_FLAG3,INTERSECT_FLAG
      CHARACTER (LEN=7) ::METHOD

      x1 = xa    ! Initial guesses
      y1 = ya
      z1 = za
      t1 = ZERO

      x2 = xb
      y2 = yb
      z2 = zb
      t2 = ONE

      CALL EVAL_F(METHOD,x1,y1,z1,Q_ID,f1,CLIP_FLAG1)
      CALL EVAL_F(METHOD,x2,y2,z2,Q_ID,f2,CLIP_FLAG2)

!======================================================================
!  The line from (x1,y1,z1) and (x2,y2,z2) is parametrized
!  from t1 = ZERO to t2 = ONE
!======================================================================

      niter = 1

      CLIP_FLAG = (CLIP_FLAG1).AND.(CLIP_FLAG2)

      if(DABS(f1)<TOL_F) then  ! ignore intersection at corner
        xc = UNDEFINED
        yc = UNDEFINED
        zc = UNDEFINED
        INTERSECT_FLAG = .FALSE.
      elseif(DABS(f2)<TOL_F) then ! ignore intersection at corner
        xc = UNDEFINED
        yc = UNDEFINED
        zc = UNDEFINED
        INTERSECT_FLAG = .FALSE.
      elseif(f1*f2 < ZERO) then
        niter = 0
        f3 = 2.0d0*TOL_F
        do while (   (abs(f3) > TOL_F)   .AND.   (niter<ITERMAX_INT)       )

          t3 = t1 - f1*(t2-t1)/(f2-f1)

          x3 = x1 + t3 * (x2 - x1)
          y3 = y1 + t3 * (y2 - y1)
          z3 = z1 + t3 * (z2 - z1)

          CALL EVAL_F(METHOD,x3,y3,z3,Q_ID,f3,CLIP_FLAG3)
          if(f1*f3<0) then
            t2 = t3
            f2 = f3
          else
            t1 = t3
            f1 = f3
          endif
          niter = niter + 1

        end do
        if (niter < ITERMAX_INT) then
           xc = x3
           yc = y3
           zc = z3
          INTERSECT_FLAG = .TRUE.
        else
           WRITE(*,*)'   Subroutine intersect_line:'
           WRITE(*,*)   'Unable to find the intersection of quadric:',Q_ID
           WRITE(*,1000)'between (x1,y1,z1)= ', xa,ya,za
           WRITE(*,1000)'   and  (x2,y2,z2)= ', xb,yb,zb
           CALL EVAL_F(METHOD,xa,ya,za,Q_ID,fa,CLIP_FLAG1)
           CALL EVAL_F(METHOD,xb,yb,zb,Q_ID,fb,CLIP_FLAG1)
           WRITE(*,1000)'f(x1,y1,z1) = ', fa
           WRITE(*,1000)'f(x2,y2,z2) = ', fb
           WRITE(*,1000)'Current Location (x3,y3,z3)= ', x3,y3,z3
           WRITE(*,1000)'Current value of abs(f) = ', DABS(f3)
           WRITE(*,1000)'Tolerance = ', TOL_F
           WRITE(*,*)'Maximum number of iterations = ', ITERMAX_INT
           WRITE(*,*)   'Please increase the intersection tolerance, '
           WRITE(*,*)   'or the maximum number of iterations, and try again.'
           WRITE(*,*)   'MFiX will exit now.'
           CALL MFIX_EXIT(myPE)
           x_intersection = UNDEFINED
           INTERSECT_FLAG = .FALSE.

        endif
      else
        xc = UNDEFINED
        yc = UNDEFINED
        zc = UNDEFINED
        INTERSECT_FLAG = .FALSE.
      endif

 1000 FORMAT(A,3(2X,G12.5))


      RETURN

      END SUBROUTINE INTERSECT_LINE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INTERSECT                                              C
!  Purpose: Intersects quadric with grid                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE INTERSECT(IJK,TYPE_OF_CELL,Xi,Yi,Zi)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE functions

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K,Q_ID,N_int_x,N_int_y,N_int_z,N_USR
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb,xc,yc,zc
      DOUBLE PRECISION :: Xi,Yi,Zi,Xc_backup,Yc_backup,Zc_backup
      LOGICAL :: INTERSECT_FLAG

      Xi = UNDEFINED
      Yi = UNDEFINED
      Zi = UNDEFINED

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

!======================================================================
!  Intersection with Edge 7 (node 7-8, Face North-Top):
!======================================================================
      xa = X_NODE(7)
      ya = Y_NODE(7)
      za = Z_NODE(7)

      xb = X_NODE(8)
      yb = Y_NODE(8)
      zb = Z_NODE(8)


      N_int_x = 0
      INTERSECT_X(IJK) = .FALSE.
      Q_ID = 1
      CALL INTERSECT_LINE('QUADRIC',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)
      IF ((INTERSECT_FLAG).AND.(xc/=Xi)) THEN
         N_int_x = N_int_x + 1
         INTERSECT_X(IJK) = .TRUE.
         xc_backup = Xi
         Xi = xc
      ENDIF

      IF(N_int_x /= 1) THEN
         Xi = UNDEFINED
         INTERSECT_X(IJK) = .FALSE.
      ENDIF


      DO Q_ID = 1, N_POLYGON
         CALL INTERSECT_LINE('POLYGON',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_X(IJK),xc,yc,zc)
         IF(INTERSECT_X(IJK)) Xi = xc
      ENDDO

      DO N_USR= 1, N_USR_DEF
         CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_X(IJK),xc,yc,zc)
         IF(INTERSECT_X(IJK)) Xi = xc
      ENDDO

!      IF(USE_STL) THEN
!         CALL INTERSECT_LINE_WITH_STL(xa,ya,za,xb,yb,zb,INTERSECT_X(IJK),xc,yc,zc)
!         IF(INTERSECT_X(IJK)) Xi = xc
!      ENDIF

      IF(TYPE_OF_CELL=='U_MOMENTUM') THEN
         IF(SNAP(IJK)) THEN
            INTERSECT_X(IJK) = .TRUE.
            I = I_OF(IJK)
            Xi = XG_E(I)
         ENDIF
      ENDIF



!======================================================================
!  Intersection with Edge 6 (node 6-8, Face East-Top):
!======================================================================
      xa = X_NODE(6)
      ya = Y_NODE(6)
      za = Z_NODE(6)

      N_int_y = 0
      INTERSECT_Y(IJK) = .FALSE.
      Q_ID = 1
         CALL INTERSECT_LINE('QUADRIC',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)
         IF ((INTERSECT_FLAG).AND.(yc/=Yi)) THEN
            N_int_y = N_int_y + 1
            INTERSECT_Y(IJK) = .TRUE.
            yc_backup = Yi
            Yi = yc
         ENDIF

      IF(N_int_y /= 1) THEN
         Yi = UNDEFINED
         INTERSECT_Y(IJK) = .FALSE.
      ENDIF

      DO Q_ID = 1, N_POLYGON
         CALL INTERSECT_LINE('POLYGON',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Y(IJK),xc,yc,zc)
         IF(INTERSECT_Y(IJK)) Yi = yc
      ENDDO

      DO N_USR= 1, N_USR_DEF
         CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Y(IJK),xc,yc,zc)
         IF(INTERSECT_Y(IJK)) Yi = yc
      ENDDO

!      IF(USE_STL) THEN
!         CALL INTERSECT_LINE_WITH_STL(xa,ya,za,xb,yb,zb,INTERSECT_Y(IJK),xc,yc,zc)
!         IF(INTERSECT_Y(IJK)) Yi = yc
!      ENDIF

      IF(TYPE_OF_CELL=='V_MOMENTUM') THEN
         IF(SNAP(IJK)) THEN
            INTERSECT_Y(IJK) = .TRUE.
            J = J_OF(IJK)
            Yi = YG_N(J)
         ENDIF
      ENDIF

      IF(DO_K) THEN
!======================================================================
!  Intersection with Edge 11 (node 4-8, Face East-North):
!======================================================================
         xa = X_NODE(4)
         ya = Y_NODE(4)
         za = Z_NODE(4)

         N_int_z = 0
         INTERSECT_Z(IJK) = .FALSE.
         Q_ID = 1
            CALL INTERSECT_LINE('QUADRIC',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)
            IF ((INTERSECT_FLAG).AND.(zc/=Zi)) THEN
               N_int_z = N_int_z + 1
               INTERSECT_Z(IJK) = .TRUE.
               zc_backup = Zi
               Zi = zc
            ENDIF

            IF(N_int_z /= 1) THEN
               Zi = UNDEFINED
               INTERSECT_Z(IJK) = .FALSE.
            ENDIF

         DO Q_ID = 1, N_POLYGON
            CALL INTERSECT_LINE('POLYGON',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Z(IJK),xc,yc,zc)
            IF(INTERSECT_Z(IJK)) Zi = zc
         ENDDO

         DO N_USR= 1, N_USR_DEF
            CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Z(IJK),xc,yc,zc)
            IF(INTERSECT_Z(IJK)) Zi = zc
         ENDDO

   !      IF(USE_STL) THEN
   !         CALL INTERSECT_LINE_WITH_STL(xa,ya,za,xb,yb,zb,INTERSECT_Z(IJK),xc,yc,zc)
   !         IF(INTERSECT_Z(IJK)) Zi = zc
   !      ENDIF

         IF(TYPE_OF_CELL=='W_MOMENTUM') THEN
            IF(SNAP(IJK)) THEN
               INTERSECT_Z(IJK) = .TRUE.
               K = K_OF(IJK)
               Zi = ZG_T(K)
            ENDIF
         ENDIF

      ENDIF


      IF(INTERSECT_X(IJK)) THEN
         POTENTIAL_CUT_CELL_AT(IJK) = .TRUE.
         POTENTIAL_CUT_CELL_AT(EAST_OF(IJK)) = .TRUE.
      ENDIF

      IF(INTERSECT_Y(IJK)) THEN
         POTENTIAL_CUT_CELL_AT(IJK) = .TRUE.
         POTENTIAL_CUT_CELL_AT(NORTH_OF(IJK)) = .TRUE.
      ENDIF


      IF(INTERSECT_Z(IJK)) THEN
         POTENTIAL_CUT_CELL_AT(IJK) = .TRUE.
         POTENTIAL_CUT_CELL_AT(TOP_OF(IJK)) = .TRUE.
      ENDIF


      RETURN

      END SUBROUTINE INTERSECT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLEAN_INTERSECT                                        C
!  Purpose: Remove Intersection flags in preparation of small cell     C
!           removal                                                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLEAN_INTERSECT(IJK,TYPE_OF_CELL,Xi,Yi,Zi)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE STL
      USE functions

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K,IM,JM,KM,IP,JP,KP
      INTEGER :: BCID
      INTEGER :: IMJK,IPJK,IJMK,IJPK,IJKM,IJKP,IMJPK,IMJKP,IPJMK,IJMKP,IPJKM,IJPKM
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb
      DOUBLE PRECISION :: Xi,Yi,Zi
      DOUBLE PRECISION :: DFC,DFC_MAX,F4,F6,F7,F8
      LOGICAL :: CLIP_FLAG,CAD,F_TEST


! When inputing geometry from CAD (STL or MSH file), the snapping procedure is
! dependent on the value of F at the cell corners
! For other gemoetry inputs (say quadrics), This is not needed, and the value
! of F_TEST is set to .TRUE. here
      CAD = USE_MSH.OR.USE_STL
      F_TEST = .TRUE.

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      IM = I - 1
      JM = J - 1
      KM = K - 1

      IP = I + 1
      JP = J + 1
      KP = K + 1

      IMJK = FUNIJK(IM,J,K)
      IPJK = FUNIJK(IP,J,K)
      IJMK = FUNIJK(I,JM,K)
      IJPK = FUNIJK(I,JP,K)
      IJKM = FUNIJK(I,J,KM)
      IJKP = FUNIJK(I,J,KP)

      IMJPK = FUNIJK(IM,JP,K)
      IMJKP = FUNIJK(IM,J,KP)

      IPJMK = FUNIJK(IP,JM,K)
      IJMKP = FUNIJK(I,JM,KP)

      IPJKM = FUNIJK(IP,J,KM)
      IJPKM = FUNIJK(I,JP,KM)


      IF(IMJK<1.OR.IMJK>DIMENSION_3) IMJK = IJK
      IF(IPJK<1.OR.IPJK>DIMENSION_3) IPJK = IJK
      IF(IJMK<1.OR.IJMK>DIMENSION_3) IJMK = IJK
      IF(IJPK<1.OR.IJPK>DIMENSION_3) IJPK = IJK
      IF(IJKM<1.OR.IJKM>DIMENSION_3) IJKM = IJK
      IF(IJKP<1.OR.IJKP>DIMENSION_3) IJKP = IJK

      IF(IMJPK<1.OR.IMJPK>DIMENSION_3) IMJPK = IJK
      IF(IMJKP<1.OR.IMJKP>DIMENSION_3) IMJKP = IJK

      IF(IPJMK<1.OR.IPJMK>DIMENSION_3) IPJMK = IJK
      IF(IJMKP<1.OR.IJMKP>DIMENSION_3) IJMKP = IJK

      IF(IPJKM<1.OR.IPJKM>DIMENSION_3) IPJKM = IJK
      IF(IJPKM<1.OR.IJPKM>DIMENSION_3) IJPKM = IJK

!======================================================================
!  Clean Intersection with Edge 7 (node 7-8, Face North-Top):
!======================================================================

      xa = X_NODE(7)
      ya = Y_NODE(7)
      za = Z_NODE(7)

      xb = X_NODE(8)
      yb = Y_NODE(8)
      zb = Z_NODE(8)

      DFC_MAX = TOL_SNAP(1) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER

      IF(INTERSECT_X(IJK)) THEN

         DFC = DABS(Xi-xa) ! DISTANCE FROM CORNER (NODE 7)

         IF(CAD) F_TEST = (F_AT(IMJK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING X-INTERSECTION ALONG EDGE 7 ONTO NODE 7'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            INTERSECT_X(IJK)  = .FALSE.
            IF(I>=IMIN1) THEN
               INTERSECT_X(IMJK) = .FALSE.
               INTERSECT_Y(IMJK)  = .FALSE.
               INTERSECT_Y(IMJPK) = .FALSE.
               IF(DO_K) INTERSECT_Z(IMJK)  = .FALSE.
               IF(DO_K) INTERSECT_Z(IMJKP) = .FALSE.

               SNAP(IMJK) = .TRUE.
            ENDIF

         ENDIF


         DFC = DABS(Xi-xb) ! DISTANCE FROM CORNER (NODE 8)

         IF(CAD) F_TEST = (F_AT(IJK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING X-INTERSECTION ALONG EDGE 7 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF


            INTERSECT_X(IJK)  = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            IF(DO_K) INTERSECT_Z(IJK)  = .FALSE.
            SNAP(IJK) = .TRUE.

            IF(I<=IMAX1) INTERSECT_X(IPJK) = .FALSE.
            IF(J<=JMAX1) INTERSECT_Y(IJPK) = .FALSE.
            IF(DO_K.AND.(K<=KMAX1)) INTERSECT_Z(IJKP) = .FALSE.


         ENDIF

         IF(USE_STL.OR.USE_MSH) THEN
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,7,F7,CLIP_FLAG,BCID)
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
            IF(F7*F8>TOL_STL**2) INTERSECT_X(IJK)  = .FALSE.
         ENDIF

      ENDIF




!======================================================================
!  Clean Intersection with Edge 6 (node 6-8, Face East-Top):
!======================================================================

      xa = X_NODE(6)
      ya = Y_NODE(6)
      za = Z_NODE(6)

      DFC_MAX = TOL_SNAP(2) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER


      IF(INTERSECT_Y(IJK)) THEN

         DFC = DABS(Yi-ya) ! DISTANCE FROM CORNER (NODE 6)

         IF(CAD) F_TEST = (F_AT(IJMK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Y-INTERSECTION ALONG EDGE 6 ONTO NODE 6'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            INTERSECT_Y(IJK)  = .FALSE.
            IF(J>=JMIN1) THEN
               INTERSECT_X(IJMK)  = .FALSE.
               INTERSECT_X(IPJMK) = .FALSE.
               INTERSECT_Y(IJMK) = .FALSE.
               IF(DO_K) INTERSECT_Z(IJMK)  = .FALSE.
               IF(DO_K) INTERSECT_Z(IJMKP) = .FALSE.

               SNAP(IJMK) = .TRUE.
            ENDIF

         ENDIF


         DFC = DABS(Yi-yb) ! DISTANCE FROM CORNER (NODE 8)

         IF(CAD) F_TEST = (F_AT(IJK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Y-INTERSECTION ALONG EDGE 6 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            INTERSECT_X(IJK)  = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            IF(DO_K) INTERSECT_Z(IJK)  = .FALSE.
            SNAP(IJK) = .TRUE.

            IF(I<=IMAX1) INTERSECT_X(IPJK) = .FALSE.
            IF(J<=JMAX1) INTERSECT_Y(IJPK) = .FALSE.
            IF(DO_K.AND.(K<=KMAX1)) INTERSECT_Z(IJKP) = .FALSE.


         ENDIF


         IF(USE_STL.OR.USE_MSH) THEN
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,6,F6,CLIP_FLAG,BCID)
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
            IF(F6*F8>TOL_STL**2) INTERSECT_Y(IJK)  = .FALSE.
         ENDIF

      ENDIF


      IF(DO_K) THEN
!======================================================================
!  Intersection with Edge 11 (node 4-8, Face East-North):
!======================================================================

         xa = X_NODE(4)
         ya = Y_NODE(4)
         za = Z_NODE(4)

         DFC_MAX = TOL_SNAP(3) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER

         IF(INTERSECT_Z(IJK)) THEN

            DFC = DABS(Zi-Za) ! DISTANCE FROM CORNER (NODE 4)

            IF(CAD) F_TEST = (F_AT(IJKM)/=ZERO)
            IF(DFC < DFC_MAX.AND.F_TEST) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 4'
                  WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
               ENDIF

               INTERSECT_Z(IJK)  = .FALSE.

               IF(K>=KMIN1) THEN
                  INTERSECT_X(IJKM)  = .FALSE.
                  INTERSECT_X(IPJKM) = .FALSE.
                  INTERSECT_Y(IJKM)  = .FALSE.
                  INTERSECT_Y(IJPKM) = .FALSE.
                  INTERSECT_Z(IJKM) = .FALSE.

                  SNAP(IJKM) = .TRUE.
               ENDIF

            ENDIF


            DFC = DABS(Zi-Zb) ! DISTANCE FROM CORNER (NODE 8)

            IF(CAD) F_TEST = (F_AT(IJK)/=ZERO)
            IF(DFC < DFC_MAX.AND.F_TEST) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 8'
                  WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
               ENDIF

               INTERSECT_X(IJK)  = .FALSE.
               INTERSECT_Y(IJK)  = .FALSE.
               INTERSECT_Z(IJK)  = .FALSE.
               SNAP(IJK) = .TRUE.
!            F_AT(IJKM) = UNDEFINED
!            IF(F_AT(IJK)/=ZERO) SNAP(IJK) = .TRUE.
!               F_AT(IJK) = ZERO

               IF(I<=IMAX1) INTERSECT_X(IPJK) = .FALSE.
               IF(J<=JMAX1) INTERSECT_Y(IJPK) = .FALSE.
               IF(K<=KMAX1) INTERSECT_Z(IJKP) = .FALSE.


            ENDIF

            IF(USE_STL.OR.USE_MSH) THEN
               CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,4,F4,CLIP_FLAG,BCID)
               CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
               IF(F4*F8>TOL_STL**2) INTERSECT_Z(IJK)  = .FALSE.
            ENDIF

         ENDIF

      ENDIF


      RETURN

      END SUBROUTINE CLEAN_INTERSECT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLEAN_INTERSECT_SCALAR                                 C
!  Purpose: Remove Intersection flags in preparation of small cell     C
!           removal                                                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-Dec-14  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLEAN_INTERSECT_SCALAR

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE STL
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK



      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_X,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Y,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Z,2)
      call send_recv(Xn_int,2)
      call send_recv(Ye_int,2)
      call send_recv(Zt_int,2)

      DO IJK = IJKSTART3, IJKEND3
         CALL SET_SNAP_FLAG(IJK,'SCALAR',Xn_int(IJK),Ye_int(IJK),Zt_int(IJK))
      END DO

      call SEND_RECEIVE_1D_LOGICAL(SNAP,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_X,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Y,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Z,2)
      call send_recv(Xn_int,2)
      call send_recv(Ye_int,2)
      call send_recv(Zt_int,2)

      DO IJK = IJKSTART3, IJKEND3
         CALL REMOVE_INTERSECT_FLAG(IJK)
      END DO

      call SEND_RECEIVE_1D_LOGICAL(SNAP,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_X,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Y,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Z,2)
      call send_recv(Xn_int,2)
      call send_recv(Ye_int,2)
      call send_recv(Zt_int,2)

      RETURN

      END SUBROUTINE CLEAN_INTERSECT_SCALAR

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_SNAP_FLAG                                          C
!  Purpose: Set SNAP flag in preparation of intersection flag removal  C
!           removal                                                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-Dec-14  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_SNAP_FLAG(IJK,TYPE_OF_CELL,Xi,Yi,Zi)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE STL
      USE functions

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K,IM,JM,KM
      INTEGER :: BCID
      INTEGER :: IMJK,IJMK,IJKM
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb
      DOUBLE PRECISION :: Xi,Yi,Zi
      DOUBLE PRECISION :: DFC,DFC_MAX,F4,F6,F7,F8
      LOGICAL :: CLIP_FLAG,CAD,F_TEST


! When inputing geometry from CAD (STL or MSH file), the snapping procedure is
! dependent on the value of F at the cell corners
! For other gemoetry inputs (say quadrics), This is not needed, and the value
! of F_TEST is set to .TRUE. here
      CAD = USE_MSH.OR.USE_STL
      F_TEST = .TRUE.

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      IM = I - 1
      JM = J - 1
      KM = K - 1

      IMJK = FUNIJK(IM,J,K)
      IJMK = FUNIJK(I,JM,K)
      IJKM = FUNIJK(I,J,KM)

      IF(IMJK<1.OR.IMJK>DIMENSION_3) IMJK = IJK
      IF(IJMK<1.OR.IJMK>DIMENSION_3) IJMK = IJK
      IF(IJKM<1.OR.IJKM>DIMENSION_3) IJKM = IJK

!======================================================================
!  Clean Intersection with Edge 7 (node 7-8, Face North-Top):
!======================================================================

      xa = X_NODE(7)
      ya = Y_NODE(7)
      za = Z_NODE(7)

      xb = X_NODE(8)
      yb = Y_NODE(8)
      zb = Z_NODE(8)


      DFC_MAX = TOL_SNAP(1) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER

      IF(INTERSECT_X(IJK)) THEN

         DFC = DABS(Xi-xa) ! DISTANCE FROM CORNER (NODE 7)

         IF(CAD) F_TEST = (F_AT(IMJK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING X-INTERSECTION ALONG EDGE 7 ONTO NODE 7'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            IF(I>=IMIN1) SNAP(IMJK) = .TRUE.

         ENDIF


         DFC = DABS(Xi-xb) ! DISTANCE FROM CORNER (NODE 8)

         IF(CAD) F_TEST = (F_AT(IJK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING X-INTERSECTION ALONG EDGE 7 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            SNAP(IJK) = .TRUE.

         ENDIF

         IF(USE_STL.OR.USE_MSH) THEN
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,7,F7,CLIP_FLAG,BCID)
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
            IF(F7*F8>TOL_STL**2) INTERSECT_X(IJK)  = .FALSE.
         ENDIF

      ENDIF


!======================================================================
!  Clean Intersection with Edge 6 (node 6-8, Face East-Top):
!======================================================================

      xa = X_NODE(6)
      ya = Y_NODE(6)
      za = Z_NODE(6)

      DFC_MAX = TOL_SNAP(2) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER


      IF(INTERSECT_Y(IJK)) THEN

         DFC = DABS(Yi-ya) ! DISTANCE FROM CORNER (NODE 6)

         IF(CAD) F_TEST = (F_AT(IJMK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Y-INTERSECTION ALONG EDGE 6 ONTO NODE 6'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            IF(J>=JMIN1) SNAP(IJMK) = .TRUE.

         ENDIF


         DFC = DABS(Yi-yb) ! DISTANCE FROM CORNER (NODE 8)

         IF(CAD) F_TEST = (F_AT(IJK)/=ZERO)
         IF(DFC < DFC_MAX.AND.F_TEST) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Y-INTERSECTION ALONG EDGE 6 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            SNAP(IJK) = .TRUE.


         ENDIF


         IF(USE_STL.OR.USE_MSH) THEN
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,6,F6,CLIP_FLAG,BCID)
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
            IF(F6*F8>TOL_STL**2) INTERSECT_Y(IJK)  = .FALSE.
         ENDIF

      ENDIF


      IF(DO_K) THEN
!======================================================================
!  Intersection with Edge 11 (node 4-8, Face East-North):
!======================================================================

         xa = X_NODE(4)
         ya = Y_NODE(4)
         za = Z_NODE(4)

         DFC_MAX = TOL_SNAP(3) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER

         IF(INTERSECT_Z(IJK)) THEN

            DFC = DABS(Zi-Za) ! DISTANCE FROM CORNER (NODE 4)

            IF(CAD) F_TEST = (F_AT(IJKM)/=ZERO)
            IF(DFC < DFC_MAX.AND.F_TEST) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 4'
                  WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
               ENDIF

               IF(K>=KMIN1) SNAP(IJKM) = .TRUE.

            ENDIF

            DFC = DABS(Zi-Zb) ! DISTANCE FROM CORNER (NODE 8)

            IF(CAD) F_TEST = (F_AT(IJK)/=ZERO)
            IF(DFC < DFC_MAX.AND.F_TEST) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 8'
                  WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
               ENDIF

               SNAP(IJK) = .TRUE.

            ENDIF


            IF(USE_STL.OR.USE_MSH) THEN
               CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,4,F4,CLIP_FLAG,BCID)
               CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
               IF(F4*F8>TOL_STL**2) INTERSECT_Z(IJK)  = .FALSE.
            ENDIF

         ENDIF

      ENDIF


      RETURN

      END SUBROUTINE SET_SNAP_FLAG

      !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REMOVE_INTERSECT_FLAG                                  C
!  Purpose: Remove Intersection flags                                  C
!           removal                                                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-Dec-14  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE REMOVE_INTERSECT_FLAG(IJK)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE STL
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,IP,JP,KP
      INTEGER :: IPJK,IJPK,IJKP


      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      IP = I + 1
      JP = J + 1
      KP = K + 1

      IPJK = FUNIJK(IP,J,K)
      IJPK = FUNIJK(I,JP,K)
      IJKP = FUNIJK(I,J,KP)

      IF(IPJK<1.OR.IPJK>DIMENSION_3) IPJK = IJK
      IF(IJPK<1.OR.IJPK>DIMENSION_3) IJPK = IJK
      IF(IJKP<1.OR.IJKP>DIMENSION_3) IJKP = IJK

       IF(SNAP(IJK)) THEN
          INTERSECT_X(IJK) = .FALSE.
          INTERSECT_Y(IJK) = .FALSE.
          INTERSECT_Z(IJK) = .FALSE.
          IF(I<=IMAX1) INTERSECT_X(IPJK) = .FALSE.
          IF(J<=JMAX1) INTERSECT_Y(IJPK) = .FALSE.
          IF(DO_K.AND.(K<=KMAX1)) INTERSECT_Z(IJKP) = .FALSE.
       ENDIF


      RETURN

      END SUBROUTINE REMOVE_INTERSECT_FLAG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_CUT_CELL_FILES                                    C
!  Purpose: Open CUT CELL related file                                 C
!           and writes headers                                         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_CUT_CELL_FILES

      USE cutcell
      USE compar

      IMPLICIT NONE

      IF(MyPE == PE_IO)  THEN
         OPEN(CONVERT='BIG_ENDIAN',UNIT = UNIT_CUT_CELL_LOG, FILE = 'CUT_CELL.LOG')
      ENDIF

      RETURN


      END SUBROUTINE OPEN_CUT_CELL_FILES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_CUT_CELL_FILES                                   C
!  Purpose: Close CUT CELL related file                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_CUT_CELL_FILES

      USE compar
      USE cutcell

      IMPLICIT NONE

      IF(MyPE == PE_IO) THEN
         CLOSE(UNIT_CUT_CELL_LOG)
      ENDIF

      RETURN


      END SUBROUTINE CLOSE_CUT_CELL_FILES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CAD_INTERSECT                                          C
!  Purpose: Intersects CAD (STL file or MSH file) geometry with grid   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CAD_INTERSECT(TYPE_OF_CELL,Xint,Yint,Zint)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE stl

      USE mpi_utility
      USE functions

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K
      INTEGER :: IM,IP,JM,JP,KM,KP,IMJK,IPJK,IJMK,IJPK,IJKM,IJKP
      INTEGER :: IJPKP,IPJKP,IPJPK
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb,xc,yc,zc
      LOGICAL :: INTERSECT_FLAG,INSIDE_FACET_a,INSIDE_FACET_b

      DOUBLE PRECISION :: X1,X2,Y1,Y2,Z1,Z2

      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: Xint,Yint,Zint

      INTEGER :: NN,I1,I2,J1,J2,K1,K2

      DOUBLE PRECISION :: X_OFFSET, Y_OFFSET, Z_OFFSET

      DOUBLE PRECISION, DIMENSION(3) :: N4,N6,N7,N8
      DOUBLE PRECISION :: CURRENT_F

      INTEGER :: N_UNDEFINED, NTOTAL_UNDEFINED,N_PROP
      INTEGER, PARAMETER :: N_PROPMAX=1000
      LOGICAL:: F_FOUND

!      CHARACTER (LEN=3) :: CAD_PROPAGATE_ORDER

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.

      Xint = UNDEFINED
      Yint = UNDEFINED
      Zint = UNDEFINED

      F_AT = UNDEFINED

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')


            X_OFFSET = ZERO
            Y_OFFSET = ZERO
            Z_OFFSET = ZERO

         CASE('U_MOMENTUM')

            X_OFFSET = HALF
            Y_OFFSET = ZERO
            Z_OFFSET = ZERO

         CASE('V_MOMENTUM')

            X_OFFSET = ZERO
            Y_OFFSET = HALF
            Z_OFFSET = ZERO

         CASE('W_MOMENTUM')

            X_OFFSET = ZERO
            Y_OFFSET = ZERO
            Z_OFFSET = HALF


         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_CELL_NODE_COORDINATES'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'SCALAR'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT


      DO NN = 1,N_FACETS


         X1 = MINVAL(VERTEX(1:3,1,NN))
         X2 = MAXVAL(VERTEX(1:3,1,NN))
         Y1 = MINVAL(VERTEX(1:3,2,NN))
         Y2 = MAXVAL(VERTEX(1:3,2,NN))
         Z1 = MINVAL(VERTEX(1:3,3,NN))
         Z2 = MAXVAL(VERTEX(1:3,3,NN))


         I1 = IEND3
         I2 = ISTART3

         IF(X2>=ZERO-DX(ISTART3).AND.X1<=XLENGTH+DX(IEND3)) THEN
            DO I = ISTART3, IEND3
               IP = I+1
               IF(XG_E(I)+X_OFFSET*DX(IP)>=X1-TOL_STL) THEN
                  I1=I
                  EXIT
               ENDIF
            ENDDO

            DO I = IEND3, ISTART3,-1
               IP = I+1
               IF(XG_E(I)-DX(I)+X_OFFSET*DX(IP)<=X2+TOL_STL) THEN
                  I2=I
                  EXIT
               ENDIF
            ENDDO
         ENDIF


         J1 = JEND3
         J2 = JSTART3

         IF(Y2>=ZERO-DY(JSTART3).AND.Y1<=YLENGTH+DY(JEND3)) THEN
            DO J = JSTART3, JEND3
               JP = J+1
               IF(YG_N(J)+Y_OFFSET*DY(JP)>=Y1-TOL_STL) THEN
                  J1=J
                  EXIT
               ENDIF
            ENDDO

            DO J = JEND3, JSTART3,-1
               JP=J+1
               IF(YG_N(J)-DY(J)+Y_OFFSET*DY(JP)<=Y2+TOL_STL) THEN
                  J2=J
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         K1 = KEND3
         K2 = KSTART3

         IF(Z2>=ZERO-DZ(KSTART3).AND.Z1<=ZLENGTH+DZ(KEND3)) THEN
            DO K = KSTART3, KEND3
               KP=K+1

               IF(ZG_T(K)+Z_OFFSET*DZ(KP)>=Z1-TOL_STL) THEN
                  K1=K
                  EXIT
               ENDIF
            ENDDO

            DO K = KEND3, KSTART3,-1
               KP = K+1
               IF(ZG_T(K)-DZ(K)+Z_OFFSET*DZ(KP)<=Z2+TOL_STL) THEN
                  K2=K
                  EXIT
               ENDIF
            ENDDO
         ENDIF




         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IJK = FUNIJK(I,J,K)


                  IM = MAX0(I - 1 ,ISTART3)
                  JM = MAX0(J - 1 ,JSTART3)
                  KM = MAX0(K - 1 ,KSTART3)

                  IP = MIN0(I + 1 ,IEND3)
                  JP = MIN0(J + 1 ,JEND3)
                  KP = MIN0(K + 1 ,KEND3)


                  IMJK = FUNIJK(IM,J,K)
                  IPJK = FUNIJK(IP,J,K)
                  IJMK = FUNIJK(I,JM,K)
                  IJPK = FUNIJK(I,JP,K)
                  IJKM = FUNIJK(I,J,KM)
                  IJKP = FUNIJK(I,J,KP)

                  IJPKP = FUNIJK(I,JP,KP)
                  IPJKP = FUNIJK(IP,J,KP)
                  IPJPK = FUNIJK(IP,JP,K)


!======================================================================
!  Get coordinates of eight nodes
!======================================================================

                  CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

!======================================================================
!  Intersection with Edge 7 (node 7-8, Face North-Top):
!======================================================================
                  xa = X_NODE(7)
                  ya = Y_NODE(7)
                  za = Z_NODE(7)

                  xb = X_NODE(8)
                  yb = Y_NODE(8)
                  zb = Z_NODE(8)

! Check if intersection occurs at corners

                  CALL IS_POINT_INSIDE_FACET(xa,ya,za,NN,INSIDE_FACET_a)

                  IF(INSIDE_FACET_a) THEN   ! corner intersection at node 7

                     F_AT(IMJK) = ZERO

                     IF(TRIM(TYPE_OF_CELL).eq.'SCALAR')  CALL ADD_FACET_AND_SET_BC_ID(IJK,NN)

                  ENDIF


                  CALL IS_POINT_INSIDE_FACET(xb,yb,zb,NN,INSIDE_FACET_b)

                  IF(INSIDE_FACET_b) THEN   ! corner intersection at node 8

                     F_AT(IJK) = ZERO

                     IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,NN)

                  ENDIF



! Check intersection within line 7-8, excluding corners

                  INTERSECT_FLAG = .FALSE.

                  IF(.NOT.(INTERSECT_X(IJK).OR.INSIDE_FACET_a.OR.INSIDE_FACET_b)) THEN
                     CALL INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,NN,INTERSECT_FLAG,xc,yc,zc)
                  ENDIF

                  IF(INTERSECT_FLAG) THEN
                     IF(INTERSECT_X(IJK)) THEN
                        IF(DABS(Xint(IJK)-xc)>TOL_STL) THEN

                           ! Ignore intersections when two intersections are detected on the same edge
                           INTERSECT_X(IJK) = .FALSE.
                        ENDIF
                     ELSE
                        INTERSECT_X(IJK) = .TRUE.
                        Xint(IJK) = xc

! Set values at corners if they are not zero

                        N7(1) = xa-xc
                        N7(2) = ya-yc
                        N7(3) = za-zc

                        IF(DABS(F_AT(IMJK))>TOL_F)   F_AT(IMJK) = -DOT_PRODUCT(N7,NORM_FACE(:,NN))

                        N8(1) = xb-xc
                        N8(2) = yb-yc
                        N8(3) = zb-zc

                        IF(DABS(F_AT(IJK))>TOL_F)   F_AT(IJK) = -DOT_PRODUCT(N8,NORM_FACE(:,NN))


                        IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,NN)
                        IF(JP<=J2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJPK,NN)
                        IF(KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJKP,NN)
                        IF(JP<=J2.AND.KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJPKP,NN)
                     ENDIF
                  ENDIF

                  IF(TYPE_OF_CELL=='U_MOMENTUM') THEN
                     IF(SNAP(IJK)) THEN
                        INTERSECT_X(IJK) = .TRUE.
                        Xn_int(IJK) = XG_E(I)
                     ENDIF
                  ENDIF

!======================================================================
!  Intersection with Edge 6 (node 6-8, Face East-Top):
!======================================================================
                  xa = X_NODE(6)
                  ya = Y_NODE(6)
                  za = Z_NODE(6)

! Check if intersection occurs at corners

                  CALL IS_POINT_INSIDE_FACET(xa,ya,za,NN,INSIDE_FACET_a)

                  IF(INSIDE_FACET_a) THEN   ! corner intersection at node 6

                     F_AT(IJMK) = ZERO

                     IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,NN)

                  ENDIF



! Check intersection within line 6-8, excluding corners



                  INTERSECT_FLAG = .FALSE.

                  IF(.NOT.(INTERSECT_Y(IJK).OR.INSIDE_FACET_a.OR.INSIDE_FACET_b)) THEN
                     CALL INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,NN,INTERSECT_FLAG,xc,yc,zc)
                  ENDIF


                  IF(INTERSECT_FLAG) THEN

                     IF(INTERSECT_Y(IJK)) THEN

                        IF(DABS(Yint(IJK)-yc)>TOL_STL) THEN

                           INTERSECT_Y(IJK) = .FALSE. ! Ignore intersections when two intersections are detected on the same edge

                        ENDIF

                     ELSE


                        INTERSECT_Y(IJK) = .TRUE.
                        Yint(IJK) = yc

! Set values at corners if they are not zero

                        N6(1) = xa-xc
                        N6(2) = ya-yc
                        N6(3) = za-zc

                        IF(DABS(F_AT(IJMK))>TOL_F)   F_AT(IJMK) = -DOT_PRODUCT(N6,NORM_FACE(:,NN))

                        N8(1) = xb-xc
                        N8(2) = yb-yc
                        N8(3) = zb-zc

                        IF(DABS(F_AT(IJK))>TOL_F)   F_AT(IJK) = -DOT_PRODUCT(N8,NORM_FACE(:,NN))


                        IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,NN)
                        IF(IP<=I2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJK,NN)
                        IF(KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJKP,NN)
                        IF(IP<=I2.AND.KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJKP,NN)

                     ENDIF

                  ENDIF

                  IF(TYPE_OF_CELL=='V_MOMENTUM') THEN
                     IF(SNAP(IJK)) THEN
                        INTERSECT_Y(IJK) = .TRUE.
                        Ye_int(IJK) = YG_N(J)
                     ENDIF
                  ENDIF

                  IF(DO_K) THEN
!======================================================================
!  Intersection with Edge 11 (node 4-8, Face East-North):
!======================================================================
                     xa = X_NODE(4)
                     ya = Y_NODE(4)
                     za = Z_NODE(4)

! Check if intersection occurs at corners

                     CALL IS_POINT_INSIDE_FACET(xa,ya,za,NN,INSIDE_FACET_a)

                     IF(INSIDE_FACET_a) THEN   ! corner intersection at node 4

                        F_AT(IJKM) = ZERO

                        IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,NN)

                     ENDIF




! Check intersection within line 4-8, excluding corners


                     INTERSECT_FLAG = .FALSE.

                     IF(.NOT.(INTERSECT_Z(IJK).OR.INSIDE_FACET_a.OR.INSIDE_FACET_b)) THEN
                        CALL INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,NN,INTERSECT_FLAG,xc,yc,zc)
                     ENDIF

                     IF(INTERSECT_FLAG) THEN

                        IF(INTERSECT_Z(IJK)) THEN

                           IF(DABS(Zint(IJK)-zc)>TOL_STL) THEN

                              INTERSECT_Z(IJK) = .FALSE. ! Ignore intersections when two intersections are detected on the same edge

                           ENDIF

                        ELSE


                           INTERSECT_Z(IJK) = .TRUE.
                           Zint(IJK) = zc


! Set values at corners if they are not zero

                           N4(1) = xa-xc
                           N4(2) = ya-yc
                           N4(3) = za-zc

                           IF(DABS(F_AT(IJKM))>TOL_F)   F_AT(IJKM) = -DOT_PRODUCT(N4,NORM_FACE(:,NN))

                           N8(1) = xb-xc
                           N8(2) = yb-yc
                           N8(3) = zb-zc

                           IF(DABS(F_AT(IJK))>TOL_F)   F_AT(IJK) = -DOT_PRODUCT(N8,NORM_FACE(:,NN))



                           IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,NN)
                           IF(IP<=I2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJK,NN)
                           IF(JP<=J2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJPK,NN)
                           IF(IP<=I2.AND.JP<=J2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJPK,NN)

                        ENDIF

                     ENDIF


                     IF(TYPE_OF_CELL=='W_MOMENTUM') THEN
                        IF(SNAP(IJK)) THEN
                           INTERSECT_Z(IJK) = .TRUE.
                           Zt_int(IJK) = ZG_T(K)
                        ENDIF
                     ENDIF

                  ENDIF


               ENDDO  ! I loop
            ENDDO  ! J loop
         ENDDO  ! K loop


      ENDDO  ! Loop over facets

      CURRENT_F = UNDEFINED

! Overwrite small values to set them to zero

      DO IJK = IJKSTART3, IJKEND3
         IF(DABS(F_AT(IJK))<TOL_STL) THEN
            F_AT(IJK)=ZERO
         ENDIF
      END DO

! Propagates node values to all interior cells
! in order defined by CAD_PROPAGATE_ORDER (outer loop)


      call send_recv(F_AT,2)
      call send_recv(Xint,2)
      call send_recv(Yint,2)
      call send_recv(Zint,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_X,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Y,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Z,2)
!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================

      DO IJK = IJKSTART3, IJKEND3

!        IF(INTERIOR_CELL_AT(IJK)) THEN

!            IF(POTENTIAL_CUT_CELL_AT(IJK))  CALL CLEAN_INTERSECT(IJK,'SCALAR',Xn_int(IJK),Ye_int(IJK),Zt_int(IJK))

            CALL CLEAN_INTERSECT(IJK,TYPE_OF_CELL,Xint(IJK),Yint(IJK),Zint(IJK))

!        ENDIF

      END DO

      call send_recv(F_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(SNAP,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_X,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Y,2)
      call SEND_RECEIVE_1D_LOGICAL(INTERSECT_Z,2)

      SELECT CASE (CAD_PROPAGATE_ORDER)

      CASE ('   ')

! After intersecting the edges of the background mesh with the STL facets,
! the end points (i.e., cell corners) are assigned a value, called F_AT, where:
! F_AT = zero if the corner point is on a facet (within some tolerance TOL_STL),
! F_AT < zero if the corner point is inside  the fluid region,
! F_AT > zero if the corner point is outside the fluid region.
! At this point F_AT is only defined across edges that intersect the STL facets,
! and it must be propagated to all cell corners to determine if uncut cells
! are inside or outside the fluid region.

! Only F_AT values that are defined and not zero are propagated to their direct
! neighbors, if it is not already defined. The propagation is repeated
! at most N_PROPMAX. The loop is exited when all F_AT values are defined.
! N_PROPMAX could be increased for very large domains.
! The propagation of F_AT will stop anytime a boundary is encountered since F_AT
! changes sign across a boundary.
!
         DO N_PROP=1,N_PROPMAX

            DO IJK = IJKSTART3, IJKEND3

! Aaron
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               IF(.NOT.IS_ON_myPE_plus1layer(I,J,K))cycle
! End aaron

               IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN

                     IMJK = IM_OF(IJK)
                     IF(F_AT(IMJK)==UNDEFINED.AND.F_AT(IMJK)/=ZERO)  F_AT(IMJK)=F_AT(IJK)

                     IPJK = IP_OF(IJK)
                     IF(F_AT(IPJK)==UNDEFINED.AND.F_AT(IPJK)/=ZERO)  F_AT(IPJK)=F_AT(IJK)

                     IJMK = JM_OF(IJK)
                     IF(F_AT(IJMK)==UNDEFINED.AND.F_AT(IJMK)/=ZERO)  F_AT(IJMK)=F_AT(IJK)

                     IJPK = JP_OF(IJK)
                     IF(F_AT(IJPK)==UNDEFINED.AND.F_AT(IJPK)/=ZERO)  F_AT(IJPK)=F_AT(IJK)

                     IJKM = KM_OF(IJK)
                     IF(F_AT(IJKM)==UNDEFINED.AND.F_AT(IJKM)/=ZERO)  F_AT(IJKM)=F_AT(IJK)

                     IJKP = KP_OF(IJK)
                     IF(F_AT(IJKP)==UNDEFINED.AND.F_AT(IJKP)/=ZERO)  F_AT(IJKP)=F_AT(IJK)

               ENDIF

            ENDDO ! IJK Loop


! Communicate F_AT accross processors for DMP runs
            call send_recv(F_AT,2)

! Count the number of undefined values of F_AT
! and exit loop if all values of F_AT have been propagated
            N_UNDEFINED = 0
            DO IJK = IJKSTART3, IJKEND3
               IF(INTERIOR_CELL_AT(IJK).AND.F_AT(IJK)==UNDEFINED) N_UNDEFINED = N_UNDEFINED + 1
            ENDDO

            call global_all_sum( N_UNDEFINED, NTOTAL_UNDEFINED )
            IF(NTOTAL_UNDEFINED==0) EXIT

         ENDDO ! N_PROP Loop


         call send_recv(F_AT,2)

! If a process still has undefined values of F_AT, this probably means
! that all cells belonging to that process are dead cells.
         IF(N_UNDEFINED>0) THEN
            WRITE(*,*)'WARNING: UNABLE TO PROPAGATE F_AT ARRAY FROM myPE=.', MyPE
            WRITE(*,*)'         THIS USUALLY INDICATE A RANK WITH NO ACTIVE CELL'
            WRITE(*,*)'         YOU MAY NEED TO ADJUST THE GRID PARTITIONNING'
            WRITE(*,*)'         TO GET BETTER LOAD_BALANCE.'
         ENDIF



      CASE ('IJK')

         DO I=ISTART3,IEND3
            DO K=KSTART3,KEND3
               F_FOUND = .FALSE.
               DO J=JSTART3,JEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO J=JSTART3,JEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO J=JSTART3,JEND3
            DO I=ISTART3,IEND3
               F_FOUND = .FALSE.
               DO K=KSTART3,KEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO K=KSTART3,KEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO K=KSTART3,KEND3
            DO J=JSTART3,JEND3
               F_FOUND = .FALSE.
               DO I=ISTART3,IEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO I=ISTART3,IEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

      CASE ('JKI')

         DO J=JSTART3,JEND3
            DO I=ISTART3,IEND3
               F_FOUND = .FALSE.
               DO K=KSTART3,KEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO K=KSTART3,KEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO K=KSTART3,KEND3
            DO J=JSTART3,JEND3
               F_FOUND = .FALSE.
               DO I=ISTART3,IEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO I=ISTART3,IEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO I=ISTART3,IEND3
            DO K=KSTART3,KEND3
               F_FOUND = .FALSE.
               DO J=JSTART3,JEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO J=JSTART3,JEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)


      CASE ('KIJ')



         DO K=KSTART3,KEND3
            DO J=JSTART3,JEND3
               F_FOUND = .FALSE.
               DO I=ISTART3,IEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO I=ISTART3,IEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO I=ISTART3,IEND3
            DO K=KSTART3,KEND3
               F_FOUND = .FALSE.
               DO J=JSTART3,JEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO J=JSTART3,JEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO


         call send_recv(F_AT,2)

         DO J=JSTART3,JEND3
            DO I=ISTART3,IEND3
               F_FOUND = .FALSE.
               DO K=KSTART3,KEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT
                  ENDIF
               ENDDO
               IF(F_FOUND) THEN
                  DO K=KSTART3,KEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO


         call send_recv(F_AT,2)

         CASE DEFAULT
            IF(myPE == PE_IO) THEN
               WRITE(*,*)'CAD_INTERSECT.'
               WRITE(*,*)'UNKNOWN CAD_PROPAGATE_ORDER:',CAD_PROPAGATE_ORDER
               WRITE(*,*)'ACCEPTABLE VALUES:'
               WRITE(*,*)'IJK'
               WRITE(*,*)'JKI'
               WRITE(*,*)'KIJ'
            ENDIF
!            CALL MFIX_EXIT(myPE)

      END SELECT

      call send_recv(F_AT,2)

      RETURN

      END SUBROUTINE CAD_INTERSECT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_FACET                                              C
!  Purpose: Add facet to list in IJK scalar cell                       C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 15-Oct-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE ADD_FACET_AND_SET_BC_ID(IJK,NN)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE stl
      USE mpi_utility

      IMPLICIT NONE
      INTEGER :: IJK,NN

      BC_ID(IJK) = BC_ID_STL_FACE(NN)             ! Set tentative BC_ID

      IF(N_FACET_AT(IJK)<DIM_FACETS_PER_CELL) THEN

         N_FACET_AT(IJK) = N_FACET_AT(IJK) + 1
         LIST_FACET_AT(IJK,N_FACET_AT(IJK)) = NN

      ELSE

         WRITE(*,*) ' FATAL ERROR: TOO MANY FACETS IN CELL: ', IJK
         CALL MFIX_EXIT(myPE)

      ENDIF

      END SUBROUTINE ADD_FACET_AND_SET_BC_ID
   END MODULE CUT_CELL_PREPROC
