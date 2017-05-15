!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_INCREMENTS                                         !
!  Author: M. Syamlal, W. Rogers                      Date: 10-DEC-91  !
!                                                                      !
!  Purpose: The purpose of this module is to create increments to be   !
!           stored in the array STORE_INCREMENT which will be added    !
!           to cell index ijk to find the effective indices of its     !
!           neighbors. These increments are found using the 'class'    !
!           of cell ijk. The class is determined based on the          !
!           neighboring cell type, i.e. wall or fluid.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_INCREMENTS

      USE compar
      USE cutcell, ONLY: CARTESIAN_GRID
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IMJK, IPJK, IJKW, IJKE  ! I+, I-, east/west
      INTEGER :: IJMK, IJPK, IJKS, IJKN  ! J+, J-, north/south
      INTEGER :: IJKM, IJKP, IJKB, IJKT  ! K+, K-, top/bottom
! DO-loop index, ranges from 1 to ICLASS
      INTEGER :: IC
! Index for the solids phase.
      INTEGER :: M
! Local DO-loop index
      INTEGER :: L
! Index denoting cell class
      INTEGER :: ICLASS
! Array of sum of increments to make the class determination faster.
      INTEGER :: DENOTE_CLASS(MAX_CLASS)
! Flags for using the 'real' I/J/K value (not cyclic.)
      LOGICAL :: SHIFT
! Used for checking iteration over core cells
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: ALREADY_VISITED
      INTEGER :: interval, j_start(2), j_end(2)
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_INCREMENTS")

! Allocate increment arrays and report an allocation errors.
      CALL ALLOCATE_ARRAYS_INCREMENTS

! Initialize the default values to Undefined_I
      IP1(:) = UNDEFINED_I
      IM1(:) = UNDEFINED_I
      JP1(:) = UNDEFINED_I
      JM1(:) = UNDEFINED_I
      KP1(:) = UNDEFINED_I
      KM1(:) = UNDEFINED_I

      DO I = ISTART3, IEND3
         SHIFT = .NOT.(I==IMIN3 .OR. I==IMIN2 .OR. &
                       I==IMAX3 .OR. I==IMAX2)

         IF(CYCLIC_X .AND. NODESI.EQ.1 .AND. DO_I .AND. SHIFT) THEN
            IP1(I) = IMAP_C(IMAP_C(I)+1)
            IM1(I) = IMAP_C(IMAP_C(I)-1)
         ELSE
            IM1(I) = MAX(ISTART3, I - 1)
            IP1(I) = MIN(IEND3,   I + 1)
         ENDIF
      ENDDO

      DO J = JSTART3, JEND3

         SHIFT = .NOT.(J==JMIN3 .OR. J==JMIN2 .OR. &
                       J==JMAX3 .OR. J==JMAX2)

         IF (CYCLIC_Y .AND. NODESJ.EQ.1 .AND. DO_J .AND. SHIFT) THEN
            JP1(J) = JMAP_C(JMAP_C(J)+1)
            JM1(J) = JMAP_C(JMAP_C(J)-1)
         ELSE
            JM1(J) = MAX(JSTART3,J - 1)
            JP1(J) = MIN(JEND3,  J + 1)
         ENDIF
      ENDDO


      DO K = KSTART3, KEND3

         SHIFT = .NOT.(K==KMIN3 .OR. K==KMIN2 .OR. &
                       K==KMAX3 .OR. K==KMAX2)

         IF(CYCLIC_Z .AND. NODESK.EQ.1 .AND. DO_K .AND. SHIFT) THEN
            KP1(K) = KMAP_C(KMAP_C(K)+1)
            KM1(K) = KMAP_C(KMAP_C(K)-1)
         ELSE
            KM1(K) = MAX(KSTART3,K - 1)
            KP1(K) = MIN(KEND3,K + 1)
         ENDIF
      ENDDO

! Loop over all cells
      DO K = KSTART3, KEND3
      DO J = JSTART3, JEND3
      DO I = ISTART3, IEND3

         IJK = FUNIJK(I,J,K)  ! Find value of IJK

         I_OF(IJK) = I
         J_OF(IJK) = J
         K_OF(IJK) = K

      ENDDO
      ENDDO
      ENDDO

      ICLASS = 0

! Loop over all cells (minus the ghost layers)
      DO K = KSTART3, KEND3
      DO J = JSTART3, JEND3
      L100: DO I = ISTART3, IEND3

         IJK = FUNIJK(I,J,K)

! Find the the effective cell-center indices for all neighbor cells
         CALL SET_INDEX1A (I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, &
            IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT)

         ICLASS = ICLASS + 1               !Increment the ICLASS counter
         IF(ICLASS > MAX_CLASS) THEN
            WRITE(ERR_MSG, 1200) trim(iVal(MAX_CLASS))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1200 FORMAT('Error 1200: The number of classes has exceeded the ',    &
         'maximum: ',A,/'Increase the MAX_CLASS parameter in param1',  &
         '_mod.f and recompile.')

         INCREMENT_FOR_N(ICLASS)  = IJKN - IJK
         INCREMENT_FOR_S(ICLASS)  = IJKS - IJK
         INCREMENT_FOR_E(ICLASS)  = IJKE - IJK
         INCREMENT_FOR_W(ICLASS)  = IJKW - IJK
         INCREMENT_FOR_T(ICLASS)  = IJKT - IJK
         INCREMENT_FOR_B(ICLASS)  = IJKB - IJK

         INCREMENT_FOR_IM(ICLASS) = IMJK - IJK
         INCREMENT_FOR_IP(ICLASS) = IPJK - IJK
         INCREMENT_FOR_JM(ICLASS) = IJMK - IJK
         INCREMENT_FOR_JP(ICLASS) = IJPK - IJK
         INCREMENT_FOR_KM(ICLASS) = IJKM - IJK
         INCREMENT_FOR_KP(ICLASS) = IJKP - IJK

         INCREMENT_FOR_NB(1,ICLASS) = INCREMENT_FOR_E(ICLASS)
         INCREMENT_FOR_NB(2,ICLASS) = INCREMENT_FOR_W(ICLASS)
         INCREMENT_FOR_NB(3,ICLASS) = INCREMENT_FOR_S(ICLASS)
         INCREMENT_FOR_NB(4,ICLASS) = INCREMENT_FOR_N(ICLASS)
         INCREMENT_FOR_NB(5,ICLASS) = INCREMENT_FOR_B(ICLASS)
         INCREMENT_FOR_NB(6,ICLASS) = INCREMENT_FOR_T(ICLASS)

         INCREMENT_FOR_MP(1,ICLASS) = INCREMENT_FOR_IM(ICLASS)
         INCREMENT_FOR_MP(2,ICLASS) = INCREMENT_FOR_IP(ICLASS)
         INCREMENT_FOR_MP(3,ICLASS) = INCREMENT_FOR_JM(ICLASS)
         INCREMENT_FOR_MP(4,ICLASS) = INCREMENT_FOR_JP(ICLASS)
         INCREMENT_FOR_MP(5,ICLASS) = INCREMENT_FOR_KM(ICLASS)
         INCREMENT_FOR_MP(6,ICLASS) = INCREMENT_FOR_KP(ICLASS)


         DENOTE_CLASS(ICLASS) = INCREMENT_FOR_N(ICLASS) + INCREMENT_FOR_S&
            (ICLASS) + INCREMENT_FOR_E(ICLASS) + INCREMENT_FOR_W(ICLASS)&
             + INCREMENT_FOR_T(ICLASS) + INCREMENT_FOR_B(ICLASS) + &
            INCREMENT_FOR_IM(ICLASS) + INCREMENT_FOR_IP(ICLASS) + &
            INCREMENT_FOR_JM(ICLASS) + INCREMENT_FOR_JP(ICLASS) + &
            INCREMENT_FOR_KM(ICLASS) + INCREMENT_FOR_KP(ICLASS)

         CELL_CLASS(IJK) = ICLASS

! Place the cell in a class based on its DENOTE_CLASS(ICLASS) value
         DO IC = 1, ICLASS - 1             !Loop over previous and present classes
!                                                !IF a possible match in cell types
            IF(DENOTE_CLASS(ICLASS) == DENOTE_CLASS(IC)) THEN
!                                                !is found, compare all increments
               IF(INCREMENT_FOR_N(ICLASS) /= INCREMENT_FOR_N(IC)) CYCLE
               IF(INCREMENT_FOR_S(ICLASS) /= INCREMENT_FOR_S(IC)) CYCLE
               IF(INCREMENT_FOR_E(ICLASS) /= INCREMENT_FOR_E(IC)) CYCLE
               IF(INCREMENT_FOR_W(ICLASS) /= INCREMENT_FOR_W(IC)) CYCLE
               IF(INCREMENT_FOR_T(ICLASS) /= INCREMENT_FOR_T(IC)) CYCLE
               IF(INCREMENT_FOR_B(ICLASS) /= INCREMENT_FOR_B(IC)) CYCLE
               IF(INCREMENT_FOR_IM(ICLASS) /= INCREMENT_FOR_IM(IC)) CYCLE
               IF(INCREMENT_FOR_IP(ICLASS) /= INCREMENT_FOR_IP(IC)) CYCLE
               IF(INCREMENT_FOR_JM(ICLASS) /= INCREMENT_FOR_JM(IC)) CYCLE
               IF(INCREMENT_FOR_JP(ICLASS) /= INCREMENT_FOR_JP(IC)) CYCLE
               IF(INCREMENT_FOR_KM(ICLASS) /= INCREMENT_FOR_KM(IC)) CYCLE
               IF(INCREMENT_FOR_KP(ICLASS) /= INCREMENT_FOR_KP(IC)) CYCLE
               CELL_CLASS(IJK) = IC        !Assign cell to a class
               ICLASS = ICLASS - 1
               CYCLE  L100                 !Go to next cell
            ENDIF
         END DO

      ENDDO L100
      ENDDO
      ENDDO

      DO M = 1, MMAX
      DO L = M, MMAX
         IF(L == M) THEN
            STORE_LM(L,M) = 0
         ELSE
            STORE_LM(L,M) = M + (L - 2)*(L - 1)/2
            STORE_LM(M,L) = M + (L - 2)*(L - 1)/2
         ENDIF
      ENDDO
      ENDDO

      USE_CORECELL_LOOP = .not.CARTESIAN_GRID

      if (USE_CORECELL_LOOP) then
         Allocate( already_visited(DIMENSION_3))
         already_visited(:) = .false.

         core_istart = istart+2
         core_iend = iend-2

         core_jstart = jstart+2
         core_jend = jend-2

         if (do_k) then
            core_kstart = kstart+2
            core_kend = kend-2
         else
            core_kstart = 1
            core_kend = 1
            kstart = 1
            kend = 1
         endif

         iclass = cell_class(funijk(core_istart,core_jstart,core_kstart))

         outer: do k = core_kstart,core_kend
            do i = core_istart,core_iend
               do j = core_jstart,core_jend
                  IJK = funijk(i,j,k)
                  ! this shouldn't happen, but we might as well check
                  if (ijk.ne. (j + c0 + i*c1 + k*c2)) then
                     USE_CORECELL_LOOP = .false.
                     exit outer
                  endif

                  ijk = (j + c0 + i*c1 + k*c2)

                  if (already_visited(ijk)) then
                     USE_CORECELL_LOOP = .false.
                     exit outer
                  endif
                  already_visited(ijk) = .true.

                  if (iclass.ne.cell_class(ijk)) then
                     USE_CORECELL_LOOP = .false.
                     exit outer
                  endif
               enddo
            enddo
         enddo outer

         j_start(1) = jstart
         j_end(1) = jend
         j_start(2) = 0 ! no iterations
         j_end(2) = -1  ! no iterations

         outer2: do k = kstart,kend
            do i = istart,iend

               if  (USE_CORECELL_LOOP) then
                  if (core_istart<= i .and. i <= core_iend .and. core_kstart <= k .and. k<=core_kend) then
                     j_start(1) = jstart
                     j_end(1) = core_jstart-1
                     j_start(2) = core_jend+1
                     j_end(2) = jend
                  else
                     j_start(1) = jstart
                     j_end(1) = jend
                     j_start(2) = 0 ! no iterations
                     j_end(2) = -1  ! no iterations
                  endif
               endif

               do interval=1,2
                  do j = j_start(interval),j_end(interval)
                     if (already_visited(funijk(i,j,k))) then
                        USE_CORECELL_LOOP = .false.
                        exit outer2
                     endif
                     already_visited(funijk(i,j,k)) = .true.
                  enddo
               enddo
            enddo
         enddo outer2

         outer3: do k = kstart,kend
            do i = istart,iend
               do j = jstart,jend
                  if (.not.already_visited(funijk(i,j,k))) then
                     USE_CORECELL_LOOP = .false.
                     exit outer3
                  endif
               enddo
            enddo
         enddo outer3

         deallocate(already_visited)

      endif

      IF(.NOT.INCREMENT_ARRAYS_ALLOCATED) THEN
         allocate(WEST_ARRAY_OF(ijkstart3:ijkend3))
         allocate(EAST_ARRAY_OF(ijkstart3:ijkend3))
         allocate(SOUTH_ARRAY_OF(ijkstart3:ijkend3))
         allocate(NORTH_ARRAY_OF(ijkstart3:ijkend3))
         allocate(BOTTOM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(TOP_ARRAY_OF(ijkstart3:ijkend3))

         allocate(IM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(IP_ARRAY_OF(ijkstart3:ijkend3))
         allocate(JM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(JP_ARRAY_OF(ijkstart3:ijkend3))
         allocate(KM_ARRAY_OF(ijkstart3:ijkend3))
         allocate(KP_ARRAY_OF(ijkstart3:ijkend3))
      ENDIF

      INCREMENT_ARRAYS_ALLOCATED = .TRUE.

      DO IJK = ijkstart3,ijkend3
         WEST_ARRAY_OF(ijk)   = WEST_OF_0(IJK)
         EAST_ARRAY_OF(ijk)   = EAST_OF_0(IJK)
         SOUTH_ARRAY_OF(ijk)  = SOUTH_OF_0(IJK)
         NORTH_ARRAY_OF(ijk)  = NORTH_OF_0(IJK)
         BOTTOM_ARRAY_OF(ijk) = BOTTOM_OF_0(IJK)
         TOP_ARRAY_OF(ijk)    = TOP_OF_0(IJK)

         IM_ARRAY_OF(ijk) = IM_OF_0(IJK)
         IP_ARRAY_OF(ijk) = IP_OF_0(IJK)
         JM_ARRAY_OF(ijk) = JM_OF_0(IJK)
         JP_ARRAY_OF(ijk) = JP_OF_0(IJK)
         KM_ARRAY_OF(ijk) = KM_OF_0(IJK)
         KP_ARRAY_OF(ijk) = KP_OF_0(IJK)
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE SET_INCREMENTS






!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RE_INDEX_ARRAYS                                        C
!  Purpose: Remove dead cells from computation.                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-MAY-11  C
!  Reviewer:                                          Date: ##-###-##  C
!                                                                      C
!  Revision Number: #                                                  C
!  Purpose: ##########                                                 C
!  Author:  ##########                                Date: ##-###-##  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE RE_INDEX_ARRAYS
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bc, only: IJK_P_G
      USE cdist
      USE compar
      USE cutcell
      USE discretelement, only: DISCRETE_ELEMENT
      USE energy
      USE exit, only: mfix_exit
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE pgcor, only :       PHASE_4_P_G
      USE physprop
      USE pscor, only :       PHASE_4_P_S
      USE run
      USE scalars
      USE sendrecv
      USE stiff_chem, only: STIFF_CHEMISTRY,notOwner
      USE stl
      USE visc_g

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
!                      Indices
      INTEGER          I, J, K, IJK, NEW_IJK,NN
!
!                      Index for the solids phase.
      INTEGER          M

      LOGICAL,DIMENSION(DIMENSION_3) :: ANY_CUT_TREATMENT, ANY_STANDARD_CELL

      LOGICAL :: ANY_GLOBAL_GHOST_CELL,NEED_TO_SKIP_CELL

      INTEGER,DIMENSION(DIMENSION_3) :: IM_COPY,IP_COPY,JM_COPY,JP_COPY,KM_COPY,KP_COPY
      INTEGER,DIMENSION(DIMENSION_3) :: WEST_COPY,EAST_COPY,SOUTH_COPY,NORTH_COPY,BOTTOM_COPY,TOP_COPY

      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: TEMP_IJK_ARRAY_OF
      INTEGER, ALLOCATABLE, DIMENSION(:)     :: TEMP_I_OF,TEMP_J_OF,TEMP_K_OF

      INTEGER, ALLOCATABLE, DIMENSION(:) :: BACKGROUND_IJKEND3_ALL,NCPP_UNIFORM_ALL

      INTEGER :: iproc,IERR

      INTEGER :: I1,I2,J1,J2,K1,K2,jj,sendsize,send_pos,recvsize,recv_pos,n_total, IC
      INTEGER :: placeholder, new_nsend1, new_nsend2,new_nrecv1,new_nrecv2
      INTEGER :: nj1,nj2

      INTEGER, DIMENSION(26) :: new_send_size, new_recv_size

      integer, pointer, dimension(:) :: new_xsend1, new_sendijk1 , new_sendproc1, new_sendtag1, &
                                        new_xsend2, new_sendijk2 , new_sendproc2, new_sendtag2, &
                                        new_xrecv1, new_recvijk1 , new_recvproc1, new_recvtag1, &
                                        new_xrecv2, new_recvijk2 , new_recvproc2, new_recvtag2

      integer :: comm

      DOUBLE PRECISION, DIMENSION(0:NumPEs-1) :: DIFF_NCPP

      INTEGER :: IJKW,IJKE,IJKS,IJKN,IJKB,IJKT
      INTEGER :: IMJK,IPJK,IJMK,IJPK,IJKM,IJKP

!                             Array index denoting a cell class, it is a
!                             column of the array STORE_INCREMENTS
      INTEGER                 ICLASS
!
!                             Array of sum of increments to make the class
!                             determination faster.
      INTEGER                 DENOTE_CLASS(MAX_CLASS)

   INTEGER :: I_SIZE,J_SIZE,K_SIZE

!======================================================================
!   Loop through useful cells and save their index
!======================================================================

      allocate(BACKGROUND_IJK_OF(DIMENSION_3))
      allocate(IJK_OF_BACKGROUND(DIMENSION_3))

      allocate(TEMP_IJK_ARRAY_OF(ISTART3-1:IEND3+1,JSTART3-1:JEND3+1,KSTART3-1:KEND3+1))
      TEMP_IJK_ARRAY_OF = IJK_ARRAY_OF

      allocate(TEMP_I_OF(DIMENSION_3))
      allocate(TEMP_J_OF(DIMENSION_3))
      allocate(TEMP_K_OF(DIMENSION_3))

      allocate(BACKGROUND_IJKEND3_ALL(0:NumPEs-1))

      TEMP_I_OF = I_OF
      TEMP_J_OF = J_OF
      TEMP_K_OF = K_OF

      TEMP_IJK_ARRAY_OF = IJK_ARRAY_OF

      DEAD_CELL_AT = .FALSE.


      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: INFO: USE_DOLOOP was set to .TRUE.'
      USE_DOLOOP = .TRUE.


!      IF(.NOT.RE_INDEXING) THEN
!         print*,'Skipping re-indexing...'
!         GOTO 999
!      ENDIF

      NEW_IJK = IJKSTART3

      IJK_OF_BACKGROUND = -999

! Step 0: Indentify dead cells

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Indentifying dead cells ...'

      ANY_CUT_TREATMENT = .FALSE.
      ANY_STANDARD_CELL = .FALSE.

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

!         IF(MyPE == PE_IO) WRITE(*,*)'IJK progress=',IJK,DBLE(IJK)/IJKEND3

         ANY_CUT_TREATMENT(IJK) =     CUT_TREATMENT_AT(IJK)&
                                  .OR.CUT_U_TREATMENT_AT(IJK)&
                                  .OR.CUT_V_TREATMENT_AT(IJK)&
                                  .OR.CUT_W_TREATMENT_AT(IJK)


         ANY_STANDARD_CELL(IJK) =     STANDARD_CELL_AT(IJK)&
                                  .OR.STANDARD_U_CELL_AT(IJK)&
                                  .OR.STANDARD_V_CELL_AT(IJK)&
                                  .OR.STANDARD_W_CELL_AT(IJK)

         ANY_GLOBAL_GHOST_CELL = (I < IMIN1).OR.(I > IMAX1)&  ! only along global ghost cells (MIN and MAX indices)
                           .OR.(J < JMIN1).OR.(J > JMAX1)&
                           .OR.(K < KMIN1).OR.(K > KMAX1)


         IF(.NOT.(ANY_CUT_TREATMENT(IJK)&
            .OR.ANY_STANDARD_CELL(IJK)&
            .OR.ANY_GLOBAL_GHOST_CELL)) THEN

            DEAD_CELL_AT(I,J,K) = .TRUE.

            IF(I==IMIN1)  DEAD_CELL_AT(IMIN3:IMIN2,J,K) = .TRUE. ! Extend dead cells to global ghost layers
            IF(I==IMAX1)  DEAD_CELL_AT(IMAX2:IMAX3,J,K) = .TRUE.

            IF(J==JMIN1)  DEAD_CELL_AT(I,JMIN3:JMIN2,K) = .TRUE.
            IF(J==JMAX1)  DEAD_CELL_AT(I,JMAX2:JMAX3,K) = .TRUE.

            IF(K==KMIN1)  DEAD_CELL_AT(I,J,KMIN3:KMIN2) = .TRUE.
            IF(K==KMAX1)  DEAD_CELL_AT(I,J,KMAX2:KMAX3) = .TRUE.

         ENDIF



      ENDDO


      IF(.NOT.MINIMIZE_SEND_RECV) THEN

         DEAD_CELL_AT(ISTART3:ISTART1,JSTART3:JEND3,KSTART3:KEND3) = .FALSE. ! Try: Keep all send/recv layers
         DEAD_CELL_AT(IEND1:IEND3,JSTART3:JEND3,KSTART3:KEND3) = .FALSE.

         DEAD_CELL_AT(ISTART3:IEND3,JSTART3:JSTART1,KSTART3:KEND3) = .FALSE.
         DEAD_CELL_AT(ISTART3:IEND3,JEND1:JEND3,KSTART3:KEND3) = .FALSE.

         DEAD_CELL_AT(ISTART3:IEND3,JSTART3:JEND3,KSTART3:KSTART1) = .FALSE.
         DEAD_CELL_AT(ISTART3:IEND3,JSTART3:JEND3,KEND1:KEND3) = .FALSE.

      ENDIF



      IF(NO_K) THEN  ! Extend dead cells to corners of ghost layers  <---------------------  SHOULD IT BE SKIPPED  ??
         DO K =  KMIN3, KMAX3,-1
            IF(DEAD_CELL_AT(IMAX1  ,JMAX1  ,K))  DEAD_CELL_AT(IMAX2:IMAX3    ,JMAX2:JMAX3    ,K) = .TRUE.
            IF(DEAD_CELL_AT(IMAX1  ,JMIN1,K))  DEAD_CELL_AT(IMAX2:IMAX3    ,JMIN3:JMIN2,K) = .TRUE.
            IF(DEAD_CELL_AT(IMIN1,JMAX1  ,K))  DEAD_CELL_AT(IMIN3:IMIN2,JMAX2:JMAX3    ,K) = .TRUE.
            IF(DEAD_CELL_AT(IMIN1,JMIN1,K))  DEAD_CELL_AT(IMIN3:IMIN2,JMIN3:JMIN2,K) = .TRUE.
         ENDDO
      ENDIF



! Step 1: Put all send and receive layer cells in a contiguous block, needed for parallel run only








      IF(NODESI>1.AND.NODESJ==1.AND.NODESK==1) THEN   ! I-DECOMPOSITION ONLY
         IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Arranging all send and receive layer cells to optimize I-decomposition ...'

! POTENTIAL RECEIVE LAYERS AT WEST
         DO I = ISTART3,ISTART2
            DO J= JSTART3,JEND3
               DO K = KSTART3, KEND3
                  IJK = FUNIJK(I,J,K)
                  IF(.NOT.DEAD_CELL_AT(I,J,K)) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

! POTENTIAL SEND LAYERS AT WEST
         DO I = ISTART1, ISTART1+1
            DO J= JSTART3,JEND3
               DO K = KSTART3, KEND3

                  IJK = FUNIJK(I,J,K)

                  IF( ANY_CUT_TREATMENT(IJK).OR.ANY_STANDARD_CELL(IJK).OR.(.NOT.DEAD_CELL_AT(I,J,K))) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     DEAD_CELL_AT(I,J,K) = .TRUE.
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

! POTENTIAL SEND LAYERS AT EAST
         DO I = IEND1-1,IEND1
            DO J= JSTART3,JEND3
               DO K = KSTART3, KEND3

                  IJK = FUNIJK(I,J,K)

                  IF( ANY_CUT_TREATMENT(IJK).OR.ANY_STANDARD_CELL(IJK).OR.(.NOT.DEAD_CELL_AT(I,J,K))) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     DEAD_CELL_AT(I,J,K) = .TRUE.
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

! POTENTIAL RECEIVE LAYERS AT EAST
         DO I = IEND2,IEND3
            DO J= JSTART3,JEND3
               DO K = KSTART3, KEND3
                  IJK = FUNIJK(I,J,K)
                  IF(.NOT.DEAD_CELL_AT(I,J,K)) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF
               ENDDO
            ENDDO
         ENDDO


         I1 = ISTART1 + 2
         I2 = IEND1   - 2

         J1 = JSTART3
         J2 = JEND3

         K1 = KSTART3
         K2 = KEND3


      ELSEIF(NODESJ>1.AND.NODESI==1.AND.NODESK==1) THEN  ! J-DECOMPOSITION ONLY
         IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Arranging all send and receive layer cells to optimize J-decomposition ...'

! POTENTIAL RECEIVE LAYERS AT SOUTH
         DO J = JSTART3,JSTART2
            DO I= ISTART3,IEND3
               DO K = KSTART3, KEND3
                  IJK = FUNIJK(I,J,K)
                  IF(.NOT.DEAD_CELL_AT(I,J,K)) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

! POTENTIAL SEND LAYERS AT SOUTH
         DO J = JSTART1, JSTART1+1
            DO I= ISTART3,IEND3
               DO K = KSTART3, KEND3

                  IJK = FUNIJK(I,J,K)

                  IF( ANY_CUT_TREATMENT(IJK).OR.ANY_STANDARD_CELL(IJK).OR.(.NOT.DEAD_CELL_AT(I,J,K))) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     DEAD_CELL_AT(I,J,K) = .TRUE.
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

! POTENTIAL SEND LAYERS AT NORTH
         DO J = JEND1-1,JEND1
            DO I= ISTART3,IEND3
               DO K = KSTART3, KEND3

                  IJK = FUNIJK(I,J,K)

                  IF( ANY_CUT_TREATMENT(IJK).OR.ANY_STANDARD_CELL(IJK).OR.(.NOT.DEAD_CELL_AT(I,J,K))) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     DEAD_CELL_AT(I,J,K) = .TRUE.
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

! POTENTIAL RECEIVE LAYERS AT NORTH
         DO J = JEND2,JEND3
            DO I= ISTART3,IEND3
               DO K = KSTART3, KEND3
                  IJK = FUNIJK(I,J,K)
                  IF(.NOT.DEAD_CELL_AT(I,J,K)) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF
               ENDDO
            ENDDO
         ENDDO


         I1 = ISTART3
         I2 = IEND3

         J1 = JSTART1 + 2
         J2 = JEND1   - 2

         K1 = KSTART3
         K2 = KEND3

      ELSEIF(NODESK>1.AND.NODESI==1.AND.NODESJ==1) THEN  ! K-DECOMPOSITION ONLY
         IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Arranging all send and receive layer cells to optimize K-decomposition ...'

! POTENTIAL RECEIVE LAYERS AT BOTTOM
         DO K = KSTART3,KSTART2
            DO J= JSTART3,JEND3
               DO I = ISTART3, IEND3
                  IJK = FUNIJK(I,J,K)
                  IF(.NOT.DEAD_CELL_AT(I,J,K)) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

! POTENTIAL SEND LAYERS AT BOTTOM
         DO K = KSTART1, KSTART1+1
            DO J= JSTART3,JEND3
               DO I = ISTART3, IEND3

                  IJK = FUNIJK(I,J,K)

                  IF( ANY_CUT_TREATMENT(IJK).OR.ANY_STANDARD_CELL(IJK).OR.(.NOT.DEAD_CELL_AT(I,J,K))) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     DEAD_CELL_AT(I,J,K) = .TRUE.
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

! POTENTIAL SEND LAYERS AT TOP
         DO K = KEND1-1,KEND1
            DO J= JSTART3,JEND3
               DO I = ISTART3, IEND3

                  IJK = FUNIJK(I,J,K)

                  IF( ANY_CUT_TREATMENT(IJK).OR.ANY_STANDARD_CELL(IJK).OR.(.NOT.DEAD_CELL_AT(I,J,K))) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     DEAD_CELL_AT(I,J,K) = .TRUE.
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

! POTENTIAL RECEIVE LAYERS AT TOP
         DO K = KEND2,KEND3
            DO J= JSTART3,JEND3
               DO I = ISTART3, IEND3
                  IJK = FUNIJK(I,J,K)
                  IF(.NOT.DEAD_CELL_AT(I,J,K)) THEN
                     CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
                  ELSE
                     IJK_OF_BACKGROUND(IJK) = -999
                  ENDIF
               ENDDO
            ENDDO
         ENDDO


         I1 = ISTART3
         I2 = IEND3

         J1 = JSTART3
         J2 = JEND3

         K1 = KSTART1 + 2
         K2 = KEND1   - 2



      ELSE                   ! SERIAL CASE OR DECOMPOSITION IN MORE THAN ONE DIRECTION



         I1 = ISTART3
         I2 = IEND3

         J1 = JSTART3
         J2 = JEND3

         K1 = KSTART3
         K2 = KEND3


      ENDIF


      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Arranging all interior cells in next contiguous block...'

! Step 2: Put all interior cells in next contiguous block

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)


         NEED_TO_SKIP_CELL =     (I < I1).OR.(I > I2)&
                             .OR.(J < J1).OR.(J > J2)

         IF(DO_K) NEED_TO_SKIP_CELL = (NEED_TO_SKIP_CELL.OR.(K < K1).OR.(K > K2))


         IF(NEED_TO_SKIP_CELL) CYCLE

         IF( .NOT.DEAD_CELL_AT(I,J,K)) THEN
            CALL  RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
         ELSE

            IJK_OF_BACKGROUND(IJK) = -999

         ENDIF

      ENDDO


      IJK_ARRAY_OF = TEMP_IJK_ARRAY_OF

!      FUNIJK = IJK_ARRAY_OF

      I_OF = TEMP_I_OF
      J_OF = TEMP_J_OF
      K_OF = TEMP_K_OF


! Save the old value of IJKEND3
      BACKGROUND_IJKEND3 = IJKEND3

      IJKEND3 = NEW_IJK - 1



      IM_COPY = IM_ARRAY_OF
      IP_COPY = IP_ARRAY_OF
      JM_COPY = JM_ARRAY_OF
      JP_COPY = JP_ARRAY_OF
      KM_COPY = KM_ARRAY_OF
      KP_COPY = KP_ARRAY_OF

      WEST_COPY   = WEST_ARRAY_OF
      EAST_COPY   = EAST_ARRAY_OF
      SOUTH_COPY  = SOUTH_ARRAY_OF
      NORTH_COPY  = NORTH_ARRAY_OF
      BOTTOM_COPY = BOTTOM_ARRAY_OF
      TOP_COPY    = TOP_ARRAY_OF


      DO IJK = IJKSTART3,IJKEND3
         IM_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(IM_COPY(BACKGROUND_IJK_OF(IJK)))
         IP_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(IP_COPY(BACKGROUND_IJK_OF(IJK)))
         JM_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(JM_COPY(BACKGROUND_IJK_OF(IJK)))
         JP_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(JP_COPY(BACKGROUND_IJK_OF(IJK)))
         KM_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(KM_COPY(BACKGROUND_IJK_OF(IJK)))
         KP_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(KP_COPY(BACKGROUND_IJK_OF(IJK)))

         WEST_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(WEST_COPY(BACKGROUND_IJK_OF(IJK)))
         EAST_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(EAST_COPY(BACKGROUND_IJK_OF(IJK)))
         SOUTH_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(SOUTH_COPY(BACKGROUND_IJK_OF(IJK)))
         NORTH_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(NORTH_COPY(BACKGROUND_IJK_OF(IJK)))
         BOTTOM_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(BOTTOM_COPY(BACKGROUND_IJK_OF(IJK)))
         TOP_ARRAY_OF(IJK) = IJK_OF_BACKGROUND(TOP_COPY(BACKGROUND_IJK_OF(IJK)))

         IF(IM_ARRAY_OF(IJK)==-999) IM_ARRAY_OF(IJK)=IJK
         IF(IP_ARRAY_OF(IJK)==-999) IP_ARRAY_OF(IJK)=IJK
         IF(JM_ARRAY_OF(IJK)==-999) JM_ARRAY_OF(IJK)=IJK
         IF(JP_ARRAY_OF(IJK)==-999) JP_ARRAY_OF(IJK)=IJK
         IF(KM_ARRAY_OF(IJK)==-999) KM_ARRAY_OF(IJK)=IJK
         IF(KP_ARRAY_OF(IJK)==-999) KP_ARRAY_OF(IJK)=IJK

         IF(WEST_ARRAY_OF(IJK)==-999)   WEST_ARRAY_OF(IJK)=IJK
         IF(EAST_ARRAY_OF(IJK)==-999)   EAST_ARRAY_OF(IJK)=IJK
         IF(SOUTH_ARRAY_OF(IJK)==-999)  SOUTH_ARRAY_OF(IJK)=IJK
         IF(NORTH_ARRAY_OF(IJK)==-999)  NORTH_ARRAY_OF(IJK)=IJK
         IF(BOTTOM_ARRAY_OF(IJK)==-999) BOTTOM_ARRAY_OF(IJK)=IJK
         IF(TOP_ARRAY_OF(IJK)==-999)    TOP_ARRAY_OF(IJK)=IJK


! Try to avoid pointing to a cell out of bound

         IF(IM_ARRAY_OF(IJK)<IJKSTART3) IM_ARRAY_OF(IJK)=IJK
         IF(IP_ARRAY_OF(IJK)<IJKSTART3) IP_ARRAY_OF(IJK)=IJK
         IF(JM_ARRAY_OF(IJK)<IJKSTART3) JM_ARRAY_OF(IJK)=IJK
         IF(JP_ARRAY_OF(IJK)<IJKSTART3) JP_ARRAY_OF(IJK)=IJK
         IF(KM_ARRAY_OF(IJK)<IJKSTART3) KM_ARRAY_OF(IJK)=IJK
         IF(KP_ARRAY_OF(IJK)<IJKSTART3) KP_ARRAY_OF(IJK)=IJK

         IF(WEST_ARRAY_OF(IJK)<IJKSTART3)   WEST_ARRAY_OF(IJK)=IJK
         IF(EAST_ARRAY_OF(IJK)<IJKSTART3)   EAST_ARRAY_OF(IJK)=IJK
         IF(SOUTH_ARRAY_OF(IJK)<IJKSTART3)  SOUTH_ARRAY_OF(IJK)=IJK
         IF(NORTH_ARRAY_OF(IJK)<IJKSTART3)  NORTH_ARRAY_OF(IJK)=IJK
         IF(BOTTOM_ARRAY_OF(IJK)<IJKSTART3) BOTTOM_ARRAY_OF(IJK)=IJK
         IF(TOP_ARRAY_OF(IJK)<IJKSTART3)    TOP_ARRAY_OF(IJK)=IJK


         IF(IM_ARRAY_OF(IJK)>IJKEND3) IM_ARRAY_OF(IJK)=IJK
         IF(IP_ARRAY_OF(IJK)>IJKEND3) IP_ARRAY_OF(IJK)=IJK
         IF(JM_ARRAY_OF(IJK)>IJKEND3) JM_ARRAY_OF(IJK)=IJK
         IF(JP_ARRAY_OF(IJK)>IJKEND3) JP_ARRAY_OF(IJK)=IJK
         IF(KM_ARRAY_OF(IJK)>IJKEND3) KM_ARRAY_OF(IJK)=IJK
         IF(KP_ARRAY_OF(IJK)>IJKEND3) KP_ARRAY_OF(IJK)=IJK

         IF(WEST_ARRAY_OF(IJK)>IJKEND3)   WEST_ARRAY_OF(IJK)=IJK
         IF(EAST_ARRAY_OF(IJK)>IJKEND3)   EAST_ARRAY_OF(IJK)=IJK
         IF(SOUTH_ARRAY_OF(IJK)>IJKEND3)  SOUTH_ARRAY_OF(IJK)=IJK
         IF(NORTH_ARRAY_OF(IJK)>IJKEND3)  NORTH_ARRAY_OF(IJK)=IJK
         IF(BOTTOM_ARRAY_OF(IJK)>IJKEND3) BOTTOM_ARRAY_OF(IJK)=IJK
         IF(TOP_ARRAY_OF(IJK)>IJKEND3)    TOP_ARRAY_OF(IJK)=IJK



      ENDDO



      IF(.NOT.ADJUST_PROC_DOMAIN_SIZE) THEN
         if(.not.allocated(NCPP_UNIFORM)) allocate( NCPP_UNIFORM(0:NumPEs-1))
         if(.not.allocated(NCPP_UNIFORM_ALL)) allocate( NCPP_UNIFORM_ALL(0:NumPEs-1))
         NCPP_UNIFORM(MyPE) = BACKGROUND_IJKEND3
      ENDIF

#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

      CALL ALLGATHER_1i (NCPP_UNIFORM(MyPE),NCPP_UNIFORM_ALL,IERR)

      CALL ALLGATHER_1i (BACKGROUND_IJKEND3,BACKGROUND_IJKEND3_ALL,IERR)

!      WRITE(*,100),'ON MyPE = ', MyPE, ' , &
!                    THE NUMBER OF ACTIVE CELLS WENT FROM ',BACKGROUND_IJKEND3, ' TO ', IJKEND3 , &
!                    ' (', DBLE(IJKEND3-BACKGROUND_IJKEND3)/DBLE(BACKGROUND_IJKEND3)*100.0D0, ' % DIFFERENCE)'

!      print*,'From set increment: MyPE,NCCP_UNIFORM=',MyPE,NCPP_UNIFORM

!       WRITE(*,*) 'set increment:',MyPE,NCPP_UNIFORM(MyPE),BACKGROUND_IJKEND3,IJKEND3

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Shifting arrays...'

      CALL SHIFT_DP_ARRAY(EP_G)
      CALL SHIFT_DP_ARRAY(EP_GO)
      CALL SHIFT_DP_ARRAY(P_G)
      CALL SHIFT_DP_ARRAY(P_GO)
      CALL SHIFT_DP_ARRAY(RO_G)
      CALL SHIFT_DP_ARRAY(RO_GO)
      CALL SHIFT_DP_ARRAY(ROP_G)
      CALL SHIFT_DP_ARRAY(ROP_GO)
      CALL SHIFT_DP_ARRAY(T_G)
      CALL SHIFT_DP_ARRAY(T_GO)
      CALL SHIFT_DP_ARRAY(GAMA_RG)
      CALL SHIFT_DP_ARRAY(T_RG)
      DO NN = 1, NMAX(0)
         CALL SHIFT_DP_ARRAY(X_G(:,NN))
         CALL SHIFT_DP_ARRAY(X_GO(:,NN))
         CALL SHIFT_DP_ARRAY(DIF_G(:,NN))
      ENDDO
      CALL SHIFT_DP_ARRAY(U_G)
      CALL SHIFT_DP_ARRAY(U_GO)
      CALL SHIFT_DP_ARRAY(V_G)
      CALL SHIFT_DP_ARRAY(V_GO)
      CALL SHIFT_DP_ARRAY(W_G)
      CALL SHIFT_DP_ARRAY(W_GO)
      CALL SHIFT_DP_ARRAY(P_s_v)
      CALL SHIFT_DP_ARRAY(P_s_f)
      CALL SHIFT_DP_ARRAY(P_s_p)
      CALL SHIFT_DP_ARRAY(P_star)
      CALL SHIFT_DP_ARRAY(P_staro)

      CALL SHIFT_DP_ARRAY(MU_G)
      CALL SHIFT_DP_ARRAY(MW_MIX_G)

      CALL SHIFT_DP_ARRAY(MU_GT)
      CALL SHIFT_DP_ARRAY(LAMBDA_GT)

      CALL SHIFT_DP_ARRAY(C_pg)
      CALL SHIFT_DP_ARRAY(K_g)

      DO M = 1, MMAX
         CALL SHIFT_DP_ARRAY(RO_S(:,M))
         CALL SHIFT_DP_ARRAY(RO_SO(:,M))
         CALL SHIFT_DP_ARRAY(ROP_S(:,M))
         CALL SHIFT_DP_ARRAY(ROP_SO(:,M))
         CALL SHIFT_DP_ARRAY(D_P(:,M))
         CALL SHIFT_DP_ARRAY(D_PO(:,M))
         CALL SHIFT_DP_ARRAY(T_S(:,M))
         CALL SHIFT_DP_ARRAY(T_SO(:,M))
         CALL SHIFT_DP_ARRAY(GAMA_RS(:,M))
         CALL SHIFT_DP_ARRAY(T_RS(:,M))
         DO NN = 1, NMAX(M)
            CALL SHIFT_DP_ARRAY(X_S(:,M,NN))
            CALL SHIFT_DP_ARRAY(X_SO(:,M,NN))
            CALL SHIFT_DP_ARRAY(DIF_S(:,M,NN))
         ENDDO
         CALL SHIFT_DP_ARRAY(U_S(:,M))
         CALL SHIFT_DP_ARRAY(U_SO(:,M))
         CALL SHIFT_DP_ARRAY(V_S(:,M))
         CALL SHIFT_DP_ARRAY(V_SO(:,M))
         CALL SHIFT_DP_ARRAY(W_S(:,M))
         CALL SHIFT_DP_ARRAY(W_SO(:,M))
         CALL SHIFT_DP_ARRAY(P_s(:,M))
         CALL SHIFT_DP_ARRAY(P_s_c(:,M))
         CALL SHIFT_DP_ARRAY(THETA_M(:,M))
         CALL SHIFT_DP_ARRAY(THETA_MO(:,M))
         CALL SHIFT_DP_ARRAY(C_ps(:,M))
         CALL SHIFT_DP_ARRAY(K_s(:,M))
      ENDDO



      DO NN=1,Nscalar
         CALL SHIFT_DP_ARRAY(SCALAR(:,NN))
         CALL SHIFT_DP_ARRAY(SCALARO(:,NN))
      ENDDO

      IF(K_Epsilon) THEN
         CALL SHIFT_DP_ARRAY(K_TURB_G)
         CALL SHIFT_DP_ARRAY(E_TURB_G)
         CALL SHIFT_DP_ARRAY(K_TURB_GO)
         CALL SHIFT_DP_ARRAY(E_TURB_GO)
      ENDIF


      CALL SHIFT_DP_ARRAY(VOL)
      CALL SHIFT_DP_ARRAY(VOL_U)
      CALL SHIFT_DP_ARRAY(VOL_V)
      CALL SHIFT_DP_ARRAY(VOL_W)

      CALL SHIFT_DP_ARRAY(AXY)
      CALL SHIFT_DP_ARRAY(AXY_U)
      CALL SHIFT_DP_ARRAY(AXY_V)
      CALL SHIFT_DP_ARRAY(AXY_W)

      CALL SHIFT_DP_ARRAY(AYZ)
      CALL SHIFT_DP_ARRAY(AYZ_U)
      CALL SHIFT_DP_ARRAY(AYZ_V)
      CALL SHIFT_DP_ARRAY(AYZ_W)

      CALL SHIFT_DP_ARRAY(AXZ)
      CALL SHIFT_DP_ARRAY(AXZ_U)
      CALL SHIFT_DP_ARRAY(AXZ_V)
      CALL SHIFT_DP_ARRAY(AXZ_W)

      CALL SHIFT_DP_ARRAY(X_U)
      CALL SHIFT_DP_ARRAY(Y_U)
      CALL SHIFT_DP_ARRAY(Z_U)

      CALL SHIFT_DP_ARRAY(X_V)
      CALL SHIFT_DP_ARRAY(Y_V)
      CALL SHIFT_DP_ARRAY(Z_V)

      CALL SHIFT_DP_ARRAY(X_W)
      CALL SHIFT_DP_ARRAY(Y_W)
      CALL SHIFT_DP_ARRAY(Z_W)

      CALL SHIFT_DP_ARRAY(NORMAL_S(:,1))
      CALL SHIFT_DP_ARRAY(NORMAL_S(:,2))
      CALL SHIFT_DP_ARRAY(NORMAL_S(:,3))

      CALL SHIFT_DP_ARRAY(REFP_S(:,1))
      CALL SHIFT_DP_ARRAY(REFP_S(:,2))
      CALL SHIFT_DP_ARRAY(REFP_S(:,3))


      CALL SHIFT_DP_ARRAY(AREA_CUT)
      CALL SHIFT_DP_ARRAY(AREA_U_CUT)
      CALL SHIFT_DP_ARRAY(AREA_V_CUT)
      CALL SHIFT_DP_ARRAY(AREA_W_CUT)

      CALL SHIFT_DP_ARRAY(DELX_Ue)
      CALL SHIFT_DP_ARRAY(DELX_Uw)
      CALL SHIFT_DP_ARRAY(DELY_Un)
      CALL SHIFT_DP_ARRAY(DELY_Us)
      CALL SHIFT_DP_ARRAY(DELZ_Ut)
      CALL SHIFT_DP_ARRAY(DELZ_Ub)

      CALL SHIFT_DP_ARRAY(DELX_Ve)
      CALL SHIFT_DP_ARRAY(DELX_Vw)
      CALL SHIFT_DP_ARRAY(DELY_Vn)
      CALL SHIFT_DP_ARRAY(DELY_Vs)
      CALL SHIFT_DP_ARRAY(DELZ_Vt)
      CALL SHIFT_DP_ARRAY(DELZ_Vb)

      CALL SHIFT_DP_ARRAY(DELX_We)
      CALL SHIFT_DP_ARRAY(DELX_Ww)
      CALL SHIFT_DP_ARRAY(DELY_Wn)
      CALL SHIFT_DP_ARRAY(DELY_Ws)
      CALL SHIFT_DP_ARRAY(DELZ_Wt)
      CALL SHIFT_DP_ARRAY(DELZ_Wb)

      CALL SHIFT_DP_ARRAY(X_U_ec)
      CALL SHIFT_DP_ARRAY(Y_U_ec)
      CALL SHIFT_DP_ARRAY(Z_U_ec)

      CALL SHIFT_DP_ARRAY(X_U_nc)
      CALL SHIFT_DP_ARRAY(Y_U_nc)
      CALL SHIFT_DP_ARRAY(Z_U_nc)

      CALL SHIFT_DP_ARRAY(X_U_tc)
      CALL SHIFT_DP_ARRAY(Y_U_tc)
      CALL SHIFT_DP_ARRAY(Z_U_tc)

      CALL SHIFT_DP_ARRAY(X_V_ec)
      CALL SHIFT_DP_ARRAY(Y_V_ec)
      CALL SHIFT_DP_ARRAY(Z_V_ec)

      CALL SHIFT_DP_ARRAY(X_V_nc)
      CALL SHIFT_DP_ARRAY(Y_V_nc)
      CALL SHIFT_DP_ARRAY(Z_V_nc)

      CALL SHIFT_DP_ARRAY(X_V_tc)
      CALL SHIFT_DP_ARRAY(Y_V_tc)
      CALL SHIFT_DP_ARRAY(Z_V_tc)

      CALL SHIFT_DP_ARRAY(X_W_ec)
      CALL SHIFT_DP_ARRAY(Y_W_ec)
      CALL SHIFT_DP_ARRAY(Z_W_ec)

      CALL SHIFT_DP_ARRAY(X_W_nc)
      CALL SHIFT_DP_ARRAY(Y_W_nc)
      CALL SHIFT_DP_ARRAY(Z_W_nc)

      CALL SHIFT_DP_ARRAY(X_W_tc)
      CALL SHIFT_DP_ARRAY(Y_W_tc)
      CALL SHIFT_DP_ARRAY(Z_W_tc)


      CALL SHIFT_DP_ARRAY(DELH_Scalar)
      CALL SHIFT_DP_ARRAY(DWALL)


      CALL SHIFT_DP_ARRAY(DELH_U)
      CALL SHIFT_DP_ARRAY(NORMAL_U(:,1))
      CALL SHIFT_DP_ARRAY(NORMAL_U(:,2))
      CALL SHIFT_DP_ARRAY(NORMAL_U(:,3))
      CALL SHIFT_DP_ARRAY(REFP_U(:,1))
      CALL SHIFT_DP_ARRAY(REFP_U(:,2))
      CALL SHIFT_DP_ARRAY(REFP_U(:,3))

      CALL SHIFT_DP_ARRAY(Theta_Ue)
      CALL SHIFT_DP_ARRAY(Theta_Ue_bar)
      CALL SHIFT_DP_ARRAY(Theta_U_ne)
      CALL SHIFT_DP_ARRAY(Theta_U_nw)
      CALL SHIFT_DP_ARRAY(Theta_U_te)
      CALL SHIFT_DP_ARRAY(Theta_U_tw)
      CALL SHIFT_DP_ARRAY(Alpha_Ue_c)
      CALL SHIFT_DP_ARRAY(NOC_U_E)
      CALL SHIFT_DP_ARRAY(Theta_Un)
      CALL SHIFT_DP_ARRAY(Theta_Un_bar)
      CALL SHIFT_DP_ARRAY(Alpha_Un_c)
      CALL SHIFT_DP_ARRAY(NOC_U_N)
      CALL SHIFT_DP_ARRAY(Theta_Ut)
      CALL SHIFT_DP_ARRAY(Theta_Ut_bar)
      CALL SHIFT_DP_ARRAY(Alpha_Ut_c)
      CALL SHIFT_DP_ARRAY(NOC_U_T)
      CALL SHIFT_DP_ARRAY(A_UPG_E)
      CALL SHIFT_DP_ARRAY(A_UPG_W)


      CALL SHIFT_DP_ARRAY(DELH_V)
      CALL SHIFT_DP_ARRAY(NORMAL_V(:,1))
      CALL SHIFT_DP_ARRAY(NORMAL_V(:,2))
      CALL SHIFT_DP_ARRAY(NORMAL_V(:,3))
      CALL SHIFT_DP_ARRAY(REFP_V(:,1))
      CALL SHIFT_DP_ARRAY(REFP_V(:,2))
      CALL SHIFT_DP_ARRAY(REFP_V(:,3))

      CALL SHIFT_DP_ARRAY(Theta_V_ne)
      CALL SHIFT_DP_ARRAY(Theta_V_se)
      CALL SHIFT_DP_ARRAY(Theta_Vn)
      CALL SHIFT_DP_ARRAY(Theta_Vn_bar)
      CALL SHIFT_DP_ARRAY(Theta_V_nt)
      CALL SHIFT_DP_ARRAY(Theta_V_st)
      CALL SHIFT_DP_ARRAY(Theta_Ve)
      CALL SHIFT_DP_ARRAY(Theta_Ve_bar)
      CALL SHIFT_DP_ARRAY(Alpha_Ve_c)
      CALL SHIFT_DP_ARRAY(NOC_V_E)
      CALL SHIFT_DP_ARRAY(Alpha_Vn_c)
      CALL SHIFT_DP_ARRAY(NOC_V_N)
      CALL SHIFT_DP_ARRAY(Theta_Vt)
      CALL SHIFT_DP_ARRAY(Theta_Vt_bar)
      CALL SHIFT_DP_ARRAY(Alpha_Vt_c)
      CALL SHIFT_DP_ARRAY(NOC_V_T)
      CALL SHIFT_DP_ARRAY(A_VPG_N)
      CALL SHIFT_DP_ARRAY(A_VPG_S)


      CALL SHIFT_DP_ARRAY(DELH_W)
      CALL SHIFT_DP_ARRAY(NORMAL_W(:,1))
      CALL SHIFT_DP_ARRAY(NORMAL_W(:,2))
      CALL SHIFT_DP_ARRAY(NORMAL_W(:,3))
      CALL SHIFT_DP_ARRAY(REFP_W(:,1))
      CALL SHIFT_DP_ARRAY(REFP_W(:,2))
      CALL SHIFT_DP_ARRAY(REFP_W(:,3))

      CALL SHIFT_DP_ARRAY(Theta_W_te)
      CALL SHIFT_DP_ARRAY(Theta_W_be)
      CALL SHIFT_DP_ARRAY(Theta_W_tn)
      CALL SHIFT_DP_ARRAY(Theta_W_bn)
      CALL SHIFT_DP_ARRAY(Theta_Wt)
      CALL SHIFT_DP_ARRAY(Theta_Wt_bar)
      CALL SHIFT_DP_ARRAY(Theta_We)
      CALL SHIFT_DP_ARRAY(Theta_We_bar)
      CALL SHIFT_DP_ARRAY(Alpha_We_c)
      CALL SHIFT_DP_ARRAY(NOC_W_E)
      CALL SHIFT_DP_ARRAY(Theta_Wn)
      CALL SHIFT_DP_ARRAY(Theta_Wn_bar)
      CALL SHIFT_DP_ARRAY(Alpha_Wn_c)
      CALL SHIFT_DP_ARRAY(NOC_W_N)
      CALL SHIFT_DP_ARRAY(Alpha_Wt_c)
      CALL SHIFT_DP_ARRAY(NOC_W_T)
      CALL SHIFT_DP_ARRAY(A_WPG_T)
      CALL SHIFT_DP_ARRAY(A_WPG_B)


      CALL SHIFT_DP_ARRAY(ONEoDX_E_U)
      CALL SHIFT_DP_ARRAY(ONEoDY_N_U)
      CALL SHIFT_DP_ARRAY(ONEoDZ_T_U)

      CALL SHIFT_DP_ARRAY(ONEoDX_E_V)
      CALL SHIFT_DP_ARRAY(ONEoDY_N_V)
      CALL SHIFT_DP_ARRAY(ONEoDZ_T_V)

      CALL SHIFT_DP_ARRAY(ONEoDX_E_W)
      CALL SHIFT_DP_ARRAY(ONEoDY_N_W)
      CALL SHIFT_DP_ARRAY(ONEoDZ_T_W)



      CALL SHIFT_INT_ARRAY(FLAG,UNDEFINED_I)
      CALL SHIFT_INT_ARRAY(FLAG_E,UNDEFINED_I)
      CALL SHIFT_INT_ARRAY(FLAG_N,UNDEFINED_I)
      CALL SHIFT_INT_ARRAY(FLAG_T,UNDEFINED_I)



      CALL SHIFT_LOG_ARRAY(INTERIOR_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(SMALL_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(BLOCKED_CELL_AT,.TRUE.)
      CALL SHIFT_LOG_ARRAY(BLOCKED_U_CELL_AT,.TRUE.)
      CALL SHIFT_LOG_ARRAY(BLOCKED_V_CELL_AT,.TRUE.)
      CALL SHIFT_LOG_ARRAY(BLOCKED_W_CELL_AT,.TRUE.)
      CALL SHIFT_LOG_ARRAY(STANDARD_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(STANDARD_U_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(STANDARD_V_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(STANDARD_W_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_U_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_V_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_W_CELL_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_TREATMENT_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_U_TREATMENT_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_V_TREATMENT_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(CUT_W_TREATMENT_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(WALL_U_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(WALL_V_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(WALL_W_AT,.FALSE.)
      CALL SHIFT_LOG_ARRAY(SNAP,.FALSE.)


      CALL SHIFT_INT_ARRAY(U_MASTER_OF,0)
      CALL SHIFT_INT_ARRAY(V_MASTER_OF,0)
      CALL SHIFT_INT_ARRAY(W_MASTER_OF,0)

      CALL SHIFT_INT_ARRAY(BC_ID,0)
      CALL SHIFT_INT_ARRAY(BC_U_ID,0)
      CALL SHIFT_INT_ARRAY(BC_V_ID,0)
      CALL SHIFT_INT_ARRAY(BC_W_ID,0)

      CALL SHIFT_INT_ARRAY(CELL_CLASS,0)

      CALL SHIFT_INT_ARRAY(PHASE_4_P_G,UNDEFINED_I)
      CALL SHIFT_INT_ARRAY(PHASE_4_P_S,UNDEFINED_I)


      CALL SHIFT_LOG_ARRAY(SCALAR_NODE_ATWALL,.FALSE.)

      CALL SHIFT_DP_ARRAY(SCALAR_NODE_XYZ)
      CALL SHIFT_DP_ARRAY(SCALAR_NODE_XYZ)
      CALL SHIFT_DP_ARRAY(SCALAR_NODE_XYZ)
      CALL SHIFT_DP_ARRAY(Ovol_around_node)



      IF (IJK_P_G /= UNDEFINED_I) IJK_P_G = IJK_OF_BACKGROUND(IJK_P_G)

      IF(STIFF_CHEMISTRY) CALL SHIFT_LOG_ARRAY(notOwner,.FALSE.)

      IF(BDIST_IO) CALL SHIFT_CONNECTIVITY_FOR_BDIST_IO


!=====================================================================
! JFD: Re-assign send and receive arrays
!=====================================================================

      IS_SERIAL = numPEs==1

      IF(.NOT.IS_SERIAL) THEN

      IF(MINIMIZE_SEND_RECV) THEN


      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Minimizing send and receive arrays for Send layer 1...'

! Send

! Layer 1

      nullify(new_xsend1, new_sendtag1, new_sendproc1, new_sendijk1)

! Get array size

      n_total = 0
      DO send_pos = lbound(sendijk1,1),ubound(sendijk1,1)
         IJK = sendijk1( send_pos )
         IF(IJK_OF_BACKGROUND(IJK)/=-999)  n_total = n_total + 1  ! count active cells
      ENDDO


      allocate( new_sendijk1( max(1,n_total) ) )
      allocate( new_xsend1(nsend1+1) )
      allocate( new_sendtag1(nsend1+1) )
      allocate( new_sendproc1(nsend1+1) )

! Fill in arrays

      new_xsend1 = 0
      send_pos = 0
      placeholder = 1
      new_nsend1 = 0

      do nn = 1,nsend1
         j1       = xsend1(nn)
         j2       = xsend1(nn+1)-1
         sendsize = j2-j1+1

         new_send_size(nn) = 0

         DO jj=j1,j2
            ijk = sendijk1( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN      ! Only keep active cells
               new_send_size(nn) = new_send_size(nn) + 1
               send_pos = send_pos + 1
               new_sendijk1(send_pos) = IJK_OF_BACKGROUND(IJK)
            ENDIF
         ENDDO

         IF(new_send_size(nn)>0) THEN
            new_nsend1 = new_nsend1 + 1
            new_xsend1(new_nsend1) = placeholder
            placeholder = placeholder + new_send_size(nn)

            new_sendtag1(new_nsend1) = sendtag1(nn)
            new_sendproc1(new_nsend1) = sendproc1(nn)

            nj1 = new_xsend1(new_nsend1)
            nj2 = nj1 + new_send_size(nn) - 1

            CALL BUBBLE_SORT_1D_INT_ARRAY(new_sendijk1(nj1:nj2),nj1,nj2)
         ENDIF

      ENDDO

      new_xsend1(new_nsend1+1)= nj2 + 1


      nsend1 = new_nsend1
      sendtag1 => new_sendtag1
      sendproc1 => new_sendproc1
      xsend1 => new_xsend1
      sendijk1 => new_sendijk1

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning send and receive arrays for Send layer 2...'

! Layer 2

      nullify(new_xsend2, new_sendtag2, new_sendproc2, new_sendijk2)

! Get array size

      n_total = 0
      DO send_pos = lbound(sendijk2,1),ubound(sendijk2,1)
         IJK = sendijk2( send_pos )
         IF(IJK_OF_BACKGROUND(IJK)/=-999)  n_total = n_total + 1  ! count active cells
      ENDDO

      allocate( new_sendijk2( max(1,n_total) ) )
      allocate( new_xsend2(nsend2+1) )
      allocate( new_sendtag2(nsend2+1) )
      allocate( new_sendproc2(nsend2+1) )

! Fill in arrays

      new_xsend2 = 0
      send_pos = 0
      placeholder = 1
      new_nsend2 = 0



      do nn = 1,nsend2
         j1       = xsend2(nn)
         j2       = xsend2(nn+1)-1
         sendsize = j2-j1+1

         new_send_size(nn) = 0

         DO jj=j1,j2
            ijk = sendijk2( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN      ! Only keep active cells
               new_send_size(nn) = new_send_size(nn) + 1
               send_pos = send_pos + 1
               new_sendijk2(send_pos) = IJK_OF_BACKGROUND(IJK)
            ENDIF
         ENDDO

         IF(new_send_size(nn)>0) THEN
            new_nsend2 = new_nsend2 + 1
            new_xsend2(new_nsend2) = placeholder
            placeholder = placeholder + new_send_size(nn)

            new_sendtag2(new_nsend2) = sendtag2(nn)
            new_sendproc2(new_nsend2) = sendproc2(nn)

            nj1 = new_xsend2(new_nsend2)
            nj2 = nj1 + new_send_size(nn) - 1

!            if (MyPE==6) print*, 'n,new_nsend2,nj1,nj2',n,new_nsend2,nj1,nj2

            CALL BUBBLE_SORT_1D_INT_ARRAY(new_sendijk2(nj1:nj2),nj1,nj2)
         ENDIF

      ENDDO

      new_xsend2(new_nsend2+1)= nj2 + 1

!      print*, 'MyPE, Laxt value of xsend2=',MyPE,new_nsend2,new_xsend2(new_nsend2+1)

      nsend2 = new_nsend2
      sendtag2 => new_sendtag2
      sendproc2 => new_sendproc2
      xsend2 => new_xsend2
      sendijk2 => new_sendijk2

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning send and receive arrays for Receive layer 1...'

! Receive

! Layer 1

      nullify(new_xrecv1, new_recvtag1, new_recvproc1, new_recvijk1)

! Get array size

      n_total = 0
      DO recv_pos = lbound(recvijk1,1),ubound(recvijk1,1)
         IJK = recvijk1( recv_pos )
         IF(IJK_OF_BACKGROUND(IJK)/=-999)  n_total = n_total + 1  ! count active cells
      ENDDO

      allocate( new_recvijk1( max(1,n_total) ) )
      allocate( new_xrecv1(nrecv1+1) )
      allocate( new_recvtag1(nrecv1+1) )
      allocate( new_recvproc1(nrecv1+1) )

! Fill in arrays

      new_xrecv1 = 0
      recv_pos = 0




      new_xrecv1 = 0
      recv_pos = 0
      placeholder = 1
      new_nrecv1 = 0

      do nn = 1,nrecv1
         j1       = xrecv1(nn)
         j2       = xrecv1(nn+1)-1
         recvsize = j2-j1+1

         new_recv_size(nn) = 0

         DO jj=j1,j2
            ijk = recvijk1( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN      ! Only keep active cells
               new_recv_size(nn) = new_recv_size(nn) + 1
               recv_pos = recv_pos + 1
               new_recvijk1(recv_pos) = IJK_OF_BACKGROUND(IJK)
            ENDIF
         ENDDO

         IF(new_recv_size(nn)>0) THEN
            new_nrecv1 = new_nrecv1 + 1
            new_xrecv1(new_nrecv1) = placeholder
            placeholder = placeholder + new_recv_size(nn)

            new_recvtag1(new_nrecv1) = recvtag1(nn)
            new_recvproc1(new_nrecv1) = recvproc1(nn)

            nj1 = new_xrecv1(new_nrecv1)
            nj2 = nj1 + new_recv_size(nn) - 1

            CALL BUBBLE_SORT_1D_INT_ARRAY(new_recvijk1(nj1:nj2),nj1,nj2)
         ENDIF

      ENDDO

      new_xrecv1(new_nrecv1+1)=nj2 + 1

      nrecv1 = new_nrecv1
      recvtag1 => new_recvtag1
      recvproc1 => new_recvproc1
      xrecv1 => new_xrecv1
      recvijk1 => new_recvijk1


      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning send and receive arrays for Receive layer 2...'
! Layer 2

      nullify(new_xrecv2, new_recvtag2, new_recvproc2, new_recvijk2)

! Get array size

      n_total = 0
      DO recv_pos = lbound(recvijk2,1),ubound(recvijk2,1)
         IJK = recvijk2( recv_pos )
         IF(IJK_OF_BACKGROUND(IJK)/=-999)  n_total = n_total + 1  ! count active cells
      ENDDO

      allocate( new_recvijk2( max(1,n_total) ) )
      allocate( new_xrecv2(nrecv2+1) )
      allocate( new_recvtag2(nrecv2+1) )
      allocate( new_recvproc2(nrecv2+1) )

! Fill in arrays

      new_xrecv2 = 0
      recv_pos = 0
      placeholder = 1
      new_nrecv2 = 0

      do nn = 1,nrecv2
         j1       = xrecv2(nn)
         j2       = xrecv2(nn+1)-1
         recvsize = j2-j1+1

         new_recv_size(nn) = 0

         DO jj=j1,j2
            ijk = recvijk2( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN      ! Only keep active cells
               new_recv_size(nn) = new_recv_size(nn) + 1
               recv_pos = recv_pos + 1
               new_recvijk2(recv_pos) = IJK_OF_BACKGROUND(IJK)
            ENDIF
         ENDDO

         IF(new_recv_size(nn)>0) THEN
            new_nrecv2 = new_nrecv2 + 1
            new_xrecv2(new_nrecv2) = placeholder
            placeholder = placeholder + new_recv_size(nn)

            new_recvtag2(new_nrecv2) = recvtag2(nn)
            new_recvproc2(new_nrecv2) = recvproc2(nn)

            nj1 = new_xrecv2(new_nrecv2)
            nj2 = nj1 + new_recv_size(nn) - 1

            CALL BUBBLE_SORT_1D_INT_ARRAY(new_recvijk2(nj1:nj2),nj1,nj2)
         ENDIF

      ENDDO


      new_xrecv2(new_nrecv2+1)=nj2 + 1

      nrecv2 = new_nrecv2
      recvtag2 => new_recvtag2
      recvproc2 => new_recvproc2
      xrecv2 => new_xrecv2
      recvijk2 => new_recvijk2



      ELSE   ! Only update IJK values

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning send and receive arrays for Send layer 1...'



! Send

! Layer 1

      nullify(new_sendijk1)

      print *, 'sendijk1=',size(sendijk1)
      allocate( new_sendijk1( size(sendijk1) ) )

! Fill in arrays


      do nn = 1,nsend1
         j1       = xsend1(nn)
         j2       = xsend1(nn+1)-1
         sendsize = j2-j1+1


         DO jj=j1,j2
            ijk = sendijk1( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN
               new_sendijk1(jj) = IJK_OF_BACKGROUND(IJK)
            ELSE
               new_sendijk1(jj) = sendijk1( j1 )
            ENDIF
         ENDDO


      ENDDO

      sendijk1 => new_sendijk1

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning send and receive arrays for Send layer 2...'

! Layer 2


      print *, 'sendijk2=',size(sendijk2)
      nullify(new_sendijk2)


      allocate( new_sendijk2( size(sendijk2) ) )

! Fill in arrays


      do nn = 1,nsend2
         j1       = xsend2(nn)
         j2       = xsend2(nn+1)-1
         sendsize = j2-j1+1

         new_send_size(nn) = 0

         DO jj=j1,j2
            ijk = sendijk2( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN      ! Only keep active cells
               new_sendijk2(jj) = IJK_OF_BACKGROUND(IJK)
            ELSE
               new_sendijk2(jj) = sendijk2( j1 )
            ENDIF
         ENDDO


      ENDDO

      sendijk2 => new_sendijk2

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning send and receive arrays for Receive layer 1...'

! Receive

! Layer 1

      print *, 'recvijk1=',size(recvijk1)
      nullify(new_recvijk1)


      allocate( new_recvijk1( size(recvijk1) ) )

! Fill in arrays


      do nn = 1,nrecv1
         j1       = xrecv1(nn)
         j2       = xrecv1(nn+1)-1
         recvsize = j2-j1+1

         new_recv_size(nn) = 0

         DO jj=j1,j2
            ijk = recvijk1( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN      ! Only keep active cells
               new_recvijk1(jj) = IJK_OF_BACKGROUND(IJK)
            ELSE
               new_recvijk1(jj) = recvijk1( j1 )
            ENDIF
         ENDDO


      ENDDO

      recvijk1 => new_recvijk1


      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning send and receive arrays for Receive layer 2...'
! Layer 2

      print *, 'secvijk2=',size(recvijk2)
      nullify(new_recvijk2)


      allocate( new_recvijk2( size(recvijk2) ) )


      do nn = 1,nrecv2
         j1       = xrecv2(nn)
         j2       = xrecv2(nn+1)-1
         recvsize = j2-j1+1

         new_recv_size(nn) = 0

         DO jj=j1,j2
            ijk = recvijk2( jj )

            IF(IJK_OF_BACKGROUND(IJK)/=-999) THEN      ! Only keep active cells
               new_recvijk2(jj) = IJK_OF_BACKGROUND(IJK)
            ELSE
               new_recvijk2(jj) = recvijk2( j1 )
            ENDIF
         ENDDO


      ENDDO


      recvijk2 => new_recvijk2

      ENDIF

#ifdef MPI
      comm = MPI_COMM_WORLD
#endif

!  INSERT NEW SEND_RECV INIT HERE

   call sendrecv_re_init_after_re_indexing(comm, 0 )

   ENDIF ! IS_SERIAL

#ifdef MPI
   call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

!goto 999

!======================================================================
!   Re-assign cell classes
!======================================================================

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Re-assigning cell classes...'

!      print*, 'before class reassignment:, iclass =',iclass


      ICLASS = 0
!
!     Loop over all cells (minus the ghost layers)
      DO K = KSTART3, KEND3
         DO J = JSTART3, JEND3
            L100: DO I = ISTART3, IEND3
               IJK = FUNIJK(I,J,K)               !Find value of IJK
!
               IF(DEAD_CELL_AT(I,J,K)) CYCLE

!          Find the the effective cell-center indices for all neighbor cells
!               CALL SET_INDEX1A (I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, &
!                  IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT)

               IJKN = NORTH_ARRAY_OF(IJK)
               IJKS = SOUTH_ARRAY_OF(IJK)
               IJKE = EAST_ARRAY_OF(IJK)
               IJKW = WEST_ARRAY_OF(IJK)
               IJKT = TOP_ARRAY_OF(IJK)
               IJKB = BOTTOM_ARRAY_OF(IJK)

               IMJK = IM_ARRAY_OF(IJK)
               IPJK = IP_ARRAY_OF(IJK)
               IJMK = JM_ARRAY_OF(IJK)
               IJPK = JP_ARRAY_OF(IJK)
               IJKM = KM_ARRAY_OF(IJK)
               IJKP = KP_ARRAY_OF(IJK)

!
               ICLASS = ICLASS + 1               !Increment the ICLASS counter
               IF (ICLASS > MAX_CLASS) THEN
                  IF(DMP_LOG)WRITE (UNIT_LOG, 2000) MAX_CLASS
                  CALL MFIX_EXIT(myPE)
               ENDIF
               INCREMENT_FOR_N(ICLASS) = IJKN - IJK
               INCREMENT_FOR_S(ICLASS) = IJKS - IJK
               INCREMENT_FOR_E(ICLASS) = IJKE - IJK
               INCREMENT_FOR_W(ICLASS) = IJKW - IJK
               INCREMENT_FOR_T(ICLASS) = IJKT - IJK
               INCREMENT_FOR_B(ICLASS) = IJKB - IJK
               INCREMENT_FOR_IM(ICLASS) = IMJK - IJK
               INCREMENT_FOR_IP(ICLASS) = IPJK - IJK
               INCREMENT_FOR_JM(ICLASS) = IJMK - IJK
               INCREMENT_FOR_JP(ICLASS) = IJPK - IJK
               INCREMENT_FOR_KM(ICLASS) = IJKM - IJK
               INCREMENT_FOR_KP(ICLASS) = IJKP - IJK


               INCREMENT_FOR_NB(1,ICLASS) = INCREMENT_FOR_E(ICLASS)
               INCREMENT_FOR_NB(2,ICLASS) = INCREMENT_FOR_W(ICLASS)
               INCREMENT_FOR_NB(3,ICLASS) = INCREMENT_FOR_S(ICLASS)
               INCREMENT_FOR_NB(4,ICLASS) = INCREMENT_FOR_N(ICLASS)
               INCREMENT_FOR_NB(5,ICLASS) = INCREMENT_FOR_B(ICLASS)
               INCREMENT_FOR_NB(6,ICLASS) = INCREMENT_FOR_T(ICLASS)


               INCREMENT_FOR_MP(1,ICLASS) = INCREMENT_FOR_IM(ICLASS)
               INCREMENT_FOR_MP(2,ICLASS) = INCREMENT_FOR_IP(ICLASS)
               INCREMENT_FOR_MP(3,ICLASS) = INCREMENT_FOR_JM(ICLASS)
               INCREMENT_FOR_MP(4,ICLASS) = INCREMENT_FOR_JP(ICLASS)
               INCREMENT_FOR_MP(5,ICLASS) = INCREMENT_FOR_KM(ICLASS)
               INCREMENT_FOR_MP(6,ICLASS) = INCREMENT_FOR_KP(ICLASS)

               DENOTE_CLASS(ICLASS) = INCREMENT_FOR_N(ICLASS) + INCREMENT_FOR_S&
                  (ICLASS) + INCREMENT_FOR_E(ICLASS) + INCREMENT_FOR_W(ICLASS)&
                   + INCREMENT_FOR_T(ICLASS) + INCREMENT_FOR_B(ICLASS) + &
                  INCREMENT_FOR_IM(ICLASS) + INCREMENT_FOR_IP(ICLASS) + &
                  INCREMENT_FOR_JM(ICLASS) + INCREMENT_FOR_JP(ICLASS) + &
                  INCREMENT_FOR_KM(ICLASS) + INCREMENT_FOR_KP(ICLASS)

               CELL_CLASS(IJK) = ICLASS

!          Place the cell in a class based on its DENOTE_CLASS(ICLASS) value
               DO IC = 1, ICLASS - 1             !Loop over previous and present classes
!                                                !IF a possible match in cell types
                  IF (DENOTE_CLASS(ICLASS) == DENOTE_CLASS(IC)) THEN
!                                                !is found, compare all increments
                     IF (INCREMENT_FOR_N(ICLASS) /= INCREMENT_FOR_N(IC)) CYCLE
                     IF (INCREMENT_FOR_S(ICLASS) /= INCREMENT_FOR_S(IC)) CYCLE
                     IF (INCREMENT_FOR_E(ICLASS) /= INCREMENT_FOR_E(IC)) CYCLE
                     IF (INCREMENT_FOR_W(ICLASS) /= INCREMENT_FOR_W(IC)) CYCLE
                     IF (INCREMENT_FOR_T(ICLASS) /= INCREMENT_FOR_T(IC)) CYCLE
                     IF (INCREMENT_FOR_B(ICLASS) /= INCREMENT_FOR_B(IC)) CYCLE
                     IF (INCREMENT_FOR_IM(ICLASS) /= INCREMENT_FOR_IM(IC)) &
                        CYCLE
                     IF (INCREMENT_FOR_IP(ICLASS) /= INCREMENT_FOR_IP(IC)) &
                        CYCLE
                     IF (INCREMENT_FOR_JM(ICLASS) /= INCREMENT_FOR_JM(IC)) &
                        CYCLE
                     IF (INCREMENT_FOR_JP(ICLASS) /= INCREMENT_FOR_JP(IC)) &
                        CYCLE
                     IF (INCREMENT_FOR_KM(ICLASS) /= INCREMENT_FOR_KM(IC)) &
                        CYCLE
                     IF (INCREMENT_FOR_KP(ICLASS) /= INCREMENT_FOR_KP(IC)) &
                        CYCLE
                     CELL_CLASS(IJK) = IC        !Assign cell to a class
                     ICLASS = ICLASS - 1
                     CYCLE  L100                 !Go to next cell
                  ENDIF
               END DO
            END DO L100
         END DO
      END DO

      IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: New number of classes = ', ICLASS

    CALL WRITE_IJK_VALUES

!      IJKEND3 = BACKGROUND_IJKEND3  ! for debugging purpose, will need to be removed

!      RETURN

      ALLOCATE( NEW_IJKSIZE3_ALL(0:NUMPES-1) )

      CALL ALLGATHER_1I (IJKEND3,NEW_IJKSIZE3_ALL,IERR)

!      print*,'MyPE, NEW_IJKSIZE3_ALL=',MyPE,NEW_IJKSIZE3_ALL

      IF(NUMPES.GT.1) THEN

#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

         IF(MyPE.EQ.0) THEN
            WRITE(*,1000)"============================================================================="
            WRITE(*,1000)"    PROCESSOR    I-SIZE     J-SIZE     K-SIZE    # CELLS    # CELLS   DIFF."
            WRITE(*,1000)"                                                 (BCKGRD) (RE-INDEXED) (%)"
            WRITE(*,1000)"============================================================================="
         ENDIF

#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

         DO IPROC = 0,NumPes-1
            IF(MyPE==IPROC) THEN
               I_SIZE = IEND1 - ISTART1 + 1
               J_SIZE = JEND1 - JSTART1 + 1
               K_SIZE = KEND1 - KSTART1 + 1
               DIFF_NCPP(IPROC) = DBLE(NEW_IJKSIZE3_ALL(IPROC)-NCPP_UNIFORM_ALL(IPROC))/DBLE(NCPP_UNIFORM_ALL(IPROC))*100.0D0
               WRITE(*,1060) IPROC,I_SIZE,J_SIZE,K_SIZE,BACKGROUND_IJKEND3_ALL(IPROC),NEW_IJKSIZE3_ALL(IPROC),DIFF_NCPP(IPROC)
            ENDIF
#ifdef MPI
            call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif
         ENDDO

         IF(MyPE.EQ.0) THEN
            WRITE(*,1000)"============================================================================="
            WRITE(*,1070)'MAX # OF CELLS (BACKGRD)    = ',MAXVAL(NCPP_UNIFORM_ALL),'     AT PROCESSOR: ',MAXLOC(NCPP_UNIFORM_ALL)-1
            WRITE(*,1070)'MAX # OF CELLS (RE-INDEXED) = ',MAXVAL(NEW_IJKSIZE3_ALL),'     AT PROCESSOR: ',MAXLOC(NEW_IJKSIZE3_ALL)-1
            WRITE(*,1080)'DIFFERENCE (%)              = ', &
                 DBLE(MAXVAL(NEW_IJKSIZE3_ALL)-MAXVAL(NCPP_UNIFORM_ALL))/DBLE(MAXVAL(NCPP_UNIFORM_ALL))*100.0
           WRITE(*,1000)"============================================================================="
         ENDIF
#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif
      ENDIF

1000  FORMAT(1x,A)
1060  FORMAT(1X,6(I10,1X),F8.1)
1070  FORMAT(1X,A,I8,A,I8)
1080  FORMAT(1X,A,F8.1)
!
!     WRITE FOLLOWING IF THERE IS AN ERROR IN MODULE
2000 FORMAT(/70('*')//'From: SET_INCREMENTS'/'Message: The number of',&
         'classes has exceeded the maximum allowed (',I8,').  Increase',&
         'MAX_CLASS in PARAM1.INC')
!

      END SUBROUTINE RE_INDEX_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RECORD_NEW_IJK_CELL                                    C
!  Purpose: Records indices for new IJK cell                           C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-MAY-11  C
!  Reviewer:                                          Date: ##-###-##  C
!                                                                      C
!  Revision Number: #                                                  C
!  Purpose: ##########                                                 C
!  Author:  ##########                                Date: ##-###-##  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE RECORD_NEW_IJK_CELL(I,J,K,IJK,NEW_IJK,TEMP_IJK_ARRAY_OF,TEMP_I_OF,TEMP_J_OF,TEMP_K_OF)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits
      USE scalars
      USE run

      USE cutcell

      USE sendrecv

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
!                      Indices
      INTEGER ::        I, J, K, IJK, NEW_IJK
      INTEGER, DIMENSION(ISTART3-1:IEND3+1,JSTART3-1:JEND3+1,KSTART3-1:KEND3+1) :: TEMP_IJK_ARRAY_OF
      INTEGER, DIMENSION(DIMENSION_3)     :: TEMP_I_OF,TEMP_J_OF,TEMP_K_OF



      BACKGROUND_IJK_OF(NEW_IJK) = IJK

      IJK_OF_BACKGROUND(IJK) = NEW_IJK

      TEMP_IJK_ARRAY_OF(I,J,K)=NEW_IJK

      TEMP_I_OF(NEW_IJK) = I
      TEMP_J_OF(NEW_IJK) = J
      TEMP_K_OF(NEW_IJK) = K

      NEW_IJK = NEW_IJK + 1


      RETURN

      END SUBROUTINE RECORD_NEW_IJK_CELL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUBBLE_SORT_1D_INT_ARRAY                               C
!  Purpose: Bubble sort a section of a 1D integer array in ascending   C
!           order. The section that is sorted out is from I1 to I2     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-MAY-11  C
!  Reviewer:                                          Date: ##-###-##  C
!                                                                      C
!  Revision Number: #                                                  C
!  Purpose: ##########                                                 C
!  Author:  ##########                                Date: ##-###-##  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!

      SUBROUTINE BUBBLE_SORT_1D_INT_ARRAY(ARRAY,I1,I2)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE indices
      USE geometry
      USE compar
      USE cutcell

      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER ::I1,I2,BUFFER,I,J
!
      INTEGER, DIMENSION(I1:I2) :: ARRAY

!-----------------------------------------------

!======================================================================
!   Bubble sort a section of a 1D integer array in ascending order
!   The section that is sorted out is from I1 to I2
!======================================================================

!     print*,'Before Bubble sorting from MyPE=',MyPE, I1,I2,ARRAY


      DO I = I1,I2-1
         DO J = I2-1,I,-1
            IF(ARRAY(J)>ARRAY(J+1)) THEN
               BUFFER     = ARRAY(J)
               ARRAY(J)   = ARRAY(J+1)
               ARRAY(J+1) = BUFFER
            ENDIF
         ENDDO
      ENDDO


!     print*,'After Bubble sorting from MyPE=',MyPE, I1,I2,ARRAY


      END SUBROUTINE BUBBLE_SORT_1D_INT_ARRAY



      SUBROUTINE SHIFT_DP_ARRAY(ARRAY)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE cutcell
      USE functions
      USE geometry
      USE indices
      USE param, only: dimension_3
      USE param1, only: undefined

      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER ::IJK
!
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: ARRAY, BUFFER

!======================================================================
!   To remove dead cells, the number of useful cells was calculated in
!   RE_INDEX_ARRAY, and is stored back in IJKEND3
!   Now, the array is shifted such that all useful values are contiguous
!   and are located between IJKSTART3 and IJKEND3
!   The array BACKGROUND_IJK_OF(IJK) points to the original cell
!======================================================================

      BUFFER = ARRAY
      ARRAY = UNDEFINED

      DO IJK = IJKSTART3, IJKEND3

          ARRAY(IJK) = BUFFER(BACKGROUND_IJK_OF(IJK))

      ENDDO


      END SUBROUTINE SHIFT_DP_ARRAY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SHIFT_INT_ARRAY                                         C
!  Purpose: Shifts an Integer array to new IJK range                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-MAY-11  C
!  Reviewer:                                          Date: ##-###-##  C
!                                                                      C
!  Revision Number: #                                                  C
!  Purpose: ##########                                                 C
!  Author:  ##########                                Date: ##-###-##  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SHIFT_INT_ARRAY(ARRAY,DEFAULT_VALUE)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE indices
      USE geometry
      USE compar
      USE cutcell
      USE functions
      USE param, only: dimension_3

      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER ::IJK
!
      INTEGER, DIMENSION(DIMENSION_3) :: ARRAY, BUFFER
      INTEGER :: DEFAULT_VALUE

!======================================================================
!   To remove dead cells, the number of useful cells was calculated in
!   RE_INDEX_ARRAY, and is stored back in IJKEND3
!   Now, the array is shifted such that all useful values are contiguous
!   and are located between IJKSTART3 and IJKEND3
!   The array BACKGROUND_IJK_OF(IJK) points to the original cell
!======================================================================

      BUFFER = ARRAY
      ARRAY = DEFAULT_VALUE

      DO IJK = IJKSTART3, IJKEND3

          ARRAY(IJK) = BUFFER(BACKGROUND_IJK_OF(IJK))

      ENDDO


      END SUBROUTINE SHIFT_INT_ARRAY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SHIFT_LOG_ARRAY                                         C
!  Purpose: Shifts an Integer array to new IJK range                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-MAY-11  C
!  Reviewer:                                          Date: ##-###-##  C
!                                                                      C
!  Revision Number: #                                                  C
!  Purpose: ##########                                                 C
!  Author:  ##########                                Date: ##-###-##  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SHIFT_LOG_ARRAY(ARRAY,DEFAULT_VALUE)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE indices
      USE geometry
      USE compar
      USE cutcell
      USE functions
      USE param, only: dimension_3

      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER ::IJK
!
      LOGICAL, DIMENSION(DIMENSION_3) :: ARRAY, BUFFER
      LOGICAL :: DEFAULT_VALUE

!======================================================================
!   To remove dead cells, the number of useful cells was calculated in
!   RE_INDEX_ARRAY, and is stored back in IJKEND3
!   Now, the array is shifted such that all useful values are contiguous
!   and are located between IJKSTART3 and IJKEND3
!   The array BACKGROUND_IJK_OF(IJK) points to the original cell
!======================================================================

      BUFFER = ARRAY
      ARRAY = DEFAULT_VALUE

      DO IJK = IJKSTART3, IJKEND3

          ARRAY(IJK) = BUFFER(BACKGROUND_IJK_OF(IJK))

      ENDDO


      END SUBROUTINE SHIFT_LOG_ARRAY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UNSHIFT_DP_ARRAY                                       C
!  Purpose: Reverts a shifted Double precision array to                C
!  original (background) IJK range                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-MAY-11  C
!  Reviewer:                                          Date: ##-###-##  C
!                                                                      C
!  Revision Number: #                                                  C
!  Purpose: ##########                                                 C
!  Author:  ##########                                Date: ##-###-##  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE UNSHIFT_DP_ARRAY(ARRAY_1,ARRAY_2)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE indices
      USE geometry
      USE compar
      USE cutcell
      USE functions
      USE param, only: dimension_3
      USE param1, only: undefined

      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER ::IJK

      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: ARRAY_1, ARRAY_2

!======================================================================

      ARRAY_2 = UNDEFINED

      DO IJK = IJKSTART3,IJKEND3

          ARRAY_2(BACKGROUND_IJK_OF(IJK)) = ARRAY_1(IJK)

      ENDDO


      END SUBROUTINE UNSHIFT_DP_ARRAY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SHIFT_CONNECTIVITY_FOR_BDIST_IO                        C
!  Purpose: Shifts connectivity for distributed IO                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 04-MAY-11  C
!  Reviewer:                                          Date: ##-###-##  C
!                                                                      C
!  Revision Number: #                                                  C
!  Purpose: ##########                                                 C
!  Author:  ##########                                Date: ##-###-##  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SHIFT_CONNECTIVITY_FOR_BDIST_IO
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE indices
      USE geometry
      USE compar
      USE cutcell
      USE functions
      USE param, only: dimension_3

      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER ::IJK, L, BCK_IJK, NN , CONN
!

      INTEGER, DIMENSION(DIMENSION_3,15) ::TEMP_CONNECTIVITY

!======================================================================
!
!   The array BACKGROUND_IJK_OF(IJK) points to the original cell
!======================================================================

      CALL SHIFT_INT_ARRAY(NUMBER_OF_NODES,0)

      TEMP_CONNECTIVITY = CONNECTIVITY

      DO IJK = 1,IJKEND3
         IF (INTERIOR_CELL_AT(IJK))      THEN
            IF (.NOT.BLOCKED_CELL_AT(IJK)) THEN

               BCK_IJK = BACKGROUND_IJK_OF(IJK)   ! Get the original IJK

               NN = NUMBER_OF_NODES(IJK)          ! Get the number of nodes for this cell, this was already shifted above

               DO L = 1, NN                       ! Loop through the connectivity list
                                                  ! and reassign each point in the list

                  CONN = TEMP_CONNECTIVITY(BCK_IJK,L)

                  IF(CONN>BACKGROUND_IJKEND3) THEN
                     CONNECTIVITY(IJK,L) = CONN - BACKGROUND_IJKEND3 + IJKEND3  ! shift new point ID
                  ELSE
                     CONNECTIVITY(IJK,L) = IJK_OF_BACKGROUND(CONN)              ! Points to the new IJK value
                  ENDIF

               ENDDO

            ENDIF
         ENDIF
      END DO



      END SUBROUTINE SHIFT_CONNECTIVITY_FOR_BDIST_IO




      SUBROUTINE WRITE_INT_TABLE(FILE_UNIT,ARRAY, ARRAY_SIZE, LSTART, LEND, NCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!

!                      FILE UNIT
      INTEGER  ::        FILE_UNIT



!                      Starting array index
      INTEGER  ::        ARRAY_SIZE


!                      Starting array index
      INTEGER  ::        LSTART
!
!                      Ending array index
      INTEGER ::         LEND
!//EFD Nov/11, avoid use of (*)
!//      DOUBLE PRECISION ARRAY(*)
      INTEGER :: ARRAY(ARRAY_SIZE)
!
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                      Number of columns in the table.  When this is changed
!                      remember to change the FORMAT statement also.
!

      INTEGER :: NCOL
!

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
!                      Number of rows
      INTEGER          NROW
!
!
!                      Local array indices
      INTEGER          L, L1, L2, L3
!-----------------------------------------------
!
      NROW = (LEND - LSTART + 1)/NCOL
!
      L2 = LSTART - 1
      DO L = 1, NROW
         L1 = L2 + 1
         L2 = L1 + NCOL - 1
         WRITE (FILE_UNIT, 1020) (ARRAY(L3),L3=L1,L2)
      END DO
      IF (NROW*NCOL < LEND - LSTART + 1) THEN
         L1 = L2 + 1
         L2 = LEND
         WRITE (FILE_UNIT, 1020) (ARRAY(L3),L3=L1,L2)
      ENDIF
      RETURN
!
 1020 FORMAT(14X,50(I12,1X))
      END SUBROUTINE WRITE_INT_TABLE







      SUBROUTINE WRITE_IJK_VALUES

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits
      USE scalars
      USE run
      USE visc_g

      USE cutcell

      USE sendrecv

      USE mpi_utility
      USE parallel

      USE cdist
      USE functions
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
!                      Indices
      INTEGER          I, J, K, IJK

      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: TEMP_IJK_ARRAY_OF
      INTEGER :: IJK_FILE_UNIT
      CHARACTER(LEN=255) :: IJK_FILE_NAME
      CHARACTER(LEN=6)  :: CHAR_MyPE

      allocate(TEMP_IJK_ARRAY_OF(ISTART3-1:IEND3+1,JSTART3-1:JEND3+1,KSTART3-1:KEND3+1))
      TEMP_IJK_ARRAY_OF = IJK_ARRAY_OF

!======================================================================
!   Write IJK value in files (for debugging or info)
!======================================================================

!      IF(NO_K) THEN

         IF(MyPE == PE_IO) WRITE(*,*)' Re-indexing: Writing IJK value in files...'

          IJK_FILE_UNIT = 1000 + MyPE
          WRITE(CHAR_MyPE,'(I6)')MyPE

          IJK_FILE_NAME = 'IJK_INFO_'//CHAR_MyPE//'.txt'

          DO I=1,LEN(TRIM(IJK_FILE_NAME))
             IF(IJK_FILE_NAME(I:I)==' ') IJK_FILE_NAME(I:I)='0'
          ENDDO

          OPEN(CONVERT='BIG_ENDIAN',UNIT=IJK_FILE_UNIT,FILE=IJK_FILE_NAME)


         WRITE(IJK_FILE_UNIT,200)'          MyPE = ',MyPE
         WRITE(IJK_FILE_UNIT,200)' ISTART1,IEND1 = ',ISTART1,IEND1
         WRITE(IJK_FILE_UNIT,200)' JSTART1,JEND1 = ',JSTART1,JEND1
         WRITE(IJK_FILE_UNIT,200)' KSTART1,KEND1 = ',KSTART1,KEND1
         WRITE(IJK_FILE_UNIT,200)' I-SIZE = ',IEND1-ISTART1+1
         WRITE(IJK_FILE_UNIT,200)' J-SIZE = ',JEND1-JSTART1+1
         WRITE(IJK_FILE_UNIT,200)' K-SIZE = ',KEND1-KSTART1+1
         WRITE(IJK_FILE_UNIT,200)' IJKSTART3 = ',IJKSTART3
         WRITE(IJK_FILE_UNIT,200)' IJKEND3 = ',IJKEND3
         WRITE(IJK_FILE_UNIT,*)''

         IF(RE_INDEXING) WRITE(IJK_FILE_UNIT,100) 'INFO: AFTER RE-INDEXING CELLS ON MyPE = ', MyPE, ' , &
                    &THE NUMBER OF ACTIVE CELLS WENT FROM ',BACKGROUND_IJKEND3, ' TO ', IJKEND3 , &
                    ' (', DBLE(IJKEND3-BACKGROUND_IJKEND3)/DBLE(BACKGROUND_IJKEND3)*100.0D0, ' % DIFFERENCE)'

         WRITE(IJK_FILE_UNIT,*)''

      IF(NO_K) THEN
          WRITE(IJK_FILE_UNIT,210) ('======',I=ISTART3,IEND3)
          K=1
          DO J=JEND3,JSTART3,-1
             DO I=ISTART3,IEND3
                IJK = funijk(I,J,K)
!                TEMP_IJK_ARRAY_OF(I,J,K) = cell_class(IJK)
                IF(DEAD_CELL_AT(I,J,K)) TEMP_IJK_ARRAY_OF(I,J,K) = 0
             ENDDO
                IF(RE_INDEXING) THEN
                   WRITE(IJK_FILE_UNIT,230) J,(TEMP_IJK_ARRAY_OF(I,J,K),I=ISTART3,IEND3)
                ELSE
                   WRITE(IJK_FILE_UNIT,230) J,(FUNIJK(I,J,K),I=ISTART3,IEND3)
                ENDIF

          ENDDO

          WRITE(IJK_FILE_UNIT,210) ('======',I=ISTART3,IEND3)
          WRITE(IJK_FILE_UNIT,220) (I,I=ISTART3,IEND3)

       ELSE
          DO IJK=IJKSTART3,IJKEND3
             WRITE(IJK_FILE_UNIT,*) IJK,I_OF(IJK),J_OF(IJK),K_OF(IJK)
          ENDDO

       ENDIF


100       FORMAT(1X,A,I6,A,I8,A,I8,A,F6.1,A)
200       FORMAT(1x,A30,2(I8))

210       FORMAT(8x,50(A))
220       FORMAT(1x,' J/I | ',50(I6))
230       FORMAT(1x,I4,' | ',50(I6))

         IF(.NOT.IS_SERIAL) THEN

         WRITE(IJK_FILE_UNIT,*)''

         WRITE(IJK_FILE_UNIT,*)' Layer    = ',1
         WRITE(IJK_FILE_UNIT,*)' nsend1    = ', nsend1
         WRITE(IJK_FILE_UNIT,*)' sendproc1 = ', sendproc1(1:nsend1)
         WRITE(IJK_FILE_UNIT,*)' sendtag1  = ', sendtag1(1:nsend1)
         WRITE(IJK_FILE_UNIT,*)' xsend1    = ', xsend1(1:nsend1)
         WRITE(IJK_FILE_UNIT,*)' size      = ', size(sendijk1)
         WRITE(IJK_FILE_UNIT,*)' sendijk1  = '
         CALL WRITE_INT_TABLE(IJK_FILE_UNIT,sendijk1, size(sendijk1), 1, size(sendijk1),5)
         WRITE(IJK_FILE_UNIT,*)''

         WRITE(IJK_FILE_UNIT,*)' nrecv1    = ', nrecv1
         WRITE(IJK_FILE_UNIT,*)' recvproc1 = ', recvproc1(1:nrecv1)
         WRITE(IJK_FILE_UNIT,*)' recvtag1  = ', recvtag1(1:nrecv1)
         WRITE(IJK_FILE_UNIT,*)' xrecv1    = ', xrecv1(1:nrecv1)
         WRITE(IJK_FILE_UNIT,*)' size      = ', size(recvijk1)
         WRITE(IJK_FILE_UNIT,*)' recvijk1  = '
         CALL WRITE_INT_TABLE(IJK_FILE_UNIT,recvijk1, size(recvijk1), 1, size(recvijk1), 5)
         WRITE(IJK_FILE_UNIT,*)''
         WRITE(IJK_FILE_UNIT,*)''

         WRITE(IJK_FILE_UNIT,*)' Layer    = ',2
         WRITE(IJK_FILE_UNIT,*)' nsend2    = ', nsend2
         WRITE(IJK_FILE_UNIT,*)' sendproc2 = ', sendproc2(1:nsend2)
         WRITE(IJK_FILE_UNIT,*)' sendtag2  = ', sendtag2(1:nsend2)
         WRITE(IJK_FILE_UNIT,*)' xsend2    = ', xsend2(1:nsend2)
         WRITE(IJK_FILE_UNIT,*)' size      = ', size(sendijk2)
         WRITE(IJK_FILE_UNIT,*)' sendijk2  = '
         CALL WRITE_INT_TABLE(IJK_FILE_UNIT,sendijk2, size(sendijk2), 1, size(sendijk2),5)
         WRITE(IJK_FILE_UNIT,*)''

         WRITE(IJK_FILE_UNIT,*)' nrecv2    = ', nrecv2
         WRITE(IJK_FILE_UNIT,*)' recvproc2 = ', recvproc2(1:nrecv2)
         WRITE(IJK_FILE_UNIT,*)' recvtag2  = ', recvtag2(1:nrecv2)
         WRITE(IJK_FILE_UNIT,*)' xrecv2    = ', xrecv2(1:nrecv2)
         WRITE(IJK_FILE_UNIT,*)' size      = ', size(recvijk2)
         WRITE(IJK_FILE_UNIT,*)' recvijk2  = '
         CALL WRITE_INT_TABLE(IJK_FILE_UNIT,recvijk2, size(recvijk2), 1, size(recvijk2), 5)
         WRITE(IJK_FILE_UNIT,*)''

         ENDIF

         CLOSE(IJK_FILE_UNIT)

!      ENDIF

#ifdef MPI
         call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

!=====================================================================
! JFD: End of Print send info
!=====================================================================
      END SUBROUTINE WRITE_IJK_VALUES
