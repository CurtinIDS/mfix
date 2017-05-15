!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_INCREMENTS3                                         !
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
      SUBROUTINE SET_INCREMENTS3

      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager
      use function3
      use functions

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
! Index denoting cell class
      INTEGER :: ICLASS
! Array of sum of increments to make the class determination faster.
      INTEGER :: DENOTE_CLASS(MAX_CLASS)
! Flag for using the 'real' I/J/K value (not cyclic.)
      LOGICAL :: SHIFT
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_INCREMENTS3")

! Initialize the default values to Undefined_I

      IP1_3(:) = UNDEFINED_I
      IM1_3(:) = UNDEFINED_I
      JP1_3(:) = UNDEFINED_I
      JM1_3(:) = UNDEFINED_I
      KP1_3(:) = UNDEFINED_I
      KM1_3(:) = UNDEFINED_I

      DO I = ISTART4, IEND4

         SHIFT = .NOT.(I==IMIN4 .OR. I==IMIN3 .OR. I==IMIN2 .OR. &
                       I==IMAX4 .OR. I==IMAX3 .OR. I==IMAX2)

         IF(CYCLIC_X .AND. NODESI.EQ.1 .AND. DO_I .AND. SHIFT) THEN
            IP1_3(I) = IMAP_C(IMAP_C(I)+1)
            IM1_3(I) = IMAP_C(IMAP_C(I)-1)
         ELSE
            IM1_3(I) = MAX(ISTART4,I - 1)
            IP1_3(I) = MIN(IEND4,I + 1)
         ENDIF
      ENDDO


      DO J = JSTART4, JEND4

         SHIFT = .NOT.(J==JMIN4 .OR. J==JMIN3 .OR. J==JMIN2 .OR. &
                       J==JMAX4 .OR. J==JMAX3 .OR. J==JMAX2)

         IF(CYCLIC_Y .AND. NODESJ.EQ.1 .AND. DO_J .AND. SHIFT) THEN
            JP1_3(J) = JMAP_C(JMAP_C(J)+1)
            JM1_3(J) = JMAP_C(JMAP_C(J)-1)
         ELSE
            JM1_3(J) = MAX(JSTART4,J - 1)
            JP1_3(J) = MIN(JEND4,J + 1)
         ENDIF
      ENDDO


      DO K = KSTART4, KEND4

         SHIFT = .NOT.(K==KMIN4 .OR. K==KMIN3 .OR. K==KMIN2 .OR. &
                       K==KMAX4 .OR. K==KMAX3 .OR. K==KMAX2)

         IF(CYCLIC_Z .AND. NODESK.EQ.1 .AND. DO_K .AND. SHIFT) THEN
            KP1_3(K) = KMAP_C(KMAP_C(K)+1)
            KM1_3(K) = KMAP_C(KMAP_C(K)-1)
         ELSE
            KM1_3(K) = MAX(KSTART4,K - 1)
            KP1_3(K) = MIN(KEND4,K + 1)
         ENDIF
      ENDDO

!     Loop over all cells
      DO K = KSTART4, KEND4
      DO J = JSTART4, JEND4
      DO I = ISTART4, IEND4
         IJK = FUNIJK3(I,J,K)

         I3_OF(IJK) = I
         J3_OF(IJK) = J
         K3_OF(IJK) = K
      ENDDO
      ENDDO
      ENDDO

! Initialize the number of cell classes
      ICLASS = 0

! Loop over all cells (minus the ghost layers)
      DO K = KSTART4, KEND4
      DO J = JSTART4, JEND4
      L100: DO I = ISTART4, IEND4

         IJK = FUNIJK3(I,J,K)

!  Find the the effective cell-center indices for all neighbor cells
         CALL SET_INDEX1A3 (I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM,&
            IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT)

! Increment the ICLASS counter
         ICLASS = ICLASS + 1
         IF(ICLASS > MAX_CLASS) THEN
            WRITE(ERR_MSG, 1200) trim(iVal(MAX_CLASS))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1200 FORMAT('Error 1200: The number of classes has exceeded the ',    &
         'maximum: ',A,/'Increase the MAX_CLASS parameter in param1',  &
         '_mod.f and recompile.')

         INCREMENT3_FOR_IM(ICLASS) = IMJK - IJK
         INCREMENT3_FOR_IP(ICLASS) = IPJK - IJK
         INCREMENT3_FOR_JM(ICLASS) = IJMK - IJK
         INCREMENT3_FOR_JP(ICLASS) = IJPK - IJK
         INCREMENT3_FOR_KM(ICLASS) = IJKM - IJK
         INCREMENT3_FOR_KP(ICLASS) = IJKP - IJK

         DENOTE_CLASS(ICLASS) =  &
            INCREMENT3_FOR_IM(ICLASS) + INCREMENT3_FOR_IP(ICLASS) + &
            INCREMENT3_FOR_JM(ICLASS) + INCREMENT3_FOR_JP(ICLASS) + &
            INCREMENT3_FOR_KM(ICLASS) + INCREMENT3_FOR_KP(ICLASS)

            CELL_CLASS3(IJK) = ICLASS

! Place the cell in a class based on its DENOTE_CLASS(ICLASS) value
! Loop over previous and present classes
         DO IC = 1, ICLASS - 1

            IF (DENOTE_CLASS(ICLASS) == DENOTE_CLASS(IC)) THEN
               IF(INCREMENT3_FOR_IM(ICLASS)/=INCREMENT3_FOR_IM(IC))CYCLE
               IF(INCREMENT3_FOR_IP(ICLASS)/=INCREMENT3_FOR_IP(IC))CYCLE
               IF(INCREMENT3_FOR_JM(ICLASS)/=INCREMENT3_FOR_JM(IC))CYCLE
               IF(INCREMENT3_FOR_JP(ICLASS)/=INCREMENT3_FOR_JP(IC))CYCLE
               IF(INCREMENT3_FOR_KM(ICLASS)/=INCREMENT3_FOR_KM(IC))CYCLE
               IF(INCREMENT3_FOR_KP(ICLASS)/=INCREMENT3_FOR_KP(IC))CYCLE
               CELL_CLASS3(IJK) = IC
               ICLASS = ICLASS - 1
               CYCLE  L100 ! Go to next cell
            ENDIF
         END DO

      END DO L100
      END DO
      END DO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_INCREMENTS3
