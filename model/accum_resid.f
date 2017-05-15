!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ACCUM_RESID                                            C
!                                                                      C
!  Purpose: Accumulate all the residuals and calculate max_resid       C
!                                                                      C
!                                                                      C
!  Author: S. Pannala                                 Date: 14-Jun-07  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE ACCUM_RESID
!
!-----------------------------------------------
!     M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE residual
      USE run
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER          M, NN

!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      INTEGER              LOCAL_INDEX

!
      IF(DEBUG_RESID) Return
!
      LOCAL_INDEX = 0

! Pack the numerators and denominators into one vector for performing single global operation

!!!$omp parallel do private( NN,M )&
!!!$omp&  REDUCTION(+:LOCAL_INDEX)
      DO NN = 2, NRESID
         DO M = 0, DIMENSION_M
            LOCAL_INDEX = LOCAL_INDEX + 1
            RESID_PACK(LOCAL_INDEX) = NUM_RESID(NN,M)
            LOCAL_INDEX = LOCAL_INDEX + 1
            RESID_PACK(LOCAL_INDEX) = DEN_RESID(NN,M)
         ENDDO
      ENDDO

      call global_all_sum(RESID_PACK)

! Unpack the numerators and denominators from the global sum vector

      LOCAL_INDEX = 0

!!!$omp parallel do private( NN,M )&
!!!$omp&  REDUCTION(+:LOCAL_INDEX)
      DO NN = 2, NRESID
         DO M = 0, DIMENSION_M
            LOCAL_INDEX = LOCAL_INDEX + 1
            NUM_RESID(NN,M) = RESID_PACK(LOCAL_INDEX)
            LOCAL_INDEX = LOCAL_INDEX + 1
            DEN_RESID(NN,M) = RESID_PACK(LOCAL_INDEX)
         ENDDO
      ENDDO

!!!$omp parallel do private( NN,M )
      DO NN = 2, NRESID
         DO M = 0, DIMENSION_M
            IF (DEN_RESID(NN,M) > ZERO) THEN
               RESID(NN,M) = NUM_RESID(NN,M)/DEN_RESID(NN,M)
            ELSE IF (NUM_RESID(NN,M) == ZERO) THEN
               RESID(NN,M) = ZERO
            ELSE
               RESID(NN,M) = UNDEFINED
!     WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!     CALL WRITE_ERROR ('ACCUM_RESID', LINE, 1)
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
