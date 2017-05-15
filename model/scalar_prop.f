!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Scalar_PROP(IER)                                       C
!  Purpose: Calculate diffusion coefficeint and sources for user-defined
!           scalars
!                                                                      C
!  Author:                                                    Date:    C
!  Reviewer:                                                  Date:    C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SCALAR_PROP()
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE scalars
      USE toleranc
      USE compar
      USE sendrecv
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
      INTEGER          L,IJK
!
!-----------------------------------------------

      IF(NScalar == 0) RETURN
!
!  ---  Remember to include all the local variables here for parallel
!  ---- processing
!!!$omp  parallel do private(ijk, L)
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
           DO L = 1, NScalar

!            d (Scalar)/dt = S
!            S is linearized as S = Scalar_c - Scalar_p * Scalar
!            Scalar_c and Scalar_p must be >= 0
!            *** Uncomment next two lines ***
              Scalar_c (IJK, L) = ZERO
              Scalar_p (IJK, L) = ZERO
!
!            Diffusion coefficient for User-defined Scalars
!            *** Uncomment next one line ***
              Dif_Scalar(IJK, L) =ZERO
           END DO
!
         ENDIF
      END DO
!\\Sendrecv operations - just to make sure all the variables computed are
!  are passed and updated locally - fool-proof approach - Sreekanth - 102199

!      call send_recv(Scalar_c,2)
!      call send_recv(Scalar_p,2)
!      call send_recv(Dif_Scalar,2)
      RETURN
      END SUBROUTINE SCALAR_PROP
