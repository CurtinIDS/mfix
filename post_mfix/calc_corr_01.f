!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_CORR_01(FINISH)                                   C
!  Purpose: Calculate correlations between different variables in a    C
!           fluidized bed hydrodynamic model                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 01-AUG-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: KMAX2, JMAX2, IMAX2, EP_g, V_g                C
!  Variables modified: STARTED, NSUM, I, J, K, IJK, SUM_EP_g, SUM_V_g  C
!                      SUM_EPxEP_g, SUM_VxV_g, AVG_EP_g, AVG_V_g       C
!                      SDV_EP_g, SDV_V_g                               C
!                                                                      C
!  Local variables: SDV_2                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_CORR_01(FINISH,INTER)
!
!
        Use param
        Use param1
        Use fldvar
        Use physprop
        Use geometry
        Use indices
        Use correl
        Use compar
        Use functions

      IMPLICIT NONE
!
!                      A flag to check whether the subroutine is called for
!                      the first time.
      INTEGER          RAND_NO
      PARAMETER (RAND_NO=946235187)
      INCLUDE 'xforms.inc'
!
!  passed variables
!
!                      Logical variable to tell the subroutine to finish
!                      calculations. i.e. To compute averages and variances
      LOGICAL          FINISH,INTER
!
!
!  local variables
!
!                       temporary stograge for SDV**2
      DOUBLE PRECISION  SDV_2
!
!                       indices
      INTEGER           I, J, K, IJK
!
!  If finish command is not given, do summations for correlations.
!
      IF(.NOT.FINISH)THEN
!
!  First check whether the subroutine is called the first time, by checking
!  whether the flag STARTED is correctly set to a specific number.  If the
!  subroutine is called the first time, do initializations.
!
        IF(STARTED .NE. RAND_NO)THEN

          Allocate(  SUM_EP_g (DIMENSION_3) )
          Allocate(  SUM_EPxEP_g (DIMENSION_3) )
          Allocate(  SUM_V_g (DIMENSION_3) )
          Allocate(  SUM_VxV_g (DIMENSION_3) )

          Allocate(  AVG_EP_g (DIMENSION_3) )
          Allocate(  SDV_EP_g (DIMENSION_3) )
          Allocate(  AVG_V_g  (DIMENSION_3) )
          Allocate(  SDV_V_g  (DIMENSION_3) )

          STARTED = RAND_NO
          NSUM = 0
          SUM_EP_g(:) = ZERO
          SUM_EPxEP_g(:) = ZERO
          SUM_V_g(:) = ZERO
          SUM_VxV_g(:) = ZERO
        ENDIF
!
!       Sum field variable values in all fluid cells
!
        NSUM = NSUM + 1
        DO K = 1, KMAX2
           DO J = 1, JMAX2
              DO I = 1, IMAX2
                 IJK = FUNIJK(I, J, K)
                 IF( FLUID_AT(IJK) ) THEN
                    SUM_EP_g(IJK)    = SUM_EP_g(IJK)    + EP_g(IJK)
                    SUM_EPxEP_g(IJK) = SUM_EPxEP_g(IJK) + EP_g(IJK)**2
                    SUM_V_g(IJK)     = SUM_V_g(IJK)     + V_g(IJK)
                    SUM_VxV_g(IJK)   = SUM_VxV_g(IJK)   + V_g(IJK)**2
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
!
!  If finish command is given, calculate averages, variances,
!  and correlations
!
      ELSE
        IF(NSUM .LE. 0)THEN
          WRITE(*,5000)
          STOP
        ENDIF
        DO K = 1, KMAX2
           DO J = 1, JMAX2
              IF (DO_XFORMS) THEN
                 CALL CHECK_INTER(INTER)
                 IF (INTER) RETURN
              END IF
              DO I = 1, IMAX2
                 IJK = FUNIJK(I, J, K)
                 IF( FLUID_AT(IJK) ) THEN
                    AVG_EP_g(IJK) = SUM_EP_g(IJK) / NSUM
                    SDV_2 = SUM_EPxEP_g(IJK) / NSUM - AVG_EP_g(IJK)**2
                    SDV_2 = MAX(ZERO,SDV_2)
                    SDV_EP_g(IJK) = SQRT( SDV_2 )
                    AVG_V_g(IJK) = SUM_V_g(IJK) / NSUM
                    SDV_2 = SUM_VxV_g(IJK) / NSUM - AVG_V_g(IJK)**2
                    SDV_2 = MAX(ZERO,SDV_2)
                    SDV_V_g(IJK) = SQRT( SDV_2 )
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     RETURN
5000 FORMAT(/1X,70('*')//' From: CALC_CORR',&
          /' Message: Attempting averaging before doing summations',&
          /1X, 70('*')/)
   END SUBROUTINE CALC_CORR_01
