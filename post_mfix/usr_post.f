!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR_POST                                               C
!                                                                      C
!  Purpose: Do user defined calculations from TIME_START to TIME_LAST  C
!                                                                      C
!  Author: M. Syamlal                               Date: 26-OCT-93    C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: V_g, EP_g, IMIN1, IMAX1, JMIN1, JMAX1         C
!                        KMIN1, KMAX1, IJKMAX2                         C
!  Variables modified: PLOT_TYPE, VAR_INDEX, LOC_X, LOC_Y, LOC_Z       C
!                      I, J, K, IJK                                    C
!                                                                      C
!  Local variables: FILE_NAME, NX, NY, NZ    L, NT                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

!     subroutine for getting average species mass fractions excluding bubbles
      SUBROUTINE USR_POST

      Use compar
      Use fldvar
      Use functions
      Use geometry
      Use indices
      Use param
      Use param1
      Use physprop
      Use post3d
      Use run

      IMPLICIT NONE
!
      REAL              TIME_START , TIME_REAL(N_SPX), TIME_FOUND
      REAL              TIME_LAST, TIME_NOW
      CHARACTER(LEN=60) :: FILE_NAME
      INTEGER           NX , NY , NZ, NSTEP_1
      INTEGER           REC_POINTER(N_SPX) , L
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)

      REAL              SUMI, SUM(JMAX2)
!
      INTEGER           NT, NI
      INTEGER           I, J, IJK
      INTEGER, EXTERNAL :: FUNIJK_LOC
!
      WRITE (*,*) ' Enter start-time and end-time'
      READ  (*,*) TIME_START, TIME_LAST
!
!
      CALL GET_FILE_NAME(FILE_NAME)
      OPEN (UNIT=40,FILE=FILE_NAME,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
!
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
!
!   Enable the required files
!
      READ_SPX(1) = .TRUE.    ! EP_g
!      READ_SPX(2) = .TRUE.    ! P_g, P_star
!      READ_SPX(3) = .TRUE.    ! U_g, V_g, W_g
!      READ_SPX(4) = .TRUE.    ! U_s, V_s, W_s
!      READ_SPX(5) = .TRUE.    ! ROP_s
!      READ_SPX(6) = .TRUE.    ! T_g, T_s1, T_s2
       READ_SPX(7) = .TRUE.    ! X_g, X_s
!      READ_SPX(8) = .TRUE.    ! Theta
!      READ_SPX(9) = .TRUE.    ! User scalar
      CALL SEEK_TIME(READ_SPX, TIME_START, REC_POINTER, TIME_FOUND)
      IF(TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        RETURN
      ENDIF
!
      NT = 0
      DO J = 1, JMAX2
        SUM(J) = ZERO
      ENDDO

!
!  Time loop -- Start
!
100   CONTINUE
      CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,&
                          TIME_NOW, TIME_REAL, NSTEP_1)
!
      IF (TIME_NOW .LT. ZERO .OR. TIME_NOW .GT. TIME_LAST) GOTO 500
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      NT = NT + 1
!
!  DO required computations
!
      DO J = 2, JMAX1

        NI = 0
        SUMI = ZERO

        DO I = 2, IMAX1
          IJK = funijk_LOC(I, J, 1)
          IF(FLUID_AT(IJK)) THEN
            IF(EP_g(IJK) < 0.5)THEN
              NI = NI + 1
              SUMI = SUMI + X_g(IJK, 1)
            ENDIF
          ENDIF
        ENDDO

        IF(NI /= 0) &
          SUM(J) = SUM(J) + SUMI/REAL(NI)
      ENDDO
!
      GOTO 100
!
!  Time loop -- End
!
500   IF (NT.EQ.0) THEN
         WRITE (*,*) ' No times found in common'
         RETURN
      END IF
!
!  Do final data processing printing
!
      DO J = 2, JMAX1
        write(40, '(G12.5, 2x, G12.5)') YDIST_SC(J), SUM(J)/REAL(NT)
      ENDDO
!
      CLOSE (UNIT=40)
!
      RETURN
    END SUBROUTINE USR_POST
!
!  The following routines are not active.  To make a routine active replace the
!  above routine with the desired routine and change its name to USR_POST
!

!     subroutine for cluster size statistics
      SUBROUTINE USR_POST1

      Use fldvar, only: d_p, ep_g, u_s, v_s, w_s, u_g, v_g, w_g
      Use functions, only: fluid_at
      Use geometry, ONLY: ijkmax2
      Use param1, only: n_spx, one, zero
      Use physprop, only: mu_g0, ro_g0

      IMPLICIT NONE
      INTEGER MAX_COUNT
      PARAMETER (MAX_COUNT=1000)
!
      REAL a1, a2, a3, Re_c, EP_c
!      PARAMETER (a1 = 250.)
      PARAMETER (a1 = 1500.)
      PARAMETER (a2 = 0.005)
      PARAMETER (a3 = 90.0)
      PARAMETER (Re_c = 5.)
      PARAMETER (EP_c = 0.92)
!
      REAL              TIME_START , TIME_REAL(N_SPX), TIME_FOUND
      REAL              TIME_LAST, TIME_NOW, SUM, FC_DIST, Re, VREL, FC
      CHARACTER(LEN=60) :: FILE_NAME
      INTEGER           NX , NY , NZ, NSTEP_1
      INTEGER           REC_POINTER(N_SPX) , L , NT
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
!
      INTEGER           FC_COUNT(MAX_COUNT)
      INTEGER           IJK
!
      WRITE (*,*) ' Enter start-time and end-time'
      READ  (*,*) TIME_START, TIME_LAST
!
      MU_g0 = 1.8e-4
      RO_g0 = 1.8e-3
      DO 10 L = 1, MAX_COUNT
        FC_COUNT(L) = 0
10    CONTINUE
!
      CALL GET_FILE_NAME(FILE_NAME)
      OPEN (UNIT=40,FILE=FILE_NAME,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
!
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
!
!   Enable the required files
!
      READ_SPX(1) = .TRUE.    ! EP_g
!      READ_SPX(2) = .TRUE.    ! P_g, P_star
      READ_SPX(3) = .TRUE.    ! U_g, V_g, W_g
      READ_SPX(4) = .TRUE.    ! U_s, V_s, W_s
!      READ_SPX(5) = .TRUE.    ! ROP_s
!      READ_SPX(6) = .TRUE.    ! T_g, T_s1, T_s2
!      READ_SPX(7) = .TRUE.    ! X_g, X_s
      CALL SEEK_TIME(READ_SPX, TIME_START, REC_POINTER, TIME_FOUND)
      IF(TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        RETURN
      ENDIF
!
      NT = 0
!
!  Time loop -- Start
!
100   CONTINUE
      CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,&
                          TIME_NOW, TIME_REAL, NSTEP_1)
!
      IF (TIME_NOW .LT. ZERO .OR. TIME_NOW .GT. TIME_LAST) GOTO 500
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      NT = NT + 1
!
!  DO required computations
!
      DO 200 IJK = 1,IJKMAX2
        IF(FLUID_AT(IJK)) THEN
          VREL = SQRT(  (U_g(IJK) - U_s(IJK,1))**2&
                      + (V_g(IJK) - V_s(IJK,1))**2&
                      + (W_g(IJK) - W_s(IJK,1))**2 )
          Re =  D_p(IJK, 1) * VREL * RO_g0 / MU_g0
          FC = (ONE + a1 * exp(-a2*(Re - Re_c)**2 &
                                - a3*(EP_g(IJK)-ep_c)**2)&
                     * Re * (1. - EP_g(IJK))               )

          L = MIN(MAX_COUNT, NINT(FC))
          FC_COUNT(L) = FC_COUNT(L) + 1
        ENDIF
200   CONTINUE
!
      GOTO 100
!
!  Time loop -- End
!
500   IF (NT.EQ.0) THEN
         WRITE (*,*) ' No times found in common'
         RETURN
      END IF
!
!  Do final data processing printing
!
      SUM = ZERO
      DO 600 L = 1, MAX_COUNT
        SUM = SUM + FC_COUNT(L)
600   CONTINUE
!
      DO 650 L = 1, MAX_COUNT
        FC_DIST = REAL(FC_COUNT(L))/SUM
        WRITE(40,*)L, FC_DIST
650   CONTINUE
!
      CLOSE (UNIT=40)
!
      RETURN
 END SUBROUTINE USR_POST1
