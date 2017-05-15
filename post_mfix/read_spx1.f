!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_SPX1(READ_SPX1,REC_POINTER,AT_EOF,                C
!                         TIME_REAL, NSTEP_1)                          c
!  Purpose: read in the time-dependent restart records (REAL)          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: TIME, NSTEP, EP_g, RO_g, P_g, P_star, U_g       C
!                        V_g, W_g, U_s, V_s, W_s, ROP_s, T_g, T_s1     C
!                        T_s2, IJKMAX2, MMAX                           C
!                                                                      C
!  Local variables:  LC, NEXT_REC, TIME_REAL                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, &
                          TIME_REAL,NSTEP_1)
!
!
      Use param, only: dimension_3
      Use param1
      Use fldvar
      Use geometry
      Use physprop
      Use run
      Use funits
      Use post3d, only: version_number, spx_open, last_rec
      Use scalars
      Use rxns
      Use machine
      Use indices
      Use compar
      Use functions

      IMPLICIT NONE
!
! passed arguments
!
!             flag whether to read a particular SPx file this time step
      LOGICAL READ_SPX(*)
!
!             pointers to next record to read in each file
      INTEGER REC_POINTER(*)
!
      LOGICAL AT_EOF(*)
      INTEGER NSTEP_1
      REAL    TIME_REAL(*)

!
!                      Dummy variable for reading T_s2
      DOUBLE PRECISION Tmp(DIMENSION_3)
!
! local variables
!
!             loop counters
      INTEGER LC,M,N,IJK
!
!             Pointer to the next record
      INTEGER NEXT_REC , num_recs

      integer :: gas_species_index , solid_species_index , solid_index
      logical :: bRead_all

      DOUBLE PRECISION :: SUM_XoRO

      common /fast_sp7/ gas_species_index , solid_species_index , &
                         solid_index , bRead_all

      num_recs = 1 + ijkmax2 / nwords_r
      if (mod(ijkmax2,nwords_r) .eq. 0) num_recs = num_recs - 1

!
! ".SP1" FILE         EP_g    [ ROP_g , RO_g must be calculated ...
!                                        not written out ]
!
      IF (READ_SPX(1).AND..NOT.AT_EOF(1)) THEN
         IF(.NOT.SPX_OPEN(1)) THEN
           WRITE(*,*)' SP1 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(1)
         IF (NEXT_REC.EQ.LAST_REC(1)) THEN
            AT_EOF(1) = .TRUE.
            RETURN
         END IF
         AT_EOF(1) = .FALSE.
         READ (UNIT_SPX+1,REC=NEXT_REC) TIME_REAL(1) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+1,EP_g,IJKMAX2,NEXT_REC)
         REC_POINTER(1) = NEXT_REC
         NSTEP_1 = NSTEP
      END IF
!
! ".SP2" FILE         P_g , P_star
!
      IF (READ_SPX(2).AND..NOT.AT_EOF(2)) THEN
         IF(.NOT.SPX_OPEN(2)) THEN
           WRITE(*,*)' SP2 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(2)
         IF (NEXT_REC.EQ.LAST_REC(2)) THEN
            AT_EOF(2) = .TRUE.
            RETURN
         END IF
         AT_EOF(2) = .FALSE.
         READ (UNIT_SPX+2,REC=NEXT_REC) TIME_REAL(2) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+2,P_g,IJKMAX2,NEXT_REC)
         CALL IN_BIN_R(UNIT_SPX+2,P_star,IJKMAX2,NEXT_REC)
         REC_POINTER(2) = NEXT_REC
      END IF
!
! ".SP3" FILE         U_g , V_g , W_g
!
      IF (READ_SPX(3).AND..NOT.AT_EOF(3)) THEN
         IF(.NOT.SPX_OPEN(3)) THEN
           WRITE(*,*)' SP3 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(3)
         IF (NEXT_REC.EQ.LAST_REC(3)) THEN
            AT_EOF(3) = .TRUE.
            RETURN
         END IF
         AT_EOF(3) = .FALSE.
         READ (UNIT_SPX+3,REC=NEXT_REC) TIME_REAL(3) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+3,U_g,IJKMAX2,NEXT_REC)
         CALL IN_BIN_R(UNIT_SPX+3,V_g,IJKMAX2,NEXT_REC)
         CALL IN_BIN_R(UNIT_SPX+3,W_g,IJKMAX2,NEXT_REC)
         REC_POINTER(3) = NEXT_REC
      END IF
!
! ".SP4" FILE         U_s , V_s , W_s
!
      IF (READ_SPX(4).AND..NOT.AT_EOF(4)) THEN
         IF(.NOT.SPX_OPEN(4)) THEN
           WRITE(*,*)' SP4 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(4)
         IF (NEXT_REC.EQ.LAST_REC(4)) THEN
            AT_EOF(4) = .TRUE.
            RETURN
         END IF
         AT_EOF(4) = .FALSE.
         READ (UNIT_SPX+4,REC=NEXT_REC) TIME_REAL(4) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 100 LC = 1,MMAX
            CALL IN_BIN_R(UNIT_SPX+4,U_s(1,LC),IJKMAX2,NEXT_REC)
            CALL IN_BIN_R(UNIT_SPX+4,V_s(1,LC),IJKMAX2,NEXT_REC)
            CALL IN_BIN_R(UNIT_SPX+4,W_s(1,LC),IJKMAX2,NEXT_REC)
100      CONTINUE
         REC_POINTER(4) = NEXT_REC
      END IF
!
! ".SP5" FILE         ROP_s
!
      IF (READ_SPX(5).AND..NOT.AT_EOF(5)) THEN
         IF(.NOT.SPX_OPEN(5)) THEN
           WRITE(*,*)' SP5 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(5)
         IF (NEXT_REC.EQ.LAST_REC(5)) THEN
            AT_EOF(5) = .TRUE.
            RETURN
         END IF
         AT_EOF(5) = .FALSE.
         READ (UNIT_SPX+5,REC=NEXT_REC) TIME_REAL(5) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 200 LC = 1,MMAX
            CALL IN_BIN_R(UNIT_SPX+5,ROP_s(1,LC),IJKMAX2,NEXT_REC)
200      CONTINUE
         REC_POINTER(5) = NEXT_REC
      END IF
!
! ".SP6" FILE         T_g  , T_s1 , T_s2
!
      IF (READ_SPX(6).AND..NOT.AT_EOF(6)) THEN
         IF(.NOT.SPX_OPEN(6)) THEN
           WRITE(*,*)' SP6 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(6)
         IF (NEXT_REC.EQ.LAST_REC(6)) THEN
            AT_EOF(6) = .TRUE.
            RETURN
         END IF
         AT_EOF(6) = .FALSE.
         READ (UNIT_SPX+6,REC=NEXT_REC) TIME_REAL(6) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+6,T_g,IJKMAX2,NEXT_REC)
         IF(VERSION_NUMBER .LE. 1.)THEN
           CALL IN_BIN_R (UNIT_SPX+6,T_s(1,1),IJKMAX2,NEXT_REC)
           IF(MMAX .GE. 2) THEN
             CALL IN_BIN_R (UNIT_SPX+6,T_s(1,2),IJKMAX2,NEXT_REC)
           ELSE
             CALL IN_BIN_R (UNIT_SPX+6,Tmp,IJKMAX2,NEXT_REC)
           ENDIF
         ELSE
           DO 220 LC = 1,MMAX
             CALL IN_BIN_R(UNIT_SPX+6,T_s(1,LC),IJKMAX2,NEXT_REC)
220        CONTINUE
         ENDIF
         REC_POINTER(6) = NEXT_REC
      END IF
!
!
! ".SP7" FILE         X_g, X_s
!
 !     IF (READ_SPX(7).AND..NOT.AT_EOF(7)) THEN
 !        IF(.NOT.SPX_OPEN(7)) THEN
 !          WRITE(*,*)' SP7 file is not open'
 !          STOP
 !        ENDIF
 !        NEXT_REC = REC_POINTER(7)
 !        IF (NEXT_REC.ge.LAST_REC(7)) THEN
 !           AT_EOF(7) = .TRUE.
 !           RETURN
 !        END IF
 !        AT_EOF(7) = .FALSE.
 !        READ (UNIT_SPX+7,REC=NEXT_REC) TIME_REAL(7) , NSTEP
  !       NEXT_REC = NEXT_REC + 1
 !        DO 250 N = 1,1 ! NMAX(0)
 !          write (*,*) n,next_rec
 !          CALL IN_BIN_R(UNIT_SPX+7,X_g(1,N),IJKMAX2,NEXT_REC)
!250      CONTINUE
!         DO 300 LC = 1,MMAX
!          DO 270 N = 1, NMAX(LC)
 !            write (*,*) lc,n,next_rec
 !            CALL IN_BIN_R(UNIT_SPX+7,X_s(1,LC, N),IJKMAX2,NEXT_REC)
!270        CONTINUE
!300      CONTINUE
!         next_rec = next_rec + 28578
!         REC_POINTER(7) = NEXT_REC
 !     END IF


      IF (READ_SPX(7).AND..NOT.AT_EOF(7)) THEN
         IF(.NOT.SPX_OPEN(7)) THEN
           WRITE(*,*)' SP7 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(7)
         IF (NEXT_REC.EQ.LAST_REC(7)) THEN
            AT_EOF(7) = .TRUE.
            RETURN
         END IF
         AT_EOF(7) = .FALSE.
         READ (UNIT_SPX+7,REC=NEXT_REC) TIME_REAL(7) , NSTEP
         NEXT_REC = NEXT_REC + 1

          DO 250 N = 1, NMAX(0)
           if (bRead_all .or. n.eq.gas_species_index) then
              CALL IN_BIN_R(UNIT_SPX+7,X_g(1,N),IJKMAX2,NEXT_REC)
           else
              next_rec = next_rec + num_recs
           end if
 250      CONTINUE
         DO 300 LC = 1,MMAX
           DO 270 N = 1, NMAX(LC)
             if (bRead_all .or. &
               (lc.eq.solid_index .and. n.eq.solid_species_index))then
                 CALL IN_BIN_R(UNIT_SPX+7,X_s(1,LC, N),IJKMAX2,NEXT_REC)
             else
                next_rec = next_rec + num_recs
             end if
270        CONTINUE
300      CONTINUE
         REC_POINTER(7) = NEXT_REC
      END IF



!
! ".SP8" FILE         THETA_m
!
      IF (READ_SPX(8).AND..NOT.AT_EOF(8)) THEN
         IF(.NOT.SPX_OPEN(8)) THEN
           WRITE(*,*)' SP8 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(8)
         IF (NEXT_REC.EQ.LAST_REC(8)) THEN
            AT_EOF(8) = .TRUE.
            RETURN
         END IF
         AT_EOF(8) = .FALSE.
         READ (UNIT_SPX+8,REC=NEXT_REC) TIME_REAL(8),NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 400 LC = 1,MMAX
            CALL IN_BIN_R(UNIT_SPX+8,THETA_m(1,LC),IJKMAX2,NEXT_REC)
400      CONTINUE
         REC_POINTER(8) = NEXT_REC
      END IF
!
! ".SP9" FILE         Scalar
!
      IF (READ_SPX(9).AND..NOT.AT_EOF(9)) THEN
         IF(.NOT.SPX_OPEN(9)) THEN
           WRITE(*,*)' SP9 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(9)
         IF (NEXT_REC.EQ.LAST_REC(9)) THEN
            AT_EOF(9) = .TRUE.
            RETURN
         END IF
         AT_EOF(9) = .FALSE.
         READ (UNIT_SPX + 9, REC=NEXT_REC) TIME_REAL(9), NSTEP
         NEXT_REC = NEXT_REC + 1
         DO LC = 1, NScalar
            CALL IN_BIN_R (UNIT_SPX + 9,Scalar(1,LC) , IJKMAX2,NEXT_REC)
         END DO
         REC_POINTER(9) = NEXT_REC
      ENDIF


!     spa : ReactionRates

      IF (READ_SPX(10).AND..NOT.AT_EOF(10)) THEN
         IF(.NOT.SPX_OPEN(10)) THEN
           WRITE(*,*)' SPA file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(10)
         IF (NEXT_REC.GE.LAST_REC(10)) THEN
            AT_EOF(10) = .TRUE.
            RETURN
         END IF
         AT_EOF(10) = .FALSE.
         READ (UNIT_SPX+10,REC=NEXT_REC) TIME_REAL(10) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO N = 1, nRR
           CALL IN_BIN_R(UNIT_SPX+10,ReactionRates(1,N),IJKMAX2, &
                                                        NEXT_REC)
         end do
         REC_POINTER(10) = NEXT_REC
      END IF
!
! ".SP11" FILE         Scalar
!
      IF (READ_SPX(11).AND..NOT.AT_EOF(11)) THEN
         IF(.NOT.SPX_OPEN(11)) THEN
           WRITE(*,*)' SP11 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(11)
         IF (NEXT_REC.EQ.LAST_REC(11)) THEN
            AT_EOF(11) = .TRUE.
            RETURN
         END IF
         AT_EOF(11) = .FALSE.
         READ (UNIT_SPX + 11, REC=NEXT_REC) TIME_REAL(11), NSTEP
         NEXT_REC = NEXT_REC + 1

         IF(K_Epsilon) THEN
           CALL IN_BIN_R (UNIT_SPX + 11,K_Turb_G , IJKMAX2,NEXT_REC)
           CALL IN_BIN_R (UNIT_SPX + 11,E_Turb_G , IJKMAX2,NEXT_REC)
         ENDIF
         REC_POINTER(11) = NEXT_REC
      ENDIF


      RETURN
      END
