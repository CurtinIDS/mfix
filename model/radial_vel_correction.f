!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: RADIAL_VEL_CORRECTION                                   !
!  Author: Akhilesh Bakshi                            Date: 28-JUL-14  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!  Reference: A. Bakshi, C. Altantzis, A.F. Ghoniem, "Towards accurate !
!  three-dimensional simulation of dense multi-phase flows using       !
!  cylindrical coordinates," Powder Technology, Volume 264, 09/2014,   !
!  Pages 242-255.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE RADIAL_VEL_CORRECTION

      use physprop
      use fldvar
      use constant
      use mpi_utility
      use functions

      IMPLICIT NONE

! phase index
      INTEGER :: M

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Wg_Temp, Wg_Temp_GL
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Ws_Temp, Ws_Temp_GL

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SUMX_G, SUMZ_G
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SUMX_S, SUMZ_S

      INTEGER :: I1, I2, J1, K1, K2
      INTEGER :: abIJK, abIJK1, abIJK2

      DOUBLE PRECISION :: Angle_Temp, WgAvg, WsAvg

!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: VAVG_U_G, VAVG_V_G, VAVG_W_G, &
                                    VAVG_U_S, VAVG_V_S, VAVG_W_S

      IF(NO_K .OR. KMAX==1) RETURN

      I1=1
      I2=2

! Fukagata Scheme

! Implement for gas
      ALLOCATE (Wg_Temp(ijkstart3:ijkend3))

      DO J1 = JSTART3, JEND3
         DO K1 = KSTART3, KEND3
         IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE

            abIJK1=FUNIJK(I2,J1,K1)
            Wg_Temp (abIJK1) = W_g(abIJK1)
         END DO ! K loop Ends
      END DO ! J Loop Ends

      ALLOCATE (Wg_Temp_GL(ijkmax3))
      CALL GATHER (Wg_Temp, Wg_Temp_GL, root)
      CALL BCAST (Wg_Temp_GL, root)

      ALLOCATE (SUMX_G(JMIN3:JMAX3))
      ALLOCATE (SUMZ_G(JMIN3:JMAX3))

      DO J1=1, JMAX3
         SUMX_G(J1)=0
         SUMZ_G(J1)=0
         DO K1=1, KMAX
            Angle_Temp = (K1-1)*2*(Pi/KMAX)
            IF (K1 .LE. 0.5*KMAX) THEN
               K2=(K1-1)+(0.5*KMAX2)
            ELSE
               K2=(K1+1)-(0.5*KMAX2)
            ENDIF

            abIJK1=FUNIJK_GL (I2,J1,K1)
            abIJK2=FUNIJK_GL (I2,J1,K2)
            WgAvg=0.5*(Wg_Temp_GL(abIJK1)-Wg_Temp_GL(abIJK2))
            SUMX_G(J1)=SUMX_G(J1)-((2*WgAvg*sin(Angle_Temp))/KMAX)
            SUMZ_G(J1)=SUMZ_G(J1)+((2*WgAvg*cos(Angle_Temp))/KMAX)
         ENDDO
      ENDDO

      DO J1= JSTART3, JEND3
         DO K1= KSTART3, KEND3
            !abIJK in local
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            abIJK = FUNIJK (I1,J1,K1)
            Angle_Temp = (-Pi/KMAX)+((K1-1)*2*(Pi/KMAX))
            U_g(abIJK)=(SUMX_G(J1)*cos(Angle_Temp)) + &
                (SUMZ_G(J1)*sin(Angle_Temp))
         ENDDO
      ENDDO

      DEALLOCATE (Wg_Temp)
      DEALLOCATE (Wg_Temp_GL)
      DEALLOCATE (SUMX_G)
      DEALLOCATE (SUMZ_G)



! Implement for Solid
      ALLOCATE (Ws_Temp(ijkstart3:ijkend3))
      ALLOCATE (Ws_Temp_GL(ijkmax3))
      ALLOCATE (SUMX_S(JMIN3:JMAX3))
      ALLOCATE (SUMZ_S(JMIN3:JMAX3))

      DO M=1, MMAX

         DO J1 = JSTART3, JEND3
            DO K1 = KSTART3, KEND3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               abIJK1=FUNIJK(I2,J1,K1)
               Ws_Temp (abIJK1) = W_s(abIJK1,M)
            ENDDO ! K loop Ends
         ENDDO ! J Loop Ends

         CALL GATHER (Ws_Temp, Ws_Temp_GL, root)
         CALL BCAST (Ws_Temp_GL, root)

         DO J1=1, JMAX3
            SUMX_S(J1)=0
            SUMZ_S(J1)=0
            DO K1=1, KMAX
               Angle_Temp = (K1-1)*2*(Pi/KMAX)
               IF (K1 .LE. 0.5*KMAX) THEN
                  K2=(K1-1)+(0.5*KMAX2)
               ELSE
                  K2=(K1+1)-(0.5*KMAX2)
               END IF
               abIJK1=FUNIJK_GL (I2,J1,K1)
               abIJK2=FUNIJK_GL (I2,J1,K2)
               WsAvg=0.5*(Ws_Temp_GL(abIJK1)-Ws_Temp_GL(abIJK2))
               SUMX_S(J1)=SUMX_S(J1)-((2*WsAvg*sin(Angle_Temp))/KMAX)
               SUMZ_S(J1)=SUMZ_S(J1)+((2*WsAvg*cos(Angle_Temp))/KMAX)
            ENDDO
         ENDDO

         DO J1= JSTART3, JEND3
            DO K1= KSTART3, KEND3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               abIJK = FUNIJK (I1,J1,K1)
               Angle_Temp = (-Pi/KMAX)+((K1-1)*2*(Pi/KMAX))
               U_s(abIJK,M)=(SUMX_S(J1)*cos(Angle_Temp)) + &
                  (SUMZ_S(J1)*sin(Angle_Temp))
            ENDDO
         ENDDO
      ENDDO ! M Loop Ends

      DEALLOCATE (Ws_Temp)
      DEALLOCATE (Ws_Temp_GL)
      DEALLOCATE (SUMX_S)
      DEALLOCATE (SUMZ_S)

      RETURN
      END SUBROUTINE RADIAL_VEL_CORRECTION
