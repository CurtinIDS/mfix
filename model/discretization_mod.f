!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DISCRETIZATION                                         C
!  Purpose: A collection of functions for higher-order discretization  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-MAR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE discretization
      IMPLICIT NONE

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION SUPERBEE (PHI_C)
      USE param1, only: one, half, zero
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C

      DOUBLE PRECISION :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         SUPERBEE = HALF*MAX(ZERO,MIN(ONE,2.d0*TH),MIN(2.d0,TH))
      ELSE IF (PHI_C == ONE) THEN
         SUPERBEE = ONE
      ELSE                                       !first order upwinding
         SUPERBEE = ZERO
      ENDIF
      RETURN
      END FUNCTION SUPERBEE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION SMART (PHI_C)
      USE param1, only: zero, half, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C

      DOUBLE PRECISION :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         SMART = HALF*MAX(ZERO,MIN(4.0D0*TH,0.75D0 + 0.25D0*TH,2.0D0))
      ELSE IF (PHI_C == ONE) THEN
         SMART = ONE
      ELSE                                       !first order upwinding
         SMART = ZERO
      ENDIF
      RETURN
      END FUNCTION SMART


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!       calculate DWF from Chi_SMART scheme                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION Chi_SMART (PHI_C, Chi)
      USE param1, only: zero, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: Chi

      if(PHI_C < ONE)then
        Chi_SMART = Chi * (3.D0/8.D0 - PHI_C/4.d0) / (ONE-PHI_C)
      elseif(Chi == zero)then ! insures that all species equations uses the same chi
        Chi_SMART = zero
      else
      Chi_SMART = ONE
      endif
      RETURN
      END FUNCTION Chi_SMART


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!       calculate CHI for SMART scheme                                 C
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION Chi4SMART (PHI_C, PHIU, PHIC, PHID)
      USE param1, only: zero, one, large_number
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: PHIU
      DOUBLE PRECISION, INTENT(IN) :: PHIC
      DOUBLE PRECISION, INTENT(IN) :: PHID

      IF(PHI_C > ZERO .AND. PHI_C <= (1.D0/6.D0))THEN
         Chi4SMART = 16.D0 * PHI_C/(3.D0 - 2.D0 * PHI_C)
      ELSEIF(PHI_C > (1.D0/6.D0) .AND. PHI_C <= (5.D0/6.D0))THEN
         Chi4SMART = ONE
      ELSEIF(PHI_C > (5.D0/6.D0) .AND. PHI_C <= ONE)THEN
         Chi4SMART = 8.D0*(ONE-PHI_C)/(3.D0 - 2.D0 * PHI_C)
      ELSEIF( (PHIU < 1d-15) .AND. (PHIC < 1d-15) .AND. (PHID < 1d-15) )THEN
! if a species Xg is less than machine precision, do not use its chi.
! this will avoid calculating chi = 0 for a vanishing species.
         Chi4SMART = LARGE_NUMBER
      ELSE
         Chi4SMART = ZERO
      ENDIF
      RETURN
      END FUNCTION Chi4SMART


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION ULTRA_QUICK (PHI_C, CF)
      USE param1, only: zero, half, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: CF

      DOUBLE PRECISION, PARAMETER :: FIVEOSIX = 5.D0/6.D0
      DOUBLE PRECISION :: TH, OCF

      OCF = MAX(ONE,ONE/MAX(1.D-2,CF))
      IF (PHI_C > ONE) THEN
         ULTRA_QUICK = HALF
      ELSE IF (PHI_C > FIVEOSIX) THEN
         ULTRA_QUICK = ONE
      ELSE IF (PHI_C > 3.D0/(8.D0*OCF - 6.D0)) THEN
         TH = PHI_C/(ONE - PHI_C)
         ULTRA_QUICK = HALF - 0.125D0*(ONE - TH)
      ELSE IF (PHI_C > ZERO) THEN
         TH = PHI_C/(ONE - PHI_C)
         ULTRA_QUICK = (OCF - ONE)*TH
      ELSE
         TH = PHI_C/(ONE - PHI_C)
         ULTRA_QUICK = HALF*TH
      ENDIF
      RETURN
      END FUNCTION ULTRA_QUICK


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION QUICKEST (PHI_C, CF, ODXC, ODXUC, ODXCD)
      USE param1, only: zero, half, one, small_number
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: CF
      DOUBLE PRECISION, INTENT(IN) :: ODXC
      DOUBLE PRECISION, INTENT(IN) :: ODXUC
      DOUBLE PRECISION, INTENT(IN) :: ODXCD

      DOUBLE PRECISION :: FCF, TH

      IF (PHI_C>ZERO .AND. PHI_C<ONE) THEN       !monotonic region
         FCF = -(ONE - CF*CF)/3.D0
         TH = PHI_C/(ONE - PHI_C + SMALL_NUMBER)
         QUICKEST = HALF*(ONE - CF) + FCF*(ODXC/ODXCD - ODXC*ODXUC*TH/ODXCD**2)
         IF(PHI_C<CF)QUICKEST=MIN(QUICKEST,(ONE/CF-ONE)*PHI_C/(ONE-PHI_C))
         QUICKEST = MAX(ZERO,MIN(ONE,QUICKEST))
      ELSE                                       !first order upwinding
         QUICKEST = ZERO
      ENDIF
      RETURN
      END FUNCTION QUICKEST


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION MUSCL (PHI_C)
      USE param1, only: zero, half, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C

      DOUBLE PRECISION :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         MUSCL = HALF*MAX(ZERO,MIN(2.0D0*TH,(ONE + TH)/2.0D0,2.0D0))
      ELSE IF (PHI_C == ONE) THEN
         MUSCL = ONE
      ELSE                                       !first order upwinding
         MUSCL = ZERO
      ENDIF
      RETURN
      END FUNCTION MUSCL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!       calculate DWF from Chi_MUSCL scheme                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION Chi_MUSCL (PHI_C, Chi)
      USE param1, only: zero, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: Chi

      if(PHI_C < ONE)then
        Chi_MUSCL = Chi /(4.D0 * (ONE - PHI_C))
      elseif(Chi == zero)then ! insures that all species equations uses the same chi
        Chi_MUSCL = zero
      else
        Chi_MUSCL = ONE
      endif
      RETURN
      END FUNCTION Chi_MUSCL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!       calculate CHI for MUSCL scheme                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION Chi4MUSCL (PHI_C, PHIU, PHIC, PHID)
      USE param1, only: zero, one, large_number
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: PHIU
      DOUBLE PRECISION, INTENT(IN) :: PHIC
      DOUBLE PRECISION, INTENT(IN) :: PHID

      IF(PHI_C > ZERO .AND. PHI_C <= 0.25d0)THEN
        Chi4MUSCL = 4.D0 * PHI_C
      ELSEIF(PHI_C > 0.25d0 .AND. PHI_C <= 0.75d0)THEN
        Chi4MUSCL = ONE
      ELSEIF(PHI_C > 0.75d0 .AND. PHI_C <= ONE)THEN
        Chi4MUSCL = 4.D0*(ONE - PHI_C)
      ELSEIF( (PHIU < 1d-15) .AND. (PHIC < 1d-15) .AND. (PHID < 1d-15) )THEN
! if a species Xg is less than machine precision, do not use its chi.
! this will avoid calculating chi = 0 for a vanishing species.
        Chi4MUSCL = LARGE_NUMBER
      ELSE
        Chi4MUSCL = ZERO
      ENDIF
      RETURN
      END FUNCTION Chi4MUSCL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION VANLEER (PHI_C)
      USE param1, only: zero, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C

      IF (PHI_C>=ZERO .AND. PHI_C<=ONE) THEN     !monotonic region
         VANLEER = PHI_C
      ELSE                                       !first order upwinding
         VANLEER = ZERO
      ENDIF
      RETURN
      END FUNCTION VANLEER


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION MINMOD (PHI_C)
      USE param1, only: zero, half, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C

      IF (PHI_C>=ZERO .AND. PHI_C<=ONE) THEN     !monotonic region
         IF (PHI_C >= HALF) THEN                 !central differencing
            MINMOD = HALF
         ELSE                                    !second order upwinding
            MINMOD = HALF*PHI_C/(ONE - PHI_C)
         ENDIF
      ELSE                                       !first order upwinding
         MINMOD = ZERO
      ENDIF

      RETURN
      END FUNCTION MINMOD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Central scheme for MMS cases.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION CENTRAL_SCHEME (PHI_C)
      USE param1, only: half
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C

      CENTRAL_SCHEME = HALF
      RETURN
      END FUNCTION CENTRAL_SCHEME


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION UMIST (PHI_C)
      USE param1, only: zero, half, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_C

      DOUBLE PRECISION :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         UMIST = HALF*MAX(ZERO,MIN(2.0D0*TH,0.75D0 + 0.25D0*TH,2.0D0))
      ELSE IF (PHI_C == ONE) THEN
         UMIST = ONE
      ELSE                                       !first order upwinding
         UMIST = ZERO
      ENDIF
      RETURN
      END FUNCTION UMIST


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION XSI (V, DWF)
      USE param1, only: zero, one
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: V
      DOUBLE PRECISION, INTENT(IN) :: DWF

      IF (V >= ZERO) THEN
         XSI = DWF
      ELSE
         XSI = ONE - DWF
      ENDIF
      RETURN
      END FUNCTION XSI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION PHI_C_OF (PHI_U, PHI_C, PHI_D)
      USE param1, only: zero
      USE toleranc, only: compare
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_U
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: PHI_D

      DOUBLE PRECISION :: DEN

      IF (COMPARE(PHI_D,PHI_U)) THEN
         PHI_C_OF = ZERO
      ELSE
         DEN = PHI_D - PHI_U
         PHI_C_OF = (PHI_C - PHI_U)/DEN
      ENDIF
      RETURN
      END FUNCTION PHI_C_OF


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION FPFOI_OF (PHI_D, PHI_C, &
                       PHI_U, PHI_UU)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_D
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: PHI_U
      DOUBLE PRECISION, INTENT(IN) :: PHI_UU

      FPFOI_OF = PHI_C + (5.0D0/16.0D0)*(PHI_D - PHI_C) + &
                         (1.0D0/4.0D0)*(PHI_C - PHI_U) - &
                         (1.0D0/16.0D0)*(PHI_U - PHI_UU)

! LIMIT THE HIGH ORDER INTERPOLATION
      FPFOI_OF = UNIV_LIMITER_OF(FPFOI_OF, PHI_D, PHI_C, PHI_U)
      RETURN
      END FUNCTION FPFOI_OF


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION UNIV_LIMITER_OF (PHI_TEMP, PHI_D,&
                        PHI_C, PHI_U)

      USE param1, only: zero
      USE run, only: c_fac
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: PHI_TEMP
      DOUBLE PRECISION, INTENT(IN) :: PHI_D
      DOUBLE PRECISION, INTENT(IN) :: PHI_C
      DOUBLE PRECISION, INTENT(IN) :: PHI_U

      DOUBLE PRECISION :: PHI_REF, DEL, CURV

      DEL = PHI_D - PHI_U
      CURV = PHI_D - 2.0*PHI_C + PHI_U
      IF (ABS (CURV) >= ABS (DEL)) THEN     ! NON-MONOTONIC
         UNIV_LIMITER_OF = PHI_C
      ELSE
         PHI_REF = PHI_U + (PHI_C-PHI_U)/C_FAC
         IF (DEL > ZERO) THEN
           IF (PHI_TEMP < PHI_C)UNIV_LIMITER_OF = PHI_C
           IF (PHI_TEMP > MIN(PHI_REF,PHI_D)) &
              UNIV_LIMITER_OF = MIN(PHI_REF,PHI_D)
           IF (PHI_TEMP >= PHI_C .AND. PHI_TEMP <= MIN(PHI_REF,PHI_D)) &
              UNIV_LIMITER_OF = PHI_TEMP
         ELSE
           IF (PHI_TEMP > PHI_C)UNIV_LIMITER_OF = PHI_C
           IF (PHI_TEMP < MAX(PHI_REF,PHI_D)) &
              UNIV_LIMITER_OF = MAX(PHI_REF,PHI_D)
           IF (PHI_TEMP >= MAX(PHI_REF,PHI_D) .AND. PHI_TEMP <= PHI_C) &
              UNIV_LIMITER_OF = PHI_TEMP
         ENDIF
      ENDIF
      RETURN
      END FUNCTION UNIV_LIMITER_OF

      END MODULE discretization
