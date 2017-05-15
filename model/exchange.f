!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: EXCHANGE                                                C
!  Purpose: Calls routines to drive calculations of the interphase     C
!           momentum, and energy exchange coefficients/terms           C
!           if directed to do so by the corresponding flags            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE EXCHANGE(IER)

! Global Variables
!---------------------------------------------------------------------//
! Flags for calculating drag coefficient.
      use coeff, only: DRAGCOEF
! Flags for calculating heat transfer coefficient.
      use coeff, only: HEAT_TR

      use param1, only: zero
      use physprop, only: smax, ro_g0
      use run, only: granular_energy
      use run, only: kt_type_enum, ia_2005

      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_CONTINUUM_HYBRID


      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(INOUT) :: IER ! Error index

! Local variables
!---------------------------------------------------------------------//
! loop counter
      INTEGER :: M, L
!---------------------------------------------------------------------//

! calculate gas-solids drag based on relatively velocity differences
      IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
         DO M = 1, SMAX
            IF (DRAGCOEF(0,M) .AND. RO_G0/=ZERO) CALL DRAG_GS(M, IER)
         ENDDO

! calculate solilds-solids drag based on relative velocity differences
         DO M = 1, SMAX
            DO L = 1, M - 1
               IF (DRAGCOEF(L,M)) CALL DRAG_SS (L, M, IER)
            ENDDO
         ENDDO
      ENDIF

! Calculate additional interphase interaction coefficients (between
! continuum solids phases)
      IF (GRANULAR_ENERGY) THEN
         SELECT CASE(KT_TYPE_ENUM)
            CASE(IA_2005)
               DO M=1,SMAX
                  DO L=1,SMAX
                     CALL COLL_MOMENTUM_COEFF_IA(L, M)
                  ENDDO
               ENDDO
            CASE DEFAULT
         END SELECT
      ENDIF

! Calculate interphase heat transfer coefficients
      DO M=1,SMAX
         IF(HEAT_TR(0,M)) CALL CALC_GAMA(M)
      ENDDO

      return
      END SUBROUTINE EXCHANGE
