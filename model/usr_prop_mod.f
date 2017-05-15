! or put the actual routine calc_usr_source here? move logic variable
!  usr_source into run or something...
      module usr_prop
      use param, only: dim_m

! If .TRUE. call user-defined physical properties
      LOGICAL :: USR_ROG
      LOGICAL :: USR_CPG 
      LOGICAL :: USR_ROS(DIM_M)
      LOGICAL :: USR_CPS(DIM_M)

! If .TRUE. call user_defined transport properties
      LOGICAL :: USR_MUG
      LOGICAL :: USR_KG
      LOGICAL :: USR_DIFG
      LOGICAL :: USR_MUS(DIM_M)
      LOGICAL :: USR_KS(DIM_M)
      LOGICAL :: USR_DIFS(DIM_M)

! Interphase momentum exchange coefficients (F_GS, F_SS)
      LOGICAL :: USR_FGS(DIM_M)
      LOGICAL :: USR_FSS( (DIM_M*(DIM_M-1)/2)+1 )
! Interphase heat transfer coefficient (GAMA)
      LOGICAL :: USR_GAMA(DIM_M)

      ENUM, BIND(C)
         ENUMERATOR :: Gas_Density, Solids_Density
         ENUMERATOR :: Gas_SpecificHeat, Solids_SpecificHeat
         ENUMERATOR :: Gas_Viscosity, Solids_Viscosity
         ENUMERATOR :: Gas_Conductivity, Solids_Conductivity
         ENUMERATOR :: Gas_Diffusivity, Solids_Diffusivity 
         ENUMERATOR :: GasSolids_Drag, SolidsSolids_Drag
         ENUMERATOR :: GasSolids_HeatTransfer
         ENUMERATOR :: BLANK
         ENUMERATOR :: DUMMY
      END ENUM


      CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_USR_PROP                                           C 
!  Purpose: Driver routine to calculate user defined physical,         C
!  transport, exchange terms that appear in the governing equations.   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_USR_PROP(lprop, lM, lL, lerr)

! Modules
!-----------------------------------------------
      use compar, only: ijkstart3, ijkend3, myPE, numPEs
      use functions, only: fluid_at, wall_at, fluidorp_flow_at
      use output, only: REPORT_NEG_DENSITY
      use output, only: REPORT_NEG_SPECIFICHEAT
      use param1, only: undefined_i, zero
      use physprop, only: nmax
      use physprop, only: k_g, dif_g, c_pg, mu_g
      use physprop, only: k_s, dif_s, c_ps
      use fldvar, only: ro_g, rop_g, ep_g
      use fldvar, only: ro_s, p_s
      use visc_g, only: mu_gt, lambda_gt
      use visc_s, only: mu_s, lambda_s
      use utilities, only: mfix_isnan

      use error_manager

      IMPLICIT NONE
! Dummy arguments
!-----------------------------------------------
! User properties available
!     ROg, ROs (gas and solids density)
!     CPg, CPs (gas and solids specific heat)
!     Mug, Mus (gas and solids viscosity)
!     Kg, Ks (gas and solids conductivity)
!     Difg, Difs (gas and solids diffusivity)
!     Fgs, Fss (gas-solids and solids-solids drag)
!     Gama (gas-solids heat transfere coefficient)
      INTEGER, INTENT(IN) :: lprop

! Phase index 
      INTEGER, OPTIONAL, INTENT(IN) :: lM

! Phase index 
      INTEGER, OPTIONAL, INTENT(IN) :: lL

! Error index
! Arrays for storing errors:
! 100 - Negative gas phase density
! 101 - Negative solids phase density
! 10x - Unclassified
      INTEGER, OPTIONAL, INTENT(INOUT) :: lerr(0:numPEs-1)

! Local variables
!-----------------------------------------------
! indices
      INTEGER :: IJK
! tmp indices
      INTEGER :: M, L, N
! Flag to write log header
      LOGICAL :: wHeader

!-----------------------------------------------

! set local values for phase index
      IF (.NOT.present(lM)) THEN
         M = UNDEFINED_I
      ELSE
         M = lM
      ENDIF

! set local values for phase index
      IF (.NOT.present(lL)) THEN
         L = UNDEFINED_I
      ELSE
         L = lL
      ENDIF
         
      N = UNDEFINED_I

      SELECT CASE(lProp)
      
! Gas physical properties: density
      CASE(Gas_Density)
         wHeader = .TRUE.
         DO IJK=IJKSTART3,IJKEND3
            IF (WALL_AT(IJK)) CYCLE
!            CALL USR_PROPERTIES(lprop, IJK, 0, N)
            CALL USR_PROP_ROG(IJK)
            ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)
! report errors
            IF(RO_G(IJK) < ZERO .or. mfix_isnan(RO_G(IJK))) THEN
               lErr(myPE) = 100  ! return appropritae error flag
               IF(REPORT_NEG_DENSITY) CALL ROgErr_LOG(IJK, wHeader)
            ENDIF

         ENDDO

! Gas physical properties: specific heat
      CASE(Gas_SpecificHeat)
         DO IJK=IJKSTART3,IJKEND3
            IF (WALL_AT(IJK)) CYCLE
!            CALL USR_PROPERTIES(lprop, IJK, 0, N)
            CALL USR_PROP_CPG(IJK)
! report errors.
            IF(C_PG(IJK) <= ZERO) THEN
!               lErr(myPE) = 103
               IF(REPORT_NEG_SPECIFICHEAT) CALL CPgErr_LOG(IJK, wHeader)
            ENDIF
         ENDDO

! Gas transport properties: viscosity
      CASE(Gas_Viscosity)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
!               CALL USR_PROPERTIES(lprop, IJK, 0, N)
               CALL USR_PROP_MUG(IJK)
               MU_GT(IJK) = MU_G(IJK)
               LAMBDA_GT(IJK) = -2./3.*MU_GT(IJK)
            ELSE ! probably unnecessary
               MU_g(IJK) = zero
               MU_GT(IJK) = ZERO
               LAMBDA_GT(IJK) = ZERO
            ENDIF
         ENDDO

! Gas transport properties: conductivity
      CASE(Gas_Conductivity)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
!               CALL USR_PROPERTIES(lprop, IJK, 0, N)
               CALL USR_PROP_KG(IJK)
            ELSE ! probably unnecessary
               K_G(IJK) = ZERO
            ENDIF
         ENDDO

! Gas transport properties: diffusivity
      CASE(Gas_Diffusivity)
         DO N=1,NMAX(0)
            DO IJK=IJKSTART3,IJKEND3
               IF (FLUID_AT(IJK)) THEN
!                  CALL USR_PROPERTIES(lprop, IJK, 0, N)
                  CALL USR_PROP_DIFG(IJK,N)
               ELSE ! probably unnecessary
                  DIF_G(IJK,N) = ZERO
               ENDIF
            ENDDO
         ENDDO
 
! Solids physical properties: density
      CASE(Solids_Density)
         wHeader = .TRUE.
         DO IJK=IJKSTART3,IJKEND3
            IF (WALL_AT(IJK)) CYCLE
!            CALL USR_PROPERTIES(lprop, IJK, M, N)
            CALL USR_PROP_ROS(IJK,M)
! report errors.
            IF(RO_S(IJK,M) <= ZERO) THEN
               lErr(myPE) = 101
               IF(REPORT_NEG_DENSITY) CALL ROsErr_LOG(IJK, M, wHeader)
            ENDIF

         ENDDO

! Solids physical properties: specific heat
      CASE(Solids_SpecificHeat)
         DO IJK=IJKSTART3,IJKEND3
            IF (WALL_AT(IJK)) CYCLE
!            CALL USR_PROPERTIES(lprop, IJK, M, N)
            CALL USR_PROP_CPS(IJK,M)
! report errors.
            IF(C_PS(IJK,M) <= ZERO) THEN
!               lErr(myPE) = 104
               IF(REPORT_NEG_SPECIFICHEAT) CALL CPsErr_LOG(IJK, M, wHeader)
            ENDIF
         ENDDO

! Solids transport properties: viscosity
      CASE(Solids_Viscosity)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
!               CALL USR_PROPERTIES(lprop, IJK, M, N)
               CALL USR_PROP_MUS(IJK,M)
            ELSE ! probably unnecessary
               MU_S(IJK,M) = ZERO
               LAMBDA_S(IJK,M) = ZERO
            ENDIF
         ENDDO

! Solids transport properties: conductivity
      CASE(Solids_Conductivity)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
!               CALL USR_PROPERTIES(lprop, IJK, M, N)
               CALL USR_PROP_KS(IJK,M)
            ELSE ! probably unnecessary
               K_S(IJK,M) = ZERO
            ENDIF
         ENDDO

! Solids transport properties: diffusivity
      CASE(Solids_Diffusivity)
         DO N=1,NMAX(M)
            DO IJK=IJKSTART3,IJKEND3
               IF (FLUID_AT(IJK)) THEN
!                  CALL USR_PROPERTIES(lprop, IJK, M, N)
                  CALL USR_PROP_DIFS(IJK,M,N)
               ELSE ! probably unnecessary
                  DIF_S(IJK,M,N) = ZERO
               ENDIF
            ENDDO
         ENDDO

! Gas-solids exchange: Interphase heat transfer coefficient (GAMA)
      CASE(GasSolids_HeatTransfer)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
!               CALL USR_PROPERTIES(lprop, IJK, M, N)
               CALL USR_PROP_GAMA(IJK,M)
            ENDIF
         ENDDO

! Gas-solids exchange: interphase drag
      CASE(GasSolids_Drag)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUIDorP_FLOW_AT(IJK)) THEN
               CALL USR_PROPERTIES(lprop, IJK, M, N)
! this hook is not accessed or functioning.
! use drag_type = usr_drag to access this entity
            ENDIF
         ENDDO

! Solids-solids exchange: interphase drag
      CASE(SolidsSolids_Drag)
         DO IJK=IJKSTART3,IJKEND3
            IF (WALL_AT(IJK)) CYCLE
!            CALL USR_PROPERTIES(lprop, IJK, M, L)
            CALL USR_PROP_FSS(IJK,L,M)
         ENDDO

! error out
      CASE DEFAULT
! should never hit this
! Initialize the error manager.
         CALL INIT_ERR_MSG("CALC_USR_PROP")
         WRITE(ERR_MSG, 1001) ival(lprop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1001 FORMAT('Error 1101: Unknown Property= ', A)
         CALL FINL_ERR_MSG
      END SELECT   ! end selection of user property 
 
      RETURN
      END SUBROUTINE CALC_USR_PROP

      end module usr_prop

