!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_K_s                                                C
!  Purpose: Calculate the effective conductivity of solids phases      C
!                                                                      C
!                                                                      C
!  Comments:                                                           C
!  This routine will not be called if k_s0(M) is defined               C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_K_S(M)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: undefined
      USE physprop, only: k_s0, k_s
      USE sendrecv, only: send_recv
! invoke user defined quantity
      USE usr_prop, only: usr_ks, calc_usr_prop
      USE usr_prop, only: solids_conductivity
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

!---------------------------------------------------------------------//


      IF (USR_Ks(M)) THEN
         CALL CALC_USR_PROP(Solids_Conductivity,lm=M)
      ELSEIF (K_s0(M) == UNDEFINED) THEN
! unncessary check but included for clarity
         CALL CALC_DEFAULT_Ks(M)
      ENDIF

      CALL send_recv(K_S, 2)

      RETURN
      END SUBROUTINE CALC_K_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DEFAULT_Ks                                         C
!  Purpose: Compute the default value for each solids phases           C
!  conductivity where a solids phase is considered to be comprised of  C
!  discrete ash particles                                              C
!                                                                      C
!  Author:M. Syamlal                                  Date: 24-APR-96  C
!                                                                      C
!  Literature/Document References:                                     C
!  Bauer & Schlunder's (1978) theory                                   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DEFAULT_Ks(M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      USE fldvar, only: ep_g, ep_s
      USE param1, only: one, zero
      USE physprop, only: k_g, k_s
      USE toleranc, only: dil_ep_s
      USE run, only: units
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids index
      INTEGER, INTENT(IN) :: M

! Local parameters
!---------------------------------------------------------------------//
! microscopic conductivity of ash in cal/s.cm.K
! (not modified by the gas phase)
      DOUBLE PRECISION :: Ks_micro
      PARAMETER (Ks_micro = 0.5258D-2)    !(2.2 J/s.m.K)
! constant in conductivity equation
      DOUBLE PRECISION :: PHI_k
      PARAMETER (PHI_k = 7.26D-3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
! Quantities in solids conductivity formula
      DOUBLE PRECISION :: BB, R_km, BoR, L_rm
!  Transform K_g(IJK) into the CGS if we work with SI
      DOUBLE PRECISION :: Kg_micro
!---------------------------------------------------------------------//

!!$omp parallel do private(IJK,B,R_km,BoR,L_rm,Kg_micro) &
!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

! All calculations are in CGS (1 cal = 4.183925J) and then converted
! to SI if needed
            IF (UNITS == 'SI') THEN
!convert k_g to CGS units (cal/s.cm.K)
               Kg_micro = K_g(IJK)/418.3925D0
            ELSE
! k_g already in CGS units (cal/s.cm.K)
               Kg_micro = K_g(IJK)
            ENDIF

            IF( EP_s(IJK,M) >  DIL_EP_s) THEN
               BB = 1.25D0 * ((ONE - EP_g(IJK))/EP_g(IJK))**(10.D0/9.D0)
               R_km = Ks_micro/Kg_micro
               BoR  = BB/R_km
               L_rm = -(2.d0/(ONE-BoR)) * &
                      ( ((R_km-ONE)/(ONE-BoR)**2)*BoR*LOG(BoR) + &
                        (BB-ONE)/(ONE-BoR) + (BB+ONE)/2.d0 )
! K_s is the macroscopic conductivity that has been modified by the
! presence of the gas phase (cal/s.cm.K)
               K_S(IJK,M) = (Phi_k*R_km + (ONE-Phi_k)*L_rm)*&
                             Kg_micro/SQRT(ONE - EP_g(IJK))
            ELSE
               K_S(IJK, M) = ZERO
            ENDIF

! An approximate average value for the solids conductivity is 2.5*K_g
!            K_S(IJK,M) = 2.5*Kg_micro            !in CGS system

         ELSE   ! else branch if(fluid_at(ijk))
            K_S(IJK,M) = ZERO
         ENDIF   ! end if/else (fluid_at(ijk))

         IF (UNITS == 'SI') K_s(IJK, M) = 418.3925D0*K_s(IJK, M)   !J/s.m.K

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE CALC_DEFAULT_KS


