! TO DO:
! 1. Check the formulation based on MCP.
! 2. The pressure correction should be based on sum of close-packed
!    solids?

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_VOL_FR                                             C
!  Purpose: Calculate volume fractions of phases used in pressure      C
!           corrections.                                               C
!                                                                      C
!  Notes: see mark_phase_4_cor for more details                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-JUL-96   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_VOL_FR(P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER)


! Modules
!---------------------------------------------------------------------//
      USE compar
      USE constant
      USE discretelement
      USE functions
      USE geometry
      USE indices
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE pgcor
      USE physprop
      USE pscor
      USE run, only: ghd_2007, kt_type_enum
      USE sendrecv
      USE solids_pressure
      USE visc_s
      use fldvar, only: RO_S, EP_S

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Solids pressure
      DOUBLE PRECISION, INTENT(IN) :: P_star(DIMENSION_3)
! Gas density
      DOUBLE PRECISION, INTENT(INOUT) :: RO_g(DIMENSION_3)
! Gas bulk density
      DOUBLE PRECISION, INTENT(INOUT) :: ROP_g(DIMENSION_3)
! Gas volume fraction
      DOUBLE PRECISION, INTENT(INOUT) :: EP_g(DIMENSION_3)
! solids bulk densities
      DOUBLE PRECISION, INTENT(INOUT) :: ROP_s(DIMENSION_3, DIMENSION_M)
! error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!---------------------------------------------------------------------//
! volume of particle type M for GHD theory
      DOUBLE PRECISION :: VOL_M
! volume fraction of close-packed region
      DOUBLE PRECISION :: EPcp
! volume fraction of solids phase
      DOUBLE PRECISION :: EPS
! sum of volume fractions
      DOUBLE PRECISION :: SUMVF
! Whichever phase is given by MF becomes the phase whose volume fraction
! is corrected based on all other phases present.  Generally MF gets
! defined as 0 so that the gas phase void fraction becomes corrected
! based on the summation of the volume fractions of all other phases
      INTEGER :: MF
! Whichever phase is given by MCPl becomes the phase whose volume
! fraction is corrected based on value of maximum close packing and all
! other solids phase that can close pack.  This is basically a local
! variable for MCP but only in cells where close packed conditions
! exist.
      INTEGER :: MCPl
! Index of solids phase
      INTEGER :: M
! Indices
      INTEGER :: IJK

! Arrays for storing errors:
! 110 - Negative volume fraciton
! 11x - Unclassified
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

      LOGICAL, PARAMETER :: REPORT_NEG_VOLFRAC = .TRUE.

! Flag to write log header
      LOGICAL :: wHeader
!---------------------------------------------------------------------//


! Initialize:
      wHeader = .TRUE.
! Initialize error flag.
      Err_l = 0

!!$omp  parallel do private(MCPl, EPCP, SUMVF, MF, M) &
!!$omp&  schedule(static)

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

! calculate the volume fraction of the solids phase based on the volume
! fractions of all other solids phases and on the value of maximum
! packing when that solids phase continuity was not solved for the
! indicated cell.
!----------------------------------------------------------------->>>
            IF (PHASE_4_P_S(IJK) /= UNDEFINED_I) THEN
! for the current cell check if a solids phase continuity was skipped.
! this value will either be undefined (no cell was skipped in any
! solids continuity equations) or be a solids phase index (indicated
! solids continuity was skipped) for a solids phase that does close pack
! (i.e., close_packed=T) and the cell exhibits close packing conditions.
! note: this branch is never entered given the existing version of
! mark_phase_4_cor.
               MCPl   = PHASE_4_P_s(IJK)

! finding the value of maximum close packing based on the expression for
! plastic pressure
               EPCP = 1. - INV_H(P_STAR(IJK),EP_g_blend_end(ijk))
! summing the solids volume fraction of all continuum solids phases
! that are marked as close_packed except for the solids phase which
! is also marked by phase_4_p_s (skipping its own continuity)
               SUMVF = ZERO
               DO M = 1, MMAX+DES_MMAX
                  IF (CLOSE_PACKED(M) .AND. M/=MCPl) SUMVF = SUMVF + EP_S(IJK,M)
               ENDDO
               ROP_S(IJK,MCPl) = (EPCP - SUMVF)*RO_S(IJK,MCPl)
            ENDIF
!-----------------------------------------------------------------<<<


! calculate the volume fraction of the 'solids' phase based on the
! volume fractions of all other phases (including gas) if the gas
! continuity was solved rather than that phases own continuity.
! if the gas continuity was not solved then calculate the void
! fraction of the gas phase based on all other solids phases.
!----------------------------------------------------------------->>>
! for the current cell check if the gas phase continuity was solved for
! the gas phase while the solids phase continuity for the indicated
! solids phase was skipped. this value will either be 0 (gas phase
! continuity was not solved) or be a solids phase index (gas continuity
! was solved while the indicated solids phase continuity was skipped)
! for a solids phase that does not close pack (i.e., close_packed=F) and
! is in greater concentration than the gas phase.
! Note: MF will always be set to 0 here given the existing version of
! mark_phase_4_cor
            MF = PHASE_4_P_G(IJK)

            SUMVF = ZERO
            IF (0 /= MF) THEN
! if gas continuity was solved rather than the solids phase, then
! include the gas phase void fraction in the summation here.
               EP_G(IJK) = ROP_G(IJK)/RO_G(IJK)
               SUMVF = SUMVF + EP_G(IJK)
            ENDIF

! evaluate mixture bulk density rop_s(mmax) for GHD theory
            IF(KT_TYPE_ENUM == GHD_2007) THEN
              ROP_S(IJK,MMAX) = ZERO  ! mixture density
              DO M = 1, SMAX
                 IF (M /= MF) THEN
                   ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + RO_S(IJK,M)*EP_S(IJK,M)
                 ENDIF
              ENDDO
            ENDIF

! summing the solids volume fraction of all continuum solids phases
! except for the continuum solids phase which was marked by phase_4_p_g
! (skipping its continuity while solving gas continuity)
            DO M = 1, DES_MMAX+MMAX
               IF(KT_TYPE_ENUM == GHD_2007 .AND. M == MMAX) CYCLE
               IF (M /= MF) SUMVF = SUMVF + EP_S(IJK,M)
            ENDDO

            IF (0 == MF) THEN
! if no gas phase continuity was solved in the current cell then correct
! the void fraction of the gas phase based on the total solids volume
! fraction of all solids phases
               EP_G(IJK) = ONE - SUMVF

! Set error flag for negative volume fraction.
               IF (EP_G(IJK) < ZERO) THEN
                  Err_l(myPE) = 110
                  IF(REPORT_NEG_VOLFRAC) CALL EPgErr_LOG(IJK, wHeader)
               ENDIF

               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)
            ELSE
! else correct the volume fraction of the solids phase that was marked
               ROP_S(IJK,MF) = (ONE - SUMVF)*RO_S(IJK,MF)
            ENDIF
!-----------------------------------------------------------------<<<

         ENDIF    ! end if (fluid_at(ijk))
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      CALL send_recv(EP_G, 2)
      CALL send_recv(ROP_G, 2)
      CALL send_recv(ROP_S, 2)

      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)


      RETURN
      END SUBROUTINE CALC_VOL_FR

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_EP_FACTORS                                          C
!  Purpose: Set multiplication factors to ep_g to control the form     C
!  of the governing equations                                          C
!                                                                      C
!  References:                                                         C
!  Anderson and Jackson, A fluid mechanical description of fluidized   C
!     beds, Ind. Eng. Chem. Fundam., 6(4) 1967, 527-539.               C
!  Ishii, Thermo-fluid dynamic theory of two-phase flow, des Etudes    C
!     et Recherches D'Electricite de France, Eyrolles, Paris, France   C
!     1975.                                                            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_EP_FACTORS

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use fldvar, only: ep_g, ep_s
      use fldvar, only: epg_ifac, eps_ifac, epg_jfac
      use functions, only: wall_at
      use param1, only: one, undefined
      use physprop, only: mmax, mu_s0
      use run, only: jackson, ishii
      IMPLICIT NONE
! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ijk, m
!---------------------------------------------------------------------//
      epg_jfac(:) = one
      epg_ifac(:) = one
      eps_ifac(:,:) = one

      if (jackson) then
         do ijk = ijkstart3, ijkend3
            if (wall_at(ijk)) cycle
            epg_jfac(ijk) = ep_g(ijk)
         enddo
      endif
      if (ishii) then
         do ijk = ijkstart3, ijkend3
            if (wall_at(ijk)) cycle
            epg_ifac(ijk) = ep_g(ijk)
            do m = 1, mmax
! This seems more appropriate for non-granular systems
! So currently only incorporate this for 'constant' viscosity cases
               IF (mu_s0(m) /= undefined) THEN
                  eps_ifac(ijk,m) = ep_s(ijk,m)
               ENDIF
            enddo
         enddo
      endif

      RETURN
      END SUBROUTINE SET_EP_FACTORS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: EPgErr_LOG                                              !
!  Purpose: Record information about the location and conditions that  !
!           resulted in a negative gas phase volume fraction.          !
!                                                                      !
!  Author: J. Musser                                  Date: 16-APR-15  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE EPgErr_LOG(IJK, tHeader)

! Simulation time
      use run, only: TIME
! Gas phase pressure.
      use fldvar, only: EP_g, EP_s, ROP_s, RO_s, U_s, V_s, W_s
      use cutcell

      use physprop, only: SMAX
      use compar
      use indices
      use functions
      use error_manager


      IMPLICIT NONE

      INTEGER, intent(in) :: IJK
      LOGICAL, intent(inout) :: tHeader

      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.

      INTEGER :: M

      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('EPgErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('EPgErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)),                         &
            status='old', position='append')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5)") TIME
         tHeader = .FALSE.
      endif

      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",ADVANCE='NO') 'EP_g:', EP_g(IJK)
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('EP_s',M)), EP_s(IJK,M)
      ENDDO
      write(lUnit,"(' ')")

      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('RO_s',M)), RO_s(IJK,M)
      ENDDO
      write(lUnit,"(' ')")

      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('ROPs',M)), ROP_s(IJK,M)
      ENDDO
      write(lUnit,"(24x,' ')")


      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('U_se',M)), U_s(IJK,M)
      ENDDO
      write(lUnit,"(24x,' ')")

      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('V_sn',M)), V_s(IJK,M)
      ENDDO
      write(lUnit,"(24x,' ')")

      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('W_st',M)), W_s(IJK,M)
      ENDDO
      write(lUnit,"(24x,' ')")

      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('U_sw',M)), U_s(WEST_OF(IJK),M)
      ENDDO
      write(lUnit,"(24x,' ')")

      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('V_ss',M)), V_s(SOUTH_OF(IJK),M)
      ENDDO
      write(lUnit,"(24x,' ')")

      write(lUnit,"(24x)", ADVANCE='NO')
      DO M=1,SMAX
         write(lUnit,"(2x,A,1X,g12.5)", ADVANCE='NO') &
            trim(iVar('W_sb',M)), W_s(BOTTOM_OF(IJK),M)
      ENDDO
      write(lUnit,"(24x,' ')")


      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT('One or more cells have reported a negative gas volume ', &
         'fraction (EP_g).',/)

 1001 FORMAT(/4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE EPgErr_LOG
