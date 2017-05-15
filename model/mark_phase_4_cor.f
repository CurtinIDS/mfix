!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: MARK_PHASE_4_COR                                   C
!  Purpose: For each cell mark the phase whose continuity equation is  C
!           being skipped/activated (see conv_rop_g/conv_rop_s for     C
!           details) and whose volume fraction is then corrected (see  C
!           calc_vol_fr for details)                                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE MARK_PHASE_4_COR(PHASE_4_P_G, PHASE_4_P_S, DO_CONT,&
            MCP, DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE fldvar
      USE physprop
      USE constant
      USE compar
      USE visc_s
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! in each cell mark the phase index of the solids phase whose continuity
! equation becomes skipped but gas phase continuity is solved and whose
! volume fraction is then corrected based on sum of gas and other solids
! volume fractions. only in cells where the solids phase does not close
! pack and whose concentration exceeds the gas phase.
      INTEGER :: PHASE_4_P_g(DIMENSION_3)
! in each cell mark the phase index of the solids phase whose continuity
! equation becomes skipped and whose volume fraction is then corrected
! based on the maximum solids packing and sum of other solids volumes
! fractions that can become close packed. only in cells where maximum
! packing occurs.
      INTEGER :: PHASE_4_P_s(DIMENSION_3)
! flag whether continuity equation needs to be solved in addition to
! pressure correction equation
      LOGICAL :: DO_CONT(0:DIMENSION_M)
! index for the solids phase that can close pack.
      INTEGER :: MCP
! flag whether solids pressure correction is needed
      LOGICAL :: DO_P_s
! flag whether different phases were used for gas pressure correction
      LOGICAL :: SWITCH_4_P_g
! flag wether different phases were used for solids pressure correction
      LOGICAL :: SWITCH_4_P_s

! Here are some notes on the outcomes of the code following the changes
! by Sof on 11/22/10:
! -currently phase_4_p_g is always set to 0 in every cell.
! -currently phase_4_p_s is always undefined in every cell.
! -do_cont for the gas phase (M=0) is always F.
!  do_cont for the solids phases (M>1) is always F.
! -mcp is assigned the index of the solids phase that has close_packed=T
!  (i.e. solids phases that reach a maximum packing).  in cases
!  involving multiple solids phases that have close_packed=T, then this
!  will be assigned the solids phase with the lowest index. if no solids
!  phase has close_packed=T, or no solids phases exist (mmax=0) then it
!  remains undefined.
! -do_p_s is set to T if mcp is assigned, otherwise it is F.
! -switch_4_p_g is always set to F (see below for logic).
! -switch_4_p_s is always set to F (see below for logic).

! If this routine is reverted to the earlier code then some of the above,
! may change, namely:
! -phase_4_p_g could also be assigned the index of the gas phase in some
!  cells, and if a solids phase has close_packed=F, then the index of
!  that solids phase in other cells.
! -phase_4_p_s could be assigned the index of a solids phase that has
!  close_packed=T (MCP) if the cell exhibits close packing. if no cell
!  exhibits close packing or no solids phase has close_packed then, this
!  will always be undefined. note this code has always been commented.
! -do_cont for gas phase (M=0) is likely to be F, however, it could be
!  set to T if there exists a solids phase that has close_packed=F (MF)
!  and that has a greater mass concentration in any given fluid cell
!  than the fluid phase. do_cont for solids phases (M>1) will still
!  always be F.
! -switch_4_p_g will only be T if there exists a solids phase that has
!  close_packed=F (MF; i.e. does not close pack) and that has a greater
!  mass concentration in any  given fluid cell than the fluid phase.

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER          IJK, M
! Index for second continuous fluid phase
      INTEGER          MF
! Count number true values
      INTEGER          True_g, True_s
! Local check for whether pressure switches are made, that is, whether
! phase_4_p_g or phase_4_p_s are changed from their default setting.
! if true for that index, then switch.
      LOGICAL          SW_g(0:DIMENSION_M), SW_s(DIMENSION_M)
!-----------------------------------------------

! Initializiations
      MF = UNDEFINED_I
      MCP = UNDEFINED_I
      SW_G(0) = .FALSE.
      DO_CONT(0) = .TRUE.

! assigning MCP to the lowest solids phase index of those solids
! phases that have close_packed=T. assigning MF to the lowest solids
! phase index of those solids phases that do not close pack (i.e.,
! close_packed=F; the 'solids' phase that can overpack)
      DO M = MMAX, 1, -1
         IF (CLOSE_PACKED(M)) THEN
            MCP = M
         ELSE
            MF = M
         ENDIF
         SW_G(M) = .FALSE.
         SW_S(M) = .FALSE.
         DO_CONT(M) = .TRUE.
      ENDDO


! Sof 11/2010. : pressure correction equation is always solved for the
! gas-phase (i.e., M=0). The user must decide which phase to designate
! as the continuous phase (i.e., M=0). This is necessary for bubble
! column simulations to work. The new code effectively marks all cells
! the same way. The code below can be uncommented to revert back to
! the original implementation.

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

! Sof: uncomment to revert back to original implementation
! if one of the solids phases does not close pack (MF=defined) then
! assign phase_4_p_g in the current cell to either that 'solids' phase
! or the fluid phase (m=0) depending on which one has the greatest
! mass concentration in that cell.
!            IF (MF /= UNDEFINED_I) THEN
!               IF (EP_G(IJK)/RO_G(IJK) > EP_S(IJK,MF)/RO_S(MF)) THEN
!                  PHASE_4_P_G(IJK) = 0
!                  SW_G(0) = .TRUE.
!               ELSE
!                  PHASE_4_P_G(IJK) = MF
!                  SW_G(MF) = .TRUE.
!               ENDIF
!            ELSE

! now always marking phase_4_p_g to 0 (the fluid phase)
               PHASE_4_P_G(IJK) = 0
               SW_G(0) = .TRUE.
!            ENDIF


! if the current cell is close packed then assign phase_4_p_s to the
! the lowest solids phase index of those solids phases that can close
! pack (i.e. MCP)
!          IF(EP_g(IJK) .LE. EP_star_array(ijk)) THEN
!            PHASE_4_P_s(IJK) = MCP
!            IF(MCP .NE. UNDEFINED_I) SW_s(MCP) = .TRUE.
!          ELSE
            PHASE_4_P_S(IJK) = UNDEFINED_I       !to indicate no need for pressure correction
!          ENDIF
         ELSE
            PHASE_4_P_G(IJK) = UNDEFINED_I       !to indicate a non-fluid cell
            PHASE_4_P_S(IJK) = UNDEFINED_I       !to indicate a non-fluid cell
         ENDIF   ! end if/else fluid_at(ijk)
      ENDDO   ! end do ijk=ijkstart3,ijkend3


! setting the local values for true_g and true_s
      TRUE_G = 0
      TRUE_S = 0
      IF (SW_G(0)) TRUE_G = TRUE_G + 1
      DO M = 1, MMAX
! sw_g(m) will be true only if one of the solids phases (M>1) does not
! close pack (MF) and if it has a greater mass concentration in any given
! cell than the fluid phase (m=0).  so true_g is likely to be 1 but may
! be 2 if such an event occurs.
         IF (SW_G(M)) TRUE_G = TRUE_G + 1
! sw_s(m) will be true only if one of the solids phaes (M>1) does close
! pack (MCP) and any given cell shows close packed conditions. hence
! true_s is likely to be 1 in any simulation where close packing may
! occur. it is unclear how true_s can exceed 1 since sw_s is only
! switched to T for one solids index (mcp).  true_s will be 0 if no
! solids phase can close pack or if no cell exhibits close pack
! conditions.
         IF (SW_S(M)) TRUE_S = TRUE_S + 1
      ENDDO


! unlikely for switch_4_p_g to be true (see above).  switch_4_p_g will
! only be T if there exists a solids phase that does not close pack (MF)
! and that has a greater mass concentration in any given fluid cell than
! the fluid phase
      IF (TRUE_G > 1) THEN
         SWITCH_4_P_G = .TRUE.
      ELSE
         SWITCH_4_P_G = .FALSE.
      ENDIF

! unlikely for switch_4_p_s to be true since it does not appear that
! true_s can ever exceed 1 (see above).
      IF (TRUE_S > 1) THEN
         SWITCH_4_P_S = .TRUE.
      ELSE
         SWITCH_4_P_S = .FALSE.
      ENDIF

! MCP will only be undefined if none of the solids phases can close
! pack.  therefore if any solids phases can close pack do_p_s will be
! set to true.
      IF (MCP == UNDEFINED_I) THEN
         DO_P_S = .FALSE.
      ELSE
         DO_P_S = .TRUE.
      ENDIF

!      DO_P_s = .FALSE.
! true_s will be 0 if none of the solids phases (M>1) can close pack or
! if none of the cells exhibit close pack conditions.  otherwise if a
! solids phase can close pack and one of the cells exhibits close pack
! conditions then true_s will be 1.
!      IF(True_s .EQ. 0)THEN
!        DO_P_s = .FALSE.
!      ELSE
!        DO_P_s = .TRUE.
!      ENDIF

! if a phase was used for pressure correction and no other phases were
! used for pressure correction, there is no need to solve its continuity
! equation.
      IF (SW_G(0) .AND. TRUE_G==1) DO_CONT(0) = .FALSE.
!      DO M = 1, MMAX
!         IF(SW_g(M) .AND. True_g .EQ. 1)DO_CONT(M) = .FALSE.
!         IF(SW_s(M) .AND. True_s .EQ. 1)DO_CONT(M) = .FALSE.
!      ENDDO


      RETURN
      END SUBROUTINE MARK_PHASE_4_COR


