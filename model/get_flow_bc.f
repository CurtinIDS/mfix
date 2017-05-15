!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GET_FLOW_BC                                             C
!  Purpose: Find and validate i, j, k locations for flow BC's. Also    C
!           set value of bc_plane for flow BC's.                       C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_X_w, BC_X_e, BC_Y_s, BC_Y_n    C
!                        BC_Z_b, BC_Z_t, DX, DY, DZ, IMAX, JMAX, KMAX  C
!                                                                      C
!  Variables modified: BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t  C
!                      ICBC_FLAG, BC_PLANE                             C
!                                                                      C
!  Local variables: BC, I, J, K, IJK, I_w, I_e, J_s, J_n, K_b, K_t     C
!                   ERROR, X_CONSTANT, Y_CONSTANT, Z_CONSTANT          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_FLOW_BC

      IMPLICIT NONE


! Stub file: TBR/JM

      RETURN
      END SUBROUTINE GET_FLOW_BC

