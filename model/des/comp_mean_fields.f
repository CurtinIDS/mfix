!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: COMP_MEAN_FIELDS                                        !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose: Driver routine for calculating continuous field variables  !
!  corresponding to discrete data (ROP_s, u_s, v_s, w_s)               !
!                                                                      !
!  o The diffusion filter is only applied to the the solids bulk       !
!    density because DEM simulations do not utilize the other field    !
!    variables within a time loop.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE COMP_MEAN_FIELDS

! Modules
!---------------------------------------------------------------------//
      use discretelement, only: DES_MMAX
      use fldvar, only: rop_s
! Flag: Diffuse DES field variables.
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG
      use physprop, only: mmax
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//

! Loop counter.
      INTEGER :: M

!......................................................................!

! Calculate field variables from particle data:
      IF(DES_INTERP_MEAN_FIELDS) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_GARG) ; CALL COMP_MEAN_FIELDS0
         CASE DEFAULT; CALL COMP_MEAN_FIELDS1
         END SELECT
      ELSE
         CALL COMP_MEAN_FIELDS1
      ENDIF

! Apply the diffusion filter.
      IF(DES_DIFFUSE_MEAN_FIELDS) THEN
         DO M=MMAX+1, MMAX+DES_MMAX
            CALL DIFFUSE_MEAN_FIELD(ROP_S(:,M),'ROP_S')
         ENDDO
      ENDIF

! Calculate the gas phase volume fraction from ROP_s.
      CALL CALC_EPG_DES

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS
