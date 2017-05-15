!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFSLIDE(V_TANG, PARTICLE_SLIDE, MU)
!  Purpose:  Check for Coulombs friction law - calculate sliding
!            friction
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer:                                          Date:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFSLIDE(V_TANG, PARTICLE_SLIDE, MU, FT_tmp, FN_tmp)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement, only: DEBUG_DES
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! tangent to the plane of contact
      DOUBLE PRECISION, INTENT(IN) :: V_TANG(3)
! logic set to T when a sliding contact occurs
      LOGICAL, INTENT(OUT) :: PARTICLE_SLIDE
! Coefficient of friction
      DOUBLE PRECISION, INTENT(IN) :: MU
! normal force
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: FN_tmp
! tangential force
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT) :: FT_tmp
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! squared magnitude of tangential and normal forces
      DOUBLE PRECISION FTMD, FNMD
!-----------------------------------------------

      FTMD = dot_product(FT_tmp(:),FT_tmp(:))
      FNMD = dot_product(FN_tmp(:),FN_tmp(:))

      IF (FTMD.GT.(MU*MU*FNMD)) THEN
! tangential force based on sliding friction
         PARTICLE_SLIDE = .TRUE.
         IF(ALL(V_TANG.EQ.0)) THEN
            FT_tmp(:) =  MU * FT_tmp(:) * SQRT(FNMD/FTMD)
         ELSE
            FT_tmp(:) = -MU * V_TANG(:) * SQRT(FNMD/dot_product(V_TANG,V_TANG))
         ENDIF

         IF(DEBUG_DES) THEN
            WRITE(*,'(7X,A)') &
                 'FROM CFSLIDE.F ---------->'
            WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
            WRITE(*,'(9X,A,2(ES15.7,1X))')&
                 'FTMD, mu*FNMD = ', FTMD, MU*FNMD
            WRITE(*,'(7X,A)') '<----------END CFSLIDE.F'
         ENDIF

      ENDIF

      RETURN
      END SUBROUTINE CFSLIDE
