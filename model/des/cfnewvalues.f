!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!
!                                                                      C
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!  pradeep : changes for parallel processing
!          1. periodic boundaries might lie in different proc. so adjust
!             particle position for periodic removed
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE derived_types, only: multisap, boxhandle
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      use geometry, only: DO_K, NO_K
      use multi_sweep_and_prune, only: aabb_t, multisap_sort, multisap_update

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: OMEGA_MAG,OMEGA_UNIT(3),ROT_ANGLE

      type(aabb_t) aabb

!-----------------------------------------------

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(IS_NONEXISTENT(L)) CYCLE                       ! Only real particles
            IF(IS_ENTERING(L).or.IS_ENTERING_GHOST(L)) CYCLE  ! Only non-entering
            IF(IS_GHOST(L)) CYCLE                             ! Skip ghost particles
            DES_ACC_OLD(L,:) = FC(L,:)/PMASS(L) + GRAV(:)
            ROT_ACC_OLD(L,:) = TOW(L,:)
         ENDDO
      ENDIF


!!$omp parallel default(none)                    &
!!$omp shared(MAX_PIP,INTG_EULER,INTG_ADAMS_BASHFORTH,fc,tow,do_nsearch,   &
!!$omp       omega_new,omega_old,pmass,grav,des_vel_new,des_pos_new,       &
!!$omp       des_vel_old,des_pos_old,dtsolid,omoi,des_acc_old,rot_acc_old, &
!!$omp       ppos,neighbor_search_rad_ratio,des_radius,DO_OLD, iGlobal_ID, &
!!$omp       particle_orientation,orientation,boxhandle,multisap,particle_state) &
!!$omp private(l,rot_angle,omega_mag,omega_unit,aabb)

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f

! Advance particle position, velocity
! first-order method
      IF (INTG_EULER) THEN
!!$omp sections
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,1) =   &
            DES_VEL_NEW(:,1) + DTSOLID*(FC(:,1)/PMASS(:) + GRAV(1))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,1) =       &
            DES_POS_NEW(:,1) + DES_VEL_NEW(:,1)*DTSOLID
         FC(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,2) =   &
            DES_VEL_NEW(:,2) + DTSOLID*(FC(:,2)/PMASS(:) + GRAV(2))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,2) =       &
            DES_POS_NEW(:,2) + DES_VEL_NEW(:,2)*DTSOLID
         FC(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,3) =   &
            DES_VEL_NEW(:,3) + DTSOLID*(FC(:,3)/PMASS(:) + GRAV(3))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,3) =       &
            DES_POS_NEW(:,3) + DES_VEL_NEW(:,3)*DTSOLID
         FC(:,3) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) OMEGA_NEW(:,1) =     &
            OMEGA_NEW(:,1) + TOW(:,1)*OMOI(:)*DTSOLID
         TOW(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE)OMEGA_NEW(:,2) =      &
            OMEGA_NEW(:,2) + TOW(:,2)*OMOI(:)*DTSOLID
         TOW(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE)OMEGA_NEW(:,3) =      &
            OMEGA_NEW(:,3) + TOW(:,3)*OMOI(:)*DTSOLID
         TOW(:,3) = ZERO
!!$omp end sections

! Second-order Adams-Bashforth/Trapezoidal scheme
      ELSEIF (INTG_ADAMS_BASHFORTH) THEN

!!$omp sections
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,1) = FC(:MAX_PIP,1)/PMASS(:MAX_PIP) + GRAV(1)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,1) = DES_VEL_OLD(:MAX_PIP,1) + 0.5d0* &
               (3.d0*FC(:MAX_PIP,1)-DES_ACC_OLD(:MAX_PIP,1) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,1) = FC(:MAX_PIP,1)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,1) = DES_POS_OLD(:MAX_PIP,1) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,1)+DES_VEL_NEW(:MAX_PIP,1))*DTSOLID
         FC(:MAX_PIP,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,2) = FC(:MAX_PIP,2)/PMASS(:MAX_PIP) + GRAV(2)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,2) = DES_VEL_OLD(:MAX_PIP,2) + 0.5d0* &
               (3.d0*FC(:MAX_PIP,2)-DES_ACC_OLD(:MAX_PIP,2) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,2) = FC(:MAX_PIP,2)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,2) = DES_POS_OLD(:MAX_PIP,2) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,2)+DES_VEL_NEW(:MAX_PIP,2))*DTSOLID
         FC(:MAX_PIP,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,3) = FC(:MAX_PIP,3)/PMASS(:MAX_PIP) + GRAV(3)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,3) = DES_VEL_OLD(:MAX_PIP,3) + 0.5d0* &
                 (3.d0*FC(:MAX_PIP,3)-DES_ACC_OLD(:MAX_PIP,3) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,3) = FC(:MAX_PIP,3)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,3) = DES_POS_OLD(:MAX_PIP,3) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,3)+DES_VEL_NEW(:MAX_PIP,3))*DTSOLID
         FC(:MAX_PIP,3) = ZERO

!!$omp section
        WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
           OMEGA_NEW(:MAX_PIP,1) = OMEGA_OLD(:MAX_PIP,1) + 0.5d0*     &
               (3.d0*TOW(:MAX_PIP,1)*OMOI(:MAX_PIP) -                  &
               ROT_ACC_OLD(:MAX_PIP,1))*DTSOLID
        WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,1) = TOW(:MAX_PIP,1)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            OMEGA_NEW(:MAX_PIP,2) = OMEGA_OLD(:MAX_PIP,2) + 0.5d0*     &
                 (3.d0*TOW(:MAX_PIP,2)*OMOI(:MAX_PIP)-                 &
                 ROT_ACC_OLD(:MAX_PIP,2) )*DTSOLID
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,2) = TOW(:MAX_PIP,2)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,2) = ZERO

!!$omp section
        WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            OMEGA_NEW(:MAX_PIP,3) = OMEGA_OLD(:MAX_PIP,3) + 0.5d0*     &
               (3.d0*TOW(:MAX_PIP,3)*OMOI(:MAX_PIP)-                   &
               ROT_ACC_OLD(:MAX_PIP,3) )*DTSOLID
        WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,3) = TOW(:MAX_PIP,3)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,3) = ZERO

!!$omp end sections
      ENDIF

#ifdef do_sap
!!$omp single
         DO L = 1, MAX_PIP
            aabb%minendpoint(:) = DES_POS_NEW(L,:)-DES_RADIUS(L)
            aabb%maxendpoint(:) = DES_POS_NEW(L,:)+DES_RADIUS(L)
            call multisap_update(multisap,aabb,boxhandle(L))
         ENDDO
!!$omp end single
#endif

! Update particle orientation - Always first order
! When omega is non-zero, compute the rotation angle, and apply the
! Rodrigues' rotation formula

      IF(PARTICLE_ORIENTATION) THEN
         DO L = 1, MAX_PIP
            OMEGA_MAG = OMEGA_NEW(L,1)**2 +OMEGA_NEW(L,2)**2 + OMEGA_NEW(L,3)**2

            IF(OMEGA_MAG>ZERO) THEN
               OMEGA_MAG=DSQRT(OMEGA_MAG)
               OMEGA_UNIT(:) = OMEGA_NEW(L,:)/OMEGA_MAG
               ROT_ANGLE = OMEGA_MAG * DTSOLID

               ORIENTATION(:,L) = ORIENTATION(:,L)*DCOS(ROT_ANGLE) + &
                  DES_CROSSPRDCT(OMEGA_UNIT,ORIENTATION(:,L))*DSIN(ROT_ANGLE) + &
                  OMEGA_UNIT(:)*DOT_PRODUCT(OMEGA_UNIT,ORIENTATION(:,L))*&
                  (ONE-DCOS(ROT_ANGLE))
            ENDIF
         ENDDO
      ENDIF

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
      IF(.NOT.DO_NSEARCH) THEN
!!$omp do reduction (.or.:do_nsearch)
         DO L = 1, MAX_PIP
            DO_NSEARCH = DO_NSEARCH .or. &
               (DES_POS_NEW(L,1) - PPOS(L,1))**2+              &
               (DES_POS_NEW(L,2) - PPOS(L,2))**2+              &
               (DES_POS_NEW(L,3) - PPOS(L,3))**2  >=           &
               (NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L))**2
         ENDDO
      ENDIF

!!$omp end parallel

#ifdef do_sap
            call multisap_sort(multisap)
#endif

      FIRST_PASS = .FALSE.

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)

      RETURN

      contains

      include 'functions.inc'

      END SUBROUTINE CFNEWVALUES
