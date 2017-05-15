!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COLLISION_WALL                                    C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!                                                                      C
!  Purpose: subroutines for particle-wall collisions when cutcell is   C
!           used. Also contains rehack of routines for cfslide and     C
!           cffctow which might be different from the stand alone      C
!           routines. Eventually all the DEM routines will be          C
!           consolidated.                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE CALC_COLLISION_WALL

      PRIVATE
      PUBLIC :: CALC_DEM_FORCE_WITH_WALL_STL,&
      & CALC_DEM_THERMO_WITH_WALL_STL

      CONTAINS

!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

      USE run
      USE param1
      USE desgrid
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE indices
      USE stl
      USE stl_functions_des
      USE functions
      Implicit none

      INTEGER :: LL
      INTEGER :: IJK, NF
      DOUBLE PRECISION ::OVERLAP_N, SQRT_OVERLAP

      DOUBLE PRECISION :: V_REL_TRANS_NORM, DISTSQ, RADSQ, CLOSEST_PT(DIMN)
! local normal and tangential forces
      DOUBLE PRECISION :: NORMAL(DIMN), VREL_T(DIMN), DIST(DIMN), DISTMOD
      DOUBLE PRECISION, DIMENSION(DIMN) :: FTAN, FNORM, OVERLAP_T

      LOGICAL :: DES_LOC_DEBUG
      INTEGER :: CELL_ID, cell_count
      INTEGER :: PHASELL

      DOUBLE PRECISION :: TANGENT(DIMN)
      DOUBLE PRECISION :: FTMD, FNMD
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      double precision :: MAG_OVERLAP_T

      double precision :: line_t
! flag to tell if the orthogonal projection of sphere center to
! extended plane detects an overlap

      DOUBLE PRECISION :: MAX_DISTSQ, DISTAPART, FORCE_COH, R_LM
      INTEGER :: MAX_NF, axis
      DOUBLE PRECISION, DIMENSION(3) :: PARTICLE_MIN, PARTICLE_MAX, POS_TMP

      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false.
      FOCUS_PARTICLE = -1

!$omp parallel default(none) private(LL,ijk,MAG_OVERLAP_T,             &
!$omp    cell_id,radsq,particle_max,particle_min,tangent,              &
!$omp    axis,nf,closest_pt,dist,r_lm,distapart,force_coh,distsq,      &
!$omp    line_t,max_distsq,max_nf,normal,distmod,overlap_n,VREL_T,     &
!$omp    v_rel_trans_norm,phaseLL,sqrt_overlap,kn_des_w,kt_des_w,      &
!$omp    etan_des_w,etat_des_w,fnorm,overlap_t,ftan,ftmd,fnmd,pos_tmp) &
!$omp shared(max_pip,focus_particle,debug_des,                         &
!$omp    pijk,dg_pijk,i_of,j_of,k_of,des_pos_new,    &
!$omp    des_radius,facets_at_dg,vertex,  &
!$omp    hert_kwn,hert_kwt,kn_w,kt_w,des_coll_model_enum,mew_w,tow,    &
!$omp    des_etan_wall,des_etat_wall,dtsolid,fc,norm_face,             &
!$omp    wall_collision_facet_id,wall_collision_PFT,use_cohesion,      &
!$omp    van_der_waals,wall_hamaker_constant,wall_vdw_outer_cutoff,    &
!$omp    wall_vdw_inner_cutoff,asperities,surface_energy)
!$omp do
      DO LL = 1, MAX_PIP

         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.

! skipping non-existent particles or ghost particles
! make sure the particle is not classified as a new 'entering' particle
! or is already marked as a potential exiting particle
         IF(.NOT.IS_NORMAL(LL)) CYCLE

         CELL_ID = DG_PIJK(LL)

! If no neighboring facet in the surrounding 27 cells, then exit
         IF(facets_at_dg(CELL_ID)%COUNT < 1) THEN
            WALL_COLLISION_FACET_ID(:,LL) = -1
            WALL_COLLISION_PFT(:,:,LL) = 0.0d0
            CYCLE
         ENDIF

! Check particle LL for wall contacts
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)

         particle_max(:) = des_pos_new( LL,:) + des_radius(LL)
         particle_min(:) = des_pos_new( LL,:) - des_radius(LL)

         DO CELL_COUNT = 1, facets_at_dg(cell_id)%count

            axis = facets_at_dg(cell_id)%dir(cell_count)

            NF = facets_at_dg(cell_id)%id(cell_count)

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN

               CALL ClosestPtPointTriangle(DES_POS_NEW(LL,:),          &
                  VERTEX(:,:,NF), CLOSEST_PT(:))

               DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
               DISTSQ = DOT_PRODUCT(DIST, DIST)
               R_LM = 2*DES_RADIUS(LL)

               IF(DISTSQ < (R_LM+WALL_VDW_OUTER_CUTOFF)**2) THEN
                  IF(DISTSQ > (WALL_VDW_INNER_CUTOFF+R_LM)**2) THEN
                     DistApart = (SQRT(DISTSQ)-R_LM)
                     FORCE_COH = WALL_HAMAKER_CONSTANT*DES_RADIUS(LL) /&
                        (12d0*DistApart**2)*(Asperities/(Asperities +  &
                        DES_RADIUS(LL)) + ONE/(ONE+Asperities/         &
                        DistApart)**2)
                  ELSE
                     FORCE_COH = 2d0*PI*SURFACE_ENERGY*DES_RADIUS(LL)* &
                        (Asperities/(Asperities+DES_RADIUS(LL)) + ONE/ &
                        (ONE+Asperities/WALL_VDW_INNER_CUTOFF)**2 )
                  ENDIF
                  FC(LL,:) = FC(LL,:) + DIST(:)*FORCE_COH/SQRT(DISTSQ)
               ENDIF
            ENDIF

            if (facets_at_dg(cell_id)%min(cell_count) >    &
               particle_max(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

            if (facets_at_dg(cell_id)%max(cell_count) <    &
               particle_min(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

! Checking all the facets is time consuming due to the expensive
! separating axis test. Remove this facet from contention based on
! a simple orthogonal projection test.

! Parametrize a line as p = p_0 + t normal and intersect with the
! triangular plane. If t>0, then point is on the non-fluid side of
! the plane. If the plane normal is assumed to point toward the fluid.

! -undefined, because non zero values will imply the sphere center
! is on the non-fluid side of the plane. Since the testing
! is with extended plane, this could very well happen even
! when the particle is well inside the domain (assuming the plane
! normal points toward the fluid). See the pic below. So check
! only when line_t is negative

!                 \   Solid  /
!                  \  Side  /
!                   \      /
!                    \    /
! Wall 1, fluid side  \  /  Wall 2, fluid side
!                      \/
!                        o particle
!
! line_t will be positive for wall 1 (incorrectly indicating center
! is outside the domain) and line_t will be negative for wall 2.
!
! Therefore, only stick with this test when line_t is negative and let
! the separating axis test take care of the other cases.

! Since this is for checking static config, line's direction is the
! same as plane's normal. For moving particles, the line's normal will
! be along the point joining new and old positions.

            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(LL,:),&
               NORM_FACE(:,NF))

! k - rad >= tol_orth, where k = -line_t, then orthogonal
! projection is false. Substituting for k
! => line_t + rad <= -tol_orth
! choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius

! Orthogonal projection will detect false positives even
! when the particle does not overlap the triangle.
! However, if the orthogonal projection shows no overlap, then
! that is a big fat negative and overlaps are not possible.
            if((line_t.le.-1.0001d0*des_radius(LL))) then  ! no overlap
               call remove_collision(LL,nf,wall_collision_facet_id)
               CYCLE
            ENDIF

            pos_tmp = DES_POS_NEW(LL,:)
            CALL ClosestPtPointTriangle(pos_tmp,             &
               VERTEX(:,:,NF), CLOSEST_PT(:))

            DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
            DISTSQ = DOT_PRODUCT(DIST, DIST)

            IF(DISTSQ .GE. RADSQ - SMALL_NUMBER) THEN !No overlap exists
               call remove_collision(LL,nf,wall_collision_facet_id)
               CYCLE
            ENDIF

            MAX_DISTSQ = DISTSQ
            MAX_NF = NF

! Assign the collision normal based on the facet with the
! largest overlap.
            NORMAL(:) = DIST(:)/sqrt(DISTSQ)

! Facet's normal is correct normal only when the intersection is with
! the face. When the intersection is with edge or vertex, then the
! normal is based on closest pt and sphere center. The definition above
! of the normal is generic enough to account for differences between
! vertex, edge, and facet.

! Calculate the particle/wall overlap.
            DISTMOD = SQRT(MAX_DISTSQ)
            OVERLAP_N = DES_RADIUS(LL) - DISTMOD

! Calculate the translational relative velocity
            CALL CFRELVEL_WALL(LL, V_REL_TRANS_NORM,VREL_T,           &
               NORMAL, DISTMOD)

! Calculate the spring model parameters.
            phaseLL = PIJK(LL,5)

! Hertz vs linear spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap
               KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt_overlap
               ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt_overlap
            ELSE
               KN_DES_W = KN_W
               KT_DES_W = KT_W
               ETAN_DES_W = DES_ETAN_WALL(phaseLL)
               ETAT_DES_W = DES_ETAT_WALL(phaseLL)
            ENDIF

! Calculate the normal contact force
            FNORM(:) = -(KN_DES_W * OVERLAP_N * NORMAL(:) + &
               ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:))

! Calculate the tangential displacement.
            OVERLAP_T(:) = DTSOLID*VREL_T(:) + GET_COLLISION(LL,       &
               NF, WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)
            MAG_OVERLAP_T = sqrt(DOT_PRODUCT(OVERLAP_T, OVERLAP_T))

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Tangential froce from spring.
               FTMD = KT_DES_W*MAG_OVERLAP_T
! Max force before the on set of frictional slip.
               FNMD = MEW_W*sqrt(DOT_PRODUCT(FNORM,FNORM))
! Direction of tangential force.
               TANGENT = OVERLAP_T/MAG_OVERLAP_T
               IF(FTMD < FNMD) THEN
                  FTAN = -FTMD * TANGENT
               ELSE
                  FTAN = -FNMD * TANGENT
                  OVERLAP_T = (FNMD/KT_DES_W) * TANGENT
               ENDIF
            ELSE
               FTAN = 0.0
            ENDIF
! Add in the tangential dashpot damping force
            FTAN = FTAN - ETAT_DES_W*VREL_T(:)

! Save the tangential displacement.
            CALL UPDATE_COLLISION(OVERLAP_T, LL, NF,                   &
               WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

! Add the collision force to the total forces acting on the particle.
            FC(LL,:) = FC(LL,:) + FNORM(:) + FTAN(:)

! Add the torque force to the toal torque acting on the particle.
            TOW(LL,:) = TOW(LL,:) + DISTMOD*DES_CROSSPRDCT(NORMAL,FTAN)

         ENDDO

      ENDDO
!$omp end do
!$omp end parallel

      RETURN

       contains

!......................................................................!
!  Function: GET_COLLISION                                             !
!                                                                      !
!  Purpose: Return the integrated (t0->t) tangential displacement.     !
!......................................................................!
      FUNCTION GET_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID,     &
          WALL_COLLISION_PFT)

      use stl_dbg_des, only: write_this_stl
      use stl_dbg_des, only: write_stls_this_dg

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION :: GET_COLLISION(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC, FREE_INDEX, LC, dgIJK

      INTEGER :: lSIZE1, lSIZE2, lSIZE3
      INTEGER, ALLOCATABLE :: tmpI2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: tmpR3(:,:,:)

      free_index = -1

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            get_collision(:) = wall_collision_PFT(:,cc,LLL)
            return
         else if (-1 == wall_collision_facet_id(cc,LLL)) then
            free_index = cc
         endif
      enddo

! Overwrite old data. This is needed because a particle moving from
! one dg cell to another may no longer 'see' an STL before it moved
! out of contact range. Therefore, the 'remove_collision' function
! does not get called to cleanup the stale data.
      if(-1 == free_index) then
         dgIJK=DG_PIJK(LLL)
         cc_lp: do cc=1, COLLISION_ARRAY_MAX
            do lc=1, facets_at_dg(dgIJK)%count
               if(wall_collision_facet_id(cc,LLL) == &
                  facets_at_dg(dgIJK)%id(LC))  cycle cc_lp
            enddo
            free_index = cc
            exit cc_lp
         enddo cc_lp
      endif

! Last resort... grow the collision array
      if(-1 == free_index) then
         free_index=COLLISION_ARRAY_MAX+1
         COLLISION_ARRAY_MAX = 2*COLLISION_ARRAY_MAX
         CALL GROW_WALL_COLLISION(COLLISION_ARRAY_MAX)
      endif

      wall_collision_facet_id(free_index,LLL) = facet_id
      wall_collision_PFT(:,free_index,LLL) = ZERO
      get_collision(:) = wall_collision_PFT(:,free_index,LLL)
      return

      END FUNCTION GET_COLLISION


!......................................................................!
!  Subroutine: GROW_WALL_COLLISION                                     !
!                                                                      !
!  Purpose: Return the integrated (t0->t) tangential displacement.     !
!......................................................................!
      SUBROUTINE GROW_WALL_COLLISION(NEW_SIZE)

      use discretelement

      use stl_dbg_des, only: write_this_stl
      use stl_dbg_des, only: write_stls_this_dg

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NEW_SIZE
      INTEGER :: lSIZE1, lSIZE2, lSIZE3
      INTEGER, ALLOCATABLE :: tmpI2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: tmpR3(:,:,:)

      lSIZE1 = size(wall_collision_facet_id,1)
      lSIZE2 = size(wall_collision_facet_id,2)

      allocate(tmpI2(NEW_SIZE, lSIZE2))
      tmpI2(1:lSIZE1,:) = WALL_COLLISION_FACET_ID(1:lSIZE1,:)
      call move_alloc(tmpI2, WALL_COLLISION_FACET_ID)
      WALL_COLLISION_FACET_ID(lSIZE1+1:NEW_SIZE,:) = -1

      lSIZE1 = size(wall_collision_pft,1)
      lSIZE2 = size(wall_collision_pft,2)
      lSIZE3 = size(wall_collision_pft,3)

      allocate(tmpR3(lSIZE1, NEW_SIZE, lSIZE3))
      tmpR3(:,1:lSIZE2,:) = WALL_COLLISION_PFT(:,1:lSIZE2,:)
      call move_alloc(tmpR3, WALL_COLLISION_PFT)

      RETURN
      END SUBROUTINE GROW_WALL_COLLISION




!......................................................................!
!  Function: UPDATE_COLLISION                                          !
!                                                                      !
!  Purpose: Update the integrated (t0->t) tangential displacement.     !
!......................................................................!
      SUBROUTINE UPDATE_COLLISION(PFT, LLL, FACET_ID,                  &
         WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

      use error_manager
      implicit none

      DOUBLE PRECISION, INTENT(IN) :: PFT(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(IN) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            wall_collision_PFT(:,cc,LLL) = PFT(:)
            return
         endif
      enddo

      CALL INIT_ERR_MSG("CALC_COLLISION_WALL_MOD: UPDATE_COLLISION")
      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error: COLLISION_ARRAY_MAX too small. ')

      END SUBROUTINE UPDATE_COLLISION

!......................................................................!
!  Function: REMOVE_COLLISION                                          !
!                                                                      !
!  Purpose: Clear the integrated (t0->t) tangential displacement once  !
!  the collision is over (contact ended).                              !
!......................................................................!
      SUBROUTINE REMOVE_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID)

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      INTEGER :: CC

      DO CC = 1, COLLISION_ARRAY_MAX
         IF (FACET_ID == WALL_COLLISION_FACET_ID(CC,LLL)) THEN
            WALL_COLLISION_FACET_ID(CC,LLL) = -1
            RETURN
         ENDIF
      ENDDO

      END SUBROUTINE REMOVE_COLLISION

      END SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CFRELVEL_WALL                                           !
!                                                                      !
!  Purpose: Calculate the normal and tangential components of the      !
!  relative velocity between a particle and wall contact.              !
!                                                                      !
!  Comments: Only the magnitude of the normal component is returned    !
!  whereas the full tangential vector is returned.                     !
!----------------------------------------------------------------------!
      SUBROUTINE CFRELVEL_WALL(LL, VRN, VRT, NORM, DIST)

! Particle translational velocity
      use discretelement, only: DES_VEL_NEW
! Particle rotational velocity
      use discretelement, only: OMEGA_NEW
! Spatial array size (parameter)
      use discretelement, only: DIMN
! Function for calculating the cross prodcut
      use discretelement, only: DES_CROSSPRDCT

      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Particle index.
      INTEGER, INTENT(IN) :: LL
! Magnitude of the total relative translational velocity.
      DOUBLE PRECISION, INTENT(OUT):: VRN
! Total relative translational velocity (vector).
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(OUT):: VRT
! Unit normal from particle center to closest point on stl (wall)
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: NORM
! Distance between particle center and stl (wall).
      DOUBLE PRECISION, INTENT(IN) :: DIST

! Local variables
!---------------------------------------------------------------------//
! Additional relative translational motion due to rotation
      DOUBLE PRECISION, DIMENSION(DIMN) :: V_ROT
! Total relative velocity at contact point
      DOUBLE PRECISION, DIMENSION(DIMN) :: VRELTRANS

! Total relative velocity + rotational contribution
      V_ROT = DIST*OMEGA_NEW(LL,:)
      VRELTRANS(:) =  DES_VEL_NEW(LL,:) + DES_CROSSPRDCT(V_ROT, NORM)

! magnitude of normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! total relative translational slip velocity at the contact point
! Equation (8) in Tsuji et al. 1992
      VRT(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL


!----------------------------------------------------------------//
! SUBROUTINE: CALC_DEM_THERMO_WITH_WALL_STL
! By: Aaron M.
! Purpose: Compute heat transfer to particles due to boundaries
!----------------------------------------------------------------//

      SUBROUTINE CALC_DEM_THERMO_WITH_WALL_STL

      USE bc
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE cutcell, only: cut_cell_at, area_cut, blocked_cell_at, small_cell_at
      USE des_thermo
      USE des_thermo_cond
      USE desgrid
      USE discretelement
      USE exit, only: mfix_exit
      USE functions
      USE geometry, only: do_k, imax1, imax2, imin1, imin2
      USE geometry, only: no_k, jmax1, jmax2, jmin1, jmin2
      USE geometry, only: kmax1, kmax2, kmin1, kmin2, zlength
      USE geometry, only: dx, dy, dz
      USE param, only: dimension_i, dimension_j, dimension_k
      USE param1
      USE physprop
      USE run, only: units, time
      USE stl
      USE stl_functions_des

      IMPLICIT NONE

      INTEGER :: LL             ! Loop index for particle
      INTEGER :: CELL_ID        ! Desgrid cell index
      INTEGER :: CELL_COUNT     ! Loop index for facets
      INTEGER :: IJK_FLUID      ! IJK index for fluid cell
      INTEGER :: I_FLUID, J_FLUID, K_FLUID
      INTEGER :: I_Facet, J_Facet, K_Facet, IJK_Facet
      INTEGER :: I1,J1,K1
      INTEGER :: phase_LL       ! Phase index for particle

      INTEGER, PARAMETER :: MAX_CONTACTS = 12
      INTEGER :: iFacet         ! loop index for facets
      INTEGER :: count_Facets   ! counter for number of facets particle contacts
      DOUBLE PRECISION, DIMENSION(3,MAX_CONTACTS) :: NORM_FAC_CONTACT
      LOGICAL :: USE_FACET
      LOGICAL :: DOMAIN_BDRY

      INTEGER :: NF             ! Facet ID
      INTEGER :: AXIS           ! Facet direction
      DOUBLE PRECISION :: RLENS_SQ ! lens radius squared
      DOUBLE PRECISION :: RLENS ! lens radius
      DOUBLE PRECISION :: RPART ! Particle radius

      ! vectors for minimum / maximum possible particle contacts
      DOUBLE PRECISION, DIMENSION(3) :: PARTPOS_MIN, PARTPOS_MAX
      DOUBLE PRECISION, DIMENSION(3) :: POS_TMP ! Temp vector
      DOUBLE PRECISION, DIMENSION(3) :: DIST, NORMAL

      DOUBLE PRECISION, DIMENSION(3) :: CLOSEST_PT

      DOUBLE PRECISION :: line_t  ! Normal projection for part/wall
      DOUBLE PRECISION :: DISTSQ ! Separation from particle and facet (squared)
      DOUBLE PRECISION :: PROJ

      INTEGER :: BC_ID ! BC ID
      INTEGER :: IBC   ! BC Loop Index

      DOUBLE PRECISION :: TWALL
      DOUBLE PRECISION :: K_gas
      DOUBLE PRECISION :: TPART
      DOUBLE PRECISION :: OVERLAP
      DOUBLE PRECISION :: QSWALL, AREA

      LOGICAL, SAVE :: OUTPUT_WARNING = .TRUE.



    !  IF(.NOT.DES_CONTINUUM_COUPLED.or.DES_EXPLICITLY_COUPLED)THEN
    !     CALL PARTICLES_IN_CELL
    !  ENDIF
      DO LL = 1, MAX_PIP
         ! Skip non-existent particles or ghost particles
         IF (.NOT.IS_NORMAL(LL)) CYCLE
         PHASE_LL = PIJK(LL,5)
         IF(.NOT.CALC_COND_DES(PHASE_LL))CYCLE
         ! Get desgrid cell index
         CELL_ID = DG_PIJK(LL)

         ! Skip cells that do not have neighboring facet
         IF (facets_at_dg(CELL_ID)%COUNT <1) CYCLE

         ! Store lens radius of particle
         RLENS = (ONE+FLPC)*DES_RADIUS(LL)
         RLENS_SQ = RLENS*RLENS

         RPART = DES_RADIUS(LL)

         ! Compute max/min particle locations
         PARTPOS_MAX(:) = des_pos_new(LL,:) + RLENS
         PARTPOS_MIN(:) = des_pos_new(LL,:) - RLENS

         ! Get fluid cell
         I_FLUID = PIJK(LL,1)
         J_FLUID = PIJK(LL,2)
         K_FLUID = PIJK(LL,3)
         IJK_FLUID = PIJK(LL,4)
         TPART = DES_T_S(LL)

         ! Sometimes PIJK is not updated every solids timestep and
         ! therefore doesn't give the correct fluid cell that
         ! contains the particles.  If so, determine actual fluid
         ! cell that contains the particle
         IF(.NOT.DES_CONTINUUM_COUPLED.or.DES_EXPLICITLY_COUPLED)THEN
            IF(I_FLUID <= ISTART3 .OR. I_FLUID >= IEND3) THEN
               CALL PIC_SEARCH(I_FLUID, DES_POS_NEW(LL,1), XE,       &
               DIMENSION_I, IMIN2, IMAX2)
            ELSE
               IF((DES_POS_NEW(LL,1) >= XE(I_FLUID-1)) .AND.         &
               (DES_POS_NEW(LL,1) <  XE(I_FLUID))) THEN
               I_FLUID = I_FLUID
               ELSEIF((DES_POS_NEW(LL,1) >= XE(I_FLUID)) .AND.       &
                  (DES_POS_NEW(LL,1) < XE(I_FLUID+1))) THEN
                  I_FLUID = I_FLUID+1
               ELSEIF((DES_POS_NEW(LL,1) >= XE(I_FLUID-2)) .AND.     &
                  (DES_POS_NEW(LL,1) < XE(I_FLUID-1))) THEN
                  I_FLUID = I_FLUID-1
               ELSE
                  CALL PIC_SEARCH(I_FLUID, DES_POS_NEW(LL,1), XE,    &
                  DIMENSION_I, IMIN2, IMAX2)
               ENDIF
            ENDIF !(I)
            IF(J_FLUID <= JSTART3 .OR. J_FLUID >= JEND3) THEN
               CALL PIC_SEARCH(J_FLUID, DES_POS_NEW(LL,2), YN,       &
               DIMENSION_J, JMIN2, JMAX2)
            ELSE
               IF((DES_POS_NEW(LL,2) >= YN(J_FLUID-1)) .AND.         &
                  (DES_POS_NEW(LL,2) <  YN(J_FLUID))) THEN
                  J_FLUID = J_FLUID
               ELSEIF((DES_POS_NEW(LL,2) >= YN(J_FLUID)) .AND.       &
                  (DES_POS_NEW(LL,2) < YN(J_FLUID+1))) THEN
                  J_FLUID = J_FLUID+1
               ELSEIF((DES_POS_NEW(LL,2) >= YN(J_FLUID-2)) .AND.     &
                  (DES_POS_NEW(LL,2) < YN(J_FLUID-1))) THEN
                  J_FLUID = J_FLUID-1
               ELSE
                  CALL PIC_SEARCH(J_FLUID, DES_POS_NEW(LL,2), YN,    &
                  DIMENSION_J, JMIN2, JMAX2)
               ENDIF
            ENDIF !(J)

            IF(NO_K) THEN
               K_FLUID = 1
            ELSE
               IF(K_FLUID <= KSTART3 .OR. K_FLUID >= KEND3) THEN
               CALL PIC_SEARCH(K_FLUID, DES_POS_NEW(LL,3), ZT,         &
                  DIMENSION_K, KMIN2, KMAX2)
               ELSE
                  IF((DES_POS_NEW(LL,3) >= ZT(K_FLUID-1)) .AND.        &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID))) THEN
                     K_FLUID = K_FLUID
                  ELSEIF((DES_POS_NEW(LL,3) >= ZT(K_FLUID)) .AND.      &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID+1))) THEN
                     K_FLUID = K_FLUID+1
                  ELSEIF((DES_POS_NEW(LL,3) >= ZT(K_FLUID-2)) .AND.    &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID-1))) THEN
                     K_FLUID = K_FLUID-1
                  ELSE
                     CALL PIC_SEARCH(K_FLUID, DES_POS_NEW(LL,3), ZT,   &
                     DIMENSION_K, KMIN2, KMAX2)
                  ENDIF
               ENDIF
            ENDIF !(K)
            IJK_FLUID = PIJK(LL,4)
         ENDIF ! (NOT CONTINUUM_COUPLED OR EXPLICIT)


! Initialize counter for number of facets
         COUNT_FACETS = 0

         ! Loop over potential facets
         DO CELL_COUNT = 1, facets_at_dg(CELL_ID)%count
            ! Get direction (axis) and facet id (NF)
            axis = facets_at_dg(cell_id)%dir(cell_count)
            NF = facets_at_dg(cell_id)%id(cell_count)

            ! Check to see if facet is out of reach of particle
            if (facets_at_dg(cell_id)%min(cell_count) >    &
               partpos_max(axis)) then
               cycle
            endif
            if (facets_at_dg(cell_id)%max(cell_count) <    &
               partpos_min(axis)) then
               cycle
            endif

            ! Compute projection (normal distance) from wall to particle
            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(LL,:),&
            &        NORM_FACE(:,NF))

            ! If normal exceeds particle radius, particle is not in contact
            if((line_t.lt.-RLENS))CYCLE

            ! Compute closest point on facet
            POS_TMP(:) = DES_POS_NEW(LL,:)
            CALL ClosestPtPointTriangle(POS_TMP, VERTEX(:,:,NF), &
            &    CLOSEST_PT(:))
            ! Compute position vector from particle to closest point on facet
            DIST(:) = CLOSEST_PT(:)-POS_TMP(:)
            DISTSQ = DOT_PRODUCT(DIST,DIST)

            ! Skip particles that are more than lens radius from facet
            IF(DISTSQ .GE. (RLENS_SQ-SMALL_NUMBER))CYCLE

            ! Only do heat transfer for particles that are directly above facet
            ! Heat transfer routines not generalized for edge contacts, but
            ! those should normally yield negligible heat transfer contributions

            ! Normalize distance vector and compare to facet normal
            NORMAL(:)=-DIST(:)/sqrt(DISTSQ)
            PROJ = sqrt(abs(DOT_PRODUCT(NORMAL, NORM_FACE(:,NF))))
            IF(ABS(ONE-PROJ).gt.1.0D-6)CYCLE

            ! Get overlap
            OVERLAP = RPART - SQRT(DISTSQ)

            ! Initialize area for facet (for post-proc. flux)
            AREA = ZERO

            ! Initialize BDRY FLAG
            DOMAIN_BDRY = .FALSE.
            ! Get BC_ID
            BC_ID = BC_ID_STL_FACE(NF)

            ! BC_ID is set to 0 in pre-proc if stl is a domain boundary
            IF(BC_ID.eq.0)then
               I1=I_FLUID
               J1=J_FLUID
               K1=K_FLUID
               DOMAIN_BDRY = .TRUE.

               ! Domain boundary, figure out real boundary ID
               IF(NORM_FACE(1,NF).ge.0.9999)THEN
                  ! WEST face
                  I1=IMIN2
                  IF(DO_K)THEN
                     AREA = DY(J_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DY(J_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(1,NF).le.-0.9999)THEN
                  ! EAST FACE
                  I1=IMAX2
                  IF(DO_K)THEN
                     AREA = DY(J_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DY(J_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(2,NF).ge.0.9999)THEN
                  ! SOUTH FACE
                  J1=JMIN2
                  IF(DO_K)THEN
                     AREA = DX(I_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DX(I_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(2,NF).le.-0.9999)THEN
                  ! NORTH FACE
                  J1=JMAX2
                  IF(DO_K)THEN
                     AREA = DX(I_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DX(I_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(3,NF).ge.0.9999)THEN
                  ! BOTTOM FACE
                  K1=KMIN2
                  AREA = DX(I_FLUID)*DY(J_FLUID)

               ELSEIF(NORM_FACE(3,NF).le.-0.9999)THEN
                  ! TOP FACE
                  K1=KMAX2
                  AREA = DX(I_FLUID)*DY(J_FLUID)
               ELSE
                  WRITE( *,*)'PROBLEM, COULD NOT FIND DOMAIN BOUNDARY'
                  WRITE(*,*)' In calc_thermo_des_wall_stl'
                  call mfix_exit(1)

               ENDIF


               ! Loop through defined BCs to see which one particle neighbors
               DO IBC = 1, DIMENSION_BC
                  IF(.NOT.BC_DEFINED(IBC))CYCLE
                  IF (I1.ge.BC_I_W(IBC).and.I1.le.BC_I_E(IBC).and.&
                      J1.ge.BC_J_S(IBC).and.J1.le.BC_J_N(IBC).and.&
                      K1.ge.BC_K_B(IBC).and.K1.le.BC_K_T(IBC))THEN
                      BC_ID = IBC
                      exit
                  ENDIF
               ENDDO
               IF(BC_ID.eq.0)then
                  IF(OUTPUT_WARNING)THEN
                     write(*,*)'WARNING: Could not find BC'
                     write(*,*)'Check input deck to make sure domain boundaries &
                     are defined'
                     write(*,*)'SUPPRESSING FUTURE OUTPUT'
                     write(*,*)'DES_POS_NEW'
                     write(*,'(3(F12.5, 3X))')(DES_POS_NEW(LL,IBC),IBC=1,3)
                     write(*,*)'I,J,K'
                     write(*,*)I1,J1,K1
                     write(*,*)'CLOSEST PT'
                     write(*,'(3(F12.5, 3X))')(CLOSEST_PT(IBC),IBC=1,3)
                     write(*,*)'NORM_FACE'
                     write(*,'(3(F12.5, 3X))')(NORM_FACE(IBC,NF),IBC=1,3)
                     OUTPUT_WARNING = .FALSE.
                  ENDIF
                  CYCLE  !
               ENDIF

            ENDIF !Domain Boundary (facet ID was 0)

            IF (BC_TYPE_ENUM(BC_ID) == NO_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == FREE_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == PAR_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_NSW .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_FSW .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_PSW) THEN



               ! CHECK TO MAKE SURE FACET IS UNIQUE
               USE_FACET=.TRUE.
               DO IFACET=1,count_facets
                  ! DO CHECK BY ENSURING NORMAL VECTOR IS NEARLY PARALLEL
                  PROJ = sqrt(abs(DOT_PRODUCT(NORMAL, NORM_FAC_CONTACT(:,IFACET))))
                  IF(ABS(ONE-PROJ).gt.1.0D-6)THEN
                     USE_FACET=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF(.NOT.USE_FACET)CYCLE

               ! FACET IS UNIQUE
               count_facets=count_facets+1
               NORM_FAC_CONTACT(:,count_facets)=NORMAL(:)

! Do heat transfer
               ! GET WALL TEMPERATURE
               TWALL = BC_TW_S(BC_ID,phase_LL)

               ! GET GAS THERMAL CONDUCTIVITY
               if(k_g0.eq.UNDEFINED)then
                  ! Compute gas conductivity as is done in calc_k_g
                  ! But use average of particle and wall temperature to be gas temperature

                  K_Gas = 6.02D-5*SQRT(HALF*(TWALL+TPART)/300.D0) ! cal/(s.cm.K)
                  ! 1 cal = 4.183925D0 J
                  IF (UNITS == 'SI') K_Gas = 418.3925D0*K_Gas !J/s.m.K
               else
                  K_Gas=k_g0
               endif
               IF(TWALL.eq.UNDEFINED)CYCLE
               QSWALL = DES_CONDUCTION_WALL(OVERLAP,K_s0(phase_LL), &
               &        K_s0(phase_LL),K_Gas,TWALL, TPART, RPART, &
               &        RLENS, phase_LL)


               Q_Source(LL) = Q_Source(LL)+QSWALL

               ! BELOW CODE IS ONLY NECESSARY FOR OUTPUTING
               ! DATA.  Need to know fluid cell that contact
               ! point resides in so that wall flux can be
               ! output correctly.

               I_FACET = I_FLUID
               J_FACET = J_FLUID
               K_FACET = K_FLUID

               ! This checks to see if the contact point was NOT on a domain boundary

               IF(.NOT.DOMAIN_BDRY)THEN
                  IF(CLOSEST_PT(1) >= XE(I_FLUID))THEN
                     I_FACET = MIN(IMAX1, I_FLUID+1)
                  ELSEIF(CLOSEST_PT(1) < XE(I_FLUID-1))THEN
                     I_FACET = MAX(IMIN1, I_FLUID-1)
                  ENDIF

                  IF(CLOSEST_PT(2) >= YN(J_FLUID))THEN
                     J_FACET = MIN(JMAX1, J_FLUID+1)
                  ELSEIF(CLOSEST_PT(2) < YN(J_FLUID-1))THEN
                     J_FACET = MAX(JMIN1, J_FLUID-1)
                  ENDIF
                  IF(DO_K)THEN
                     IF(CLOSEST_PT(3) >= ZT(K_FLUID))THEN
                        K_FACET = MIN(KMAX1, K_FLUID+1)
                     ELSEIF(CLOSEST_PT(3) < ZT(K_FLUID-1))THEN
                        K_FACET = MAX(KMIN1, K_FLUID-1)
                     ENDIF
                  ENDIF
                  IJK_FACET=funijk(I_facet,J_facet,K_facet)
                  AREA=AREA_CUT(IJK_FACET)

               ! AREA is left undefined if contact was with cut-cell surface
                  IF (AREA.eq.ZERO)then
                     I_FACET = I_FLUID
                     J_FACET = J_FLUID
                     K_FACET = K_FLUID
                     IJK_FACET=funijk(I_facet,J_facet,K_facet)
                     IF(NORM_FACE(1,NF).ge.0.9999)THEN
                     ! WEST face
                        IF(DO_K)THEN
                           AREA = DY(J_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DY(J_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(1,NF).le.-0.9999)THEN
                     ! EAST FACE
                        IF(DO_K)THEN
                           AREA = DY(J_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DY(J_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(2,NF).ge.0.9999)THEN
                     ! SOUTH FACE
                        IF(DO_K)THEN
                           AREA = DX(I_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DX(I_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(2,NF).le.-0.9999)THEN
                     ! NORTH FACE
                        IF(DO_K)THEN
                           AREA = DX(I_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DX(I_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(3,NF).ge.0.9999)THEN
                     ! BOTTOM FACE
                        AREA = DX(I_FACET)*DY(J_FACET)
                     ELSEIF(NORM_FACE(3,NF).le.-0.9999)THEN
                     ! TOP FACE
                        AREA = DX(I_FACET)*DY(J_FACET)
                     ENDIF
                  ENDIF ! Area==0 (because cut-cell facet exists in non cut-cell)
               ENDIF ! NOT DOMAIN_BDRY
               IJK_FACET = FUNIJK(I_FACET,J_FACET,K_FACET)
               ! AM error check
               IF(.NOT.FLUID_AT(IJK_FACET).AND. &
               &  .NOT.BLOCKED_CELL_AT(IJK_FACET))THEN
                  write(*,*)'ERROR: Cell containing facet is not a fluid &
                  &  cell or a blocked cell'
                  write(*,*)FLUID_AT(IJK_FACET), BLOCKED_CELL_AT(IJK_FACET)
                  write(*,*)'PART POS',(DES_POS_NEW(LL,IBC),IBC=1,3)
                  write(*,*)'FACET NORM',(NORM_FACE(IBC,NF),IBC=1,3)
                  write(*,*)'BC_ID', BC_ID
                  write(*,*)'I,J,K (Facet)', I_FACET,J_FACET,K_FACET

                  call mfix_exit(1)
               ENDIF

               IF(FLUID_AT(IJK_FACET))THEN
                  DES_QW_Cond(IJK_FACET,phase_LL) = &
                     DES_QW_Cond(IJK_FACET, phase_LL) + QSWALL/AREA
               ENDIF

            ENDIF ! WALL BDRY
         ENDDO                  ! CELL_COUNT (facets)
      ENDDO  ! LL
      RETURN

      END SUBROUTINE CALC_DEM_THERMO_WITH_WALL_STL

      end module CALC_COLLISION_WALL
