!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Calculate contact force and torque on particle from        !
!           particle-particle and particle-wall collisions. Treats     !
!           wall interaction also as a two-particle interaction but    !
!           accounting for the wall properties                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
#include "version.inc"
      SUBROUTINE CALC_FORCE_DEM

! Modules
!---------------------------------------------------------------------//
      USE calc_collision_wall
      USE constant, ONLY: Pi
      USE derived_types, only: multisap, boxhandle
      USE des_thermo
      USE des_thermo_cond
      USE discretelement
      USE run
      use pair_manager
      use param1, only: one, small_number, zero

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0
! particle no. indices
      INTEGER :: I, LL, cc, CC_START, CC_END
! the overlap occuring between particle-particle or particle-wall
! collision in the normal direction
      DOUBLE PRECISION :: OVERLAP_N, OVERLAP_T(3)
! square root of the overlap
      DOUBLE PRECISION :: SQRT_OVERLAP
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      DOUBLE PRECISION :: R_LM,DIST_CI,DIST_CL
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM, rad
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(3), NORMAL(3), DIST_MAG, POS(3)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: VREL_T(3)
! normal and tangential forces
      DOUBLE PRECISION :: FN(3), FT(3)
! temporary storage of force
      DOUBLE PRECISION :: FC_TMP(3)
! temporary storage of force for torque
      DOUBLE PRECISION :: TOW_FORCE(3)
! temporary storage of torque
      DOUBLE PRECISION :: TOW_TMP(3,2)
! temporary storage of conduction/radiation
      DOUBLE PRECISION :: QQ_TMP

! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER :: PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION :: ETAN_DES, ETAT_DES
      DOUBLE PRECISION :: KN_DES, KT_DES
! local values used for calculating cohesive forces
      DOUBLE PRECISION :: FORCE_COH, EQ_RADIUS, DistApart

      LOGICAL, PARAMETER :: report_excess_overlap = .FALSE.

      DOUBLE PRECISION :: FNMD, FTMD, MAG_OVERLAP_T, TANGENT(3)

      integer :: nn, mm, box_id, box_id2
      logical :: found
      integer :: pair(2)

!......................................................................!

! Initialize cohesive forces
      IF(USE_COHESION) PostCohesive(:) = ZERO

      CALL CALC_DEM_FORCE_WITH_WALL_STL

#ifdef do_sap
      ! do nn=0, size(multisap%saps)-1
      !    !print *,"nn = ",nn
      !    if (.not.check_boxes(multisap%saps(nn))) ERROR_STOP __LINE__
      !    if (.not.check_sort(multisap%saps(nn))) ERROR_STOP __LINE__
      ! enddo

!print *,"CALC_FORCE_DEM =================================================================================="

print *," TOTAL NUM OF NEIGHBORS IS ",NEIGHBOR_INDEX(MAX_PIP-1)

open (unit=123,file="neighbors.txt",action="write",status="replace")

DO LL = 1, MAX_PIP
   CC_START = 1
   IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
   CC_END   = NEIGHBOR_INDEX(LL)

   DO CC = CC_START, CC_END-1
      I  = NEIGHBORS(CC)
      write (123,*) ll,i
   enddo
enddo

close (unit=123)

             call reset_pairs(multisap%hashtable)
             do
                call get_pair(multisap%hashtable,pair)
                if (pair(1).eq.0 .and. pair(2).eq.0) exit

                if ( des_radius(pair(1))+des_radius(pair(2))> sqrt(dot_product(DES_POS_NEW(pair(1),:)-DES_POS_NEW(pair(2),:),DES_POS_NEW(pair(1),:)-DES_POS_NEW(pair(2),:)))) then
                   print *,"invalid pair: ",pair(1),pair(2)
                   stop __LINE__
                endif
             enddo
#endif

! Check particle LL neighbor contacts
!---------------------------------------------------------------------//

!!$omp parallel default(none) private(pos,rad,cc,cc_start,cc_end,ll,i,  &
!!$omp    overlap_n,vrel_t,v_rel_trans_norm,sqrt_overlap,dist,r_lm,     &
!!$omp    kn_des,kt_des,hert_kn,hert_kt,phasell,phasei,etan_des,        &
!!$omp    etat_des,fn,ft,overlap_t,tangent,mag_overlap_t,               &
!!$omp    eq_radius,distapart,force_coh,dist_mag,NORMAL,ftmd,fnmd,      &
!!$omp    dist_cl, dist_ci, fc_tmp, tow_tmp, tow_force, qq_tmp, box_id, box_id2, found)         &
!!$omp shared(max_pip,neighbors,neighbor_index,des_pos_new,des_radius,  &
!!$omp    des_coll_model_enum,kn,kt,pft_neighbor,pijk,                  &
!!$omp    des_etan,des_etat,mew,use_cohesion, calc_cond_des, dtsolid,   &
!!$omp    van_der_waals,vdw_outer_cutoff,vdw_inner_cutoff,              &
!!$omp    hamaker_constant,asperities,surface_energy,                   &
!!$omp    tow, fc, energy_eq, grav_mag, postcohesive, pmass, q_source, multisap, boxhandle)

!!$omp do

      DO LL = 1, MAX_PIP
         IF(IS_NONEXISTENT(LL)) CYCLE
         pos = DES_POS_NEW(LL,:)
         rad = DES_RADIUS(LL)

         CC_START = 1
         IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
         CC_END   = NEIGHBOR_INDEX(LL)

         DO CC = CC_START, CC_END-1
            I  = NEIGHBORS(CC)
            IF(IS_NONEXISTENT(I)) CYCLE

            R_LM = rad + DES_RADIUS(I)
            DIST(:) = DES_POS_NEW(I,:) - POS(:)
            DIST_MAG = dot_product(DIST,DIST)

            FC_TMP(:) = ZERO

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
               IF(DIST_MAG < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
                  EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) /     &
                    (DES_RADIUS(LL)+DES_RADIUS(I))
                  IF(DIST_MAG > (VDW_INNER_CUTOFF+R_LM)**2) THEN
                     DistApart = (SQRT(DIST_MAG)-R_LM)
                     FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS /           &
                        (12d0*DistApart**2) * (Asperities/(Asperities+    &
                        EQ_RADIUS) + ONE/(ONE+Asperities/DistApart)**2 )
                  ELSE
                     FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS *  &
                       (Asperities/(Asperities+EQ_RADIUS) + ONE/          &
                       (ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                  ENDIF
                  FC_TMP(:) = DIST(:)*FORCE_COH/SQRT(DIST_MAG)
                  TOW_TMP(:,:) = ZERO

! just for post-processing mag. of cohesive forces on each particle
                  PostCohesive(LL) = PostCohesive(LL) + FORCE_COH / PMASS(LL)
               ENDIF
            ENDIF

            IF (ENERGY_EQ) THEN
            ! Calculate conduction and radiation for thermodynamic neighbors
               IF(CALC_COND_DES(PIJK(LL,5))) THEN
                  QQ_TMP = DES_CONDUCTION(LL, I, sqrt(DIST_MAG), PIJK(LL,5), PIJK(LL,4))

!!$omp atomic
                  Q_Source(LL) = Q_Source(LL) + QQ_TMP

!!$omp atomic
                  Q_Source(I) = Q_Source(I) - QQ_TMP
               ENDIF
            ENDIF

            IF(DIST_MAG > (R_LM - SMALL_NUMBER)**2) THEN
               PFT_NEIGHBOR(:,CC) = 0.0
               CYCLE
            ENDIF

#ifdef do_sap
               if (.not.is_pair(multisap%hashtable,ll,i)) then

                  print *,"SAP DIDNT FIND PAIR: ",ll,i
                  print *,"PARTICLE (",ll,"):  ",des_pos_new(:,ll), " WITH RADIUS: ",des_radius(ll)
                  print *,"PARTICLE (",i,"):  ",des_pos_new(:,i), " WITH RADIUS: ",des_radius(i)

                  print *,""
                  print *," ******   ",sqrt(dot_product(des_pos_new(:,ll)-des_pos_new(:,i),des_pos_new(:,ll)-des_pos_new(:,i))),"     *********"
                  print *,""

                  ! print *,"LLLLLLLLL ",boxhandle(ll)%list(:)
                  ! print *,"IIIIIIIII ",boxhandle(i)%list(:)

                  do mm=1,size(boxhandle(ll)%list)
                     if (boxhandle(ll)%list(mm)%sap_id < 0 ) cycle
                     print *," PARTICLE ",ll," IS IN ",boxhandle(ll)%list(mm)%sap_id
                     box_id = boxhandle(ll)%list(mm)%box_id

                     found = .false.
                     do nn=1,size(boxhandle(i)%list)
                        if (boxhandle(i)%list(nn)%sap_id .eq. boxhandle(ll)%list(mm)%sap_id) then
                           ! print *," PARTICLE ",i," IS ALSO IN ",boxhandle(i)%list(nn)
                           box_id2 = boxhandle(i)%list(nn)%box_id
                           found = .true.
                        endif
                     enddo

                     if (.not.found) cycle

                     ! print *,"BOTH ",ll,i," ARE IN ",boxhandle(ll)%list(mm)

                  enddo

                  ERROR_STOP __LINE__
               endif
#endif

            IF(DIST_MAG == 0) THEN
               WRITE(*,8550) LL, I
               ERROR_STOP "division by zero"
 8550 FORMAT('distance between particles is zero:',2(2x,I10))
            ENDIF

            DIST_MAG = SQRT(DIST_MAG)
            NORMAL(:)= DIST(:)/DIST_MAG

! Calcuate the normal overlap
            OVERLAP_N = R_LM-DIST_MAG
            IF(REPORT_EXCESS_OVERLAP) CALL PRINT_EXCESS_OVERLAP

! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
            CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, VREL_T,            &
               NORMAL(:), DIST_MAG)

            phaseLL = PIJK(LL,5)
            phaseI = PIJK(I,5)

! Hertz spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
               KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
               ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap

! Linear spring-dashpot contact model
            ELSE
               KN_DES = KN
               KT_DES = KT
               ETAN_DES = DES_ETAN(phaseLL,phaseI)
               ETAT_DES = DES_ETAT(phaseLL,phaseI)
            ENDIF

! Calculate the normal contact force
            FN(:) =  -(KN_DES * OVERLAP_N * NORMAL(:) + &
               ETAN_DES * V_REL_TRANS_NORM * NORMAL(:))

! Calcuate the tangential overlap
            OVERLAP_T(:) = DTSOLID*VREL_T(:) + PFT_NEIGHBOR(:,CC)
            MAG_OVERLAP_T = sqrt(dot_product(OVERLAP_T,OVERLAP_T))

! Calculate the tangential contact force.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Tangential froce from spring.
               FTMD = KT_DES*MAG_OVERLAP_T
! Max force before the on set of frictional slip.
               FNMD = MEW*sqrt(dot_product(FN,FN))
! Direction of tangential force.
               TANGENT = OVERLAP_T/MAG_OVERLAP_T
! Frictional slip
               IF(FTMD > FNMD) THEN
                  FT = -FNMD * TANGENT
                  OVERLAP_T = (FNMD/KT_DES) * TANGENT
               ELSE
                  FT = -FTMD * TANGENT
               ENDIF
            ELSE
               FT = 0.0
            ENDIF

! Add in the tangential dashpot damping force
            FT = FT - ETAT_DES * VREL_T(:)

! Save tangential displacement history
            PFT_NEIGHBOR(:,CC) = OVERLAP_T(:)

! calculate the distance from the particles' centers to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
            DIST_CL = DIST_MAG/2.d0 + (DES_RADIUS(LL)**2 - &
               DES_RADIUS(I)**2)/(2.d0*DIST_MAG)

            DIST_CI = DIST_MAG - DIST_CL

            TOW_force(:) = DES_CROSSPRDCT(NORMAL(:), FT(:))
            TOW_TMP(:,1) = DIST_CL*TOW_force(:)
            TOW_TMP(:,2) = DIST_CI*TOW_force(:)

! Calculate the total force FC of a collision pair
! total contact force ( FC_TMP may already include cohesive force)
            FC_TMP(:) = FC_TMP(:) + FN(:) + FT(:)

            FC(LL,:) = FC(LL,:) + FC_TMP(:)

!!$omp atomic
            FC(I,1) = FC(I,1) - FC_TMP(1)
!!$omp atomic
            FC(I,2) = FC(I,2) - FC_TMP(2)
!!$omp atomic
            FC(I,3) = FC(I,3) - FC_TMP(3)

! for each particle the signs of norm and ft both flip, so add the same torque
            TOW(LL,:) = TOW(LL,:) + TOW_TMP(:,1)

!!$omp atomic
            TOW(I,1)  = TOW(I,1)  + TOW_TMP(1,2)
!!$omp atomic
            TOW(I,2)  = TOW(I,2)  + TOW_TMP(2,2)
!!$omp atomic
            TOW(I,3)  = TOW(I,3)  + TOW_TMP(3,2)

         ENDDO
      ENDDO
!!$omp end do

!!$omp end parallel

! just for post-processing mag. of cohesive forces on each particle
      IF(USE_COHESION .AND. VAN_DER_WAALS .AND. GRAV_MAG > ZERO) THEN
         PostCohesive(:) = PostCohesive(:)/GRAV_MAG
      ENDIF

      RETURN

      contains

        include 'functions.inc'

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: print_excess_overlap                                    !
!                                                                      !
!  Purpose: Print overlap warning messages.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PRINT_EXCESS_OVERLAP

      use error_manager

      IF(OVERLAP_N > flag_overlap*DES_RADIUS(LL) .OR.                  &
         OVERLAP_N > flag_overlap*DES_RADIUS(I)) THEN

         WRITE(ERR_MSG,1000) trim(iVAL(LL)), trim(iVAL(I)), S_TIME,    &
            DES_RADIUS(LL), DES_RADIUS(I), OVERLAP_N

         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1000 FORMAT('WARNING: Excessive overplay detected between ',          &
         'particles ',A,' and ',/A,' at time ',g11.4,'.',/             &
         'RADII:  ',g11.4,' and ',g11.4,4x,'OVERLAP: ',g11.4)

      END SUBROUTINE PRINT_EXCESS_OVERLAP

    END SUBROUTINE CALC_FORCE_DEM
