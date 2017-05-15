!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: MPPIC_COMP_EULERIAN_VELS_NON_CG                         !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MPPIC_COMP_EULERIAN_VELS_NON_CG
        USE compar
        USE constant
        use desmpi
        USE discretelement
        use fldvar, only: u_s, v_s, w_s
        USE functions
        USE geometry
        USE indices
        USE mfix_pic
        USE parallel
        USE param
        USE param1
        USE physprop, only: mmax
        IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP

! index of solid phase that particle NP belongs to
      INTEGER :: M
      INTEGER :: MMAX_TOT
!-----------------------------------------------

      MMAX_TOT = MMAX+DES_MMAX
      DO IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            !U_so(IJK, :) = U_s(IJK, :)
            !V_so(IJK, :) = V_s(IJK, :)
            !W_so(IJK, :) = W_s(IJK, :)
! could these be replaced by u_s, w_s, v_s?
            PIC_U_S(IJK, :) = ZERO
            PIC_V_S(IJK, :) = ZERO
            PIC_W_S(IJK, :) = ZERO
            IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE
            IF(.NOT.FLUID_AT(IJK)) CYCLE

            IF(I.GE.IMIN1.AND.I.LT.IMAX1) then !because we don't want to
               !calculate solids velocity at the wall cells.
               IPJK = IP_OF(IJK)
               DO M = MMAX+1, MMAX_TOT
                  PIC_U_S(IJK,M) = 0.5d0*(U_S(IJK, M) + U_S(IPJK,M))
               ENDDO
            ENDIF

            if(J.GE.JMIN1.AND.J.LT.JMAX1) then !because we don't want to
               !calculate solids velocity at the wall cells.
               IJPK = JP_OF(IJK)
               DO M = MMAX+1, MMAX_TOT
                  PIC_V_S(IJK,M) = 0.5d0*(V_S(IJK, M) + V_S(IJPK,M))
               ENDDO
            ENDIF


            if(K.GE.KMIN1.AND.K.LT.KMAX1.AND.DO_K) then !because we don't want to
               !calculate solids velocity at the wall cells.
               IJKP = KP_OF(IJK)
               DO M = MMAX+1, MMAX_TOT
                  PIC_W_S(IJK,M) = 0.5d0*(W_S(IJK, M) + W_S(IJKP,M))
               ENDDO
            ENDIF
         ENDDO
         !CALL SET_WALL_BC(IER)
         !the above routine will apply noslip or free slip BC as per the mfix  convention.
         !currently, this implies NSW or FSW wall BC's will be re-applied to gas-phase
         !field as well. This can be changed later on to be more specific to MPPIC case
         !CALL WRITE_MPPIC_VEL_S
         CALL MPPIC_BC_U_S
         CALL MPPIC_BC_V_S
         IF(DO_K) CALL MPPIC_BC_W_S

      END SUBROUTINE MPPIC_COMP_EULERIAN_VELS_NON_CG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: MPPIC_COMP+EULERIAN_VELS_CG                             !
!  Author: R. Garg                                                     !
!                                                                      !
!  Puryyse:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MPPIC_COMP_EULERIAN_VELS_CG
        USE compar
        USE constant
        USE cutcell
        USE desmpi
        USE discretelement
        use fldvar, only: u_s, v_s, w_s
        USE functions
        USE geometry
        USE indices
        USE mfix_pic
        USE parallel
        USE param
        USE param1
        use physprop, only: mmax
        IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP

! index of solid phase that particle NP belongs to
      INTEGER :: M
      INTEGER :: MMAX_TOT
!-----------------------------------------------

      MMAX_TOT = MMAX+DES_MMAX
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         !U_so(IJK, :) = U_s(IJK, :)
         !V_so(IJK, :) = V_s(IJK, :)
         !W_so(IJK, :) = W_s(IJK, :)
         PIC_U_S(IJK, :) = ZERO
         PIC_V_S(IJK, :) = ZERO
         PIC_W_S(IJK, :) = ZERO

         IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE

         IF(WALL_U_AT(IJK)) THEN
            PIC_U_S(IJK, :) = ZERO
            !currently only No slip BC is being set on this mean
            !solid's velocity field. Later this wall part can be
            !treated separately and U_S set only for scalar cells
            !where FLUID_AT(IJK) is true.
         ELSE
            if(.not.FLUID_AT(IJK)) cycle

            IPJK = IP_OF(IJK)
            IF(FLUID_AT(IPJK)) THEN
               DO M = MMAX+1, MMAX_TOT
                  PIC_U_S(IJK,M) = 0.5d0*(U_S(IJK, M) + U_S(IPJK,M))
               ENDDO
            ELSE
               PIC_U_S(IJK,:) = U_S(IJK, :)
            ENDIF
         ENDIF

         IF(WALL_V_AT(IJK)) THEN
            PIC_V_S(IJK, :) = ZERO
         ELSE
            if(.not.FLUID_AT(IJK)) cycle
            IJPK = JP_OF(IJK)
            IF(FLUID_AT(IJPK)) THEN
               DO M = MMAX+1, MMAX_TOT
                  PIC_V_S(IJK,M) = 0.5d0*(V_S(IJK, M) + V_S(IJPK,M))
               ENDDO
            ELSE
               PIC_V_S(IJK,:) = V_S(IJK, :)
            ENDIF
         ENDIF

         IF(DO_K) THEN
            IF(WALL_W_AT(IJK)) THEN
               PIC_W_S(IJK, :) = ZERO
            ELSE
               if(.not.FLUID_AT(IJK)) cycle
               IJKP = KP_OF(IJK)
               IF(FLUID_AT(IJKP)) THEN
                  DO M = MMAX+1, MMAX_TOT
                     PIC_W_S(IJK,M) = 0.5d0*(W_S(IJK, M) + W_S(IJKP,M))
                  ENDDO
               ELSE
                  PIC_W_S(IJK,:) = W_S(IJK, :)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      END SUBROUTINE MPPIC_COMP_EULERIAN_VELS_CG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: MPPIC_APPLY_PS_GRAD_PART                                !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MPPIC_APPLY_PS_GRAD_PART(L)

      USE param
      USE param1
      USE parallel
      USE constant
      USE discretelement
      USE mpi_utility
      USE mfix_pic
      USE cutcell
      USE fldvar, only: ep_g, u_s, v_s, w_s
      USE fun_avg
      USE functions
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER M, IDIM
      INTEGER IJK, IJK_C

      DOUBLE PRECISION COEFF_EN, COEFF_EN2

      DOUBLE PRECISION DELUP(DIMN), PS_FORCE(DIMN), VEL_ORIG(DIMN)

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION :: MEANUS(DIMN, DIMENSION_M)
      DOUBLE PRECISION :: RELVEL(DIMN)
      DOUBLE PRECISION :: MEANVEL(DIMN)
      DOUBLE PRECISION :: VEL_NEW(DIMN)
!      INTEGER :: TOT_CASE, case1_count, case2_count, case3_count, case4_count

      LOGICAL :: INSIDE_DOMAIN
!-----------------------------------------------

      M = PIJK(L,5)
      IJK = PIJK(L,4)
      COEFF_EN  = MPPIC_COEFF_EN1
      COEFF_EN2 = MPPIC_COEFF_EN2

      VEL_ORIG(1:DIMN) = DES_VEL_NEW(L,:)
      VEL_NEW (1:DIMN) = DES_VEL_NEW(L,:)
      !IF(L.eq.1) WRITE(*,*) 'MPPIC COEFFS = ', COEFF_EN, COEFF_EN2
      IF(L.EQ.FOCUS_PARTICLE) THEN

         WRITE(*,'(A20,2x,3(2x,i4))') 'CELL ID = ', PIJK(L,1:3)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'EPS = ', 1.d0 - EP_g(PIJK(L,4))
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', DES_VEL_NEW(L,:)

         WRITE(*,'(A20,2x,3(2x,g17.8))') 'FC = ', FC(L,:)
      ENDIF

      MEANVEL(1) = U_S(IJK,M)
      MEANVEL(2) = V_S(IJK,M)
      IF(DO_K) MEANVEL(3) = W_S(IJK,M)

      PS_FORCE(:) = PS_GRAD(:,L)
      !IF(ABS(PS_FORCE(2)).GT.ZERO)  WRITE(*,*) 'PS_FORCE = ', PS_FORCE
      DELUP(:) = -PS_FORCE(:)

      MEANUS(:,M) =  AVGSOLVEL_P (:,L)
      !MEANUS(:,M) = MEANVEL(:)
      RELVEL(:) = DES_VEL_NEW(L,:) - MEANUS(:,M)

      !IF(EPg_P(L).gt.1.2d0*ep_star) RETURN

      DO IDIM = 1, DIMN

         IF(ABS(PS_FORCE(IDIM)).eq.zero) cycle

         IF(VEL_ORIG(IDIM)*MEANUS(IDIM,M).GT.ZERO) THEN

            IF(VEL_ORIG(IDIM)*DELUP(IDIM).GT.ZERO) THEN

               IF(ABS(MEANUS(IDIM,M)) .GT. ABS(VEL_ORIG(IDIM))) THEN
                  IJK_C = IJK
                  IF(IDIM.eq.1) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                  ELSEIF(IDIM.eq.2) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                  ELSEIF(IDIM.eq.3) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                  ENDIF
                  INSIDE_DOMAIN = .false.
                  INSIDE_DOMAIN = fluid_at(IJK_C)!and.(.not.cut_cell_at(IJK_C))

                  if(INSIDE_DOMAIN) then
                     VEL_NEW(IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
                  endif
!                  case4_count = case4_count + 1
               ELSE
                  !do nothing
               ENDIF
            ELSE
               IF(ABS(VEL_ORIG(IDIM)).GT.ABS(MEANUS(IDIM,M))) then
                  !VEL_NEW(IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
                  IJK_C = IJK
                  IF(IDIM.eq.1) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                  ELSEIF(IDIM.eq.2) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                  ELSEIF(IDIM.eq.3) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                  ENDIF
                  !VEL_NEW(IDIM) = (MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM))

                  INSIDE_DOMAIN = .false.
                  INSIDE_DOMAIN = fluid_at(IJK_C)!.and.(.not.cut_cell_at(IJK_C))

                  if(INSIDE_DOMAIN) then
                     VEL_NEW(IDIM) = (MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM))
                  else
                     VEL_NEW(IDIM) = -COEFF_EN*VEL_ORIG(IDIM)
                  endif

                  IJK_C = IJK
                  IF(IDIM.eq.1) then
                     if(VEL_NEW(IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                     if(VEL_NEW(IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                  ELSEIF(IDIM.eq.2) then
                     if(VEL_NEW(IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                     if(VEL_NEW(IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                  ELSEIF(IDIM.eq.3) then
                     if(VEL_NEW(IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                     if(VEL_NEW(IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                  ENDIF
                  INSIDE_DOMAIN = .false.
                  INSIDE_DOMAIN = fluid_at(IJK_C) !.and.(.not.cut_cell_at(IJK_C))

                  !if(.not.INSIDE_DOMAIN) then
                  !   VEL_NEW(IDIM) = -COEFF_EN*VEL_ORIG(IDIM)
                  !ENDIF

!                  case1_count = case1_count + 1
               ELSE
                  !do nothing
                  VEL_NEW(IDIM) = COEFF_EN2 * VEL_ORIG(IDIM)
                  !turning on the above would make the model uncondtionally stable
!                  case1_count = case1_count + 1

               ENDIF
            ENDIF
         ELSE
            IF(MEANUS(IDIM,M)*DELUP(IDIM).GT.ZERO) THEN
               VEL_NEW(IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
!               case2_count = case2_count + 1
            ELSE
!               case3_count = case3_count + 1
               !DO NOTHING
            ENDIF
         ENDIF

         IF(MPPIC_GRAV_TREATMENT) THEN
            IF(DELUP(IDIM)*GRAV(IDIM).LT.ZERO.AND.VEL_ORIG(IDIM)*GRAV(IDIM).GT.ZERO) THEN
               VEL_NEW(IDIM) = -COEFF_EN*VEL_ORIG(IDIM)
            ENDIF
         ENDIF
         DES_VEL_NEW( L,IDIM) = VEL_NEW(IDIM)
      ENDDO

         !
      if(L.eq.FOCUS_PARTICLE) THEN
         !iF((IJK.eq.epg_min_loc(1).or.IJK_OLD.eq.epg_min_loc(1)).and.epg_min2.lt.0.38) then
         !if(j.ne.2) cycle

         WRITE(*,'(A20,2x,3(2x,i5))') 'PIJK I, J, K =', I_OF(IJK),J_OF(IJK),K_OF(IJK)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', VEL_ORIG(:)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL NEW = ', DES_VEL_NEW(L,:)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'MEANUS = ', MEANUS(:,1)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_POS_NEW = ', DES_POS_NEW(L,:)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'GRAD PS = ', PS_FORCE(:)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DELUP =  ', DELUP(:)

         !WRITE(*,'(A20,2x,3(2x,g17.8))') 'UPRIMETAU_INT = ', UPRIMETAU_INT(:)
         !WRITE(*,'(A20,2x,3(2x,g17.8))') 'UPRIMETAU = ', UPRIMETAU(:)
         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) WRITE(*,'(A20,2x,3(2x,g17.8))') 'U*DT, MFP =', UPRIMEMOD*DTSOLID, MEAN_FREE_PATH
         read(*,*)
      ENDIF

      !TOT_CASE = case1_count + case2_count + case3_count + case4_count
      !IF(TOT_CASE.GT.0) THEN
      !WRITE(*,'(A, 4(2x,i10))') 'CASE COUNT NUMBERS  = ', case1_count ,case2_count ,case3_count ,case4_count
      !WRITE(*,'(A, 4(2x,g12.7))') 'CASE COUNT %AGE = ', real(case1_count)*100./real(tot_case), &
      !     real(case2_count)*100./real(tot_case), real(case3_count)*100./real(tot_case), real(case4_count)*100./real(tot_case)
      !ENDIF
      RETURN

      !MEANUS(:,M) = MEANVEL(:)
      !RELVEL(:) = DES_VEL_NEW(L,:) - MEANUS(:,M)
      !DO IDIM = 1, DIMN
      !    IF(ABS(PS_FORCE(IDIM)).eq.zero) cycle

      !   IF(RELVEL(IDIM)*DELUP(IDIM).GT.ZERO) THEN
      !do nothing
      !    ELSE
      !       DES_VEL_NEW(L,IDIM) = MEANUS(IDIM,M) - 0.4d0*RELVEL(IDIM)

      !IF(DES_VEL_NEW(L,IDIM)*DELUP(IDIM).LT.ZERO) DES_VEL_NEW(L,IDIM) = -0.5d0*DES_VEL_NEW(L,IDIM)
      !    ENDIF
      ! ENDDO
      ! CYCLE


      END SUBROUTINE MPPIC_APPLY_PS_GRAD_PART


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MPPIC_BC_U_S                                            !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MPPIC_BC_U_S

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE compar
      USE mfix_pic
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I1, J1, K1, IJK,&
                 IJK_WALL
!-----------------------------------------------

! Set the default boundary conditions
      IF (DO_K) THEN
         K1 = 1
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1+1)
               PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
            END DO
         END DO

         K1 = KMAX2
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1-1)
               PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
            END DO
         END DO
      ENDIF

      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1,J1+1,K1)
            PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
         END DO
      END DO
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1,J1-1,K1)
            PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
         END DO
      END DO

      RETURN
      END SUBROUTINE MPPIC_BC_U_S



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MPPIC_BC_V_S                                            !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MPPIC_BC_V_S


!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE compar
      USE mfix_pic
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER          I1, J1, K1, IJK,&
                       IJK_WALL
!-----------------------------------------------

! Set the default boundary conditions
      IF (DO_K) THEN
         K1 = 1
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1+1)
               PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)
            END DO
         END DO
         K1 = KMAX2
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1-1)
               PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)
            END DO
         END DO
      ENDIF

      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1+1,J1,K1)
            PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)

         END DO
      END DO
      I1 = IMAX2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1-1,J1,K1)
            PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)

         END DO
      END DO
      END SUBROUTINE MPPIC_BC_V_S



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MPPIC_BC_W_S                                            !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MPPIC_BC_W_S

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE compar
      USE mfix_pic
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I1, J1, K1, IJK,&
                 IJK_WALL
!-----------------------------------------------

! Set the default boundary conditions
      J1 = 1
      DO K1 = kmin3,kmax3
         DO I1 = imin3,imax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1+1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
           IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1-1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO
      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1+1,J1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO
      I1 = IMAX2
      DO K1 = kmin3,kmax3
         DO J1 = jmin3,jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1-1,J1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO


      END SUBROUTINE MPPIC_BC_W_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: WRITE_NODEDATA                                          !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_NODEDATA(funit)
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE geometry
      USE indices
      USE compar
      USE discretelement
      use desmpi
      USE functions
      USE fun_avg

      IMPLICIT NONE
      integer, intent(in) :: funit

      integer :: ijk, i, j,k

      write(funit,*)'VARIABLES= ',' "I" ',' "J" ',' "K" ',&
         ' "DES_ROPS_NODE" '

      write(funit, *)'ZONE F=POINT, I=', (IEND1-ISTART2)+1,  ', J=', &
         JEND1-JSTART2+1, ', K=', KEND1-KSTART2 + 1

      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IJK = funijk(I, J, K)
               !WRITE(*,*) 'IJK = ', IJK, I, J, K , SIZE(BUFIN,1)

               WRITE(funit, '(3(2x, i10), 3x, g17.8)') I, J, K,&
                  DES_ROPS_NODE(IJK,1)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(funit, status = 'keep')
      end SUBROUTINE WRITE_NODEDATA

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: WRITE_MPPIC_VEL_S                                       !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_MPPIC_VEL_S
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE geometry
      USE indices
      USE compar
      USE fldvar, only : ep_g, u_s, v_s, w_s
      use physprop, only: mmax
      USE discretelement
      USE mfix_pic
      USE functions
      implicit none
      integer :: i, j, k, ijk, fluid_ind, LL, PC, IDIM
      double precision :: zcor
      character(LEN=255) :: filename

      WRITE(filename,'(A,"_",I5.5,".dat")') TRIM(RUN_NAME)//'_U_S_',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted', &
           status='unknown',CONVERT='BIG_ENDIAN')
      IF(DIMN.eq.2) then
         write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ', &
              ' "EP_s " ', ' "pU_S" ', ' "pV_S" ',' "dU_s" ',&
              ' "dV_s" '!, ' "P_S_FOR1" ', ' "P_S_FOR2" '
         write(1000,*)'ZONE F=POINT, I=', (IEND3-ISTART3)+1,  ', J=',&
            JEND3-JSTART3+1, ', K=', KEND3-KSTART3 + 1
      else
         write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ', &
              ' "EP_s " ', ' "pU_S" ', ' "pV_S" ', ' "pW_S" ',&
              ' "dU_s" ', ' "dV_s" ', ' "dW_s" '!, &
        ! & ' "P_S_FOR1" ', ' "P_S_FOR2" ', ' "P_S_FOR3" '
         write(1000,*)'ZONE F=POINT, I=', (IEND3-ISTART3)+1,  ', J=',&
            JEND3-JSTART3+1, ', K=', KEND3-KSTART3 + 1
      ENDIF
      DO K=KSTART3, KEND3
         DO J=JSTART3, JEND3
            DO I=ISTART3, IEND3
               IJK  = FUNIJK(I,J,K)
               IF(FLUID_AT(IJK)) THEN
                  FLUID_IND = 1
               ELSE
                  FLUID_IND = 0
               END IF
               IF(DIMN.EQ.2) THEN
                  ZCOR = ZT(K)
                  WRITE(1000,'(3(2X,G17.8),4( 2X, G17.8))') &
                     XE(I-1)+DX(I), YN(J-1)+DY(J), ZCOR, 1.D0-EP_G(IJK),&
                     PIC_U_S(IJK,MMAX+1), PIC_V_S(IJK,MMAX+1),&
                     U_S(IJK,MMAX+1), V_S(IJK,MMAX+1) !, PS_FORCE_PIC(1,IJK), PS_FORCE_PIC(2,IJK)
               ELSE
                  ZCOR = ZT(K-1) + DZ(K)
                  WRITE(1000,'(3(2X,G17.8),4( 2X, G17.8))') XE(I-1)+DX(I),&
                     YN(J-1)+DY(J), ZCOR, 1.D0-EP_G(IJK), &
                     PIC_U_S(IJK,MMAX+1), PIC_V_S(IJK,MMAX+1), PIC_W_S(IJK,MMAX+1),&
                     U_S(IJK,MMAX+1), V_S(IJK,MMAX+1), W_S(IJK,MMAX+1)!, PS_FORCE_PIC(1,IJK), PS_FORCE_PIC(2,IJK),  PS_FORCE_PIC(3,IJK)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      CLOSE(1000, STATUS='KEEP')

      return

      WRITE(FILENAME,'(A,"_",I5.5,".DAT")') &
         TRIM(RUN_NAME)//'_PS_FORCE_',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted',&
           status='unknown',CONVERT='BIG_ENDIAN')

      IF(DIMN.eq.3) write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',&
           ' "DELPX" ', '"DELPY"', '"DELPZ" ',&
           ' "US_part" ', '"VS_part"' , '"WS_part"', '"EP_s_part"'
      IF(DIMN.eq.2) write(1000,*)'VARIABLES= ',' "X" ',' "Y" ', &
           ' "DELPX" ', '"DELPY"', ' "US_part" ', '"VS_part"' , &
           '"EP_S_part"'

      PC = 1
      DO LL = 1, MAX_PIP

         IF(PC .GT. PIP) EXIT
         IF(IS_NONEXISTENT(LL)) CYCLE
         pc = pc+1
         IF(IS_GHOST(LL) .OR. IS_ENTERING_GHOST(LL) .OR. &
            IS_EXITING_GHOST(LL)) CYCLE

         WRITE(1000,'(10( 2x, g17.8))') &
            (DES_POS_NEW(LL,IDIM), IDIM=1,DIMN), &
            (PS_GRAD(IDIM,LL), IDIM=1,DIMN), &
            (AVGSOLVEL_P(IDIM,LL),IDIM=1,DIMN), 1-EPg_P(LL)
      ENDDO
      close(1000, status='keep')

      !write(*,*) 'want to quit ?', LL, mAX_PIP, PIP
      !read(*,*) finish
      !if(finish) STOp
      END SUBROUTINE WRITE_MPPIC_VEL_S
