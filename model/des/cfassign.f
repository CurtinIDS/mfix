!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CFASSIGN                                                C
!                                                                      C
!  Purpose: Assign the necessary values for DEM computation. For       C
!           example:                                                   C
!     - assigning DEM boundaries from the values entered for           C
!       MFIX input in mfix.datat                                       C
!     - assigning DEM gravity vector from MFIX input.                  C
!     - calculating DTSOLID based on particle properties: spring       C
!       coefficient, damping factor & mass                             C
!                                                                      C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!  Reviewer: Rahul Garg                               Date: 25-Mar-14  C
!  Comments: Breaking this subroutine into several subroutines for DEM C
!            and PIC models                                            C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE constant, only : GRAVITY_X, GRAVITY_Y, GRAVITY_Z
      USE discretelement
      USE mfix_pic
      use error_manager
! Flag: DEM solids present.
      use run, only: DEM_SOLIDS
! Runtime flag specifying MPPIC solids
      use run, only: PIC_SOLIDS

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

      CALL INIT_ERR_MSG("CFASSIGN")

! The common assignments are done in this routine.
! The model spcific assignmets are moved to the specific subroutines

      GRAV(1) = GRAVITY_X
      GRAV(2) = GRAVITY_Y
      GRAV(3) = GRAVITY_Z
      GRAV_MAG = sqrt(dot_product(GRAV,GRAV))

      IF(DEM_SOLIDS) CALL CFASSIGN_DEM
      IF(PIC_SOLIDS) CALL CFASSIGN_PIC

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CFASSIGN



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CFASSIGN_PIC                                            C
!                                                                      C
!  Purpose: PIC related cfassign source code moved here                C
!           example:                                                   C
!     - calculating DTSOLID based on particle response time            C
!                                                                      C
!                                                                      C
!  Reviewer: Rahul Garg                               Date: 25-Mar-14  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN_PIC


!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement, only: des_mmax, dtsolid
      use param1, only: large_number
      USE physprop, only: MU_g0, mmax, d_p0, ro_s0
      USE mfix_pic, only : dtpic_taup, des_tau_p
      USE mfix_pic, only: MPPIC_PDRAG_IMPLICIT
      use error_manager
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: M, MMAX_TOT
!-----------------------------------------------

      CALL INIT_ERR_MSG("CFASSIGN_PIC")

      MMAX_TOT = DES_MMAX+MMAX
      DO M = MMAX+1, MMAX_TOT
! account for a possible offset index when using d_p0 and ro_s.
         DES_TAU_P(M) = RO_S0(M)*(D_P0(M)**2.d0)/(18.d0*MU_g0)
         WRITE(err_msg,'(/A,I2,A,G17.8)') &
         'TAU_P FOR ', M,'th SOLID PHASE= ', DES_TAU_P(M)
         CALL FLUSH_ERR_MSG (Header = .false., Footer = .false.)

      ENDDO

      DTSOLID = MINVAL(DES_TAU_P(MMAX+1:MMAX_TOT))
      DTPIC_TAUP = DTSOLID      !maximum dt for point-particles based on taup
      if(MPPIC_PDRAG_IMPLICIT) DTPIC_TAUP = LARGE_NUMBER

      WRITE(ERR_MSG, 1000) DTSolid
      CALL FLUSH_ERR_MSG(Header = .false.)

 1000 format('MPPIC: Point-particle ',&
      'approximation for particle-particle and ', /, &
      'particle-wall collisions', /, &
      'DTSOLID based on particle time response Taup', /, &
      'DTSOLID = ', 2x, E17.10)

! Finalize the error manager.
      CALL FINL_ERR_MSG

      END SUBROUTINE CFASSIGN_PIC




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CFASSIGN_DEM                                            C
!                                                                      C
!  Purpose: DEM related cfassign source code moved here                C
!           example:                                                   C
!     - calculating DTSOLID based on particle properties: spring       C
!       coefficient, damping factor & mass                             C
!                                                                      C
!                                                                      C
!  Reviewer: Rahul Garg                               Date: 25-Mar-14  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN_DEM


!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE discretelement
      use error_manager
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

      CALL INIT_ERR_MSG("CFASSIGN_DEM")

!      WRITE(ERR_MSG,'(A)') 'Setting collision model parameters for DEM'
!      CALL FLUSH_ERR_MSG (Footer = .false.)

! Finalize the error manager.
      CALL FINL_ERR_MSG
      END SUBROUTINE CFASSIGN_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! subroutine: compute_volume_of_nodes                                      C
! Author: Rahul Garg                                                       C
! Dare: Dec 20, 2011                                                       C
! Purpose:  Due to the staggered grid, the interpolaiton of mean fields    C
! is always first done at the nodes (of the scalar cell) for any quantity. C
! For a quantity like drag force or ep_s, one needs to assign a geometric  C
! volume to a node. In the past this was done on-the-fly in drag_fgs.f.    C
! VCELL was the variable that was used and it used to be incorrecty set to C
! the volume of the scalar cell that the particle belongs to.              C
! This will be ok for uniform grid and will not hold for non-uniform grids.C
! In addition, at the edges (in 2-D, the node volume is only half of       C
! the standard cell volume. And for corner (in 2-D), it is only one fourth.C
! This was also done on-the-fly earlier                                    C
! But again the volume of the cell was used, which was not correct         C
! also not extendable to cut-cell. So this routine computes the geoemetric C
! volume of the nodes.                                                     C
!                                                                          C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE compute_volume_of_nodes

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE cutcell
      USE functions
      implicit none
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! ijk index of fluid grid and corresponding i, j, k indices
      integer :: ijk, iraw, jraw, kraw
! i, j, k and (i+1, j+1, k+1) indices corrected for any
! cyclic ghost boundaries on fluid grid
      INTEGER :: I, J, K, ip, jp, kp
      integer :: ipjk, ijpk, ipjpk, ijkp, ipjkp, ijpkp, ipjpkp
! volume of indicated grid
      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp
      double precision :: vol_node_count, vol_node_actual_count
! weighting factor
      double precision :: avg_factor
! not used?
      double precision :: vol_node_uncorr
!-----------------------------------------------

      avg_factor = merge(0.25d0, 0.125d0, NO_K)

! compute the volume at the grid nodes
! grid nodes start from istart2 to iend1
      vol_node_count = merge(4., 8., NO_K)

      DO ijk = ijkstart3,ijkend3
         des_vol_node(ijk) = zero
         iraw = i_of(ijk)
         jraw = j_of(ijk)
         kraw = k_of(ijk)



! start at 1 (ghost cell) and go to last fluid cell. why start at a
! ghost cell and not a fluid cell?
! See below

! Since we are interested in nodes making up the interpolation stencil,
! their numbering goes from 1 to iend1.
! Think of a case with IMAX = 3. Here the nodes where the interpolation will be
! done will run from 1 (=istart2) to 4 (=iend1).
         IF(iraw.LT.istart2 .OR. iraw.GT.iend1) CYCLE
         IF(jraw.LT.jstart2 .OR. jraw.GT.jend1) CYCLE
         IF(kraw.LT.kstart2 .OR. kraw.GT.kend1) CYCLE

! this will force indices of ghost cells on cyclic borders back to
! the corresponding fluid cells. since we are using i, j, k indices and
! not just a composite ijk index we need these to be shifted to account
! for periodicity
         I = imap_c(iraw)
         J = jmap_c(jraw)
         K = kmap_c(kraw)
         IP = imap_c(iraw+1)
         JP = jmap_c(jraw+1)

! using a function like ip_of(ijk) should work the same as getting funijk
! of the shifted i, j, k indices.  however, small differences will
! occur on the 'edges/corners'. so retaining the latter method at this
! time. see j. galvin for discussion
         ipjk = funijk(IP,J,K)
         ijpk = funijk(I,JP,K)
         ipjpk = funijk(IP,JP,K)

! the existing variable vol(ijk) is not used here for cut-cell reasons
! see r. garg for discussion
         vol_ijk   = dx(I) *dy(J) *dz(K)
         vol_ipjk  = dx(IP)*dy(J) *dz(K)
         vol_ijpk  = dx(I) *dy(JP)*dz(K)
         vol_ipjpk = dx(IP)*dy(JP)*dz(K)

         vol_node_uncorr = avg_factor*(vol_ijk + vol_ipjk + vol_ijpk + vol_ipjpk)
         vol_node_actual_count = vol_node_count

         IF(.NOT.FLUID_AT(ijk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijk  = zero
         ENDIF

         IF(.NOT.FLUID_AT(ipjk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjk  = zero
         ENDIF

         IF(.NOT.FLUID_AT(ijpk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijpk  = zero
         ENDIF

         IF(.NOT.FLUID_AT(ipjpk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjpk = zero
         ENDIF

! this will have non-zero values for non-fluid cells at the
! west/south/bottom borders but not for east/north/top borders?
         des_vol_node(ijk) = avg_factor*(vol_ijk + vol_ipjk + &
            vol_ijpk + vol_ipjpk)

         IF (DO_K) THEN
            KP     = kmap_c(kraw+1)
            ijkp   = funijk(I, J, KP)
            ipjkp  = funijk(IP,J, KP)
            ijpkp  = funijk(I, JP,KP)
            ipjpkp = funijk(IP,JP,KP)

            vol_ijkp   = dx(I) *dy(J) *dz(KP)
            vol_ipjkp  = dx(IP)*dy(J) *dz(KP)
            vol_ijpkp  = dx(I) *dy(JP)*dz(KP)
            vol_ipjpkp = dx(IP)*dy(JP)*dz(KP)

            vol_node_uncorr = avg_factor*(vol_node_uncorr + vol_ijkp + &
               vol_ipjkp + vol_ijpkp + vol_ipjpkp)

            IF(.NOT.FLUID_AT(ijkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ijkp   = zero
            ENDIF

            IF(.NOT.FLUID_AT(ipjkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ipjkp  = zero
            ENDIF

            IF(.NOT.FLUID_AT(ijpkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ijpkp  = zero
            ENDIF

            IF(.NOT.FLUID_AT(ipjpkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ipjpkp = zero
            ENDIF

            des_vol_node(ijk) = des_vol_node(ijk) + avg_factor*&
               (vol_ijkp + vol_ipjpkp + vol_ijpkp + vol_ipjkp)

         ENDIF


      ENDDO   ! end do ijk=ijkstart3,ijkend3

      RETURN
      RETURN
      END SUBROUTINE compute_volume_of_nodes

