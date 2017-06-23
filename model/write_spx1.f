!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SPX1                                             C
!  Purpose: write out the time-dependent restart records (REAL)        C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: TIME, NSTEP, EP_g, P_g, P_star, U_g           C
!                        V_g, W_g, U_s, V_s, W_s, ROP_s, T_g, T_s      C
!                        IJKMAX2, MMAX                                 C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables:  LC, N, NEXT_REC, NUM_REC                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_SPX1(L, unit_add)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE cdist
      USE compar
      USE cutcell
      USE exit, only: mfix_exit
      USE fldvar
      USE funits
      USE geometry
      USE machine
      USE mpi_utility
      USE output
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE scalars
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DISCRETE_ELEMENT
      use discretelement, only: PRINT_DES_DATA

!//       USE tmp_array
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!             flag whether to write a particular SPx file
      INTEGER L

!              offset for use in post_mfix
      INTEGER  unit_add
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
! local variables
!
!//
      double precision, allocatable :: array1(:)     !//
      double precision, allocatable :: array2(:)     !//

!             loop counters
      INTEGER LC, NN
!
!             Pointer to the next record
      INTEGER NEXT_REC
!
!              Number of records written each time step
      INTEGER  NUM_REC

      INTEGER  uspx   ! UNIT_SPX + offset from post_mfix
      CHARACTER(LEN=50), DIMENSION(1) :: LINE   !error message
      double precision, dimension(:), allocatable :: TMP_VAR

      allocate(TMP_VAR(DIMENSION_3))

!-----------------------------------------------
      uspx = UNIT_SPX + unit_add

!
      if (myPE .eq.PE_IO) then
         allocate (array1(ijkmax2))   !//
         allocate (array2(ijkmax3))   !//
      else
         allocate (array1(1))   !//
         allocate (array2(1))   !//
      end if
!
! ".SP1" FILE         EP_g    [ ROP_g, RO_g  must be calculated ...
!                                        not written out ]
!
      SELECT CASE (L)
      CASE (1)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if

! Explicitly coupled simulations do not need to rebin particles to
! the fluid grid every time step. However, this implies that the
! fluid cell information and interpolation weights become stale.
         IF(DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_COUPLED) THEN
! Bin particles to fluid grid.
            CALL PARTICLES_IN_CELL
! Calculate interpolation weights
            CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
            CALL COMP_MEAN_FIELDS
         ENDIF

         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               CALL UNSHIFT_DP_ARRAY(EP_g,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
            ELSE
               call OUT_BIN_R(uspx+L,EP_g,size(EP_g),NEXT_REC)
            ENDIF
!           call OUT_BIN_R(uspx+L,EP_g,size(EP_g),NEXT_REC)
         else
            call gatherWriteSpx (EP_g,array2, array1, uspx+L, NEXT_REC)   !//
         end if
         if (myPE .eq. PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if

! The call made in make_arrays captures the initial state of the system
! as the input and RES files for DES runs are read afte the the first
! call to this routine.
!         IF(DISCRETE_ELEMENT.AND.PRINT_DES_DATA) THEN
!            IF(TIME /= ZERO .OR. TRIM(RUN_TYPE)=='RESTART_1') &
!               CALL WRITE_DES_DATA
!         ENDIF

!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP2" FILE         P_g , P_star
!
      CASE (2)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               CALL UNSHIFT_DP_ARRAY(P_g,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               CALL UNSHIFT_DP_ARRAY(P_star,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
            ELSE
              call OUT_BIN_R(uspx+L,P_g,size(P_g),NEXT_REC)
              call OUT_BIN_R(uspx+L,P_star,size(P_star),NEXT_REC)
            ENDIF
!           call OUT_BIN_R(uspx+L,P_g,size(P_g),NEXT_REC)
!           call OUT_BIN_R(uspx+L,P_star,size(P_star),NEXT_REC)
         else
           call gatherWriteSpx (P_g,array2, array1, uspx+L, NEXT_REC)   !//
           call gatherWriteSpx (P_star,array2, array1, uspx+L, NEXT_REC)   !//
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP3" FILE         U_g , V_g , W_g
!
      CASE (3)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               CALL UNSHIFT_DP_ARRAY(U_g,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               CALL UNSHIFT_DP_ARRAY(V_g,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               CALL UNSHIFT_DP_ARRAY(W_g,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
            ELSE
               call OUT_BIN_R(uspx+L,U_g,size(U_g),NEXT_REC)
               call OUT_BIN_R(uspx+L,V_g,size(V_g),NEXT_REC)
               call OUT_BIN_R(uspx+L,W_g,size(W_g),NEXT_REC)
            ENDIF
!           call OUT_BIN_R(uspx+L,U_g,size(U_g),NEXT_REC)
!           call OUT_BIN_R(uspx+L,V_g,size(V_g),NEXT_REC)
!           call OUT_BIN_R(uspx+L,W_g,size(W_g),NEXT_REC)
         else
           call gatherWriteSpx (U_g,array2, array1, uspx+L, NEXT_REC)   !//
           call gatherWriteSpx (V_g,array2, array1, uspx+L, NEXT_REC)   !//
           call gatherWriteSpx (W_g,array2, array1, uspx+L, NEXT_REC)   !//
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP4" FILE         U_s , V_s , W_s
!
      CASE (4)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            DO LC = 1, MMAX
               IF(RE_INDEXING) THEN
                  CALL UNSHIFT_DP_ARRAY(U_s(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
                  CALL UNSHIFT_DP_ARRAY(V_s(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
                  CALL UNSHIFT_DP_ARRAY(W_s(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               ELSE
                  call OUT_BIN_R(uspx+L,U_s(:,LC),Size(U_s(:,LC)),NEXT_REC)
                  call OUT_BIN_R(uspx+L,V_s(:,LC),Size(V_s(:,LC)),NEXT_REC)
                  call OUT_BIN_R(uspx+L,W_s(:,LC),Size(W_s(:,LC)),NEXT_REC)
               ENDIF
            ENDDO
!        DO LC = 1, MMAX
!          call OUT_BIN_R(uspx+L,U_s(:,LC),Size(U_s(:,LC)),NEXT_REC)
!          call OUT_BIN_R(uspx+L,V_s(:,LC),Size(V_s(:,LC)),NEXT_REC)
!          call OUT_BIN_R(uspx+L,W_s(:,LC),Size(W_s(:,LC)),NEXT_REC)
!        END DO
         else
         DO LC = 1, MMAX
            call gatherWriteSpx (U_s(:,LC),array2, array1, uspx+L, NEXT_REC)
            call gatherWriteSpx (V_s(:,LC),array2, array1, uspx+L, NEXT_REC)
            call gatherWriteSpx (W_s(:,LC),array2, array1, uspx+L, NEXT_REC)
         END DO
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP5" FILE         ROP_s
!
      CASE (5)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            DO LC = 1, MMAX
               IF(RE_INDEXING) THEN
                  CALL UNSHIFT_DP_ARRAY(ROP_s(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               ELSE
                  call OUT_BIN_R(uspx+L,ROP_s(:,LC),size(ROP_s(:,LC)), NEXT_REC)
               ENDIF
            ENDDO
!          DO LC = 1, MMAX
!            call OUT_BIN_R(uspx+L,ROP_s(:,LC),size(ROP_s(:,LC)), NEXT_REC)
!         END DO
         else
         DO LC = 1, MMAX
            call gatherWriteSpx (ROP_s(:,LC),array2, array1, uspx+L, NEXT_REC)
         END DO
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP6" FILE         T_g  , T_s
!
      CASE (6)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               CALL UNSHIFT_DP_ARRAY(T_g,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               DO LC = 1, MMAX
                  CALL UNSHIFT_DP_ARRAY(T_s(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR), NEXT_REC)
               END DO
            ELSE
               call OUT_BIN_R(uspx+L,T_g,size(T_g), NEXT_REC)
               DO LC = 1, MMAX
                  call OUT_BIN_R(uspx+L,T_s(:,LC),size(T_s(:,LC)), NEXT_REC)
               END DO
            ENDIF
!           call OUT_BIN_R(uspx+L,T_g,size(T_g), NEXT_REC)
!          DO LC = 1, MMAX
!            call OUT_BIN_R(uspx+L,T_s(:,LC),size(T_s(:,LC)), NEXT_REC)
!          END DO
         else
         call gatherWriteSpx (T_g,array2, array1, uspx+L, NEXT_REC)   !//
          DO LC = 1, MMAX
            call gatherWriteSpx (T_s(:,LC),array2, array1, uspx+L, NEXT_REC)
          END DO
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!
! ".SP7" FILE         X_g, X_s
!
      CASE (7)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               DO NN = 1, NMAX(0)
                  CALL UNSHIFT_DP_ARRAY(X_G(:,NN),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               END DO
               DO LC = 1, MMAX
                  DO NN = 1, NMAX(LC)
                     CALL UNSHIFT_DP_ARRAY(X_s(:,LC,NN),TMP_VAR)
                     call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR), NEXT_REC)
                  ENDDO
               END DO
            ELSE
               DO NN = 1, NMAX(0)
                  call OUT_BIN_R(uspx+L,X_G(:,nn),size(X_G(:,nn)), NEXT_REC)
               END DO
               DO LC = 1, MMAX
                  DO NN = 1, NMAX(LC)
                     call OUT_BIN_R(uspx+L,X_s(:,LC,nn),size(X_s(:,LC,nn)), NEXT_REC)
                  END DO
               END DO
            ENDIF

!           DO NN = 1, NMAX(0)
!             call OUT_BIN_R(uspx+L,X_G(:,nn),size(X_G(:,nn)), NEXT_REC)
!           END DO
!           DO LC = 1, MMAX
!            DO NN = 1, NMAX(LC)
!               call OUT_BIN_R(uspx+L,X_s(:,LC,nn),size(X_s(:,LC,nn)), NEXT_REC)
!            END DO
!           END DO

         else
          DO NN = 1, NMAX(0)
            call gatherWriteSpx (X_G(:,nn),array2, array1, uspx+L, NEXT_REC)
          END DO
          DO LC = 1, MMAX
            DO NN = 1, NMAX(LC)
               call gatherWriteSpx (X_s(:,LC,nn),array2, array1, uspx+L, NEXT_REC)
            END DO
          END DO
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP8" FILE         THETA_m
!
      CASE (8)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               DO LC = 1, MMAX
                  CALL UNSHIFT_DP_ARRAY(THETA_m(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               ENDDO
            ELSE
               DO LC = 1, MMAX
                  call OUT_BIN_R(uspx+L,THETA_m(:,LC),size(THETA_m(:,LC)), NEXT_REC)
               END DO
            ENDIF
!          DO LC = 1, MMAX
!            call OUT_BIN_R(uspx+L,THETA_m(:,LC),size(THETA_m(:,LC)), NEXT_REC)
!         END DO
         else
          DO LC = 1, MMAX
            call gatherWriteSpx (THETA_m(:,LC),array2, array1, uspx+L, NEXT_REC)
          END DO
         end if
         if (myPE.eq.PE_IO .or. bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP9" FILE         Scalar
!
      CASE (9)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               DO LC = 1, Nscalar
                  CALL UNSHIFT_DP_ARRAY(Scalar(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               ENDDO
            ELSE
               DO LC = 1, Nscalar
                 call OUT_BIN_R(uspx+L,Scalar(:,LC),size(Scalar(:,LC)), NEXT_REC)
               END DO
            ENDIF
!          DO LC = 1, Nscalar
!            call OUT_BIN_R(uspx+L,Scalar(:,LC),size(Scalar(:,LC)), NEXT_REC)
!          END DO
         else
          DO LC = 1, Nscalar
            call gatherWriteSpx (Scalar(:,LC),array2, array1, uspx+L, NEXT_REC)
          END DO
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
      CASE (10)  ! Reaction rates

         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (bDist_IO) then
            IF(RE_INDEXING) THEN
               DO LC = 1, nRR
                  CALL UNSHIFT_DP_ARRAY(ReactionRates(:,LC),TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               ENDDO
            ELSE
               DO LC = 1, nRR
                  call OUT_BIN_R(uspx+L,ReactionRates(:,LC),size(ReactionRates(:,LC)), NEXT_REC)
               END DO
            ENDIF
!          DO LC = 1, nRR
!            call OUT_BIN_R(uspx+L,ReactionRates(:,LC),size(ReactionRates(:,LC)), NEXT_REC)
!          END DO
         else
          DO LC = 1, nRR
            call gatherWriteSpx (ReactionRates(:,LC),array2, array1, uspx+L, NEXT_REC)
          END DO
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!
! ".SP11" FILE         turbulence
!
      CASE (11)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         if (K_Epsilon) then
            if (bDist_IO) then
            IF(RE_INDEXING) THEN
               CALL UNSHIFT_DP_ARRAY(K_Turb_G,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               CALL UNSHIFT_DP_ARRAY(E_Turb_G,TMP_VAR)
               call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
            ELSE
               call OUT_BIN_R(uspx+L,K_Turb_G,size(K_Turb_G), NEXT_REC)
               call OUT_BIN_R(uspx+L,E_Turb_G,size(E_Turb_G), NEXT_REC)
            ENDIF
!             call OUT_BIN_R(uspx+L,K_Turb_G,size(K_Turb_G), NEXT_REC)
!             call OUT_BIN_R(uspx+L,E_Turb_G,size(E_Turb_G), NEXT_REC)

            else
             call gatherWriteSpx (K_Turb_G,array2, array1, uspx+L, NEXT_REC)
             call gatherWriteSpx (E_Turb_G,array2, array1, uspx+L, NEXT_REC)
            end if
          end if

           if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
           end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
!
      CASE DEFAULT
            LINE(1) = 'Unknown SPx file index'
            CALL WRITE_ERROR ('WRITE_SPX1', LINE, 1)
            CALL MFIX_EXIT(myPE)
      END SELECT

!//      call unlock_tmp_array
!
      deallocate (array1)
      deallocate (array2)
      deallocate (TMP_VAR)
!
      RETURN
      END SUBROUTINE WRITE_SPX1

      subroutine gatherWriteSpx(VAR, array2, array1, uspxL, NEXT_REC)
        USE geometry
        USE compar           !//
        USE mpi_utility      !//d pnicol : for gatherWriteSpx
        USE sendrecv         !//d pnicol : for gatherWriteSpx
        USE cutcell
        USE in_binary_512
        USE param, only: dimension_3
        USE param1, only: undefined
        IMPLICIT NONE
        integer uspxL, NEXT_REC
        double precision, dimension(ijkmax2) :: array1
        double precision, dimension(ijkmax3) :: array2
        double precision, dimension(DIMENSION_3) :: VAR
        double precision, dimension(:), allocatable :: TMP_VAR

        allocate(TMP_VAR(DIMENSION_3))

!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        IF(RE_INDEXING) THEN
           TMP_VAR = UNDEFINED
           CALL UNSHIFT_DP_ARRAY(VAR,TMP_VAR)
           CALL gather (TMP_VAR,array2,root)
        ELSE
           CALL gather (VAR,array2,root)
        ENDIF
!        call gather (VAR,array2,root)  !//d pnicol

!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        if (myPE.eq.PE_IO) then
           call convert_to_io_dp(array2,array1,ijkmax2)
           CALL OUT_BIN_R (uspxL, array1, IJKMAX2, NEXT_REC)
        end if

        deallocate(TMP_VAR)

      End subroutine gatherWriteSpx



      subroutine gatherWriteSpx_netcdf(VAR, arr1, arr2 , arr4d, ncid, varid , &
                nx,ny,nz,ijkmax2_use , ijkmax3_use)


     USE geometry
     use param, only: dimension_3
     USE compar           !//
     USE mpi_utility      !//d pnicol : for gatherWriteSpx
     USE sendrecv         !//d pnicol : for gatherWriteSpx
     USE MFIX_netcdf
     USE in_binary_512

     IMPLICIT NONE

     integer          :: ncid , varid , nx,ny,nz , ijkmax2_use , ijkmax3_use
     integer          :: ii , jj , kk , ijk

     double precision ::  arr1(ijkmax2_use)
     double precision ::  arr2(ijkmax3_use)
     double precision ::  arr4d(nx,ny,nz,1)
     double precision ::  var(dimension_3)

     call gather(var,arr2,root)
     if (myPE .eq. PE_IO) then
        call convert_to_io_dp(arr2,arr1,ijkmax2_use)

        ijk = 0
        do kk = 1,nz
           do jj = 1,ny
              do ii = 1,nx
                 ijk = ijk + 1
                 arr4d(ii,jj,kk,1) = arr1(ijk)
              end do
           end do
        end do

        call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid,arr4d) )

     end if


    End subroutine gatherWriteSpx_netcdf



      subroutine gatherWriteSpx_netcdf_int(VAR, arr1, arr2 , arr4d, ncid, &
                varid , nx,ny,nz,ijkmax2_use , ijkmax3_use)


     USE geometry
     use param, only: dimension_3
     USE compar           !//
     USE mpi_utility      !//d pnicol : for gatherWriteSpx
     USE sendrecv         !//d pnicol : for gatherWriteSpx
     USE MFIX_netcdf
     USE in_binary_512i

     IMPLICIT NONE

     integer          :: ncid , varid , nx,ny,nz , ijkmax2_use , ijkmax3_use
     integer          :: ii , jj , kk , ijk

     integer ::  arr1(ijkmax2_use)
     integer ::  arr2(ijkmax3_use)
     integer :: arr4d(nx,ny,nz,1)
     integer ::   var(dimension_3)

     call gather(var,arr2,root)
     if (myPE .eq. PE_IO) then
        call convert_to_io_i(arr2,arr1,ijkmax2_use)

        ijk = 0
        do kk = 1,nz
           do jj = 1,ny
              do ii = 1,nx
                 ijk = ijk + 1
                 arr4d(ii,jj,kk,1) = arr1(ijk)
              end do
           end do
        end do

        call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid,arr4d) )

     end if


    End subroutine gatherWriteSpx_netcdf_int

        subroutine copy_d_to_r(darr,rarr,nx,ny,nz)
        implicit none

        integer          :: nx , ny , nz
        double precision :: darr(*)
        real             :: rarr(nx,ny,*)
        integer          :: i , j , k , ijk


        ijk = 0

        do i = 1,nx
           do j = 1,ny
              do k = 1,nz
                 ijk = ijk + 1
                 rarr(i,j,k) = real(darr(ijk))
              end do
           end do
        end do

        return
        end subroutine copy_d_to_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                          write_mesh_netcdf
!
!
        SUBROUTINE write_mesh_netcdf

        USE param
        USE param1
        USE fldvar
        USE geometry
        USE physprop
        USE run
!       USE funits
        USE scalars
!       USE output
        USE rxns
        USE cdist
        USE compar
        USE mpi_utility
        USE MFIX_netcdf
!       USE tmp_array

        implicit none

        integer   :: ncid    , x_dimid , y_dimid , z_dimid  , t_dimid
        integer   :: varid_x , varid_y , varid_z , L , dimids(4)
        integer   :: varid_flag , coords_dimid , varid_coords , coords

        character(LEN=80) :: fname

        double precision, dimension(:) , allocatable :: xloc
        double precision, dimension(:) , allocatable :: yloc
        double precision, dimension(:) , allocatable :: zloc

        integer, dimension(:) , allocatable :: arr1
        integer, dimension(:) , allocatable :: arr2
        integer, dimension(:,:,:,:) , allocatable :: arr4d


        if (.not. MFIX_usingNETCDF()) return

        if (.not. bGlobalNetcdf) return  ! no netCDF writes asked for

        if (.not. bFirst_netcdf_write) return

        if (myPE.eq.PE_IO) then
                allocate ( arr1(ijkmax2))
                allocate ( arr2(ijkmax3))
                allocate ( arr4d(imax2,jmax2,kmax2,1))
                allocate ( xloc(imax2) )
                allocate ( yloc(jmax2) )
                allocate ( zloc(kmax2) )
        else
                allocate ( arr1(1))
                allocate ( arr2(1))
                allocate ( arr4d(1,1,1,1))
                allocate ( xloc(1) )
                allocate ( yloc(1) )
                allocate ( zloc(1) )
        end if

        if (myPE.eq.PE_IO) then
                xloc(1) = -dx(1)
                do L = 2,imax2
                        xloc(L) = xloc(L-1) + dx(L)
                end do

                yloc(1) = -dy(1)
                do L = 2,jmax2
                        yloc(L) = yloc(L-1) + dy(L)
                end do

                zloc(1) = -dz(1)
                do L = 2,kmax2
                        zloc(L) = zloc(L-1) + dz(L)
                end do

                fname = trim(run_name) // "_MESH.nc"
                call MFIX_check_netcdf( MFIX_nf90_create(fname, NF90_64BIT_OFFSET, ncid) )

                call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "x"   , imax2   , x_dimid  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "y"   , jmax2   , y_dimid  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "z"   , kmax2   , z_dimid  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "coordinates"   , 1   , coords_dimid  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "t"   , 1       , t_dimid  ) )  ! 4

                dimids =  (/ x_dimid , y_dimid, z_dimid , t_dimid/)


                ! The dimids array is used to pass the IDs of the dimensions of
                ! the variables. Note that in fortran arrays are stored in
                ! column-major format.
                call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "x"  , NF90_DOUBLE, x_dimid, varid_x  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "y"  , NF90_DOUBLE, y_dimid, varid_y  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "z"  , NF90_DOUBLE, z_dimid, varid_z  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "coordinates"  , NF90_INT, coords_dimid, varid_coords  ) )
                call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "flag"  , NF90_INT, dimids, varid_flag  ) )  ! 9


                call MFIX_check_netcdf( MFIX_nf90_enddef(ncid) )

                coords = 0
                if (COORDINATES .eq. 'CYLINDRICAL') coords = 1

                call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_coords,coords) )
                call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_x,xloc) )
                call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_y,yloc) )
                call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_z,zloc) )
        end if


        ! needs to be called by all processes

        call gatherWriteSpx_netcdf_int(flag, arr1, arr2 , arr4d, ncid, &
               varid_flag , imax2,jmax2,kmax2,ijkmax2 , ijkmax3)



        if (myPE.eq.PE_IO) then
                call MFIX_check_netcdf( MFIX_nf90_close(ncid) )
        end if

        deallocate ( arr1  )
        deallocate ( arr2  )
        deallocate ( arr4d )
        deallocate ( xloc  )
        deallocate ( yloc  )
        deallocate ( zloc  )

        return
        end subroutine write_mesh_netcdf

        SUBROUTINE write_netcdf(L, unit_add, the_time)

        USE param
        USE param1
        USE fldvar
        USE geometry
        USE physprop
        USE run
!       USE funits
        USE scalars
!       USE output
        USE rxns
        USE cdist
        USE compar
        USE mpi_utility
        USE MFIX_netcdf
!       USE tmp_array


        implicit none

        integer :: L , unit_add , I , nn , ii

        integer   :: ncid , x_dimid , y_dimid , z_dimid
        integer   :: t_dimid
        integer   :: dimids(4) , varid_epg , varid_pg

        integer   :: varid_pstar  , varid_ug , varid_vg , varid_wg
        integer   :: varid_tg , varid_x , varid_y , varid_z , varid_t
        integer   :: varid_coords , coords_dimid , coords

        integer   :: varid_us(20) , varid_vs(20) , varid_ws(20)  !! MMAX
        integer   :: varid_rops(20)  , varid_ts(20) !! mmax
        integer   :: varid_thetam(20) !! mmax

        integer   :: varid_xg(20)  ! nmax(0)
        integer   :: varid_xs(20,20)  ! mmax , MAX(nmax(1:mmax))

        integer   :: varid_scalar(20)  ! nscalar
        integer   :: varid_rr(20)      ! nRR

        integer   :: varid_kturbg , varid_eturbg


        character(LEN=80) :: fname, var_name
        character(LEN=9) :: fname_index

        double precision, dimension(:) , allocatable :: arr1
        double precision, dimension(:) , allocatable :: arr2

        double precision, dimension(:,:,:,:) , allocatable :: arr4d


        double precision, dimension(:) , allocatable :: xloc
        double precision, dimension(:) , allocatable :: yloc
        double precision, dimension(:) , allocatable :: zloc

        double precision :: the_time
        logical          :: file_exists


! bWrite_netcdf(1)  : EP_g
! bWrite_netcdf(2)  : P_g
! bWrite_netcdf(3)  : P_star
! bWrite_netcdf(4)  : U_g / V_g / W_g
! bWrite_netcdf(5)  : U_s / V_s / W_s
! bWrite_netcdf(6)  : ROP_s
! bWrite_netcdf(7)  : T_g
! bWrite_netcdf(8)  : T_s
! bWrite_netcdf(9)  : X_g
! bWrite_netcdf(10) : X_s
! bWrite_netcdf(11) : Theta_m
! bWrite_netCDF(12) : Scalar
! bWrite_netCDF(13) : ReactionRates
! bWrite_netCDF(14) : k_turb_g , e_turb_g

        if (.not. MFIX_usingNETCDF()) return
        if (.not. bGlobalNetcdf) return

        call write_mesh_netcdf

        if (myPE.eq.PE_IO .and. .not.bDist_IO) then
           allocate (arr1(ijkmax2))
           allocate (arr2(ijkmax3))
           allocate (arr4d(imax2,jmax2,kmax2,1))
           allocate ( xloc(imax2) )
           allocate ( yloc(jmax2) )
           allocate ( zloc(kmax2) )




        else
           allocate (arr1(1))
           allocate (arr2(1))
           allocate (arr4d(1,1,1,1))
           allocate ( xloc(1) )
           allocate ( yloc(1) )
           allocate ( zloc(1) )
       end if

!       CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
        if (myPE .ne. PE_IO) goto 1234


        xloc(1) = -dx(1)
        do II = 2,imax2
           xloc(II) = xloc(II-1) + dx(II)
        end do


        yloc(1) = -dy(1)
        do II = 2,jmax2
           yloc(II) = yloc(II-1) + dy(II)
        end do

        zloc(1) = -dz(1)
        do II = 2,kmax2
           zloc(II) = zloc(II-1) + dz(II)
        end do




        if (bFirst_netcdf_write .and. MFIX_usingNETCDF()) then
           bFirst_netcdf_write = .false.
           fname = trim(run_name) // '_netcdf_index.txt'
           inquire (file=fname,exist=file_exists)

           if (file_exists) then
               open (unit=11,file=fname,status='old')
               read (11,*) netCDF_file_index
               close (unit=11)
           else
               netCDF_file_index = 0
           end if
        end if

        fname_index = '_xxxxx.nc'
        write (fname_index(2:6),'(i5.5)') netCDF_file_index
        fname = trim(run_name)// fname_index
        call MFIX_check_netcdf( MFIX_nf90_create(fname, NF90_64BIT_OFFSET, ncid) )
        netCDF_file_index = netCDF_file_index + 1

        if (MFIX_usingNETCDF()) then
                fname = trim(run_name) // '_netcdf_index.txt'
                open (unit=11,file=fname,status='unknown')
                write (11,*) netCDF_file_index
                close (unit=11)
        end if


        call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "x"   , imax2   , x_dimid  ) )  ! 1
        call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "y"   , jmax2   , y_dimid  ) )  ! 2
        call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "z"   , kmax2   , z_dimid  ) )  ! 3
        call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "coordinates"   , 1   , coords_dimid  ) )  ! 3
        call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "t"   , 1       , t_dimid  ) )  ! 4


        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "x"  , NF90_DOUBLE, x_dimid, varid_x  ) )  ! 5
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "y"  , NF90_DOUBLE, y_dimid, varid_y  ) )  ! 6
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "z"  , NF90_DOUBLE, z_dimid, varid_z  ) )  ! 7
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "coordinates"  , NF90_INT, coords_dimid, varid_coords  ) )  ! 7
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "t"  , NF90_DOUBLE, t_dimid, varid_t  ) )  ! 8

        dimids =  (/ x_dimid , y_dimid, z_dimid , t_dimid/)

        if (bWrite_netcdf(1)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "EP_g"  , NF90_DOUBLE, dimids, varid_epg  ) )  ! 9
        if (bWrite_netcdf(2)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "P_g"   , NF90_DOUBLE, dimids, varid_pg  ) )   ! 10

        if (bWrite_netcdf(3)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "P_star", NF90_DOUBLE, dimids, varid_pstar) )
        if (bWrite_netcdf(4)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "U_g"   , NF90_DOUBLE, dimids, varid_ug   ) )
        if (bWrite_netcdf(4)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "V_g"   , NF90_DOUBLE, dimids, varid_vg   ) )
        if (bWrite_netcdf(4)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "W_g"   , NF90_DOUBLE, dimids, varid_wg   ) )
        if (bWrite_netcdf(7)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "T_g"   , NF90_DOUBLE, dimids, varid_tg   ) )
        do i = 1,1   ! mmax
           var_name = 'U_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           if (bWrite_netcdf(5)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_us(I)) )

           var_name = 'V_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           if (bWrite_netcdf(5)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_vs(I)) )

           var_name = 'W_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           if (bWrite_netcdf(5)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_ws(I)) )

           var_name = 'ROP_s_xxx'
           write (var_name(7:10),'(i3.3)') I
           if (bWrite_netcdf(6)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_rops(I)) )

           var_name = 'T_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           if (bWrite_netcdf(8)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_ts(I)) )

           var_name = 'Theta_m_xxx'
           write (var_name(9:11),'(i3.3)') I
           if (bWrite_netcdf(11)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_thetam(I)) )

           DO NN = 1, NMAX(i)
              var_name = 'X_s_xxx_xxx'
              write (var_name(5:7) ,'(i3.3)') I
              write (var_name(9:11),'(i3.3)') nn
              if (bWrite_netcdf(10)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_xs(I,nn)) )
           END DO


        end do

        do i = 1,nmax(0)
           var_name = 'X_g_xxx'
           write (var_name(5:7),'(i3.3)') I
           if (bWrite_netcdf(9)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_xg(I)) )
        end do

        do i = 1,nscalar
           var_name = 'Scalar_xxx'
           write (var_name(8:10),'(i3.3)') I
           if (bWrite_netCDF(12)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_scalar(I)) )
        end do

        do i = 1,nRR
           var_name = 'RRates_xxx'
           write (var_name(8:10),'(i3.3)') I
           if (bWrite_netCDF(13)) call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_rr(I)) )
        end do


        if (bWrite_netcdf(14) .and. k_Epsilon) then
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, 'k_turb_g', NF90_DOUBLE, dimids, varid_kturbg) )
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, 'e_turb_g', NF90_DOUBLE, dimids, varid_eturbg) )
        end if



        call MFIX_check_netcdf( MFIX_nf90_enddef(ncid) ) ! 11

 1234   continue
           bFirst_netcdf_write = .false.

!       CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)


        if (myPE .eq. PE_IO) then
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_t,the_time) )   ! 12
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_x,xloc) )
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_y,yloc) )
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_z,zloc) )
           coords = 0
           if (COORDINATES .eq. 'CYLINDRICAL') coords = 1
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_coords,coords) )
        end if

        if (bWrite_netcdf(1)) then

            call gatherWriteSpx_netcdf(EP_g, arr1, arr2 , arr4d, ncid, varid_epg , &
                imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

        end if

        if (bWrite_netcdf(2)) then

            call gatherWriteSpx_netcdf(P_g, arr1, arr2 , arr4d, ncid, varid_pg , &
                imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

        end if




        if (bWrite_netcdf(3)) then

            call gatherWriteSpx_netcdf(P_star, arr1, arr2 , arr4d, ncid, varid_pstar , &
                imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

        end if

        if (bWrite_netcdf(4)) then

            call gatherWriteSpx_netcdf(U_g, arr1, arr2 , arr4d, ncid, varid_ug , &
                imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

            call gatherWriteSpx_netcdf(V_g, arr1, arr2 , arr4d, ncid, varid_vg , &
                imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

            call gatherWriteSpx_netcdf(W_g, arr1, arr2 , arr4d, ncid, varid_wg , &
                imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

        end if

        if (bWrite_netcdf(7)) then

            call gatherWriteSpx_netcdf(t_g, arr1, arr2 , arr4d, ncid, varid_tg , &
                imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

        end if


        do i = 1,1   ! mmax

           if (bWrite_netcdf(5)) then

              call gatherWriteSpx_netcdf(u_s(:,i) , arr1, arr2 , arr4d, ncid, varid_us(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

              call gatherWriteSpx_netcdf(v_s(:,i) , arr1, arr2 , arr4d, ncid, varid_vs(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

              call gatherWriteSpx_netcdf(w_s(:,i) , arr1, arr2 , arr4d, ncid, varid_ws(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

           end if

           if (bWrite_netcdf(6)) then

              call gatherWriteSpx_netcdf(ROP_s(:,i) , arr1, arr2 , arr4d, ncid, varid_rops(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

           end if

           if (bWrite_netcdf(8)) then

              call gatherWriteSpx_netcdf(T_s(:,i) , arr1, arr2 , arr4d, ncid, varid_ts(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

           end if

           if (bWrite_netcdf(11)) then

              call gatherWriteSpx_netcdf(Theta_m(:,i) , arr1, arr2 , arr4d, ncid, varid_thetam(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

           end if

           if (bWrite_netcdf(10)) then
              do nn = 1,nmax(i)

                 call gatherWriteSpx_netcdf(X_s(:,i,NN) , arr1, arr2 , arr4d, ncid, varid_xs(i,nn) , &
                      imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

              end do
           end if


        end do

        if (bWrite_netcdf(9)) then
           do i = 1,nmax(0)

               call gatherWriteSpx_netcdf(X_g(:,i) , arr1, arr2 , arr4d, ncid, varid_xg(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

           end do
        end if

       if (bWrite_netcdf(12)) then
          do i = 1,nscalar

             call gatherWriteSpx_netcdf(Scalar(:,i) , arr1, arr2 , arr4d, ncid, varid_scalar(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

           end do
        end if


        if (bWrite_netcdf(13)) then
           do i = 1,nRR

             call gatherWriteSpx_netcdf(ReactionRates(:,i) , arr1, arr2 , arr4d, ncid, varid_rr(i) , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

          end do
        end if

        if (bWrite_netcdf(14) .and. k_Epsilon) then

           call gatherWriteSpx_netcdf(k_turb_g , arr1, arr2 , arr4d, ncid, varid_kturbg , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

           call gatherWriteSpx_netcdf(e_turb_g , arr1, arr2 , arr4d, ncid, varid_eturbg , &
                   imax2,jmax2,kmax2,ijkmax2 , ijkmax3)

        end if

        ! Close the file. This frees up any internal netCDF resources
        ! associated with the file, and flushes any buffers.
!       CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
        if (myPE .eq. PE_IO) then
           call MFIX_check_netcdf( MFIX_nf90_close(ncid) )
        end if

        deallocate (arr1)
        deallocate (arr2)
        deallocate (arr4d)
        deallocate (xloc)
        deallocate (yloc)
        deallocate (zloc)

        return
        end subroutine write_netcdf
