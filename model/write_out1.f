!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_OUT1                                             C
!  Purpose: write out the field variables to standard output           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: TIME, P_g, EP_g, RO_g, ROP_g, MMAX, ROP_s     C
!                        T_g, T_s, U_g, V_g, W_g, U_s, V_s, W_s C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, N                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_OUT1
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE fldvar
      USE run
      USE scalars
      USE funits
      USE rxns
      USE compar             !//d
      USE mpi_utility        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LC, NN
!-----------------------------------------------
!
      double precision, allocatable :: array1(:)    !//d
!


      if (myPE == PE_IO) then
         allocate (array1(ijkmax3))     !//d
      else
         allocate (array1(1))           !//d
      end if
!
!             form feed character = CHAR(12)
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1000) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (P_g,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'P_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1050) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (P_star,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'P_star')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1100) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (EP_g,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'EP_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1200) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (RO_g,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'RO_g')
!
      DO LC = 1, MMAX
         if (myPE == PE_IO) WRITE (UNIT_OUT, 1400) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (ROP_s(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'ROP_s')
         if (LC == 2) then
            if (myPE == PE_IO) WRITE (UNIT_OUT, *) 'solid phase ', LC, 'at t= ', TIME
            call gather (RO_S(:,LC),array1,root)
            if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'RO_S')
            call gather (ERR_ARRAY(:,LC),array1,root)
            if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'ERR_EPS')
         endif
      END DO
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1500) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (T_g,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'T_g')
!
      DO LC = 1, MMAX
         if (myPE == PE_IO) WRITE (UNIT_OUT, 1600) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (T_s(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'T_s')
      END DO
      IF (SPECIES_EQ(0)) THEN
         DO NN = 1, NMAX(0)
            if (myPE == PE_IO) WRITE (UNIT_OUT, 1710) CHAR(12), NN, TIME
!           call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
            call gather (X_g(:,NN),array1,root)    !//
!           call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
            if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'X_g')
         END DO
      ENDIF
!
      DO LC = 1, MMAX
         IF (SPECIES_EQ(LC)) THEN
            DO NN = 1, NMAX(LC)
               if (myPE == PE_IO) WRITE (UNIT_OUT, 1720) CHAR(12), LC, NN, TIME
!              call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
               call gather (X_s(:,LC,NN),array1,root)    !//
!              call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
               if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'X_s')
            END DO
         ENDIF
      END DO
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1800) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (U_g,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!//SP - Changed U_G to array1
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'U_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1900) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (V_g,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!//SP - Changed V_G to array1
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'V_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 2000) CHAR(12), TIME
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      call gather (W_g,array1,root)    !//
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!//SP - Changed W_G to array1
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'W_g')
!
      DO LC = 1, MMAX
         if (myPE == PE_IO) WRITE (UNIT_OUT, 2100) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (U_s(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'U_s')
!
         if (myPE == PE_IO) WRITE (UNIT_OUT, 2200) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (V_s(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'V_s')
!
         if (myPE == PE_IO) WRITE (UNIT_OUT, 2300) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (W_s(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'W_s')
!
!         IF(GRANULAR_ENERGY)THEN
         if (myPE == PE_IO) WRITE (UNIT_OUT, 2400) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (THETA_m(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'Theta_m')
      END DO
      if (myPE == PE_IO) WRITE (UNIT_OUT, '(/1X,1A1)') CHAR(12)
      IF (CALL_USR) CALL USR_WRITE_OUT1
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here

      IF(NScalar /= 0) THEN
        DO LC = 1, NScalar
         if (myPE == PE_IO) WRITE (UNIT_OUT, 2500) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (Scalar(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'Scalar')
        END DO
      ENDIF

      IF(K_Epsilon) THEN
         if (myPE == PE_IO) WRITE (UNIT_OUT, 2600) CHAR(12), TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (K_Turb_G,array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'K_Turb_G')
         call gather (E_Turb_G,array1,root)    !//
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'E_Turb_G')
      ENDIF
      IF(nRR /= 0) THEN
        DO LC = 1, nRR
         if (myPE == PE_IO) WRITE (UNIT_OUT, 2500) CHAR(12), LC, TIME
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         call gather (ReactionRates(:,LC),array1,root)    !//
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'RRates')
        END DO
      ENDIF

      deallocate(array1)  !//

!
!             form feed character = CHAR(12)
      if (myPE == PE_IO) WRITE (UNIT_OUT, '(/1X,1A1)') CHAR(12)
      IF (CALL_USR) CALL USR_WRITE_OUT1
      RETURN
 1000 FORMAT(1X,A1,/5X,'--- Gas pressure (P_g) at time ',G12.5,' ---',2/)
 1050 FORMAT(1X,A1,/5X,'--- Solids pressure (P_star) at time ',G12.5,' ---',2/)
 1100 FORMAT(1X,A1,/5X,'--- Void fraction (EP_g) at time ',G12.5,' ---',2/)
 1200 FORMAT(1X,A1,/5X,'--- Gas density (RO_g) at time ',G12.5,' ---',2/)
 1400 FORMAT(1X,A1,/5X,'--- Solids Phase-',I1,' density x volume',&
         ' fraction (ROP_s) at time ',G12.5,' ---',2/)
 1500 FORMAT(1X,A1,/5X,'--- Gas temperature (T_g) at time ',G12.5,' ---',2/)
 1600 FORMAT(1X,A1,/5X,'--- Solids Phase-',I1,' temperature (T_s)',' at time ',&
         G12.5,' ---',2/)
 1710 FORMAT(1X,A1,/5X,'--- Mass fraction of gas species (X_g) ',&
         I2,' at time ',&
         G12.5,' ---',2/)
 1720 FORMAT(1X,A1,/5X,'--- Mass fraction of solids-',I1,' species (X_s)',I2,&
         ' at time ',G12.5,' ---',2/)
 1800 FORMAT(1X,A1,/5X,'--- X-component of gas velocity (U_g) at time ',G12.5,&
         ' ---',2/)
 1900 FORMAT(1X,A1,/5X,'--- Y-component of gas velocity (V_g) at time ',G12.5,&
         ' ---',2/)
 2000 FORMAT(1X,A1,/5X,'--- Z-component of gas velocity (W_g) at time ',G12.5,&
         ' ---',2/)
 2100 FORMAT(1X,A1,/5X,'--- X-component of Solids Phase-',I1,&
         ' velocity (U_s) at time ',G12.5,' ---',2/)
 2200 FORMAT(1X,A1,/5X,'--- Y-component of Solids Phase-',I1,&
         ' velocity (V_s) at time ',G12.5,' ---',2/)
 2300 FORMAT(1X,A1,/5X,'--- Z-component of Solids Phase-',I1,&
         ' velocity (W_s) at time ',G12.5,' ---',2/)
 2400 FORMAT(1X,A1,/5X,'--- Granular temperature of Solids Phase-',I1,&
         ' (Theta_m) at time ',G12.5,' ---',2/)
 2500 FORMAT(1X,A1,/5X,'--- Scalar Field-',I2, ' (Scalar) at time ',G12.5,' ---',2/)
 2600 FORMAT(1X,A1,/5X,'--- Turbulence Field-', ' (K-Epsilon) at time ',G12.5,' ---',2/)
      END SUBROUTINE WRITE_OUT1
