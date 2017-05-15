!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_RES1                                             C
!  Purpose: write out the time-dependent restart records               C
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
!  Variables referenced: TIME, NSTEP, EP_g, P_g, P_star, RO_g, ROP_g   C
!                        T_g, T_s, U_g, V_g, W_g, ROP_s, U_s    C
!                        V_s, W_s, IJKMAX2                             C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, N, NEXT_REC                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_RES1
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE cdist
      USE compar           !//
      USE energy
      USE fldvar
      USE funits
      USE geometry
      USE machine, only: flush_res
      USE mpi_utility      !//d pnicol : for gather
      USE output
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE scalars
      USE sendrecv         !//d pnicol : for gather
!//d pnicol  ... not needed    USE tmp_array
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
!
!
!
!
!//d pnicol : allocate arrays for gather/convert_to_io_dp
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)


!             loop counter
      INTEGER :: LC, NN
!
!             pointer to first time-dependent record in restart file
      INTEGER :: NEXT_REC
!-----------------------------------------------
!

!//d pnicol
!      if (myPE.eq.PE_IO .and. .not.distio) then
         allocate (array1(ijkmax2))
         allocate (array2(ijkmax3))
!      else
!         allocate (array1(1))
!         allocate (array2(1))
!      end if


      if (myPE.eq.PE_IO .or. bDist_IO) then
         READ (UNIT_RES, REC=3) NEXT_REC
         WRITE (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP
         NEXT_REC = NEXT_REC + 1
      end if
!
!\\SP Local Send Receive - need to be moved to source later!!

      if (.not. bDist_IO) then

          call send_recv(EP_g,2)
          call send_recv(P_g,2)
          call send_recv(P_star,2)
          call send_recv(RO_g,2)
          call send_recv(ROP_g,2)
          call send_recv(X_g,2)
          call send_recv(T_g,2)
          call send_recv(U_g,2)
          call send_recv(V_g,2)
          call send_recv(W_g,2)
          call send_recv(ROP_S,2)
          call send_recv(T_S,2)
          call send_recv(U_S,2)
          call send_recv(V_S,2)
          call send_recv(W_S,2)
          call send_recv(THETA_M,2)
          call send_recv(X_S,2)
          if(NScalar > 0)call send_recv(Scalar,2)
          if(K_Epsilon) THEN
              call send_recv(K_Turb_G,2)
              call send_recv(E_Turb_G,2)
          endif
          call send_recv(GAMA_RG,2)
          call send_recv(T_RG,2)
          call send_recv(GAMA_RS,2)
          call send_recv(T_RS,2)
          if(nRR > 0)call send_recv(ReactionRates,2)

      end if


      call gatherWriteRes (EP_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (P_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (P_star,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (RO_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (ROP_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (T_g,array2, array1, NEXT_REC)  !//d pnicol
!
      DO NN = 1, NMAX(0)
            call gatherWriteRes (X_g(:,nn),array2, array1, NEXT_REC)  !//d pnicol
      END DO
!
      call gatherWriteRes (U_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (V_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (W_g,array2, array1, NEXT_REC)  !//d pnicol
!
      DO LC = 1, MMAX
!
        call gatherWriteRes (ROP_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol

        IF(ANY(SOLVE_ROs)) &
            call gatherWriteRes (RO_S(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (T_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (U_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (V_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (W_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (THETA_M(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
         DO NN = 1, NMAX(LC)
            call gatherWriteRes (X_s(:,LC,NN),array2, array1, NEXT_REC)  !//d pnicol
         END DO
      END DO
!
!     Version 1.3

      DO LC = 1, NScalar
            call gatherWriteRes (Scalar(:,LC),array2, array1, NEXT_REC)  !//d pnicol
      END DO
!
!     Version 1.4 -- write radiation variables in write_res1
      call gatherWriteRes (GAMA_RG,array2, array1, NEXT_REC)  !//d pnicol

      call gatherWriteRes (T_RG,array2, array1, NEXT_REC)  !//d pnicol

      DO LC = 1, MMAX
        call gatherWriteRes (GAMA_RS(1,LC),array2, array1, NEXT_REC)  !//d pnicol

        call gatherWriteRes (T_RS(1,LC),array2, array1, NEXT_REC)  !//d pnicol
      ENDDO

!
!     Version 1.5
      DO LC = 1, nRR
            call gatherWriteRes (ReactionRates(:,LC),array2, array1, NEXT_REC)  !//d pnicol
      END DO
!
!     Version 1.6

      if (K_epsilon) then
            call gatherWriteRes (K_turb_G,array2, array1, NEXT_REC)  !//d pnicol
          call gatherWriteRes (E_turb_G,array2, array1, NEXT_REC)
      endif
!---------------------------------------------------------------------

      if ( (myPE.eq.PE_IO .and. .not.bDist_IO) .or. bDist_IO) then
           CALL FLUSH_res (UNIT_RES)
        end if

!      call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
      deallocate (array1)  !//d pnicol
      deallocate (array2)  !//d pnicol
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here

      call write_res1_netcdf
!
      RETURN
      END SUBROUTINE WRITE_RES1

      subroutine gatherWriteRes(VAR, array2, array1, NEXT_REC)

      USE geometry
      USE funits
      USE cdist
      USE compar           !//
      USE mpi_utility      !//d pnicol : for gather
      USE sendrecv         !//d pnicol : for gather

      USE cutcell
      USE in_binary_512
      USE param, only: dimension_3

      IMPLICIT NONE

      double precision, dimension(ijkmax2) :: array1
      double precision, dimension(ijkmax3) :: array2
      double precision, dimension(DIMENSION_3) :: VAR,TMP_VAR

      INTEGER :: NEXT_REC

!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (.not.bDist_IO) then


         IF(RE_INDEXING) THEN
            CALL UNSHIFT_DP_ARRAY(VAR,TMP_VAR)
            CALL gather (TMP_VAR,array2,root)
         ELSE
            CALL gather (VAR,array2,root)
         ENDIF
!         call gather (VAR,array2,root)  !//d pnicol

!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
         if (myPE.eq.PE_IO) then
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC)
         end if

      else

         IF(RE_INDEXING) THEN
            CALL UNSHIFT_DP_ARRAY(VAR,TMP_VAR)
            CALL OUT_BIN_512 (UNIT_RES, TMP_VAR, size(TMP_VAR), NEXT_REC)
         ELSE
            CALL OUT_BIN_512 (UNIT_RES, var, size(var), NEXT_REC)
         ENDIF
!         CALL OUT_BIN_512 (UNIT_RES, var, size(var), NEXT_REC)

      end if

      End subroutine gatherWriteRes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!                                      write_res1_netcdf                       !
!                                                                              !
      subroutine write_res1_netcdf

        USE param
        USE param1
        USE fldvar
        USE geometry
        USE physprop
        USE run
        USE scalars
        USE rxns
        USE cdist
        USE compar
        USE mpi_utility
        USE MFIX_netcdf
        USE energy
        USE in_binary_512
!       USE tmp_array


        implicit none

        integer :: I , nn

        integer   :: ncid , xyz_dimid
        integer   :: varid_time , t_dimid
        integer   :: dimids(1) , varid_epg , varid_pg,dims_time(1)
        integer   :: varid_pstar  , varid_ug , varid_vg , varid_wg
        integer   :: varid_tg  , varid_ropg
        integer   :: varid_rog , varid_gamaRG , varid_TRG
        integer   :: varid_gamaRS(20) , varid_TRS(20)

        integer   :: varid_us(20) , varid_vs(20) , varid_ws(20)  !! MMAX
        integer   :: varid_rops(20)  , varid_ts(20) !! mmax
        integer   :: varid_thetam(20) !! mmax

        integer   :: varid_xg(20)  ! nmax(0)
        integer   :: varid_xs(20,20)  ! mmax , MAX(nmax(1:mmax))

        integer   :: varid_scalar(20)  ! nscalar
        integer   :: varid_rr(20)      ! nRR

        integer   :: varid_kturbg , varid_eturbg


        character(LEN=80) :: fname, var_name

        double precision, allocatable :: arr1(:)
        double precision, allocatable :: arr2(:)

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

        return ! for now ... do not write restart files using netCDF

        if (myPE .eq. PE_IO) then
           allocate (arr1(ijkmax2))
           allocate (arr2(ijkmax3))
        else
           allocate (arr1(1))
           allocate (arr2(1))
        end if

!       call mpi$barrier(MPI_COMM_WORLD,mpierr)
        if (myPE .ne. PE_IO) goto 1234

        fname = trim(run_name) // "_RES1.nc"
        call MFIX_check_netcdf( MFIX_nf90_create(fname, NF90_CLOBBER, ncid) )

        call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "xyz", ijkmax2, xyz_dimid) )
        call MFIX_check_netcdf( MFIX_nf90_def_dim(ncid, "t"  ,         1        ,   t_dimid) )

        ! The dimids array is used to pass the IDs of the dimensions of
        ! the variables. Note that in fortran arrays are stored in
        ! column-major format.
        dims_time(1) = t_dimid
        dimids =  (/ xyz_dimid /)

        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid ,"time"  , NF90_DOUBLE ,dims_time , varid_time ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "EP_g"  , NF90_DOUBLE, dimids    , varid_epg  ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "P_g"   , NF90_DOUBLE, dimids    , varid_pg   ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "P_star", NF90_DOUBLE, dimids    , varid_pstar) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "U_g"   , NF90_DOUBLE, dimids    , varid_ug   ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "V_g"   , NF90_DOUBLE, dimids    , varid_vg   ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "W_g"   , NF90_DOUBLE, dimids    , varid_wg   ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "T_g"   , NF90_DOUBLE, dimids    , varid_tg   ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "ROP_g" , NF90_DOUBLE, dimids    , varid_ropg ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "RO_g"  , NF90_DOUBLE, dimids    , varid_rog  ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "gamaRG"  , NF90_DOUBLE, dimids    , varid_gamaRG  ) )
        call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, "TRG"  , NF90_DOUBLE, dimids       , varid_TRG  ) )

        do i = 1,1   ! mmax
           var_name = 'U_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_us(I)) )

           var_name = 'V_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_vs(I)) )

           var_name = 'W_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_ws(I)) )

           var_name = 'ROP_s_xxx'
           write (var_name(7:10),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_rops(I)) )

           var_name = 'T_s_xxx'
           write (var_name(5:7),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_ts(I)) )

           var_name = 'Theta_m_xxx'
           write (var_name(9:11),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_thetam(I)) )

           var_name = 'gamaRS_xxx'
           write (var_name(8:10),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_gamaRS(I)) )

           var_name = 'TRS_xxx'
           write (var_name(5:7),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_trs(I)) )

           DO NN = 1, NMAX(i)
              var_name = 'X_s_xxx_xxx'
              write (var_name(5:7) ,'(i3.3)') I
              write (var_name(9:11),'(i3.3)') nn
              call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_xs(I,nn)) )
           END DO


        end do

        do i = 1,nmax(0)
           var_name = 'X_g_xxx'
           write (var_name(5:7),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_xg(I)) )
        end do

        do i = 1,nscalar
           var_name = 'Scalar_xxx'
           write (var_name(8:10),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_scalar(I)) )
        end do

        do i = 1,nRR
           var_name = 'RRates_xxx'
           write (var_name(8:10),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, var_name, NF90_DOUBLE, dimids, varid_rr(I)) )
        end do


        if (k_Epsilon) then
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, 'k_turb_g', NF90_DOUBLE, dimids, varid_kturbg) )
           call MFIX_check_netcdf( MFIX_nf90_def_var(ncid, 'e_turb_g', NF90_DOUBLE, dimids, varid_eturbg) )
        end if


        call MFIX_check_netcdf( MFIX_nf90_enddef(ncid) )

 1234   continue
!       call mpi$barrier(MPI_COMM_WORLD,mpierr)

        if (myPE .eq. PE_IO) then
!           call check_netcdf( nf90_put_var(ncid,varid_time,the_time) )   ! ???
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid,varid_time,time) )       ! ???
        end if

        call gather(EP_g,arr2,root)

        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_epg, arr1) )
        end if

        call gather(P_g,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_pg, arr1) )
        end if

        call gather(P_Star,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_pstar, arr1) )
        end if

        call gather(U_g,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_ug, arr1) )
        end if

        call gather(V_g,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_vg, arr1) )
        end if

        call gather(W_g,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_wg, arr1) )
        end if

        call gather(T_g,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_tg, arr1) )
        end if

        call gather(gama_rg,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_gamarg, arr1) )
        end if

        call gather(ro_g,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_rog, arr1) )
        end if

        call gather(rop_g,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_ropg, arr1) )
        end if

        call gather(t_rg,arr2,root)
        if (myPE .eq. PE_IO) then
           call convert_to_io_dp(arr2,arr1,ijkmax2)
           call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_trg, arr1) )
        end if

        do i = 1,1   ! mmax

           call gather(U_s(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_us(i), arr1) )
           end if

           call gather(V_s(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_vs(i), arr1) )
           end if


           call gather(W_s(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_ws(i), arr1) )
           end if

           call gather(ROP_s(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_rops(i), arr1 ) )
           end if

           call gather(T_s(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_ts(i), arr1 ) )
           end if

           call gather(Theta_m(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_thetam(i), arr1 ) )
           end if

           call gather(gama_rs(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_gamaRS(i), arr1 ) )
           end if

           call gather(T_rs(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_TRS(i) , arr1 ) )
           end if

           do nn = 1,nmax(i)
              call gather(X_s(:,i,NN),arr2,root)
              if (myPE .eq. PE_IO) then
                 call convert_to_io_dp(arr2,arr1,ijkmax2)
                 call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_xs(i,NN), arr1 ) )
              end if
           end do

        end do

        do i = 1,nmax(0)
           call gather(X_g(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_xg(i), arr1 ) )
           end if
        end do

        do i = 1,nscalar
           call gather(Scalar(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_scalar(i), arr1 ) )
           end if
        end do

        do i = 1,nRR
           call gather(ReactionRates(:,i),arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_rr(i), arr1 ) )
           end if
        end do

        if (k_epsilon) then
           call gather(k_turb_g,arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_kturbg, arr1) )
           end if

           call gather(e_turb_g,arr2,root)
           if (myPE .eq. PE_IO) then
              call convert_to_io_dp(arr2,arr1,ijkmax2)
              call MFIX_check_netcdf( MFIX_nf90_put_var(ncid, varid_eturbg, arr1) )
           end if
        end if

 2222   continue

        ! Close the file. This frees up any internal netCDF resources
        ! associated with the file, and flushes any buffers.
        if (myPE .eq. PE_IO) then
           call MFIX_check_netcdf( MFIX_nf90_close(ncid) )
        end if

        deallocate (arr1)
        deallocate (arr2)

        return

      end subroutine write_res1_netcdf



