!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_RES1                                              C
!  Purpose: read in the time-dependent restart records                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2, MMAX, DT                             C
!  Variables modified: TIME, NSTEP, EP_g, P_g, P_star, RO_g            C
!                      ROP_g, T_g, T_s,  U_g, V_g, W_g, ROP_s    C
!                      U_s, V_s, W_s                                   C
!                                                                      C
!  Local variables: TIME_READ, LC, NEXT_REC                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE READ_RES1
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE physprop
      USE run
      USE rxns
      USE scalars
      USE funits
      USE energy
      USE compar
      USE cdist
      USE mpi_utility
      USE sendrecv
      USE in_binary_512

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
!             loop counter
      INTEGER LC
!
!                      Local species index
      INTEGER          NN
!
!             pointer to the next record
      INTEGER NEXT_REC
!
!                file version id
      CHARACTER(LEN=512) :: VERSION
!
!                version number
      REAL       VERSION_NUMBER
!
!                      Dummy array
      DOUBLE PRECISION DT_SAVE

!//PAR_I/O declare global scratch arrays
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)
!-----------------------------------------------
!
      if (myPE .eq. PE_IO .or. .not.bStart_with_one_res) then
         allocate (array1(ijkmax2))
         allocate (array2(ijkmax3))
      else
         allocate (array1(1))
         allocate (array2(1))
      end if

!      call MPI_barrier(MPI_COMM_WORLD,mpierr)
!
!     Use DT from data file if DT_FAC is set to 1.0
      IF (DT_FAC == ONE) DT_SAVE = DT
!

!//PAR_I/O only PE_IO reads the restart file
      if (myPE == PE_IO .or. (bDist_IO .and. .not.bStart_with_one_RES)) then
         READ (UNIT_RES, REC=1) VERSION
         READ (VERSION(6:512), *) VERSION_NUMBER

         READ (UNIT_RES, REC=3) NEXT_REC
         IF (VERSION_NUMBER >= 1.12) THEN
            READ (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP
         ELSE
            READ (UNIT_RES, REC=NEXT_REC) TIME, NSTEP
         ENDIF
         NEXT_REC = NEXT_REC + 1
      end if



      if (.not.bDist_IO  .or. bStart_with_one_RES) then
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)

         call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c
         call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r
         call bcast(TIME, PE_IO)           !//PAR_I/O BCAST0d
         call bcast(NSTEP, PE_IO)          !//PAR_I/O BCAST0i
         if (VERSION_NUMBER >= 1.12) call bcast(DT, PE_IO)   !//PAR_I/O BCAST0d
        end if
!      call MPI_barrier(MPI_COMM_WORLD,mpierr)

!AE TIME 091501 Store the timestep counter level at the begin of RESTART run
        NSTEPRST = NSTEP

! for now ... do not do RES1 file in netCDF format
!       call read_res1_netcdf
!        goto 999
!


!
      call readScatterRes(EP_G,array2, array1, 0, NEXT_REC)

      call readScatterRes(P_G,array2, array1, 0, NEXT_REC)

      call readScatterRes(P_STAR,array2, array1, 1, NEXT_REC)

      call readScatterRes(RO_G,array2, array1, 0, NEXT_REC)

      call readScatterRes(ROP_G,array2, array1, 0, NEXT_REC)

      call readScatterRes(T_G,array2, array1, 1, NEXT_REC)
!

      IF (VERSION_NUMBER < 1.15) THEN
         call readScatterRes (T_s(:,1),array2, array1, 1, NEXT_REC)

         IF (MMAX >= 2) THEN
            call readScatterRes (T_s(:,2),array2, array1, 1, NEXT_REC)
         ELSE
            if (myPE == PE_IO) &
                   CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC)
         ENDIF
      ENDIF

      IF (VERSION_NUMBER >= 1.05) THEN
         DO NN = 1, NMAX(0)
            call readScatterRes (X_g(:,nn),array2, array1, 1, NEXT_REC)
         END DO
      ENDIF

      call readScatterRes(U_G, array2, array1, 0, NEXT_REC)
      call readScatterRes(V_G, array2, array1, 0, NEXT_REC)
      call readScatterRes(W_G, array2, array1, 0, NEXT_REC)
!
      DO LC = 1, MMAX
         call readScatterRes(ROP_S(:,LC), array2, array1, 0, NEXT_REC)

         IF(ANY(SOLVE_ROs)) &
            CALL readScatterRes(RO_S(:,LC), array2, array1, 0, NEXT_REC)

         IF (VERSION_NUMBER >= 1.15) THEN
            call readScatterRes(T_S(:,LC), array2, array1, 1, NEXT_REC)
         END IF
         call readScatterRes(U_S(:,LC), array2, array1, 0, NEXT_REC)
         call readScatterRes(V_S(:,LC), array2, array1, 0, NEXT_REC)
         call readScatterRes(W_S(:,LC), array2, array1, 0, NEXT_REC)

         IF (VERSION_NUMBER >= 1.2) then
            call readScatterRes(THETA_M(:,LC), array2, array1, 1, NEXT_REC)
         end if
         IF (VERSION_NUMBER >= 1.05) THEN
            DO NN = 1, NMAX(LC)
               call readScatterRes(X_S(:,LC,NN), array2, array1, 1, NEXT_REC)
            END DO
         ENDIF
      END DO

      IF (VERSION_NUMBER >= 1.3) THEN
        DO NN = 1, NScalar
          call readScatterRes(Scalar(:,NN), array2, array1, 1, NEXT_REC)
        END DO
      ENDIF

      IF (VERSION_NUMBER >= 1.4) THEN
        call readScatterRes(GAMA_RG, array2, array1, 1, NEXT_REC)

        call readScatterRes(T_RG, array2, array1, 0, NEXT_REC)

        DO LC = 1, MMAX
          call readScatterRes(GAMA_RS(1,LC), array2, array1, 1, NEXT_REC)

          call readScatterRes(T_RS(1,LC), array2, array1, 0, NEXT_REC)

        ENDDO
      ELSE
        GAMA_RG(:)   = ZERO
        T_RG (:)     = ZERO
        GAMA_RS(:,:) = ZERO
        T_RS(:,:)    = ZERO
      ENDIF


      IF (VERSION_NUMBER >= 1.5) THEN
        DO NN = 1, nRR
          call readScatterRes(ReactionRates(:,NN), array2, array1, 1, NEXT_REC)
        END DO
      ENDIF

      IF (VERSION_NUMBER >= 1.6 .AND. K_Epsilon) THEN
          call readScatterRes(K_Turb_G, array2, array1, 1, NEXT_REC)
          call readScatterRes(E_Turb_G, array2, array1, 1, NEXT_REC)
      ENDIF
!------------------------------------------------------------------------
!

!      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      deallocate( array1 )
      deallocate( array2 )
!      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (.not.bDist_IO .or. bStart_with_one_RES) then
         call send_recv(rop_g)
         call send_recv(ro_g)
         call send_recv(rop_s)
      end if


      IF (DT_FAC == ONE) DT = DT_SAVE
!

!     We may no longer need PATCH_AFTER_RESTART
!     CALL PATCH_AFTER_RESTART

      RETURN
      END SUBROUTINE READ_RES1

      subroutine readScatterRes(VAR, array2, array1, init, NEXT_REC)
        use param1, only: zero, undefined
        use param, only: dimension_3
        USE geometry
        USE funits
        USE compar
        USE cdist
        USE mpi_utility
        USE sendrecv
        USE in_binary_512
        IMPLICIT NONE
        double precision, dimension(ijkmax2) :: array1
        double precision, dimension(ijkmax3) :: array2
        double precision, dimension(DIMENSION_3) :: VAR
        INTEGER :: init  ! define VAR initialization, 0: undefin, 1: zero
        INTEGER :: NEXT_REC

!// Reset global scratch arrays
        if( init==0 ) then
          array1(:) = Undefined
          array2(:) = Undefined
        else
          array1(:) = zero
          array2(:) = zero
        endif

        if (.not.bDist_IO .or. bStart_with_one_RES) then
         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC)
            CALL convert_from_io_dp(array1, array2, IJKMAX2)
         end if
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)
         call scatter(VAR, array2, PE_IO)
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)
      else
         CALL IN_BIN_512 (UNIT_RES, var, size(var) , NEXT_REC)
        end if

      End subroutine readScatterRes


      subroutine readScatterRes_netcdf(VAR, array2, array1, ncid , varid)
      USE param, only: dimension_3
      USE geometry
      USE compar
      USE cdist
      USE mpi_utility
      USE sendrecv
      USE MFIX_netcdf
      USE in_binary_512

      IMPLICIT NONE

      double precision, dimension(ijkmax2)     :: array1
      double precision, dimension(ijkmax3)     :: array2
      double precision, dimension(DIMENSION_3) :: VAR

      integer :: ncid , varid



      if (myPE .eq. PE_IO) then
         call MFIX_check_netcdf( MFIX_nf90_get_var(ncid , varid   , array1   ) )
         call convert_from_io_dp(array1,array2,ijkmax2)
      end if

!     call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(var,array2,PE_IO)
!     call MPI_barrier(MPI_COMM_WORLD,mpierr)

      End subroutine readScatterRes_netcdf





!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization: array1, array2
!// 400 Added sendrecv module and send_recv calls for COMMunication
!// 400 Added mpi_utility module and other global reduction (bcast) calls



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!                                      read_res1_netcdf                       !
!                                                                              !
      subroutine read_res1_netcdf

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
        USE energy


        implicit none

        integer :: I , nn

        integer   :: ncid
        integer   :: varid_time
        integer   :: varid_epg , varid_pg
        integer   :: varid_pstar  , varid_ug , varid_vg , varid_wg
        integer   :: varid_tg
        integer   :: varid_rog , varid_gamaRG , varid_TRG
        integer   :: varid_gamaRS(20) , varid_TRS(20) , varid_ropg

        integer   :: varid_us(20) , varid_vs(20) , varid_ws(20)  !! MMAX
        integer   :: varid_rops(20)  , varid_ts(20) !! mmax
        integer   :: varid_thetam(20) !! mmax

        integer   :: varid_xg(20)  ! nmax(0)
        integer   :: varid_xs(20,20)  ! mmax , MAX(nmax(1:mmax))

        integer   :: varid_scalar(20)  ! nscalar
        integer   :: varid_rr(20)      ! nRR

        integer   :: varid_kturbg , varid_eturbg


        character(LEN=80) :: fname, var_name

        integer nDim , nVars , nAttr , unID , formatNUM

        integer xyz_id , xyz_dim

        character(LEN=80) :: varname
        integer vartype,nvdims,vdims(10),nvatts,rcode

        double precision, allocatable :: array1(:)
        double precision, allocatable :: array2(:)


!        integer :: MFIX_nf90_create
!        integer :: MFIX_nf90_def_dim
!        integer :: MFIX_nf90_close
!        integer :: MFIX_nf90_open
!       integer :: MFIX_nf90_inquire
!       integer :: MFIX_nf90_inq_dimid
!       integer :: MFIX_nf90_inquire_dimension
!       integer :: MFIX_nf90_inq_varid








!
!


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

  !    return

        if (.not. MFIX_usingNETCDF()) return

      if (myPE .eq. PE_IO) then
         allocate (array1(ijkmax2))
         allocate (array2(ijkmax3))
      else
         allocate (array1(1))
         allocate (array2(1))
      end if



!        call MPI_barrier(MPI_COMM_WORLD,mpierr)
      if (myPE .eq. PE_IO) then
         fname = trim(run_name) // "_RES1.nc"

         call MFIX_check_netcdf( MFIX_nf90_open(fname, NF90_NOWRITE, ncid) )

         call MFIX_check_netcdf( MFIX_nf90_inquire(ncid, nDim , nVars , nAttr , unID , formatNUM) )

!         write (*,*) ' nDim      = ' , nDim
!         write (*,*) ' nVars     = ' , nVars
!         write (*,*) ' nAttr     = ' , nAttr
!         write (*,*) ' unID      = ' , unID
!         write (*,*) ' formatNum = ' , formatNum

         call MFIX_check_netcdf( MFIX_nf90_inq_dimid(ncid,"xyz",xyz_id) )
         call MFIX_check_netcdf( MFIX_nf90_inquire_dimension(ncid,xyz_id,len=xyz_dim) )


         do i = 1,nVars
            call MFIX_ncvinq(ncid,i,varname,vartype,nvdims,vdims,nvatts,rcode)
         end do

         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid ,"time"  , varid_time  ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "EP_g"  , varid_epg   ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "P_g"   , varid_pg    ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "P_star", varid_pstar ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "U_g"   , varid_ug    ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "V_g"   , varid_vg    ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "W_g"   , varid_wg    ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "T_g"   , varid_tg    ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "RO_g"   , varid_rog  ) )
         call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, "ROP_g"  , varid_ropg ) )
         call MFIX_check_netcdf( MFIX_nf90_get_var(ncid,varid_time,time) )

      end if


!      if (myPE .eq. PE_IO) then
!        call check_netcdf( nf90_get_var(ncid , varid_epg   , array1   ) )
!         call convert_from_io_dp(array1,array2,ijkmax2)
!      end if
!      call scatter(EP_g,array2,PE_IO)

      call readScatterRes_netcdf(EP_g   , array2, array1, ncid , varid_epg)
      call readScatterRes_netcdf(P_g    , array2, array1, ncid , varid_pg)
      call readScatterRes_netcdf(P_star , array2, array1, ncid , varid_pstar)
      call readScatterRes_netcdf(ro_g   , array2, array1, ncid , varid_rog)
      call readScatterRes_netcdf(rop_g  , array2, array1, ncid , varid_ropg)
      call readScatterRes_netcdf(u_g    , array2, array1, ncid , varid_ug)
      call readScatterRes_netcdf(v_g    , array2, array1, ncid , varid_vg)
      call readScatterRes_netcdf(w_g    , array2, array1, ncid , varid_wg)
      call readScatterRes_netcdf(t_g    , array2, array1, ncid , varid_tg)

      do i = 1,1   ! mmax

         if (myPe .eq. PE_IO) then
            var_name = 'U_s_xxx'
            write (var_name(5:7),'(i3.3)') I
           call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_us(I)) )
         end if
         call readScatterRes_netcdf(U_s(:,i) , array2, array1, ncid , varid_us(i))

         if (myPe .eq. PE_IO) then
            var_name = 'V_s_xxx'
            write (var_name(5:7),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_vs(I)) )
         end if
         call readScatterRes_netcdf(V_s(:,i) , array2, array1, ncid , varid_vs(i))

         if (myPe .eq. PE_IO) then
            var_name = 'W_s_xxx'
            write (var_name(5:7),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_ws(I)) )
         end if
         call readScatterRes_netcdf(W_s(:,i) , array2, array1, ncid , varid_ws(i))

         if (myPe .eq. PE_IO) then
            var_name = 'ROP_s_xxx'
            write (var_name(7:10),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_rops(I)) )
         end if
         call readScatterRes_netcdf(ROP_s(:,i) , array2, array1, ncid , varid_rops(i))

         if (myPe .eq. PE_IO) then
            var_name = 'T_s_xxx'
            write (var_name(5:7),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_ts(I)) )
         end if
         call readScatterRes_netcdf(T_s(:,i) , array2, array1, ncid , varid_ts(i))

         if (myPe .eq. PE_IO) then
            var_name = 'Theta_m_xxx'
            write (var_name(9:11),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_thetam(I)) )
         end if
         call readScatterRes_netcdf(theta_m(:,i) , array2, array1, ncid , varid_thetam(i))

         if (myPe .eq. PE_IO) then
            var_name = 'gamaRS_xxx'
            write (var_name(8:10),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_gamaRS(I)) )
         end if
         call readScatterRes_netcdf(gama_rs(:,i) , array2, array1, ncid , varid_gamaRS(i))

         if (myPe .eq. PE_IO) then
            var_name = 'TRS_xxx'
            write (var_name(5:7),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_TRS(I)) )
         end if
         call readScatterRes_netcdf(T_rs(:,i) , array2, array1, ncid , varid_TRS(i))

         DO NN = 1, NMAX(i)
            if (myPe .eq. PE_IO) then
               var_name = 'X_s_xxx_xxx'
               write (var_name(5:7) ,'(i3.3)') I
               write (var_name(9:11),'(i3.3)') nn
               call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_xs(I,nn)) )
            end if
            call readScatterRes_netcdf(X_s(:,i,nn) , array2, array1, ncid , varid_xs(i,nn))
         END DO

      end do

      do i = 1,nmax(0)
         if (myPe .eq. PE_IO) then
            var_name = 'X_g_xxx'
            write (var_name(5:7),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_xg(I)) )
         end if
         call readScatterRes_netcdf(X_g(:,i) , array2, array1, ncid , varid_xg(i))
      end do

      do i = 1,nscalar
         if (myPe .eq. PE_IO) then
            var_name = 'Scalar_xxx'
            write (var_name(8:10),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_scalar(I)) )
         end if
         call readScatterRes_netcdf(scalar(:,i) , array2, array1, ncid , varid_scalar(i))
      end do

      do i = 1,nRR
         if (myPe .eq. PE_IO) then
            var_name = 'RRates_xxx'
            write (var_name(8:10),'(i3.3)') I
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, var_name, varid_rr(I)) )
         end if
         call readScatterRes_netcdf(ReactionRates(:,i) , array2, array1, ncid , varid_rr(i))
      end do

      if (k_Epsilon) then
         if (myPe .eq. PE_IO) then
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, 'k_turb_g', varid_kturbg) )
         end if
         call readScatterRes_netcdf(k_turb_g , array2, array1, ncid , varid_kturbg)

         if (myPe .eq. PE_IO) then
            call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, 'e_turb_g', varid_eturbg) )
         end if
         call readScatterRes_netcdf(e_turb_g , array2, array1, ncid , varid_eturbg)
      end if

     if (myPe .eq. PE_IO) then
        call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, 'gamaRG', varid_gamaRG) )
     end if
     call readScatterRes_netcdf(gama_rg , array2, array1, ncid , varid_gamaRG)

     if (myPe .eq. PE_IO) then
        call MFIX_check_netcdf( MFIX_nf90_inq_varid(ncid, 'TRG', varid_TRG) )
     end if
     call readScatterRes_netcdf(T_rg , array2, array1, ncid , varid_TRG)

        ! Close the file. This frees up any internal netCDF resources
        ! associated with the file, and flushes any buffers.
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)
        if (myPE .eq. PE_IO) then
           call MFIX_check_netcdf( MFIX_nf90_close(ncid) )
       end if
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)

  !      stop


        return

      end subroutine read_res1_netcdf



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PATCH_AFTER_RESTART                                    C
!  Purpose: Patch new fluid cells after a restart                      C
!           This could occur when restarting with a different          C
!           grid partition when the RESTART file was generated         C
!           prior to the Dec. 4th 2014 bug fix.                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 14-APR-15  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PATCH_AFTER_RESTART
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE physprop
      USE run
      USE rxns
      USE scalars
      USE funits
      USE energy
      USE compar
      USE cdist
      USE mpi_utility
      USE sendrecv
      USE cutcell
      use functions

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
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I,J,K, IJK, IJKNB
      INTEGER :: M,NN
      INTEGER :: NB
      INTEGER, DIMENSION(6) :: NBCELL
      LOGICAL :: NB_FOUND


!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK).AND.EP_G(IJK)==UNDEFINED) THEN

! Detects new fluid cells that used to be blocked cells with undefined
! values. When a fluid cell has undefined void fraction, this means all
! variables need to be patched. Typically, this will be a fairly small
! cut cell that was flagged as BLOCKED cell with a different partition.
! If a valid fluid cell is found next to this undefined cell, all field
! variables will be copied over. If no valid fluid cell is found, the
! code will continue and will likely stop during the check_data_30
! (zero species mass fractions will yield a zero specific heat).
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            NBCELL(1) = IM_OF(IJK)
            NBCELL(2) = JM_OF(IJK)
            NBCELL(3) = KM_OF(IJK)
            NBCELL(4) = IP_OF(IJK)
            NBCELL(5) = JP_OF(IJK)
            NBCELL(6) = KP_OF(IJK)

            NB_FOUND = .FALSE.

            DO NB = 1,6

               IJKNB = NBCELL(NB)

               IF(FLUID_AT(IJKNB).AND.EP_G(IJKNB)/=UNDEFINED) THEN
                  NB_FOUND = .TRUE.
                  WRITE (*, 1010) MyPE, I,J,K

                  EP_G(IJK)   = EP_G(IJKNB)
                  P_G(IJK)    = P_G(IJKNB)
                  P_STAR(IJK) = P_STAR(IJKNB)
                  RO_G(IJK)   = RO_G(IJKNB)
                  ROP_G(IJK)  = ROP_G(IJKNB)
                  T_G(IJK)    = T_G(IJKNB)

                  T_s(IJK,1:MMAX) = T_s(IJKNB,1:MMAX)

                  U_s(IJK,1:MMAX) = U_s(IJKNB,1:MMAX)
                  V_s(IJK,1:MMAX) = V_s(IJKNB,1:MMAX)
                  W_s(IJK,1:MMAX) = W_s(IJKNB,1:MMAX)

                  X_g(IJK,1:NMAX(0)) = X_g(IJKNB,1:NMAX(0))

                  U_G(IJK) = U_G(IJKNB)
                  V_G(IJK) = V_G(IJKNB)
                  W_G(IJK) = W_G(IJKNB)

                  ROP_S(IJK,1:MMAX) = ROP_S(IJKNB,1:MMAX)

                  IF(ANY(SOLVE_ROs)) RO_S(IJK,1:MMAX) = RO_S(IJKNB,1:MMAX)


                  THETA_M(IJK,1:MMAX) = THETA_M(IJKNB,1:MMAX)

                  DO M = 1,MMAX
                     DO NN = 1, NMAX(M)
                        X_S(IJK,M,NN)= X_S(IJKNB,M,NN)
                     ENDDO
                  ENDDO


                  DO NN = 1, NScalar
                     Scalar(IJK,NN) = Scalar(IJKNB,NN)
                  END DO

                  GAMA_RG(IJK) = GAMA_RG(IJKNB)
                  T_RG(IJK)    = T_RG(IJKNB)

                  GAMA_RS(IJK,1:MMAX) = GAMA_RS(IJKNB,1:MMAX)
                  T_RS(IJK,1:MMAX)    = T_RS(IJKNB,1:MMAX)


                  DO NN = 1, nRR
                     ReactionRates(IJK,NN) = ReactionRates(IJKNB,NN)
                  END DO

                  IF (K_Epsilon) THEN
                     K_Turb_G(IJK) = K_Turb_G(IJKNB)
                     E_Turb_G(IJK) = E_Turb_G(IJKNB)
                  ENDIF

                  EXIT ! Exit as soon as first valid neighbor cell is found
               ENDIF  ! NB is a fluid cell

            ENDDO ! NB Loop

            IF(.NOT.NB_FOUND) WRITE (*, 1020) MyPE, I,J,K   ! NO FLUID CELL AMONG NEIGBHORS

         ENDIF ! New fuid cell

      ENDDO ! IJK loop

1010  FORMAT(1X,'PATCHING NEW FLUID CELL UPON RESTART: MyPE,I,J,K =' ,I6,I6,I6,I6)
1020  FORMAT(1X,'UNABLE TO PATCH NEW FLUID CELL UPON RESTART: MyPE,I,J,K =' ,I6,I6,I6,I6)
      END SUBROUTINE PATCH_AFTER_RESTART
