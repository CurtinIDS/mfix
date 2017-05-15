!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: POST_MFIX                                              C
!  Purpose: main routine for postprocessing MFIX results               C
!                                                                      C
!  Author: P. Nicoletti                               Date: 14-FEB-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE F_INIT
!
      USE cdist
      USE constant
      USE fldvar
      USE funits
      USE geometry
      USE gridmap
      USE indices
      USE machine
      USE parallel_mpi
      USE param
      USE param1
      USE physprop
      USE post3d
      USE read_input
      USE run
!
      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      CHARACTER(LEN=30)  FILE_NAME
      INTEGER L
      INTEGER REC_POINTER(N_SPX)
      LOGICAL READ_SPX(N_SPX)
      LOGICAL AT_EOF(N_SPX)
      REAL    TIME_REAL(N_SPX), TIME_LAST
      LOGICAL OPEN_FILEP
!
      INTEGER M_PASS,N_PASS
      integer :: gas_species_index , solid_species_index , solid_index
      logical :: bRead_all

      logical :: bAllocateAll
      common /ballocate/ bAllocateAll

      common /fast_sp7/ gas_species_index , solid_species_index , &
                         solid_index , bRead_all

      solid_species_index = 0
      solid_index    = 0
      gas_species_index = 0
      bRead_all = .true.
      bAllocateAll = .false.

      bDoing_postmfix = .true.

      nodesi = 1
      nodesj = 1
      nodesk = 1
      call parallel_init()
!// Partition the domain and set indices
!
!
! set up machine constants
!
      CALL MACHINE_CONS
!
! get the RUN_NAME from the user
!
      IF (DO_XFORMS) THEN
         L = INDEX(DIR_NAME,'.RES')
         IF (L.EQ.0) GOTO 10
         RUN_NAME = DIR_NAME(1:L-1)
         GOTO 20
      END IF
!
10    WRITE (*,'(A)',ADVANCE='NO') ' Enter the RUN_NAME to post_process > '
      READ  (*,'(A)') RUN_NAME
      CALL MAKE_UPPER_CASE(RUN_NAME,60)
 20   IF(.NOT. OPEN_FILEP(RUN_NAME,'RESTART_1',N_SPX))GOTO 10
      CALL INIT_NAMELIST
      nodesi = 1
      nodesj = 1
      nodesk = 1


      imin2 = 1
      jmin2 = 1
      kmin2 = 1
      CALL READ_RES0
! this shouldn't be necessary and can create issues:
!      if(no_k .and. dz(1) .eq. UNDEFINED)DZ(1)=ONE
      call SET_MAX2
      call GRIDMAP_INIT
      CALL SET_INCREMENTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     CALL READ_RES1
!
!  Initializations
!
      IF(COORDINATES .EQ. 'CYLINDRICAL') THEN
        CYLINDRICAL = .TRUE.
      ELSE
        CYLINDRICAL = .FALSE.
      ENDIF
!
      IF(IMAX2 .GE. 3) THEN
        DO_I = .TRUE.
        NO_I = .FALSE.
      ELSE
        DO_I = .FALSE.
        NO_I = .TRUE.
      ENDIF
!
      IF(JMAX2 .GE. 3) THEN
        DO_J = .TRUE.
        NO_J = .FALSE.
      ELSE
        DO_J = .FALSE.
        NO_J = .TRUE.
      ENDIF
!
      IF(KMAX2 .GE. 3) THEN
        DO_K = .TRUE.
        NO_K = .FALSE.
      ELSE
        DO_K = .FALSE.
        NO_K = .TRUE.
      ENDIF
!
      CALL SET_MAX2
      CALL SET_GEOMETRY
      CALL SET_CONSTANTS
      CALL SET_INCREMENTS
      TIME_LAST = TIME
!
! READ INITIAL RECORDS OF THE .SPx files ...
!
      DO L = 1,N_SPX
         READ_SPX(L)    = .TRUE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
      CALL READ_SPX0(READ_SPX)
      DO L = 1,N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
!
! CALCULATE DISTANCES AND VOLUMES
!
      CALL CALC_DISTANCE (XMIN,DX,IMAX2,XDIST_SC,XDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DY,JMAX2,YDIST_SC,YDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DZ,KMAX2,ZDIST_SC,ZDIST_VEC)
      write (*,*) ' after call to calc distance'
      CALL CALC_VOL
      write (*,*) ' after call to calc_vol'
!
      IF (DO_XFORMS) RETURN
50    CALL HEADER_MAIN
      IF (SELECTION.EQ.0) THEN
        STOP
      ELSEIF (SELECTION.EQ.1) THEN
        CALL EXAMINE_DATA
      ELSEIF (SELECTION.EQ.2) THEN
        CALL RES_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)
      ELSEIF (SELECTION.EQ.3) THEN
        CALL INTERP_RES
      ELSEIF (SELECTION.EQ.4) THEN
        CALL CALC_QUANTITIES
      ELSEIF (SELECTION.EQ.5) THEN
        CALL PRINT_OUT
      ELSEIF (SELECTION.EQ.6) THEN
        CALL USR_POST
      ELSEIF (SELECTION.EQ.7) THEN
        CALL SELECT_SPX_REC
      ELSEIF (SELECTION.EQ.8) THEN
        CALL TIME_AVG
      ELSEIF (SELECTION.EQ.9) THEN
        CALL ornl_routines
      ELSEIF (SELECTION.EQ.10) THEN
        CALL scavenger
      ENDIF
      GOTO 50
      END
!
      SUBROUTINE CLOSE_OLD_RUN()
!
      Use param
      Use param1
      Use funits
!
      IMPLICIT NONE
!
      INTEGER L
!
      CLOSE (UNIT=UNIT_RES)
      DO L = 1,N_SPX
         CLOSE (UNIT=UNIT_SPX+L)
      END DO
!
      RETURN
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  scavenger modifications !!!!!!!!!!!!!!!!!



! *************************  subroutine scavenger *********************

      subroutine scavenger
        use param, only: dimension_scalar
        Use param1
        Use run
        Use geometry
        Use fldvar
        Use ParallelData
        Use physprop
        Use scalars
        Use rxns
        Use drag
        Use energy, only: hor_g, hor_s, gama_gs, gama_rg, gama_rs, t_rs, t_rg

      implicit none

      integer   :: L , nb , n , i , nArrays , kfile

      ! deallocate some of MFIX variables so we can
      ! allocate in ProcessSpxFile without worrying
      ! about virtual memory

      deAllocate(  F_gs )
      deAllocate(  F_ss )


!energy
      deAllocate(  HOR_g  )
      deAllocate(  HOR_s )
      deAllocate(  GAMA_gs )
      deAllocate(  GAMA_Rg )
      deAllocate(  GAMA_Rs  )
      deAllocate(  T_Rg  )
      deAllocate(  T_Rs )

!fldvar
      deAllocate(  EP_g  )
      deAllocate(  EP_go  )
      deAllocate(  P_g  )
      deAllocate(  P_go  )
      deAllocate(  RO_g  )
      deAllocate(  RO_go)
      deAllocate(  ROP_g  )
      deAllocate(  ROP_go )
      deAllocate(  ROP_s  )
      deAllocate(  ROP_so  )
      deAllocate(  T_g )
      deAllocate(  T_s  )
      deAllocate(  T_go  )
      deAllocate(  T_so  )
      deAllocate(  X_g  )
      deAllocate(  X_s )
      deAllocate(  X_go  )
      deAllocate(  X_so  )
      deAllocate(  U_g )
      deAllocate(  U_go )
      deAllocate(  U_s  )
      deAllocate(  U_so  )
      deAllocate(  V_g )
      deAllocate(  V_go )
      deAllocate(  V_s  )
      deAllocate(  V_so )
      deAllocate(  W_g )
      deAllocate(  W_go )
      deAllocate(  W_s  )
      deAllocate(  W_so  )
      deAllocate(  P_s )
      deAllocate(  P_s_c  )
      deAllocate(  P_s_v )
      deAllocate(  P_s_f )
      deAllocate(  P_s_p )
      deAllocate(  P_star )
      deAllocate(  P_staro )
      deAllocate(  THETA_m  )
      deAllocate(  THETA_mo  )



      IF(K_Epsilon)then
        deAllocate(  K_Turb_G  )
        deAllocate(  K_Turb_Go )
        deAllocate(  E_Turb_G )
        deAllocate(  E_Turb_Go )
      ENDIF

      IF(DIMENSION_Scalar /= 0)then
        deAllocate(  Scalar  )
        deAllocate(  Scalaro  )

      ENDIF

      ext = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      write (*,*) 'enter number of processors'
      read (*,*) np

      write (*,*) ' '
      write (*,*) 'what files should be stitched together ?'
      write (*,*) ' '
      write (*,*) '    -1 : RES and all SPX'
      write (*,*) '     0 : RES only'
      write (*,*) '   k>0 : The kth SPx file (k=1 to 11)'
      write (*,*) ' '
      read  (*,*) kfile
      write (*,*) ' '




      allocate ( cell_map(ijkmax2) )
      allocate ( r_tmp(ijkmax2) )

      allocate ( is3(np) )
      allocate ( ie3(np) )
      allocate ( js3(np) )
      allocate ( je3(np) )
      allocate ( ks3(np) )
      allocate ( ke3(np) )
      allocate ( n_cells(np) )
      allocate ( cr(np) )

! DETERMINE THE FIRST BLANK CHARCATER IN RUN_NAME

      DO L = 1,LEN(RUN_NAME)
         IF (RUN_NAME(L:L).EQ.' ') THEN
            NB = L
            GOTO 100
         END IF
      END DO

100   CONTINUE

      call ReadProcessorInfo       ! istart3/iend3 etc for each processor
      call OpenScavengerFiles(NB,kfile)


      ! process the files


      do L = 1,N_SPX

         nArrays = 0
         if (L .eq.  1) nArrays = 1
         if (L .eq.  2) nArrays = 2
         if (L .eq.  3) nArrays = 3
         if (L .eq.  4) nArrays = 3 * MMAX
         if (L .eq.  5) nArrays = MMAX
         if (L .eq.  6) nArrays = MMAX + 1
         if (L .eq.  7) then
             nArrays = NMAX(0)
             do i = 1,MMAX
                nArrays = nArrays + NMAX(i)
             end do
         end if
         if (L .eq.  8) nArrays = MMAX
         if (L .eq.  9) nArrays = nScalar
         if (L .eq. 10) nArrays = nRR
         if (L .eq. 11 .and. k_Epsilon) nArrays = 2

         if ( kfile.lt.0 .or. (kfile.gt.0 .and. kfile.eq.l) )  &
                     call ProcessSpxFile(L,nb,nArrays,kfile)

      end do



      ! now let's look at the RES files

      call CreateOpen_RES_files(NB,kfile)

      deallocate( cell_map )
      deallocate ( is3 )
      deallocate ( ie3 )
      deallocate ( js3 )
      deallocate ( je3 )
      deallocate ( ks3 )
      deallocate ( ke3 )
      deallocate ( cr )
      deallocate ( n_cells )
      deallocate ( r_tmp )

      stop
      end

! *************************  subroutine ProcessSpxFile *********************

      subroutine ProcessSpxFile(L,nb,nArrays,kfile)

      Use param1
      Use run
      Use geometry
      Use ParallelData

      implicit none

      integer   :: L , n , nb , last_r , ijk , k , nArrays , NA
      integer   :: next_rec , tstep , num_rec , kfile
      real      :: time_r
      real, dimension(:,:) , allocatable :: r_array

      write (fname_scav(NB+8:NB+8),'(a1)') ext(L:L)

      if (nArrays.eq.0) then
         write (*,*) fname_scav(1:nb+8) , ' : no time records added'
         if (kfile .gt. 0) write (*,*) ' '
         return
      end if

      write (*,*)  fname_scav(1:nb+8) , ' : processing'

      allocate( r_array(ijkmax2,nArrays) )

      open (unit=10,file=fname_scav,status='old',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')

      read(unit=10,rec=3) next_rec

      cr = 4


      do while (.true.)

         do n=1,np

           fname_dist = run_name(1:NB-1) // '_xxxxx.SPx'
           write (fname_dist(nb+1:nb+5),'(i5.5)') n-1
           write (fname_dist(nb+9:nb+9),'(a1)') ext(L:L)

           open(unit=20,file=fname_dist,status='old',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')

           read (20,rec=3) last_r

           if (cr(n).ge.last_r .and. n.eq.np) then
              close(unit=20)
              close(unit=10)
              deallocate( r_array )
              if (kfile .gt. 0) write (*,*) ' '
              return
           end if

           tstep = -1

           if (cr(n) .lt. last_r) then
              read (20,rec=cr(n)) time_r , tstep
              cr(n) = cr(n) + 1

              do NA = 1,nArrays
                 call in_bin_512r(20,r_tmp,n_cells(n),cr(n))


                 do k = 1,cellcount(n)
                    ijk = cell_map_v2(n,k)   ! ijk is global IJK
                    r_array(ijk,NA) = r_tmp( cell_map(ijk)%ijk)
                 end do


!                 do ijk = 1,ijkmax2
!                    if ( cell_map(ijk)%proc .eq. n ) then
!                       k = cell_map(ijk)%ijk
!                           r_array( ijk , NA ) = r_tmp ( k )
!                    end if
!                 end do
              end do
           end if

           close (unit=20)

         end do

         if (tstep .ne. -1) then
            write (10,rec=next_rec) time_r , tstep
            num_rec = next_rec
            next_rec = next_rec + 1
            do NA=1,nArrays
              call out_bin_512r(10,r_array(1,NA),ijkmax2,next_rec)
            end do
            num_rec = next_rec - num_rec
            write (unit=10,rec=3) next_rec,num_rec
         end if

      end do

      deallocate(r_array)

      if (kfile .gt. 0) write (*,*) ' '

      close (unit=10)

      return
      end

! *************************  subroutine ReadProcessInfo *********************

      subroutine ReadProcessorInfo

      Use param
      Use param1
      Use run
      Use geometry
      Use ParallelData

      implicit none

      integer   :: L , idummy , ijk_proc , ijk_io , i , cellcount_max

      integer , allocatable :: cell_index(:)

      character(len=80) :: fname

      ! read in the cell info for each processor (using p_info_xxxxx.txt'
      ! (xxxx = processor number)


      allocate ( cellcount(np) )
      allocate ( cell_index(np) )
      do l = 1,np
         cellcount(l)  = 0;
         cell_index(l) = 0
      end do

      do L = 1,np
         fname = 'p_info_xxxxx.txt'
         write (fname(8:12),'(i5.5)') L-1

         open (unit=10,file=fname,status='old',err=100,convert='big_endian')
         goto 101
 100     continue
         fname = 'p_info_xxxx.txt'
         write (fname(8:11),'(i4.4)') L-1

         open (unit=10,file=fname,status='old',convert='big_endian')

 101     continue
         read (10,*)
         read (10,*) idummy , is3(L) , ie3(L)
         read (10,*) idummy , js3(L) , je3(L)
         read (10,*) idummy , ks3(L) , ke3(L)
         read (10,*)
         read (10,*)
         read (10,*)

         n_cells(L) = (ie3(L)-is3(L)+1) * (je3(L)-js3(L)+1) * &
                                        (ke3(L)-ks3(L)+1)

 !        if (n_cells(L) .gt. ncells_max) ncells_max = n_cells(L)

         cr(L) = 4

         do i = 1,n_cells(L)
            read (10,*) idummy,idummy,idummy,idummy,idummy, &
                      ijk_proc , ijk_io

            if (ijk_io .gt. 0  .and. ijk_io .le. ijkmax2) then
               cell_map(ijk_io)%proc = L
               cell_map(ijk_io)%ijk  = ijk_proc
            end if
         end do

         close (unit=10)
      end do

      do L = 1,ijkmax2
         i = cell_map(L)%proc
         cellcount(i) = cellcount(i) + 1
      end do

      cellcount_max = 0
      do l = 1,np
         if (cellcount(l) .gt. cellcount_max) cellcount_max = cellcount(l)
      end do

      allocate(cell_map_v2(np,cellcount_max))

      do L = 1,ijkmax2
         i = cell_map(L)%proc
         cell_index(i) = cell_index(i) + 1
         cell_map_v2(i,cell_index(i)) = L  !cell_map(L)%ijk
      end do

      return
      end

! *************************  subroutine OpenScavengerFiles *********************

      subroutine OpenScavengerFiles(nb,kfile)

      Use param
      Use param1
      Use run
      Use ParallelData

      implicit none

      integer   :: L , NB , kfile

      if (kfile .eq. 0) return  ! only doing RES file

      fname_scav = run_name(1:NB-1) // '_SCAV.SPx'

      ! create the scavenger files ... write out the first 3 records
      do L = 1,N_SPX

         if (kfile.gt.0 .and. l.ne.kfile) goto 100

         write (fname_scav(NB+8:NB+8),'(a1)') ext(L:L)

         open (unit=10,file=fname_scav,status='unknown',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')

         fname_dist = run_name(1:NB-1) // '_00000.SPx'
         write (fname_dist(nb+9:nb+9),'(a1)') ext(L:L)

         open(unit=20,file=fname_dist,status='old',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')

         read(20,rec=1)  pbuffer
         write(10,rec=1) pbuffer
         read(20,rec=2)  pbuffer
         write(10,rec=2) pbuffer

         write(10,rec=3) 4,-1
         close (unit=20)
         close (unit=10)

 100     continue

      end do

      return
      end




! *************************  subroutine CreateOpen_RES_files *********************

      subroutine CreateOpen_RES_files(nb,kfile)

      Use param1
      Use run
      Use geometry
      Use fldvar
      Use ParallelData
      Use physprop
      Use scalars
      Use rxns

      implicit none

      integer :: L , NB , next_rec , ntot
      integer :: max_dim , n , ijk , k , kfile

      double precision, dimension(:) , allocatable :: d_tmp
      double precision, dimension(:) , allocatable :: array

      character(len=80) :: fname

      if (kfile .gt. 0) return

      call close_old_run()

      fname_scav = run_name(1:NB-1) // '_SCAV.RES'
      fname_dist = run_name(1:NB-1) // '.RES'        ! the "00000" file (but currently
                                                      ! still has original name

      write (*,*) fname_scav(1:nb+8) , ' : processing'


      open (unit=10,file=fname_scav(1:nb+8),status='unknown',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')

      open(unit=20,file=fname_dist,status='old',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')

!      open(unit=11,file='deb.txt',status='unknown')

      ! the following writes all of the RES0 header data
      read(20,rec=3) next_rec

      do L = 1,next_rec

         read(20,rec=L)  pbuffer
         write(10,rec=L) pbuffer

      end do
      write (10,rec=3) next_rec
      next_rec = next_rec + 1
      close(unit=20)

      ! calculate the maximum dimesnion needed for all processors
      max_dim = 0
      do L = 1,np

         if (n_cells(L) .gt. max_dim) max_dim = n_cells(L)

         call set_res_name(run_name,nb,L,fname)

         open(unit=21,file=fname,status='old',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')

         read(21,rec=3) cr(L)

         cr(L) = cr(L) + 1   ! TIME, DT, NSTEP
!         write (*,*) ' l,cr(l) = ' , L,CR(L)

         close(unit=21)

      end do

      allocate( d_tmp(max_dim) )
      allocate( array(ijkmax2) )


      ! calculate the number of arrays in the RES1 section

      ntot = 6                ! EP_g , P_g , P_star , RO_g , ROP_g , T_g
      ntot = ntot + NMAX(0)   ! X_g
      ntot = ntot + 3         ! U_g , V_g , W_g

      ntot = ntot + mmax*6    ! ROP_s , T_s , U_s , V_s , W_s , THETA_m

      do L=1,MMAX
         ntot = ntot + NMAX(L) ! X_s
      end do

      ntot = ntot + nScalar
      ntot = ntot + nRR
      if (k_Epsilon) ntot = ntot + 2

      ntot = ntot + 2 ! gama_rg , T_rg

      ntot = ntot + 2*mmax  ! Gama_rs , T_rs

      do L = 1,ntot
         do n = 1,np

            call set_res_name(run_name,nb,N,fname)

            open(unit=20,file=fname,status='old',recl=512, &
                         access='direct',form='unformatted',convert='big_endian')
            call in_bin_512(20,d_tmp,n_cells(N),cr(N))
            close (unit=20)

!
! newer, faster method
!
            do k = 1,cellcount(n)                            !change ... was cellcount(NP)
               ijk = cell_map_v2(n,k)   ! ijk is global IJK
               array(ijk) = d_tmp( cell_map(ijk)%ijk)
            end do

!
! previous method
!            do ijk = 1,ijkmax2
!
!                if ( cell_map(ijk)%proc .eq. n ) then
!                   k = cell_map(ijk)%ijk
!                   array( ijk ) = d_tmp ( k )
!               end if
!            end do

         end do

         call out_bin_512(10,array,ijkmax2,next_rec)

      end do

      close (unit=10)
!      close(unit=20)

      deallocate (array)
      deallocate (d_tmp)

      write (*,*) ' '

      return
      end

      subroutine set_res_name(run_name,nb,L,fname)
      implicit none

      integer :: nb,L
      character(len=*) :: run_name , fname

      if (l .eq. 1) then
         fname = run_name(1:nb-1) // '.RES'
      else
         fname = run_name(1:nb-1) // '_xxxxx.RES'
         write (fname(nb+1:nb+5),'(i5.5)') L-1
      end if

      return
      end
