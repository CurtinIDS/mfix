!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TIME_AVG                                               C
!  Purpose: Select records from SPx files interactively and write to a C
!           new SPx file                                               C
!           NOTE : USES N_SPX                                          C
!                                                                      C
!  Author: P. Nicoeltti                               Date:13-JUN-2002 C
!  Reviewer:                                                           C
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
      SUBROUTINE TIME_AVG
!
!
      Use param
      Use param1
      Use geometry
      Use indices
      Use run
      Use machine
      Use funits
      Use post3d
      Use physprop
      Use fldvar
      Use scalars
      Use rxns
!
      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      REAL    :: TIME_FOR_RES, TIME_FOUND
      LOGICAL :: AT_EOF(N_SPX), READ_SPX(N_SPX),SELECT
      INTEGER :: REC_POINTER(N_SPX), REC_POINTER_t(N_SPX)
      INTEGER :: NSTEP_1 , ERROR_CODE , solmax
      CHARACTER(LEN=1)   :: IANS
      CHARACTER(LEN=13)  :: line
      REAL    :: TIME_REAL(N_SPX)
      LOGICAL :: ERROR

!
! variables for time averaging
!
      double precision,allocatable :: tavg(:,:,:)
      double precision,allocatable :: tavgs(:,:,:)
      integer :: tcount
      real    :: tstart , tend
!
      INTEGER L, L_SPX , LL , M , i , NB
      integer :: unit_add = 10
!
      CHARACTER(LEN=35) :: EXT_END
!-----------------------------------------------

      ext_end = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
      ERROR  = .FALSE.
      SELECT = .TRUE.
      L_SPX  = SPX_NUM

      write (*,*) ' enter starting time for time averaging'
      read  (*,*) tstart
      write (*,*) ' enter ending time for time averaging'
      read  (*,*) tend

      WRITE (*,'(A)',ADVANCE='NO') 'Enter the RUN_NAME for time averaged data > '
      READ  (*,'(A60)') TEMP_FILE

      do i = 1,len(temp_file)
          if (temp_file(i:i) .eq. ' ') then
             nb = i
             goto 10
          end if
      end do

 10   call make_upper_case(temp_file,60)
      temp_file(nb:nb+3) = '.SPx'
      unit_add = n_spx + 1

      do l_spx = 1,n_spx  ! for now

         IF (.NOT.DO_XFORMS) THEN
            write (*,*) ' '
            write (*,*) ' process old SPx file : ' , l_spx
            temp_file(nb+3:nb+3) = ext_end(l_spx:l_spx)
            OPEN (UNIT=UNIT_SPX+L_SPX+unit_add,FILE=TEMP_FILE, &
                 STATUS='NEW', &
                 RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED', &
                 ERR=101,CONVERT='BIG_ENDIAN')
         END IF
!
! allocate variables as needed  (what if no scalars or nRR ??)
!
         if (l_spx .eq. 1) then          ! ep_g
            allocate (tavg(ijkmax2,1,1))
         else if (l_spx .eq. 2) then     ! p_g and P_star
            allocate (tavg(ijkmax2,2,1))
         else if (l_spx .eq. 3) then     ! U_g , V_g , W_g
            allocate (tavg(ijkmax2,3,1))
         else if (l_spx .eq. 4) then     ! U_s , V_s , W_s
            allocate (tavg(ijkmax2,3,MMAX))
         else if (l_spx .eq. 5) then     ! ROP_s
            allocate (tavg(ijkmax2,MMAX,1))
         else if (l_spx .eq. 6) then     ! T_g , T_s
            allocate (tavg(ijkmax2,1,1))
            allocate (tavgs(ijkmax2,MMAX,1))
         else if (l_spx .eq. 7) then     ! X_g , X_s
            allocate (tavg(ijkmax2,nmax(0),1))
            solmax = 0
            do i = 1,mmax
               solmax = max(solmax,nmax(i))
            end do
            allocate (tavgs(ijkmax2,mmax,solmax))
         else if (l_spx .eq. 8) then     ! THETA_m
            allocate (tavg(ijkmax2,MMAX,1))
         else if (l_spx .eq. 9 .and. nscalar.gt.0) then      ! Scalars
            allocate (tavg(ijkmax2,nscalar,1))
         else if (l_spx .eq. 10 .and. nRR.gt.0) then     ! ReactionRate
            allocate (tavg(ijkmax2,nRR,1))
         else if (l_spx .eq. 9 .and. nscalar.eq.0) then      ! Scalars
            allocate (tavg(ijkmax2,1,1))
         else if (l_spx .eq. 10 .and. nRR.eq.0) then     ! ReactionRate
            allocate (tavg(ijkmax2,1,1))
         else if (l_spx .eq. 11 .and. k_Epsilon) then        ! turbulence
            allocate (tavg(ijkmax2,2,1))
         else
            allocate (tavg(ijkmax2,1,1))
         end if

         tavg(:,:,:) = 0.0
         if ( (l_spx.eq.6) .or. (l_spx.eq.7) ) tavgs(:,:,:) = 0.0
!
!
         CALL WRITE_SPX0(L_SPX,unit_add)
!
         DO L = 1, N_SPX
            READ_SPX(L) = .FALSE.
            REC_POINTER(L) = 4
            AT_EOF(L) = .FALSE.
         END DO

         READ_SPX(L_SPX) = .TRUE.
         tcount = 0
!
         L = 0
100      continue
         CALL READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, TIME_REAL,NSTEP_1)
         IF (.NOT.AT_EOF(L_SPX)) THEN
            L = L + 1
            IF (.NOT.DO_XFORMS) THEN
               if (time_real(l_spx) .gt. tend) goto 99
               if (time_real(l_spx).ge.tstart) then
                  tcount = tcount + 1
                  TIME = DBLE(TIME_REAL(L_SPX))

               if (l_spx .eq. 1) then
                  tavg(:,1,1) = tavg(:,1,1) + ep_g(:)
               else if (l_spx .eq. 2) then
                  tavg(:,1,1) = tavg(:,1,1) + P_g(:)
                  tavg(:,2,1) = tavg(:,2,1) + P_star(:)
               else if (l_spx .eq. 3) then
                  tavg(:,1,1) = tavg(:,1,1) + U_g(:)
                  tavg(:,2,1) = tavg(:,2,1) + V_g(:)
                  tavg(:,3,1) = tavg(:,3,1) + W_g(:)
               else if (l_spx .eq. 4) then
                  do m = 1,MMAX
                     tavg(:,1,m) = tavg(:,1,m) + U_s(:,m)
                     tavg(:,2,m) = tavg(:,2,m) + V_s(:,m)
                     tavg(:,3,m) = tavg(:,3,m) + W_s(:,m)
                  end do
               else if (l_spx .eq. 5) then
                  do m = 1,MMAX
                     tavg(:,m,1) = tavg(:,m,1) + ROP_s(:,m)
                  end do
               else if (l_spx .eq. 6) then
                  tavg(:,1,1) = tavg(:,1,1) + T_g(:)
                  do m = 1,MMAX
                     tavgs(:,m,1) = tavgs(:,m,1) + T_s(:,m)
                  end do
               else if (l_spx .eq. 7) then
                  do i = 1,nmax(0)
                     tavg(:,i,1) = tavg(:,i,1) + X_g(:,i)
                  end do
                  do m = 1,MMAX
                     do i = 1,nmax(m)
                        tavgs(:,m,i) = tavgs(:,m,i) + X_s(:,m,i)
                     end do
                  end do
               else if (l_spx .eq. 8) then
                  do m = 1,MMAX
                     tavg(:,m,1) = tavg(:,m,1) + THETA_m(:,m)
                  end do
               else if (l_spx .eq. 9) then
                  do m = 1,nscalar
                     tavg(:,m,1) = tavg(:,m,1) + Scalar(:,m)
                  end do
               else if (l_spx .eq. 10) then
                  do m = 1,nRR
                     tavg(:,m,1) = tavg(:,m,1) + ReactionRates(:,m)
                  end do
               else if (l_spx .eq. 11) then
                  if (k_Epsilon) then
                     tavg(:,1,1) = tavg(:,1,1) + k_turb_g(:)
                     tavg(:,2,1) = tavg(:,1,1) + e_turb_g(:)
                  end if
               else
                  tavg = 0
               end if

              !  write (*,*) 'processed data ... time = ' , time_real(l_spx)

            END IF
         ENDIF
         GOTO 100
      ENDIF

 99   continue

      time = 0.0

      if (tcount .eq. 0) tcount = 1

      if (l_spx .eq. 1) then
          EP_g(:) = tavg(:,1,1) / real(tcount)
      else if (l_spx .eq. 2) then
          P_g(:) = tavg(:,1,1) / real(tcount)
          p_star(:) = tavg(:,2,1) / real(tcount)
      else if (l_spx .eq. 3) then
          U_g(:) = tavg(:,1,1) / real(tcount)
          V_g(:) = tavg(:,2,1) / real(tcount)
          W_g(:) = tavg(:,3,1) / real(tcount)
      else if (l_spx .eq. 4) then
          do m = 1,mmax
             U_s(:,m) = tavg(:,1,m) / real(tcount)
             V_s(:,m) = tavg(:,2,m) / real(tcount)
             W_s(:,m) = tavg(:,3,m) / real(tcount)
          end do
      else if (l_spx .eq. 5) then
          do m = 1,mmax
             ROP_s(:,m) = tavg(:,m,1) / real(tcount)
          end do
      else if (l_spx .eq. 6) then
          T_g(:) = tavg(:,1,1) / real(tcount)
          do m = 1,mmax
             T_s(:,m) = tavgs(:,m,1) / real(tcount)
          end do
      else if (l_spx .eq. 7) then
          do i = 1,nmax(0)
             X_g(:,i) = tavg(:,i,1) / real(tcount)
          end do
          do m = 1,mmax
             do i = 1,nmax(m)
                X_s(:,m,i) = tavgs(:,m,i) / real(tcount)
             end do
          end do
      else if (l_spx .eq. 8) then
          do m = 1,mmax
             THETA_m(:,m) = tavg(:,m,1) / real(tcount)
          end do
      else if (l_spx .eq. 9) then
          do m = 1,nscalar
             Scalar(:,m) = tavg(:,m,1) / real(tcount)
          end do
      else if (l_spx .eq. 10) then
          do m = 1,nRR
             ReactionRates(:,m) = tavg(:,m,1) / real(tcount)
          end do
      else if (l_spx .eq. 11) then
          if (k_Epsilon) then
             K_TURB_G(:) = tavg(:,1,1) / real(tcount)
             E_TURB_G(:) = tavg(:,2,1) / real(tcount)
          endif
      end if

      CALL WRITE_SPX1(L_SPX,unit_add)

      deallocate(tavg)
      if ( (l_spx.eq.6) .or. (l_spx.eq.7) ) deallocate(tavgs)

      CLOSE(UNIT_SPX+L_SPX+unit_add)

      WRITE (*,*) ' number of records processed = ' , tcount
      WRITE(*,*) ' time averaged file written >' , temp_file(1:nb+3)

 101  continue

      end do ! loop of l_spx

      write (*,*) ' '
      write (*,*) ' '
      RETURN
      END
