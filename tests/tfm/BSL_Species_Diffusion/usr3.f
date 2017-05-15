!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is     C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.              C
!           This routine is not called from an IJK loop, hence         C
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE USR3
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      use usr, only         : tecplot_output
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!


        CALL calculate_de_norms

        if(tecplot_output) CALL write_tecplot_data

        CALL deallocate_usr_variables


      RETURN

      END SUBROUTINE USR3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_de_norms	                               !
!  Purpose: Calculates discretization error norms for the problem      !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: May 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_de_norms
      Use param1, only    : zero
      Use compar, only    : ijkstart3, ijkend3, myPE
      Use indices, only   : i_of, j_of, k_of
      Use functions, only : IS_ON_myPE_owns
      Use usr, only       : de_x_g_1, de_x_g_2
      Use usr, only       : x_g_1_ex, x_g_2_ex
      Use usr, only       : lnorms_x_g_1, lnorms_x_g_2
      Use fldvar, only    : X_g
      Use geometry, only  : imax, imin1, imax1
      Use geometry, only  : jmax, jmin1, jmax1
      Use funits, only    : newunit
      Use compar, only    : myPE, PE_IO
      IMPLICIT NONE

! number of data points for norm calculation
        integer   :: var_size

! looping variables
        integer   :: ijk, i, j, k

! file unit
        integer   :: f1


        call calculate_exact_solution

! scalar variable (temperature)
        de_x_g_1 = zero
        de_x_g_2 = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle

          de_x_g_1(ijk) = X_g(ijk,1) - x_g_1_ex(ijk)
          de_x_g_2(ijk) = X_g(ijk,2) - x_g_2_ex(ijk)
        end do

        var_size = imax*jmax
        call calculate_lnorms(var_size, de_x_g_1, lnorms_x_g_1)
        call calculate_lnorms(var_size, de_x_g_2, lnorms_x_g_2)

! Output DE norms data to a file
        if(myPE == PE_IO) then
          open(unit=newunit(f1), file="de_norms.dat", status='unknown')
          write(f1,*) "# DE Norms for 1D species diffusion:"
          write(f1,*) "# imax= ",imax, " jmax=", jmax
          write(f1,*) "# 1st line: L1 Norms, 2nd line: L2 Norms, &
                       &3rd line: Linf Norms"
          write(f1,*) "# Columns: X_g_1 X_g_2"
          write(f1,*) lnorms_x_g_1(1), lnorms_x_g_2(1)
          write(f1,*) lnorms_x_g_1(2), lnorms_x_g_2(2)
          write(f1,*) lnorms_x_g_1(3), lnorms_x_g_2(3)
          close(f1)
        end if


      RETURN

      END SUBROUTINE calculate_de_norms


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_lnorms	                               !
!  Purpose: Calculates L1, L2 and Lifinity norms for any input         !
!  variable of defined size                                            !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: May 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_lnorms(var_size, var, lnorms)
      Use compar, only      : ijkstart3, ijkend3
      Use param, only       : dimension_3
      Use param1, only      : zero
      Use mpi_utility, only : global_sum, global_max
      IMPLICIT NONE

! total number of data points for norm calculation
        integer, intent(in)                             :: var_size

! variable for which norms are to be calculated
        double precision, dimension(dimension_3), intent(in)   :: var

! the three norms: L1, L2 and Linfinity
        double precision, dimension(3), intent(out)     :: lnorms

! sum or max of lnorms (over all processors
        double precision, dimension(3)                  :: lnorms_all

! looping indices
        integer         :: ijk


! calculate L1, L2 and Linfinity norms
        lnorms(1:3) = zero
        do ijk = ijkstart3, ijkend3
            lnorms(1) = lnorms(1) + abs(var(ijk))
            lnorms(2) = lnorms(2) + abs(var(ijk))**2
            lnorms(3) = max(lnorms(3), abs(var(ijk)))
        end do

! save global sum in lnorms_all variables
        call global_sum(lnorms(1), lnorms_all(1))
        call global_sum(lnorms(2), lnorms_all(2))
        call global_max(lnorms(3), lnorms_all(3))

! put final result in lnorms
        lnorms(1) = lnorms_all(1)/dble(var_size)
        lnorms(2) = sqrt(lnorms_all(2)/dble(var_size))
        lnorms(3) = lnorms_all(3)


      RETURN

      END SUBROUTINE calculate_lnorms


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_exact_solution                               !
!  Purpose: Calculates exact solution for the 1D heat conduction       !
!  through a plane slab and specificed temperatures at boundaries.     !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Aug 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_exact_solution
      Use compar, only    : ijkstart3, ijkend3
      Use indices, only   : i_of, j_of, k_of
      Use geometry, only  : dx, dy, dz, imax1
      Use usr, only       : xtr, ytr, ztr
      Use usr, only       : x_g_1_ex, x_g_2_ex
      Use param1, only    : zero, half
      IMPLICIT NONE

! user-defined 1D steady heat conduction parameters
! (must match the input in mfix.dat)
! (better to set these up via an input file in future)
        double precision, parameter :: xh2O0  = 0.9d0   ! X_g at x=0
        double precision, parameter :: xh2OL  = 0.1d0   ! X_g at x=L
        double precision, parameter :: L      = 0.5d0     ! Length

! looping indices
        integer :: ijk, i, j, k

        integer :: ii, jj, kk

! temporary variables
        double precision  :: xt, yt, zt


! generate grid locations for exact solution calculation
        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          xt = zero - dx(1)
          yt = zero - dy(1)
          zt = zero - dz(1)

          do ii = 1, i
            xt = xt + dx(ii)
          end do

          do jj = 1, j
            yt = yt + dy(jj)
          end do

          do kk = 1, k
            zt = zt + dz(kk)
          end do

          xtr(ijk) = xt
          ytr(ijk) = yt
          ztr(ijk) = zt
        end do

! set exact solution
        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)

! scalar variable (species mass fraction)
          xt = xtr(ijk) - dx(i)*half

          x_g_1_ex(ijk) = xH2OL*(xH2O0/xH2OL)**(xt/L)
          x_g_2_ex(ijk) = 1.0d0 - xH2OL*(xH2O0/xH2OL)**(xt/L)

        end do


      RETURN

      END SUBROUTINE calculate_exact_solution


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: deallocate_usr_variables                               !
!  Purpose: Deallocate allocatable variables defined in usr_mod.f      !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Aug 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE deallocate_usr_variables
      Use usr, only     : x_g_1_ex, x_g_2_ex
      Use usr, only     : lnorms_x_g_1, lnorms_x_g_2
      Use usr, only     : xtr, ytr, ztr
      Use usr, only     : de_x_g_1, de_x_g_2
      IMPLICIT NONE

        ! allocate variables defined in usr_mod.f
        deallocate(x_g_1_ex)
        deallocate(x_g_2_ex)
        deallocate(lnorms_x_g_1)
        deallocate(lnorms_x_g_2)
        deallocate(de_x_g_1)
        deallocate(de_x_g_2)

        deallocate(xtr)
        deallocate(ytr)
        deallocate(ztr)


      RETURN

      END SUBROUTINE deallocate_usr_variables


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: write_tecplot_data                                     !
!  Purpose: Write data for visualization in Tecplot.                   !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Aug 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE write_tecplot_data
      use geometry, only    : ijkmax3
      use geometry, only    : imin2, imin1, imax1, imax2
      use geometry, only    : jmin2, jmin1, jmax1, jmax2
      use geometry, only    : kmin2, kmin1, kmax1, kmax2
      use geometry, only    : imax, jmax, kmax
      use geometry, only    : dx, dy, dz
      use compar, only      : myPE, PE_IO
      use funits, only      : newunit
      use functions, only   : funijk_gl
      use fldvar, only      : X_g
      use usr, only         : x_g_1_ex, x_g_2_ex
      use usr, only         : xtr, ytr, ztr
      use usr, only         : tec_output_block, tec_no_k, &
                              x_g_profile
      use param1, only      : half
      use mpi_utility, only : gather
      IMPLICIT NONE

! indices
        integer             :: i, j, k, ijk
        integer             :: imjk, ijmk, ijkm

! file units
        integer             :: fxcc, fxgp

! arrays to gather variables
        double precision, allocatable :: arr_Xg1(:)
        double precision, allocatable :: arr_Xg2(:)
        double precision, allocatable :: arr_Xg1ex(:)
        double precision, allocatable :: arr_Xg2ex(:)
        double precision, allocatable :: arr_xtr(:)
        double precision, allocatable :: arr_ytr(:)
        double precision, allocatable :: arr_ztr(:)

! temporary variables for x,y,z coordinates
        double precision              :: xt, yt, zt

! x,y,z coordinate variables (3D arrays)
        double precision, &
          dimension(imin2:imax1,jmin2:jmax1,kmin2:kmax1)  :: &
                                          x_tmp, y_tmp, z_tmp

! output variables (3D arrays)
        double precision, &
          dimension(imin2:imax2,jmin2:jmax2,kmin2:kmax2)  :: &
                                          Xg1_tmp, Xg1ex_tmp, &
                                          Xg2_tmp, Xg2ex_tmp


! allocate variables as ijkmax3 size 1D arrays
        allocate(arr_Xg1(ijkmax3))
        allocate(arr_Xg1ex(ijkmax3))
        allocate(arr_Xg2(ijkmax3))
        allocate(arr_Xg2ex(ijkmax3))
        allocate(arr_xtr(ijkmax3))
        allocate(arr_ytr(ijkmax3))
        allocate(arr_ztr(ijkmax3))

! gather the MFIX global variables into ijkmax3 sized variables
        call gather(X_g(:,1),arr_Xg1,PE_IO)
        call gather(X_g(:,2),arr_Xg2,PE_IO)
        call gather(X_g_1_ex,arr_Xg1ex,PE_IO)
        call gather(X_g_2_ex,arr_Xg2ex,PE_IO)
        call gather(xtr,arr_xtr,PE_IO)
        call gather(ytr,arr_ytr,PE_IO)
        call gather(ztr,arr_ztr,PE_IO)

! write out data from PE_IO processor
        if(myPE==PE_IO) then

! create 3D arrays for cell-centered data

! create 3D arrays for mesh nodes
          do k = kmin2, kmax1
          do j = jmin2, jmax1
          do i = imin2, imax1
             ijk = funijk_gl(i,j,k)
             x_tmp(I,J,K) = arr_xtr(IJK)
             y_tmp(I,J,K) = arr_ytr(IJK)
             z_tmp(I,J,K) = arr_ztr(IJK)
          end do
          end do
          end do

! create 3D arrays for solution variables at cell-centers
          do k = kmin1, kmax1
          do j = jmin1, jmax1
          do i = imin1, imax1

            ijk = funijk_gl(i,j,k)
            imjk = funijk_gl(i-1,j,k)
            ijmk = funijk_gl(i,j-1,k)
            if(tec_no_k) then
              ijkm = funijk_gl(i,j,k)
            else
              ijkm = funijk_gl(i,j,k-1)
            end if

            Xg1_tmp(i,j,k) = arr_Xg1(ijk)
            Xg2_tmp(i,j,k) = arr_Xg2(ijk)
            Xg1ex_tmp(i,j,k) = arr_Xg1ex(ijk)
            Xg2ex_tmp(i,j,k) = arr_Xg2ex(ijk)

          end do
          end do
          end do

! write true tecplot cell-centered data (block format)
! This is the traditional tecplot cell-centered data. Be wary while
! using contour plots which could do inaccurate interpolation near the
! boundaries. It's better to use tecplot's "primary flood" with this
! option.

          if(tec_output_block) then

            open(unit=newunit(fxcc), file="solution_tec_block.dat", &
                                      status='unknown')

            write(fxcc,"(9a)") 'variables = "x""y""z"&
                &"Xg1""Xg2""Xg1ex""Xg2ex"'
            write(fxcc,*) 'zone T="',0,'" '
            write(fxcc,*) 'I=',IMAX1,' J=',JMAX1,' K=',KMAX1

            write(fxcc,*) 'DATAPACKING=BLOCK'
            write(fxcc,*) "VARLOCATION=([4-7]=CELLCENTERED)"

            ! x coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(fxcc,*) (x_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! y coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(fxcc,*) (y_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! z coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(fxcc,*) (z_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! solution variables
            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(fxcc,*) (Xg1_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            ! solution variables
            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(fxcc,*) (Xg2_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(fxcc,*) (Xg1ex_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(fxcc,*) (Xg2ex_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            close(fxcc)

          end if ! end of if(tec_output_block)

! write tecplot data (point format) for the species mass fraction profile at 
! horizontal centerline:
! - numerical solution for x_g_1 and x_g_2          
! - exact solution for x_g_1 and x_g_2
! - error for x_g_1 and x_g_2

          if(x_g_profile) then

            open(unit=newunit(fxgp), &
             file="solution_xg_profile.dat", status='unknown')

            write(fxgp,"(6a)") 'variables = "x""y""z""Xg1""Xg2"&
                                "Xg1ex""Xg2ex"&
                                &"Xg1Err""Xg2Err"'
            write(fxgp,*) 'zone T="',0,'" '
            write(fxgp,*) 'I= ',IMAX

            j = jmax/2 + 1
            k = kmin1
            do i = imin1, imax1

              ijk = funijk_gl(i,j,k)

              xt = arr_xtr(ijk) - dx(i)*half
              yt = arr_ytr(ijk) - dy(j)*half
              zt = arr_ztr(ijk) - dz(k)*half

              Xg1_tmp(I,J,K) = arr_Xg1(IJK)
              Xg2_tmp(I,J,K) = arr_Xg2(IJK)
              Xg1ex_tmp(I,J,K) = arr_Xg1ex(IJK)
              Xg2ex_tmp(I,J,K) = arr_Xg2ex(IJK)

              write(fxgp,*) xt, yt, zt, &
                Xg1_tmp(I,J,K), &
                Xg2_tmp(I,J,K), &
                Xg1ex_tmp(I,J,K), &
                Xg2ex_tmp(I,J,K), &
                Xg1_tmp(I,J,K) - Xg1ex_tmp(I,J,K), &
                Xg2_tmp(I,J,K) - Xg2ex_tmp(I,J,K)

            end do

            close(fxgp)

          endif ! end of if(x_g_profile)


        end if ! end of if(myPE==PE_IO)



      RETURN

      END SUBROUTINE write_tecplot_data
