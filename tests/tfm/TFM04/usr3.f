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


        if(tecplot_output) then
!          CALL generate_grid_locations
          CALL write_tecplot_data
        endif
        
        CALL deallocate_usr_variables


      RETURN

      END SUBROUTINE USR3


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: write_tecplot_data                                     !
!  Purpose: Write data for visualization in Tecplot.                   !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jul 2015   !
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
      use fldvar, only      : p_g, u_g, v_g
      use usr, only         : p_g_ex, u_g_ex, v_g_ex
      use usr, only         : xtr, ytr, ztr
      use usr, only         : tec_output_block, tec_no_k, &
                              error_summary
      use param1, only      : half, zero
      use mpi_utility, only : gather
      use run, only         : discretize
      IMPLICIT NONE

! indices
        integer             :: i, j, k, ijk
        integer             :: imjk, ijmk, ijkm
        integer             :: ipjk, ijpk
        integer             :: ipjmk, imjmk, imjpk

! file units
        integer             :: ftcc, fes

! solution functionals
        double precision    :: TKE, TKEex, PgErrL1

! arrays to gather variables
        double precision, allocatable :: arr_Pg(:)
        double precision, allocatable :: arr_Ug(:)
        double precision, allocatable :: arr_Vg(:)
        double precision, allocatable :: arr_Pgex(:)
        double precision, allocatable :: arr_Ugex(:)
        double precision, allocatable :: arr_Vgex(:)
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
                                          Pg_tmp, Pg_ex_tmp, &
                                          Ug_tmp, Ug_ex_tmp, &
                                          Vg_tmp, Vg_ex_tmp, &
                                          Og_tmp, Og_ex_tmp
                                          ! Og = Vorticity    

! temporary gradients
        double precision              :: DV_DX, DU_DY      

        character(len=12), dimension(0:9) :: DISCR_NAME

        data DISCR_NAME/'FOUP', 'FOUP', 'Superbee', 'Smart', &
          'Ultra-Quick', &
          'QUICKEST', 'Muscl', 'VanLeer', 'Minmod', 'Central'/

! allocate variables as ijkmax3 size 1D arrays
        allocate(arr_Pg(ijkmax3))
        allocate(arr_Ug(ijkmax3))
        allocate(arr_Vg(ijkmax3))
        allocate(arr_Pgex(ijkmax3))
        allocate(arr_Ugex(ijkmax3))
        allocate(arr_Vgex(ijkmax3))
        allocate(arr_xtr(ijkmax3))
        allocate(arr_ytr(ijkmax3))
        allocate(arr_ztr(ijkmax3))

! gather the MFIX global variables into ijkmax3 sized variables
        call gather(P_G,arr_Pg,PE_IO)
        call gather(U_G,arr_Ug,PE_IO)
        call gather(V_G,arr_Vg,PE_IO)
        call gather(P_G_ex,arr_Pgex,PE_IO)
        call gather(U_G_ex,arr_Ugex,PE_IO)
        call gather(V_G_ex,arr_Vgex,PE_IO)
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
            ipjk = funijk_gl(i-1,j,k)
            ijpk = funijk_gl(i,j+1,k)
            ipjmk = funijk_gl(i+1,j-1,k)
            imjmk = funijk_gl(i-1,j-1,k)
            imjpk = funijk_gl(i-1,j+1,k) 


            Pg_tmp(i,j,k) = arr_Pg(ijk)

            Ug_tmp(i,j,k) = half*(arr_Ug(ijk)+arr_Ug(imjk))

            Vg_tmp(i,j,k) = half*(arr_Vg(ijk)+arr_Vg(ijmk))      

            Pg_ex_tmp(i,j,k) = arr_Pgex(ijk)

            Ug_ex_tmp(i,j,k) = half*(arr_Ugex(ijk)+arr_Ugex(imjk))

            Vg_ex_tmp(i,j,k) = half*(arr_Vgex(ijk)+arr_Vgex(ijmk))

! Vorticity variables
! uniform grid spacing assumed
            DV_DX = ( half*(arr_Vg(ipjk)+arr_Vg(ipjmk)) - &
                      half*(arr_Vg(imjk)+arr_Vg(imjmk)) )/&
                      (2.0d0*(arr_xtr(ijk)-arr_xtr(imjk)))

            DU_DY = ( half*(arr_Ug(ijpk)+arr_Ug(imjpk)) - &
                      half*(arr_Ug(ijmk)+arr_Ug(imjmk)) )/&
                      (2.0d0*(arr_ytr(ijk)-arr_ytr(ijmk)))

            Og_tmp(i,j,k) = DV_DX - DU_DY


            DV_DX = ( half*(arr_Vgex(ipjk)+arr_Vgex(ipjmk)) - &
                      half*(arr_Vgex(imjk)+arr_Vgex(imjmk)) )/&
                      (2.0d0*(arr_xtr(ijk)-arr_xtr(imjk)))

            DU_DY = ( half*(arr_Ugex(ijpk)+arr_Ugex(imjpk)) - &
                      half*(arr_Ugex(ijmk)+arr_Ugex(imjmk)) )/&
                      (2.0d0*(arr_ytr(ijk)-arr_ytr(ijmk)))

            Og_ex_tmp(i,j,k) = DV_DX - DU_DY

          end do
          end do
          end do

! write true tecplot cell-centered data (block format)
! This is the traditional tecplot cell-centered data. Be wary while
! using contour plots which could do inaccurate interpolation near the
! boundaries. It's better to use tecplot's "primary flood" with this
! option.

          if(tec_output_block) then

            open(unit=newunit(ftcc), file="solution_tec_block.dat", &
                  position='append', status='unknown')

            write(ftcc,"(11a)") 'variables = "x""y""z"&
                &"Pg""Ug""Vg""Og""Pgex""Ugex""Vgex""Ogex"'
            write(ftcc,*) 'zone T="',DISCR_NAME(DISCRETIZE(1)),'" '
            write(ftcc,*) 'I=',IMAX1,' J=',JMAX1,' K=',KMAX1

            write(ftcc,*) 'DATAPACKING=BLOCK'
            write(ftcc,*) "VARLOCATION=([4-11]=CELLCENTERED)"

            ! x coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(ftcc,*) (x_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! y coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(ftcc,*) (y_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! z coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(ftcc,*) (z_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! solution variables
            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Pg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Ug_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Vg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Og_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Pg_ex_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Ug_ex_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Vg_ex_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Og_ex_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            close(ftcc)

          end if ! end of if(tec_output_block)

! write summary for some solution functionals
! (1) Total KE = Sum of KE at all locations in the domain 
! (2) L1 norm of error in P_G          

          if(error_summary) then

            open(unit=newunit(fes), &
             file="error_summary.dat", status='unknown', & 
             position='append')
!            write(fes,"(4a)") 'variables = "TKE""TKEex""TKEerrPerc"&
!                              "PgErrL1"'
!            write(ftcc,*) 'zone T="',DISCRETIZE(1),'" '
!            write(ftcc,*) 'I=',1

            TKE = zero
            TKEex = zero
            PgErrL1 = zero

            do k=kmin1,kmax1
            do j=jmin1,jmax1
            do i=imin1,imax1
              TKE = TKE + &
                    half*(Ug_tmp(i,j,k)**2+Vg_tmp(i,j,k)**2)
              TKEex = TKEex + &
                      half*(Ug_ex_tmp(i,j,k)**2+Vg_ex_tmp(i,j,k)**2)
              PgErrL1 = PgErrL1 + abs((Pg_tmp(i,j,k)-Pg_ex_tmp(i,j,k))/&
                                  Pg_ex_tmp(i,j,k))*100.0d0
            end do
            end do
            end do

            PgErrL1 = PgErrL1/float(imax*jmax)

            write(fes,*) DISCR_NAME(DISCRETIZE(1)), TKE, TKEex, &
                        abs(TKE-TKEex)/TKEex*100.0d0, PgErrL1

            close(fes)

          endif ! end of if(error_summary)


        end if ! end of if(myPE==PE_IO)

! deallocate local arrays        
        deallocate(arr_Pg)
        deallocate(arr_Ug)
        deallocate(arr_Vg)
        deallocate(arr_Pgex)
        deallocate(arr_Ugex)
        deallocate(arr_Vgex)        
        deallocate(arr_xtr)
        deallocate(arr_ytr)
        deallocate(arr_ztr)


      RETURN

      END SUBROUTINE write_tecplot_data


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: deallocate_usr_variables                               !
!  Purpose: Deallocate allocatable variables defined in usr_mod.f      !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE deallocate_usr_variables
      Use usr, only     : xtr, ytr, ztr
      Use usr, only     : p_g_ex, u_g_ex, v_g_ex
      IMPLICIT NONE

! deallocate variables defined in usr_mod.f
        deallocate(xtr)
        deallocate(ytr)
        deallocate(ztr)

        deallocate(p_g_ex)
        deallocate(u_g_ex)
        deallocate(v_g_ex)

      RETURN

      END SUBROUTINE deallocate_usr_variables
