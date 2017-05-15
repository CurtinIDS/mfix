!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR3                                                   !
!  Purpose: This routine is called after the time loop ends and is     !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.              ! 
!           This routine is not called from an IJK loop, hence         !
!           all indices are undefined.                                 !
!                                                                      !
!  Author:                                            Date: dd-mmm-yy  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!  Revision Number:                                                    !
!  Purpose:                                                            !
!  Author:                                            Date: dd-mmm-yy  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!                                                                      !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
      SUBROUTINE USR3
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      use mms, only         : deallocate_mms_vars
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


      Call calculate_de_norms

      if(tecplot_output) Call write_tecplot_data

      Call deallocate_mms_vars


      RETURN
      END SUBROUTINE USR3



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_de_norms	                               !
!  Purpose: Calculates discretization error norms for the problem      !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Feb 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_de_norms
      use param, only     : dimension_3
      use param1, only    : zero
      use compar, only    : ijkstart3, ijkend3
      use functions, only : funijk_gl
      use functions, only : IS_ON_myPE_owns
      use indices, only   : i_of, j_of, k_of
      use usr             ! all variables
      use mms             ! all variables
      use fldvar, only    : p_g, u_g, v_g, w_g, u_s, v_s, w_s, &
                            t_g, t_s, ep_g, rop_s, theta_m
      use geometry, only  : imax, imin1, imax1
      use geometry, only  : jmax, jmin1, jmax1
      use geometry, only  : kmax, kmin1, kmax1
      use funits, only    : newunit
      use compar, only    : myPE, PE_IO
      use mpi_utility, only : global_all_sum
      IMPLICIT NONE

! number of data points for norm calculation
        integer   :: var_size

! looping variables
        integer   :: ijk, i, j, k

! file unit
        integer   :: f1

! global ijk
        integer   :: ijk_gl


! allocate de_ variables and lnorms_ variables
        allocate(de_ep_g(dimension_3))
        allocate(de_p_g(dimension_3))
        allocate(de_u_g(dimension_3))
        allocate(de_v_g(dimension_3))
        allocate(de_w_g(dimension_3))
        allocate(de_t_g(dimension_3))
        allocate(de_rop_s(dimension_3))
        allocate(de_u_s(dimension_3))
        allocate(de_v_s(dimension_3))
        allocate(de_w_s(dimension_3))
        allocate(de_t_s(dimension_3))
        allocate(de_theta_m(dimension_3))

        allocate(lnorms_ep_g(3))
        allocate(lnorms_p_g(3))
        allocate(lnorms_u_g(3))
        allocate(lnorms_v_g(3))
        allocate(lnorms_w_g(3))
        allocate(lnorms_t_g(3))
        allocate(lnorms_rop_s(3))
        allocate(lnorms_u_s(3))
        allocate(lnorms_v_s(3))
        allocate(lnorms_w_s(3))
        allocate(lnorms_t_s(3))
        allocate(lnorms_theta_m(3))

! determine pressure shift value
        delta_p_g = zero
        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.not.(IS_ON_myPE_owns(i,j,k))) cycle

          ijk_gl = funijk_gl(i,j,k)
          if(ijk_sh.eq.ijk_gl) then
            delta_p_g = P_G(ijk) - MMS_P_G(ijk)
            write(*,*) "delta_p_g= ", delta_p_g
          endif

        end do

        call global_all_sum(delta_p_g)

! scalar variables
        de_ep_g = zero
        de_p_g = zero
        de_t_g = zero
        de_rop_s = zero
        de_t_s = zero
        de_theta_m = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle

! select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle
          if((k.lt.kmin1).or.(k.gt.kmax1)) cycle

          de_ep_g(ijk)    = ep_g(ijk)     - mms_ep_g(ijk)
          de_p_g(ijk)     = p_g(ijk)      - mms_p_g(ijk) - delta_p_g
          de_t_g(ijk)     = t_g(ijk)      - mms_t_g(ijk)
          de_rop_s(ijk)   = rop_s(ijk,1)  - mms_rop_s(ijk)
          de_t_s(ijk)     = t_s(ijk,1)    - mms_t_s(ijk)
          de_theta_m(ijk) = theta_m(ijk,1)- mms_theta_m(ijk)

        end do

        var_size = imax*jmax*kmax
        call calculate_lnorms(var_size, de_ep_g, lnorms_ep_g)
        call calculate_lnorms(var_size, de_p_g, lnorms_p_g)
        call calculate_lnorms(var_size, de_t_g, lnorms_t_g)
        call calculate_lnorms(var_size, de_rop_s, lnorms_rop_s)
        call calculate_lnorms(var_size, de_t_s, lnorms_t_s)
        call calculate_lnorms(var_size, de_theta_m, lnorms_theta_m)

! vector variables (x)
        de_u_g = zero
        de_u_s = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.ge.imax1)) cycle  ! note: >=
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle
          if((k.lt.kmin1).or.(k.gt.kmax1)) cycle

          de_u_g(ijk) = u_g(ijk) - mms_u_g(ijk)
          de_u_s(ijk) = u_s(ijk,1) - mms_u_s(ijk)
        end do

        var_size = (imax-1)*jmax*kmax
        call calculate_lnorms(var_size, de_u_g, lnorms_u_g)
        call calculate_lnorms(var_size, de_u_s, lnorms_u_s)

! vector variables (y)
        de_v_g = zero
        de_v_s = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.ge.jmax1)) cycle  ! note: >=
          if((k.lt.kmin1).or.(k.gt.kmax1)) cycle

          de_v_g(ijk) = v_g(ijk) - mms_v_g(ijk)
          de_v_s(ijk) = v_s(ijk,1) - mms_v_s(ijk)
        end do

        var_size = imax*(jmax-1)*kmax
        call calculate_lnorms(var_size, de_v_g, lnorms_v_g)
        call calculate_lnorms(var_size, de_v_s, lnorms_v_s)

! vector variables (z)
        de_w_g = zero
        de_w_s = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle
          if((k.lt.kmin1).or.(k.ge.kmax1)) cycle  ! note: >=

          de_w_g(ijk) = w_g(ijk) - mms_w_g(ijk)
          de_w_s(ijk) = w_s(ijk,1) - mms_w_s(ijk)
        end do

        var_size = imax*jmax*(kmax-1)
        call calculate_lnorms(var_size, de_w_g, lnorms_w_g)
        call calculate_lnorms(var_size, de_w_s, lnorms_w_s)

! Output DE norms data to a file
        if(myPE == PE_IO) then
          open(unit=newunit(f1), file="de_norms.dat", status='unknown')
          write(f1,*) "# DE Norms for MMS test cases:"
          write(f1,*) "# imax= ",imax, " jmax=", jmax, " kmax=", kmax
          write(f1,*) "# 1st line: L1 Norms, 2nd line: L2 Norms, &
                       &3rd line: Linf Norms"
          write(f1,*) "# variables='p_g''u_g''v_g''w_g''u_s''v_s''w_s'&
                      &'t_g''t_s''ep_g''rop_s''theta_m'"
          write(f1,*) lnorms_p_g(1), &
            lnorms_u_g(1), lnorms_v_g(1), lnorms_w_g(1), &
            lnorms_u_s(1), lnorms_v_s(1), lnorms_w_s(1), &
            lnorms_t_g(1), lnorms_t_s(1), &
            lnorms_ep_g(1), lnorms_rop_s(1), &
            lnorms_theta_m(1)
          write(f1,*) lnorms_p_g(2), &
            lnorms_u_g(2), lnorms_v_g(2), lnorms_w_g(2), &
            lnorms_u_s(2), lnorms_v_s(2), lnorms_w_s(2), &
            lnorms_t_g(2), lnorms_t_s(2), &
            lnorms_ep_g(2), lnorms_rop_s(2), &
            lnorms_theta_m(2)
          write(f1,*) lnorms_p_g(3), &
            lnorms_u_g(3), lnorms_v_g(3), lnorms_w_g(3), &
            lnorms_u_s(3), lnorms_v_s(3), lnorms_w_s(3), &
            lnorms_t_g(3), lnorms_t_s(3), &
            lnorms_ep_g(3), lnorms_rop_s(3), &
            lnorms_theta_m(3)
          close(f1)
        end if

! de allocate de_ variables and lnorms_ variables
        deallocate(de_ep_g)
        deallocate(de_p_g)
        deallocate(de_u_g)
        deallocate(de_v_g)
        deallocate(de_w_g)
        deallocate(de_t_g)
        deallocate(de_rop_s)
        deallocate(de_u_s)
        deallocate(de_v_s)
        deallocate(de_w_s)
        deallocate(de_t_s)
        deallocate(de_theta_m)

        deallocate(lnorms_ep_g)
        deallocate(lnorms_p_g)
        deallocate(lnorms_u_g)
        deallocate(lnorms_v_g)
        deallocate(lnorms_w_g)
        deallocate(lnorms_t_g)
        deallocate(lnorms_rop_s)
        deallocate(lnorms_u_s)
        deallocate(lnorms_v_s)
        deallocate(lnorms_w_s)
        deallocate(lnorms_t_s)
        deallocate(lnorms_theta_m)

      RETURN

      END SUBROUTINE calculate_de_norms


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_lnorms	                               !
!  Purpose: Calculates L1, L2 and Lifinity norms for any input         !
!  variable of defined size                                            !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jan 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_lnorms(var_size, var, lnorms)
      use compar, only      : ijkstart3, ijkend3
      use param, only       : dimension_3
      use param1, only      : zero
      use mpi_utility, only : global_sum, global_max
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


        if(var_size.eq.0) then
          lnorms(:) = zero
          return
        end if

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
!  Module name: write_tecplot_data                                     !
!  Purpose: Write data for visualization in Tecplot. Only variables    !
!  that are relevant to MMS setup have been coded, but similar code    !
!  can be added for other variables, if needed. This module uses a     !
!  simple strategy of gathering all data to root for output.           !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Feb 2015   !
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
      use geometry, only    : imin2, imin1, imax, imax1, imax2
      use geometry, only    : jmin2, jmin1, jmax, jmax1, jmax2
      use geometry, only    : kmin2, kmin1, kmax, kmax1, kmax2
      use geometry, only    : dx, dy, dz
      use compar, only      : myPE, PE_IO
      use funits, only      : newunit
      use functions, only   : funijk_gl
      use fldvar, only      : p_g, u_g, v_g, w_g, u_s, v_s, w_s, &
                              t_g, t_s, ep_g, rop_s, theta_m
      use mms               ! all variables
      use usr, only         : tec_output_point, tec_output_block, &
                              tec_output_raw, tec_no_k
      use usr, only         : delta_p_g
      use param1, only      : half
      use mpi_utility, only : gather
      IMPLICIT NONE

! indices
        integer             :: i, j, k, ijk
        integer             :: imjk, ijmk, ijkm

! file units
        integer             :: fcc, ftcc, fs, fvx, fvy, fvz

! arrays to gather variables
        double precision, allocatable :: arr_Pg(:)
        double precision, allocatable :: arr_MMSPg(:)
!        double precision, allocatable :: arr_PgSh(:)
        double precision, allocatable :: arr_Epg(:)
        double precision, allocatable :: arr_MMSEpg(:)
        double precision, allocatable :: arr_Ug(:)
        double precision, allocatable :: arr_MMSUg(:)
        double precision, allocatable :: arr_Vg(:)
        double precision, allocatable :: arr_MMSVg(:)
        double precision, allocatable :: arr_Wg(:)
        double precision, allocatable :: arr_MMSWg(:)
        double precision, allocatable :: arr_ROPs(:)
        double precision, allocatable :: arr_MMSROPs(:)
        double precision, allocatable :: arr_Us(:)
        double precision, allocatable :: arr_MMSUs(:)
        double precision, allocatable :: arr_Vs(:)
        double precision, allocatable :: arr_MMSVs(:)
        double precision, allocatable :: arr_Ws(:)
        double precision, allocatable :: arr_MMSWs(:)
        double precision, allocatable :: arr_Tg(:)
        double precision, allocatable :: arr_MMSTg(:)
        double precision, allocatable :: arr_Ts(:)
        double precision, allocatable :: arr_MMSTs(:)
        double precision, allocatable :: arr_Ths(:)
        double precision, allocatable :: arr_MMSThs(:)
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
                                          Pg_tmp, MMSPg_tmp, &
                                          Epg_tmp, MMSEPg_tmp, &
                                          Ug_tmp, MMSUg_tmp, &
                                          Vg_tmp, MMSVg_tmp, &
                                          Wg_tmp, MMSWg_tmp, &
                                          Rops_tmp, MMSRops_tmp, &
                                          Us_tmp, MMSUs_tmp, &
                                          Vs_tmp, MMSVs_tmp, &
                                          Ws_tmp, MMSWs_tmp, &
                                          Tg_tmp, MMSTg_tmp, &
                                          Ts_tmp, MMSTs_tmp, &
                                          Ths_tmp, MMSThs_tmp


! allocate variables as ijkmax3 size 1D arrays
        allocate(arr_Pg(ijkmax3))
        allocate(arr_MMSPg(ijkmax3))
!        allocate(arr_PgSh(ijkmax3))
        allocate(arr_Epg(ijkmax3))
        allocate(arr_MMSEpg(ijkmax3))
        allocate(arr_Ug(ijkmax3))
        allocate(arr_MMSUg(ijkmax3))
        allocate(arr_Vg(ijkmax3))
        allocate(arr_MMSVg(ijkmax3))
        allocate(arr_Wg(ijkmax3))
        allocate(arr_MMSWg(ijkmax3))
        allocate(arr_ROPs(ijkmax3))
        allocate(arr_MMSROPs(ijkmax3))
        allocate(arr_Us(ijkmax3))
        allocate(arr_MMSUs(ijkmax3))
        allocate(arr_Vs(ijkmax3))
        allocate(arr_MMSVs(ijkmax3))
        allocate(arr_Ws(ijkmax3))
        allocate(arr_MMSWs(ijkmax3))
        allocate(arr_Tg(ijkmax3))
        allocate(arr_MMSTg(ijkmax3))
        allocate(arr_Ts(ijkmax3))
        allocate(arr_MMSTs(ijkmax3))
        allocate(arr_Ths(ijkmax3))
        allocate(arr_MMSThs(ijkmax3))
        allocate(arr_xtr(ijkmax3))
        allocate(arr_ytr(ijkmax3))
        allocate(arr_ztr(ijkmax3))

! gather the MFIX global variables into ijkmax3 sized variables
        call gather(P_G,arr_Pg,PE_IO)
        call gather(MMS_P_G,arr_MMSPg,PE_IO)
!        call gather(P_G_Sh,arr_PgSh,PE_IO) ! P_G_Sh not used globally
        call gather(EP_g,arr_Epg,PE_IO)
        call gather(MMS_EP_G,arr_MMSEpg,PE_IO)
        call gather(U_G,arr_Ug,PE_IO)
        call gather(MMS_U_G,arr_MMSUg,PE_IO)
        call gather(V_G,arr_Vg,PE_IO)
        call gather(MMS_V_G,arr_MMSVg,PE_IO)
        call gather(W_G,arr_Wg,PE_IO)
        call gather(MMS_W_G,arr_MMSWg,PE_IO)
        call gather(ROP_s(:,1),arr_ROPs,PE_IO)
        call gather(MMS_ROP_s,arr_MMSROPs,PE_IO)
        call gather(U_s(:,1),arr_Us,PE_IO)
        call gather(MMS_U_s,arr_MMSUs,PE_IO)
        call gather(V_s(:,1),arr_Vs,PE_IO)
        call gather(MMS_V_s,arr_MMSVs,PE_IO)
        call gather(W_s(:,1),arr_Ws,PE_IO)
        call gather(MMS_W_s,arr_MMSWs,PE_IO)
        call gather(T_g,arr_Tg,PE_IO)
        call gather(MMS_T_g,arr_MMSTg,PE_IO)
        call gather(T_s(:,1),arr_Ts,PE_IO)
        call gather(MMS_T_s,arr_MMSTs,PE_IO)
        call gather(Theta_m(:,1),arr_Ths,PE_IO)
        call gather(MMS_Theta_m,arr_MMSThs,PE_IO)
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

            Pg_tmp(i,j,k) = arr_Pg(IJK) - delta_p_g
            MMSPg_tmp(i,j,k) = arr_MMSPg(IJK)

            Epg_tmp(i,j,k) = arr_Epg(IJK)
            MMSEpg_tmp(i,j,k) = arr_MMSEpg(IJK)

            Ug_tmp(i,j,k) = HALF*(arr_Ug(ijk)+arr_Ug(imjk))
            MMSUg_tmp(i,j,k) = HALF*(arr_MMSUg(ijk)+arr_MMSUg(imjk))

            Vg_tmp(i,j,k) = HALF*(arr_Vg(ijk)+arr_Vg(ijmk))
            MMSVg_tmp(i,j,k) = HALF*(arr_MMSVg(ijk)+arr_MMSVg(ijmk))

            Wg_tmp(i,j,k) = HALF*(arr_Wg(ijk)+arr_Wg(ijkm))
            MMSWg_tmp(i,j,k) = HALF*(arr_MMSWg(ijk)+arr_MMSWg(ijkm))

            Rops_tmp(i,j,k) = arr_ROPs(IJK)
            MMSRops_tmp(i,j,k) = arr_MMSROPs(IJK)

            Us_tmp(i,j,k) = HALF*(arr_Us(ijk)+arr_Us(imjk))
            MMSUs_tmp(i,j,k) = HALF*(arr_MMSUs(ijk)+arr_MMSUs(imjk))

            Vs_tmp(i,j,k) = HALF*(arr_Vs(ijk)+arr_Vs(ijmk))
            MMSVs_tmp(i,j,k) = HALF*(arr_MMSVs(ijk)+arr_MMSVs(ijmk))

            Ws_tmp(i,j,k) = HALF*(arr_Ws(ijk)+arr_Ws(ijkm))
            MMSWs_tmp(i,j,k) = HALF*(arr_MMSWs(ijk)+arr_MMSWs(ijkm))

            Tg_tmp(i,j,k) = arr_Tg(IJK)
            MMSTg_tmp(i,j,k) = arr_MMSTg(IJK)

            Ts_tmp(i,j,k) = arr_Ts(IJK)
            MMSTs_tmp(i,j,k) = arr_MMSTs(IJK)

            Ths_tmp(i,j,k) = arr_Ths(IJK)
            MMSThs_tmp(i,j,k) = arr_MMSThs(IJK)

          end do
          end do
          end do

! write the cell centered data in point format
! Here, we calculate and output cell-centered data at the cell-centers
! locations.  Otheriwse, when cell-centered data is outputted with
! VARLOCATION=CELLCENTERED option, tecplot does inaccurate interpolation
! near the boundary cells.  So this is not the traditional
! tecplot cell-centered data file.
          
          if(tec_output_point) then

            open(unit=newunit(fcc), file="solution_tec_point.dat", &
                                    status='unknown')
            write(fcc,"(27a)") 'variables = "x""y""z"&
                  &"Pg""Ug""Vg""Wg""Us""Vs""Ws"&
                  &"Tg""Ts""Epg""Rops""Ths"&
                  &"MMSPg""MMSUg""MMSVg""MMSWg"&
                  &"MMSUs""MMSVs""MMSWs"&
                  &"MMSTg""MMSTs""MMSEpg""MMSRops""MMSThs"'
            write(fcc,*) 'zone T="',0,'" '
            write(fcc,*) 'I=',IMAX,' J=',JMAX,' K=',KMAX


! write all data at cell-centers to output file
            do k = kmin1, kmax1
            do j = jmin1, jmax1
            do i = imin1, imax1

              ijk = funijk_gl(i,j,k)

              xt = arr_xtr(ijk) - dx(i)*half
              yt = arr_ytr(ijk) - dy(j)*half
              zt = arr_ztr(ijk) - dz(k)*half

              write(fcc,*) xt, yt, zt, &
                Pg_tmp(i,j,k), &
                Ug_tmp(i,j,k), Vg_tmp(i,j,k), Wg_tmp(i,j,k),&
                Us_tmp(i,j,k), Vs_tmp(i,j,k), Ws_tmp(i,j,k),&
                Tg_tmp(i,j,k), Ts_tmp(i,j,k), &
                Epg_tmp(i,j,k), Rops_tmp(i,j,k),&
                Ths_tmp(i,j,k), &
                MMSPg_tmp(i,j,k), &
                MMSUg_tmp(i,j,k), MMSVg_tmp(i,j,k), MMSWg_tmp(i,j,k),&
                MMSUs_tmp(i,j,k), MMSVs_tmp(i,j,k), MMSWs_tmp(i,j,k),&
                MMSTg_tmp(i,j,k), MMSTs_tmp(i,j,k), &
                MMSEpg_tmp(i,j,k), MMSRops_tmp(i,j,k),&
                MMSThs_tmp(i,j,k)

            end do
            end do
            end do

            close(fcc)

          end if ! end of if(tec_output_point)

! write true tecplot cell-centered data (block format)
! This is the traditional tecplot cell-centered data. Be wary while
! using contour plots which could do inaccurate interpolation near the
! boundaries. It's better to use tecplot's "primary flood" with this
! option.

          if(tec_output_block) then

            open(unit=newunit(ftcc), file="solution_tec_block.dat", &
                                      status='unknown')

            write(ftcc,"(27a)") 'variables = "x""y""z"&
                &"Pg""Ug""Vg""Wg""Us""Vs""Ws"&
                &"Tg""Ts""Epg""Rops""Ths"&
                &"MMSPg""MMSUg""MMSVg""MMSWg"&
                &"MMSUs""MMSVs""MMSWs"&
                &"MMSTg""MMSTs""MMSEpg""MMSRops""MMSThs"'
            write(ftcc,*) 'zone T="',0,'" '
            write(ftcc,*) 'I=',IMAX1,' J=',JMAX1,' K=',KMAX1

            write(ftcc,*) 'DATAPACKING=BLOCK'
            write(ftcc,*) "VARLOCATION=([4-27]=CELLCENTERED)"

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
              write(ftcc,*) (Wg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Us_tmp(i,j,k),i=imin1,imax1)
            end do
            end do


            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Vs_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Ws_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Tg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Ts_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Epg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Rops_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Ths_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSPg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSUg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSVg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSWg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSUs_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSVs_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSWs_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSTg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSTs_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSEpg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSRops_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (MMSThs_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

          end if  ! end of if(tec_output_block)

! Output data where calculated:
! i.e., at cell centers for scalar variables, and
!       at face centers for vector variables
! This is useful for debugging

          if(tec_output_raw) then

! output scalar variables
            open(unit=newunit(fs), file="solution_scalar.dat", &
                                   status='unknown')

            write(fs,"(21a)") 'variables = "x""y""z"&
              &"Pg""Tg""Ts""Epg""Rops""Ths"&
              &"MMSPg""MMSTg""MMSTs""MMSEpg""MMSRops""MMSThs"&
              &"PgErr""TgErr""TsErr""EpgErr""RopsErr""ThsErr"'
            write(fs,*) 'zone T="',0,'" '
            write(fs,*) 'I=',IMAX2,' J=',JMAX2,' K=',KMAX2

            do k = kmin2, kmax2
            do j = jmin2, jmax2
            do i = imin2, imax2

              ijk = funijk_gl(i,j,k)

              xt = arr_xtr(ijk) - dx(i)*half
              yt = arr_ytr(ijk) - dy(j)*half
              zt = arr_ztr(ijk) - dz(k)*half

              Pg_tmp(i,j,k) = arr_Pg(IJK) - delta_p_g
              MMSPg_tmp(i,j,k) = arr_MMSPg(IJK)

              Tg_tmp(i,j,k) = arr_Tg(IJK)
              MMSTg_tmp(i,j,k) = arr_MMSTg(IJK)

              Ts_tmp(i,j,k) = arr_Ts(IJK)
              MMSTs_tmp(i,j,k) = arr_MMSTs(IJK)

              Epg_tmp(i,j,k) = arr_Epg(IJK)
              MMSEpg_tmp(i,j,k) = arr_MMSEpg(IJK)

              Rops_tmp(i,j,k) = arr_ROPs(IJK)
              MMSRops_tmp(i,j,k) = arr_MMSROPs(IJK)

              Ths_tmp(i,j,k) = arr_Ths(IJK)
              MMSThs_tmp(i,j,k) = arr_MMSThs(IJK)

              write(fs,*) xt, yt, zt, &
                Pg_tmp(i,j,k), &
                Tg_tmp(i,j,k), Ts_tmp(i,j,k), &
                Epg_tmp(i,j,k), Rops_tmp(i,j,k),&
                Ths_tmp(i,j,k), &
                MMSPg_tmp(i,j,k), &
                MMSTg_tmp(i,j,k), MMSTs_tmp(i,j,k), &
                MMSEpg_tmp(i,j,k), MMSRops_tmp(i,j,k),&
                MMSThs_tmp(i,j,k), &
                Pg_tmp(i,j,k) - MMSPg_tmp(i,j,k), &
                Tg_tmp(i,j,k) - MMSTg_tmp(i,j,k), &
                Ts_tmp(i,j,k) - MMSTs_tmp(i,j,k), &
                Epg_tmp(i,j,k) - MMSEpg_tmp(i,j,k), &
                Rops_tmp(i,j,k) - MMSRops_tmp(i,j,k), &
                Ths_tmp(i,j,k) - MMSThs_tmp(i,j,k)

            end do
            end do
            end do

            close(fs)

! output vector variables (x-direction)
            open(unit=newunit(fvx), file="solution_vectorx.dat", &
                                   status='unknown')

            write(fvx,"(9a)") 'variables = "x""y""z""Ug""Us"&
              &"MMSUg""MMSUs""UgErr""UsErr"'
            write(fvx,*) 'zone T="',0,'" '
            write(fvx,*) 'I=',IMAX2,' J=',JMAX2,' K=',KMAX2

            do k = kmin2, kmax2
            do j = jmin2, jmax2
            do i = imin2, imax2

              ijk = funijk_gl(i,j,k)

              xt = arr_xtr(ijk)
              yt = arr_ytr(ijk) - dy(j)*half
              zt = arr_ztr(ijk) - dz(k)*half

              Ug_tmp(I,J,K) = arr_Ug(IJK)
              MMSUg_tmp(I,J,K) = arr_MMSUg(IJK)

              Us_tmp(I,J,K) = arr_Us(IJK)
              MMSUs_tmp(I,J,K) = arr_MMSUs(IJK)

              write(fvx,*) xt, yt, zt, &
                Ug_tmp(i,j,k), &
                Us_tmp(i,j,k), &
                MMSUg_tmp(i,j,k), &
                MMSUs_tmp(i,j,k), &
                Ug_tmp(i,j,k) - MMSUg_tmp(i,j,k), &
                Us_tmp(i,j,k) - MMSUs_tmp(i,j,k)

            end do
            end do
            end do

            close(fvx)

! output vector variables (y-direction)
            open(unit=newunit(fvy), file="solution_vectory.dat", &
                                   status='unknown')

            write(fvy,"(9a)") 'variables = "x""y""z""Vg""Vs"&
              &"MMSVg""MMSVs""VgErr""VsErr"'
            write(fvx,*) 'zone T="',0,'" '
            write(fvx,*) 'I=',IMAX2,' J=',JMAX2,' K=',KMAX2

            do k = kmin2, kmax2
            do j = jmin2, jmax2
            do i = imin2, imax2

              ijk = funijk_gl(i,j,k)

              xt = arr_xtr(ijk) - dx(i)*half
              yt = arr_ytr(ijk)
              zt = arr_ztr(ijk) - dz(k)*half

              Vg_tmp(I,J,K) = arr_Vg(IJK)
              MMSVg_tmp(I,J,K) = arr_MMSVg(IJK)

              Vs_tmp(I,J,K) = arr_Vs(IJK)
              MMSVs_tmp(I,J,K) = arr_MMSVs(IJK)

              write(fvx,*) xt, yt, zt, &
                Vg_tmp(i,j,k), &
                Vs_tmp(i,j,k), &
                MMSVg_tmp(i,j,k), &
                MMSVs_tmp(i,j,k), &
                Vg_tmp(i,j,k) - MMSVg_tmp(i,j,k), &
                Vs_tmp(i,j,k) - MMSVs_tmp(i,j,k)

            end do
            end do
            end do

            close(fvx)

! output vector variables (z-direction)
            open(unit=newunit(fvz), file="solution_vectorz.dat", &
                                   status='unknown')

            write(fvz,"(9a)") 'variables = "x""y""z""Wg""Ws"&
              &"MMSWg""MMSWs""WgErr""WsErr"'
            write(fvz,*) 'zone T="',0,'" '
            write(fvz,*) 'I=',IMAX2,' J=',JMAX2,' K=',KMAX2

            do k = kmin2, kmax2
            do j = jmin2, jmax2
            do i = imin2, imax2

              ijk = funijk_gl(i,j,k)

              xt = arr_xtr(ijk) - dx(i)*half
              yt = arr_ytr(ijk) - dy(j)*half
              zt = arr_ztr(ijk)

              Wg_tmp(I,J,K) = arr_Wg(IJK)
              MMSWg_tmp(I,J,K) = arr_MMSWg(IJK)

              Ws_tmp(I,J,K) = arr_Ws(IJK)
              MMSWs_tmp(I,J,K) = arr_MMSWs(IJK)

              write(fvz,*) xt, yt, zt, &
                Wg_tmp(i,j,k), &
                Ws_tmp(i,j,k), &
                MMSWg_tmp(i,j,k), &
                MMSWs_tmp(i,j,k), &
                Wg_tmp(i,j,k) - MMSWg_tmp(i,j,k), &
                Ws_tmp(i,j,k) - MMSWs_tmp(i,j,k)

            end do
            end do
            end do

            close(fvz)


          end if ! end if(raw_output)

        end if ! end of if((myPE==PE_IO)


        RETURN

      END SUBROUTINE write_tecplot_data
