!------------------------------------------------------------------------
! Module           : desmpi
! Purpose          : Contains wrapper class for mpi communications- send,recv
!
! Author           : Pradeep.G
!
! Purpose          : Module contains subroutines and variables related to
!                    des mpi communication.
!
! Comments         : do_nsearch flag should be set to true before calling
!                    des_par_exchange; when do_nsearch is true ghost particles of the
!                    system will be updated, which will be later used to generate
!                    neighbour list.
!
! Contains following subroutines:
!    des_addnodevalues, des_addnodevalues2, des_addnodevalues_mean_fields
!------------------------------------------------------------------------
      module mpi_node_des

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use parallel_mpi
      use mpi_utility
      use discretelement
      use desgrid
      use compar
      use physprop
      use sendrecv
      use des_bc
      use desmpi_wrapper
      use sendrecvnode
      use mfix_pic
      use des_thermo
      use run, only: ENERGY_EQ,ANY_SPECIES_EQ
      use param, only: DIMENSION_N_s
      use des_rxns
      use desmpi

      use mpi_comm_des, only: desmpi_sendrecv_init
      use mpi_comm_des, only: desmpi_sendrecv_wait

      contains


!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues_mean_fields
! Purpose          : This routine is specially used for computing mean
!                    fields by backward interpolation.
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues_mean_fields()

      use functions
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: li, lj, lk
      integer :: lm,lijkmin,lijkmax
!-----------------------------------------------

! fill the temporary buffer
      DO LM = 1,DES_MMAX+MMAX
         CALL DES_EXCHANGENODE(DES_ROPS_NODE(:,LM),PADD=.TRUE.)
         DO LI =1,DIMN
            CALL DES_EXCHANGENODE(DES_VEL_NODE(:,LI,LM),PADD=.TRUE.)
         END DO
      END DO

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            lijkmin = funijk(1,lj,lk)
            lijkmax = funijk(imax1,lj,lk)
            des_rops_node(lijkmin,:)  = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_vel_node(lijkmin,:,:) = des_vel_node(lijkmin,:,:)+des_vel_node(lijkmax,:,:)
            des_rops_node(lijkmax,:)  = des_rops_node(lijkmin,:)
            des_vel_node(lijkmax,:,:) = des_vel_node(lijkmin,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            lijkmin = funijk(li,1,lk)
            lijkmax = funijk(li,jmax1,lk)
            des_rops_node(lijkmin,:)  = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_vel_node(lijkmin,:,:) = des_vel_node(lijkmin,:,:)+des_vel_node(lijkmax,:,:)
            des_rops_node(lijkmax,:)  = des_rops_node(lijkmin,:)
            des_vel_node(lijkmax,:,:) = des_vel_node(lijkmin,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1 .and. do_K) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            lijkmin = funijk(li,lj,1)
            lijkmax = funijk(li,lj,kmax1)
            des_rops_node(lijkmin,:)  = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_vel_node(lijkmin,:,:) = des_vel_node(lijkmin,:,:)+des_vel_node(lijkmax,:,:)
            des_rops_node(lijkmax,:)  = des_rops_node(lijkmin,:)
            des_vel_node(lijkmax,:,:) = des_vel_node(lijkmin,:,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues_mean_fields



!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues
! Purpose          : This routine is specially used for des_drag_gs
!                    The backward interpolation in des_drag_gs computes
!                    the grid node values of drag_am and drag_bm
!                    node values are from istart2 to iend1;
!                    hence a separate module is created to exchange
!                    node values
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues()

      use functions
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: li, lj, lk
      integer :: lijkmin,lijkmax
!-----------------------------------------------

! fill the temporary buffer
      call des_exchangenode(drag_am, padd=.true.)
      do li =1,dimn
         call des_exchangenode(drag_bm(:,li), padd=.true.)
      end do

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            lijkmin = funijk(1,lj,lk)
            lijkmax = funijk(imax1,lj,lk)
            drag_am(lijkmin) = drag_am(lijkmin)+drag_am(lijkmax)
            drag_bm(lijkmin,:) = drag_bm(lijkmin,:)+drag_bm(lijkmax,:)
            drag_am(lijkmax) = drag_am(lijkmin)
            drag_bm(lijkmax,:) = drag_bm(lijkmin,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            lijkmin = funijk(li,1,lk)
            lijkmax = funijk(li,jmax1,lk)
            drag_am(lijkmin) = drag_am(lijkmin)+drag_am(lijkmax)
            drag_bm(lijkmin,:) = drag_bm(lijkmin,:)+drag_bm(lijkmax,:)
            drag_am(lijkmax) = drag_am(lijkmin)
            drag_bm(lijkmax,:) = drag_bm(lijkmin,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1 .and. do_K) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            lijkmin = funijk(li,lj,1)
            lijkmax = funijk(li,lj,kmax1)
            drag_am(lijkmin) = drag_am(lijkmin)+drag_am(lijkmax)
            drag_bm(lijkmin,:) = drag_bm(lijkmin,:)+drag_bm(lijkmax,:)
            drag_am(lijkmax) = drag_am(lijkmin)
            drag_bm(lijkmax,:) = drag_bm(lijkmin,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues


!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues2
! Purpose          : This routine is specially used for calc_des_rop_s
!                    The backward interpolation in calc_des_rop_s computes
!                    the grid node values of des_rops_node
!                    node values are from istart2 to iend1;
!                    hence a separate module is created to exchange
!                    node values
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues2()

      use functions
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: li, lj, lk
      integer :: lm,lijkmin,lijkmax
!-----------------------------------------------

! fill the temporary buffer
      do lm = 1,DES_MMAX+MMAX
         call des_exchangenode(des_rops_node(:,lm),padd=.true.)
      end do

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            lijkmin = funijk(1,lj,lk)
            lijkmax = funijk(imax1,lj,lk)
            des_rops_node(lijkmin,:) = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_rops_node(lijkmax,:) = des_rops_node(lijkmin,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            lijkmin = funijk(li,1,lk)
            lijkmax = funijk(li,jmax1,lk)
            des_rops_node(lijkmin,:) = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_rops_node(lijkmax,:) = des_rops_node(lijkmin,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1 .and. do_K) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            lijkmin = funijk(li,lj,1)
            lijkmax = funijk(li,lj,kmax1)
            des_rops_node(lijkmin,:) = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_rops_node(lijkmax,:) = des_rops_node(lijkmin,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues2



      end module mpi_node_des
