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

!------------------------------------------------------------------------
      module mpi_comm_des

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

!-----------------------------------------------

! generic interface definition
      interface des_gather
         module procedure des_gather_l,des_gather_i,des_gather_d
      end interface

      contains

!------------------------------------------------------------------------
! Subroutine       : desmpi_sendrecv_init
! Purpose          : posts asynchronous send and recv and updates the request id
!
! Parameters       : pface - face number (1to6)
!                    debug - for printing debug statments
!------------------------------------------------------------------------
      subroutine desmpi_sendrecv_init(pface,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer,intent(in) :: pface
      integer,intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character(len=80), parameter :: name = 'desmpi_sendrecv_init'
      integer :: ldebug,ltag,lerr,lrecvface
!-----------------------------------------------

! set the debug flag
      ldebug = 0
      if (present(pdebug)) then
        ldebug = pdebug
      endif



!direct copy in case of single processor
      lrecvface = pface+mod(pface,2)-mod(pface+1,2)




      if (ineighproc(pface).eq.mype) then
         drecvbuf(1+mod(lrecvface,2))%facebuf(1:isendcnt(pface)) = &
            dsendbuf(1+mod(pface,2))%facebuf(1:isendcnt(pface))
      else
         ltag = message_tag(ineighproc(pface),mype,pface)
         call des_mpi_irecv(drecvbuf(1+mod(pface,2))%facebuf(:),imaxbuf, &
                            ineighproc(pface),ltag,irecvreq(pface),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )

         ltag = message_tag(mype,ineighproc(pface),lrecvface)
         call des_mpi_isend(dsendbuf(1+mod(pface,2))%facebuf(:),isendcnt(pface), &
                        ineighproc(pface),ltag,isendreq(pface),lerr)
         call mpi_check( name //':mpi_isend ', lerr )

      end if
      return

    contains

      integer function message_tag(lsource,ldest,lrecvface)
        implicit none
        integer, intent(in) :: lsource,ldest,lrecvface
        message_tag = lsource+numpes*ldest+numpes*numpes*lrecvface+100
      end function message_tag

    end subroutine desmpi_sendrecv_init

!------------------------------------------------------------------------
! Subroutine       : desmpi_sendrecv_wait
! Purpose          : waits for the communication for the specified interface
!
! Parameters       : pface - face number (1to6)
!                    debug - for printing debug statments
!------------------------------------------------------------------------
      subroutine desmpi_sendrecv_wait(pface,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer,intent(in) :: pface
      integer,intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character(len=80), parameter :: name = 'desmpi_sendrecv_wait'
      integer :: ldebug,lerr
!-----------------------------------------------

! set the debug flag
      ldebug = 0
      if (present(pdebug)) then
        ldebug = pdebug
      endif

! wait for both send and recv request completes
      if (ineighproc(pface).ne.mype) then
         call des_mpi_wait(isendreq(pface),lerr)
         call mpi_check( name //':mpi_wait-send', lerr )
         call des_mpi_wait(irecvreq(pface),lerr)
         call mpi_check( name //':mpi_wait-recv', lerr )
      end if
      return
      end subroutine desmpi_sendrecv_wait


!------------------------------------------------------------------------
! Subroutine       : desmpi_scatterv
! Purpose          : scatters the particle from PE_IO
! Parameters       : ptype - flag for datatype integer (1) or double precision (2)
!                    pdebug - optional flag for debugging
!------------------------------------------------------------------------
      subroutine desmpi_scatterv(ptype,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: ptype
      integer, intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lroot,lidebug,lerr
      character(len=80), parameter :: name = 'desmpi_scatterv'
!-----------------------------------------------

      lroot = pe_io
      if (.not. present(pdebug)) then
         lidebug = 0
      else
         lidebug = pdebug
      endif

      if (ptype .eq. 1) then
         call des_MPI_Scatterv(irootbuf,iscattercnts,idispls, &
                               iprocbuf,iscr_recvcnt,lroot,lerr)
      else
         call des_MPI_Scatterv(drootbuf,iscattercnts,idispls, &
                               dprocbuf,iscr_recvcnt,lroot,lerr)
      end if
      call MPI_Check( name //':MPI_Scatterv', lerr )

      return
      end subroutine desmpi_scatterv


!------------------------------------------------------------------------
! Subroutine       : desmpi_gatherv
! Purpose          : gathers the particle from local proc to root proc
! Parameters       : ptype - flag for datatype integer (1) or double precision (2)
!                    pdebug - optional flag for debugging
!------------------------------------------------------------------------
      subroutine desmpi_gatherv(ptype,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: ptype
      integer, intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lroot,lidebug,lerr
      character(len=80), parameter :: name = 'des_gather'
!-----------------------------------------------

      lroot = pe_io
      if (.not. present(pdebug)) then
         lidebug = 0
      else
         lidebug = pdebug
      endif
      if(ptype.eq.1) then
         call des_MPI_Gatherv(iprocbuf,igath_sendcnt,irootbuf, &
                              igathercnts,idispls,lroot,lerr)
      else
         call des_MPI_Gatherv(dprocbuf,igath_sendcnt,drootbuf, &
                              igathercnts,idispls,lroot,lerr)
      end if
      call MPI_Check( name //':MPI_Gatherv', lerr )

      return
      end subroutine desmpi_gatherv


!------------------------------------------------------------------------
! Subroutine       : des_gather_d
! Purpose          : gathers double precision array from local to root
! Parameters       :
!                    parray - array to be writen
!------------------------------------------------------------------------
      subroutine des_gather_d(parray)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      double precision, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar,lparcount,lcount
!-----------------------------------------------

! pack the variables in case of
      lparcount = 1
      lcount = 0
      do lcurpar = 1, max_pip
         if (lparcount.gt.pip) exit
         if (is_nonexistent(lcurpar)) cycle
         lparcount = lparcount +1
         if(is_ghost(lcurpar) .or. is_entering_ghost(lcurpar) .or. is_exiting_ghost(lcurpar)) cycle
         lcount = lcount + 1
         dprocbuf(lcount) = parray(lcurpar)
      end do
      call desmpi_gatherv(ptype=2)
      end subroutine des_gather_d

!------------------------------------------------------------------------
! Subroutine       : des_gather_l
! Purpose          : gathers logical array from local to root
! Parameters       :
!                    parray - array to be writen
!------------------------------------------------------------------------
      subroutine des_gather_l(parray)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      logical, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar,lparcount,lcount
!-----------------------------------------------

! pack the variables in proc buffer
      lparcount = 1
      lcount = 0
      do lcurpar = 1, max_pip
         if (lparcount.gt.pip) exit
         if (is_nonexistent(lcurpar)) cycle
         lparcount = lparcount +1
         if(is_ghost(lcurpar) .or. is_entering_ghost(lcurpar) .or. is_exiting_ghost(lcurpar)) cycle
         lcount = lcount + 1
         if(parray(lcurpar)) then
            iprocbuf(lcount) = 1
         else
            iprocbuf(lcount) = 0
         end if
      end do
      call desmpi_gatherv(ptype=1)

      end subroutine des_gather_l

!------------------------------------------------------------------------
! Subroutine       : des_gather_i
! Purpose          : gathers integer array from local to root
! Parameters       :
!                    parray - array to be writen
!                    ploc2glb - this flag is used to conver local particle
!                    number into global particle number (used for history
!                    and neighbour terms)
!------------------------------------------------------------------------
      subroutine des_gather_i(parray,ploc2glb)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, dimension(:) :: parray
      logical,optional :: ploc2glb
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar,lparcount,lcount
      logical :: lloc2glb
!-----------------------------------------------

      if (present(ploc2glb)) then
         lloc2glb = ploc2glb
      else
         lloc2glb = .false.
      end if
! pack the variables in proc buffer
      lparcount = 1
      lcount = 0
      if (lloc2glb) then
         do lcurpar = 1, max_pip
            if (lparcount.gt.pip) exit
            if (is_nonexistent(lcurpar)) cycle
            lparcount = lparcount +1
            if(is_ghost(lcurpar) .or. is_entering_ghost(lcurpar) .or. is_exiting_ghost(lcurpar)) cycle
            lcount = lcount + 1
            if(parray(lcurpar).gt.0) then
               iprocbuf(lcount) = iglobal_id(parray(lcurpar))
            else
               iprocbuf(lcount) = 0
            end if
         end do
      else
         do lcurpar = 1, max_pip
            if (lparcount.gt.pip) exit
            if (is_nonexistent(lcurpar)) cycle
            lparcount = lparcount +1
            if(is_ghost(lcurpar) .or. is_entering_ghost(lcurpar) .or. is_exiting_ghost(lcurpar)) cycle
            lcount = lcount + 1
            iprocbuf(lcount) = parray(lcurpar)
         end do
      end if
      call desmpi_gatherv(ptype=1)

      end subroutine des_gather_i

      end module mpi_comm_des
