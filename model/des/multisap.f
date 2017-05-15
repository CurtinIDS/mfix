!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: multi_sweep_and_prune                                       !
!                                                                      !
!  Purpose: Divides space into a grid of sap_t instances and adds      !
!           AABB's to the corresponding sap_t instance(s)              !
!                                                                      !
!  Reference: http://www.codercorner.com/SAP.pdf                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

#include "version.inc"

module multi_sweep_and_prune

  use sweep_and_prune

  type multisap_t
     ! grid particle, e.g. 20x20x20
     integer :: grid(3)

     ! bounds of the grid (grid saps on the boundary extend to infinity)
     real :: minbounds(3), maxbounds(3)

     ! grid(:)/(maxbounds(:)-minbounds(:))
     real :: one_over_cell_length(3)

     ! list saps of length grid(1)*grid(2)*grid(3)
     type(sap_t), dimension(:), allocatable :: saps

     ! union of all hashtables for all saps in this multisap
     type(hashtable_t) :: hashtable
  end type multisap_t

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type: boxhandle_t                                                   !
  !                                                                      !
  !  Purpose: Represents a reference to a particle box in a particular   !
  !           sap for this multisap.                                     !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  type boxhandle_t
     integer :: sap_id
     integer :: box_id
  end type boxhandle_t

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type: boxhandlelist_t                                               !
  !                                                                      !
  !  Purpose: List of boxhandle_t instances that corresponding to a      !
  !           particular particle with id particle_id.                   !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  ! The maximum number of saps an aabb can belong to. For 2d this could be 4,
  ! but no harm in just setting it to 8.
  integer, parameter :: MAX_SAPS = 8

  type boxhandlelist_t
     integer :: particle_id
     type(boxhandle_t) :: list(MAX_SAPS)
  end type boxhandlelist_t

 public :: init_multisap, multisap_add, multisap_del, multisap_update, multisap_sort, multisap_quicksort, multisap_sweep, boxhandle_grow
  private :: multisap_raster

contains
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: init_sap                                                !
  !                                                                      !
  !  Purpose: multisap_t constructor                                     !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine init_multisap(this,x_grid,y_grid,z_grid2,minbounds,maxbounds)
    use geometry
    implicit none
    type(multisap_t), intent(inout) :: this
    integer, intent(in) :: x_grid,y_grid, z_grid2
    real, dimension(3), intent(in) :: minbounds, maxbounds
    integer :: ii,jj,kk, z_grid, sap_id

    if (NO_K) then
       z_grid = 1
    else
       z_grid = z_grid2
    endif

    this%grid(1) = max(x_grid,1)
    this%grid(2) = max(y_grid,1)
    this%grid(3) = max(z_grid,1)

    allocate(this%saps(0:this%grid(1)*this%grid(2)*this%grid(3)-1))
    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = (ii*this%grid(2)+jj)*this%grid(3)+kk
             call init_sap(this%saps(sap_id),sap_id)
          enddo
       enddo
    enddo

    this%one_over_cell_length(:) = this%grid(:)/(maxbounds(:)-minbounds(:))

    this%minbounds(:) = minbounds(:)
    this%maxbounds(:) = maxbounds(:)

  end subroutine init_multisap

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: multisap_raster                                         !
  !                                                                      !
  !  Purpose: For a given aabb, return the list sap_ids of ids of the    !
  !           saps that the aabb intersects according to the grid.       !
  !           Used by multisap_add and multisap_del.                     !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine multisap_raster(this,aabb,sap_ids)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(out) :: sap_ids(MAX_SAPS)
    integer :: max_grid(3),min_grid(3)

    sap_ids(:) = -1

    ! bin to grid
    min_grid(:) = floor((aabb%minendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:))
    max_grid(:) = floor((aabb%maxendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:))

    min_grid(1) = min(max(min_grid(1),0),this%grid(1)-1)
    min_grid(2) = min(max(min_grid(2),0),this%grid(2)-1)
    min_grid(3) = min(max(min_grid(3),0),this%grid(3)-1)
    max_grid(1) = min(max(max_grid(1),0),this%grid(1)-1)
    max_grid(2) = min(max(max_grid(2),0),this%grid(2)-1)
    max_grid(3) = min(max(max_grid(3),0),this%grid(3)-1)

    call add_to_set((min_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + min_grid(3))
    call add_to_set((min_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + max_grid(3))
    call add_to_set((min_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + min_grid(3))
    call add_to_set((min_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + max_grid(3))
    call add_to_set((max_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + min_grid(3))
    call add_to_set((max_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + max_grid(3))
    call add_to_set((max_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + min_grid(3))
    call add_to_set((max_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + max_grid(3))

  contains

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !                                                                      !
    !  Subroutine: add_to_set                                              !
    !                                                                      !
    !  Purpose: Add sap_id to sap_ids if it's not already in there.        !
    !                                                                      !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    subroutine add_to_set(sap_id)
      integer, intent(in) :: sap_id
      integer :: mm

      do mm=1, size(sap_ids)
         if (sap_ids(mm).eq.sap_id) then
            ! already in the list
            return
         endif
         if (sap_ids(mm).eq.-1) then
            sap_ids(mm) = sap_id
            return
         endif
      enddo

    end subroutine add_to_set

  end subroutine multisap_raster

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: multisap_add                                            !
  !                                                                      !
  !  Purpose: Add boxes in appropriate saps for aabb representing        !
  !           particle_id and return list of (sap_id,box_id) handles     !
  !           in handlelist.                                             !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine multisap_add(this,aabb,particle_id,handlelist)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(in) :: particle_id
    type(boxhandlelist_t), intent(out) :: handlelist
    integer :: sap_ids(MAX_SAPS)
    integer :: nn

    call multisap_raster(this,aabb,sap_ids)

    handlelist%particle_id = particle_id

    ! add to each individual SAP
    do nn=1, size(sap_ids)
       handlelist%list(nn)%sap_id = sap_ids(nn)
       if (sap_ids(nn) >= 0) then
          call add_box(this%saps(sap_ids(nn)),aabb,particle_id,handlelist%list(nn)%box_id)
       endif
    enddo

  end subroutine multisap_add

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: multisap_del                                            !
  !                                                                      !
  !  Purpose: Delete all sap entries corresponding to handlelist.        !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine multisap_del(this,handlelist)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(boxhandlelist_t), intent(in) :: handlelist
    integer :: nn

    ! remove from each SAP listed for this id
    do nn=1, size(handlelist%list)
       if (handlelist%list(nn)%sap_id >= 0) then
          call del_box(this%saps(handlelist%list(nn)%sap_id),handlelist%list(nn)%box_id)
       endif
    enddo

  end subroutine multisap_del

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: multisap_update                                         !
  !                                                                      !
  !  Purpose: Update the sap entries corresponding to handlelist with    !
  !           new location aabb.  Returns an updated handlelist with     !
  !           possibly different sap entries.                            !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine multisap_update(this,aabb,handlelist)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    type(boxhandlelist_t), intent(inout) :: handlelist
    integer, DIMENSION(MAX_SAPS) :: new_sap_ids
    logical :: found
    integer :: mm,nn,first_blank, box_id

    call multisap_raster(this,aabb,new_sap_ids)

    ! update for each SAP listed for this id
    do nn=1, size(new_sap_ids)
       if (new_sap_ids(nn) < 0) exit

       found = .false.
       first_blank = -1
       do mm=1, size(handlelist%list)
          if (handlelist%list(mm)%sap_id < 0) first_blank = mm

          if (handlelist%list(mm)%sap_id .eq. new_sap_ids(nn)) then
             ! update existing SAPs
             box_id = handlelist%list(mm)%box_id
             call update_box(this%saps(new_sap_ids(nn)),handlelist%list(mm)%box_id,aabb)
             found = .true.
             exit
          endif
       enddo
       if (.not.found) then
          ! add SAPs yet not listed for this id

          if (first_blank .eq. -1) then
             print *,"AABB is in more than 8 saps? Something is wrong:  ",first_blank,handlelist%list
             ERROR_STOP __LINE__
          endif
          call add_box(this%saps(new_sap_ids(nn)),aabb,handlelist%particle_id,handlelist%list(first_blank)%box_id)
          handlelist%list(first_blank)%sap_id = new_sap_ids(nn)
       endif
    enddo

    do mm=1, size(handlelist%list)
       if (handlelist%list(mm)%sap_id < 0) cycle

       found = .false.
       do nn=1, size(new_sap_ids)
          if (new_sap_ids(nn) < 0) exit

          if (handlelist%list(mm)%sap_id .eq. new_sap_ids(nn)) then
             found = .true.
             exit
          endif
       enddo

       if (.not.found) then
          ! remove SAP not found in handle%sap_id
          call del_box(this%saps(handlelist%list(mm)%sap_id),handlelist%list(mm)%box_id)
          handlelist%list(mm)%sap_id = -1
       endif

    enddo

  end subroutine multisap_update

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: multisap_sort                                           !
  !                                                                      !
  !  Purpose: Call sort on each of the saps in this multisap,            !
  !           then set the multisap's hashtable to the union of the      !
  !           hashtables of all the saps.                                !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine multisap_sort(this)
    ! use discretelement
    ! use functions, only: is_nonexistent
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk, sap_id
    integer :: pair(2)
    integer :: sighs, cc,ll,i, cc_start, cc_end

    call init_pairs(this%hashtable)

!$omp parallel default(none) private(ii,jj,kk,sap_id,pair) shared(this)
!$omp do collapse(3)
    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             call sort(this%saps(sap_id))
          enddo
       enddo
    enddo
!$omp end parallel

    print *," NUMBER OF DESGRIDS IS: ",this%grid(1)*this%grid(2)*this%grid(3)

 sighs = 0

 open (unit=235,file="pairs.txt",action="write",status="replace")

   do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             call reset_pairs(this%saps(sap_id)%hashtable)
             sighs = sighs + this%saps(sap_id)%hashtable%table_size
             do
                call get_pair(this%saps(sap_id)%hashtable,pair)
                if (pair(1).eq.0 .and. pair(2).eq.0) exit
                write (235,*) pair(1),pair(2)
                call add_pair(this%hashtable,pair(1),pair(2))

             enddo
          enddo
       enddo
    enddo

    close (unit=235)

    ! call reset_pairs(this%hashtable)

    ! do
    !    call get_pair(this%hashtable,pair)
    !    if (pair(1).eq.0 .and. pair(2).eq.0) exit
    !    write (235,*) pair(1),pair(2)
    ! enddo

    ! print *,"multisize ",this%hashtable%table_size
    ! print *,"sighs ",sighs
    ! stop __LINE__

  end subroutine multisap_sort

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: multisap_add                                            !
  !                                                                      !
  !  Purpose: Call quicksort on each of the saps in this multisap.       !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine multisap_quicksort(this)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    integer :: sap_id, boxcount

    boxcount = 0

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             call quicksort(this%saps(sap_id))
             boxcount = boxcount + this%saps(sap_id)%boxes_len
          enddo
       enddo
    enddo

  end subroutine multisap_quicksort

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: multisap_sweep                                          !
  !                                                                      !
  !  Purpose: Call sweep on each of the saps in this multisap.           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine multisap_sweep(this)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    integer :: sap_id

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = (ii*this%grid(2)+jj)*this%grid(3)+kk
             call sweep(this%saps(sap_id))
          enddo
       enddo
    enddo

  end subroutine multisap_sweep

  subroutine multisap_check(this)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    integer :: sap_id

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = (ii*this%grid(2)+jj)*this%grid(3)+kk

             call check(this%saps(sap_id))
          enddo
       enddo
    enddo

  end subroutine multisap_check

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: boxhandle_GROW                                          !
  !                                                                      !
  !  Purpose: resize array of box_t                                      !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine boxhandle_GROW(boxhandles,new_size)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: new_size
    type(boxhandlelist_t), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: boxhandles
    type(boxhandlelist_t), DIMENSION(:), ALLOCATABLE :: boxhandle_tmp
    INTEGER lSIZE

    lSIZE = size(boxhandles,1)
    allocate(boxhandle_tmp(new_size))
    boxhandle_tmp(1:lSIZE) = boxhandles(1:lSIZE)
    call move_alloc(boxhandle_tmp,boxhandles)

  end subroutine boxhandle_GROW

end module multi_sweep_and_prune
