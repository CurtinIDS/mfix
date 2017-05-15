!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: sweep_and_prune                                             !
!                                                                      !
!  Purpose: maintains a sorted list of axis-aligned bounding boxes     !
!                                                                      !
!           main type is sap_t                                         !
!                                                                      !
!  Reference: http://www.codercorner.com/SAP.pdf                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

#include "version.inc"

module sweep_and_prune

  use pair_manager
  use list

  ! endpoint data structure should be 64 bits in size
  type endpoint_t
     ! owner of the endpoint
     ! if Min endpoint, id = -box_id
     ! if Max endpoint, id = box_id
     integer(kind=4) :: box_id

     ! actual value of the endpoint; single precision to save space
     real :: value
  end type endpoint_t

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type: aabb_t                                                        !
  !                                                                      !
  !  Purpose: Represents an axis-aligned bounding box                    !
  !           Used to add new boxes to sap_t and multisap_t              !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  type aabb_t
     real :: maxendpoint(3)
     real :: minendpoint(3)
  end type aabb_t

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type: box_t                                                         !
  !                                                                      !
  !  Purpose: Represents the AABB of a particle                          !
  !           maxendpoint(1) and minendpoint_id(1) point to x_endpoints  !
  !           maxendpoint(2) and minendpoint_id(2) point to y_endpoints  !
  !           maxendpoint(3) and minendpoint_id(3) point to z_endpoints  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  type box_t
     integer :: maxendpoint_id(3)
     integer :: minendpoint_id(3)
     integer :: particle_id
  end type box_t

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type: sap_t                                                         !
  !                                                                      !
  !  Purpose: Sorts the endpoints array                                  !
  !           axis is (1,2,3) to represent x,y,z                         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  type sap_t
     ! id to identify within a multisap
     integer :: id

     ! list of endpoints of boxes
     type(endpoint_t), dimension(:), allocatable :: x_endpoints
     type(endpoint_t), dimension(:), allocatable :: y_endpoints
     type(endpoint_t), dimension(:), allocatable :: z_endpoints
     integer :: x_endpoints_len
     integer :: y_endpoints_len
     integer :: z_endpoints_len

     ! list of boxes
     type(box_t), dimension(:), allocatable :: boxes
     integer :: boxes_len

     ! list of indices in boxes that have value 0 (deleted boxes)
     type(list_t) :: deleted_boxes

     ! list of pairs of particle ids that have intersecting boxes
     type(hashtable_t) :: hashtable
  end type sap_t

  public :: init_sap, add_box, del_box, update_box, sort, sweep
  private :: partition, fullcheck, boxes_grow, endpoints_grow, sort_endpoints, quicksort_endpoints

contains

  subroutine check(this)
    use discretelement
    implicit none
    type(sap_t), intent(inout) :: this
    integer :: ii, pid
    do ii = 1, this%x_endpoints_len

       if (0 > this%x_endpoints(ii)%box_id) then

       if (ii .ne. this%boxes(-this%x_endpoints(ii)%box_id)%minendpoint_id(1)  ) then
          print *,"ii is ",ii
          print *,"this%x_endpoints(ii)%box_id is ",this%x_endpoints(ii)%box_id
          print *,"this%boxes(this%x_endpoints(ii)%box_id)%minendpoint_id(1)",this%boxes(-this%x_endpoints(ii)%box_id)%minendpoint_id(1)
          stop __LINE__
       endif
else
       if (ii .ne. this%boxes(this%x_endpoints(ii)%box_id)%maxendpoint_id(1)  ) then
          print *,"ii is ",ii
          print *,"this%x_endpoints(ii)%box_id is ",this%x_endpoints(ii)%box_id
          print *,"this%boxes(this%x_endpoints(ii)%box_id)%minendpoint_id(1)",this%boxes(this%x_endpoints(ii)%box_id)%maxendpoint_id(1)

          stop __LINE__
       endif
endif
    enddo

    do ii = 1, this%y_endpoints_len
       if (0 > this%y_endpoints(ii)%box_id) then
          if (ii .ne. this%boxes(-this%y_endpoints(ii)%box_id)%minendpoint_id(2)  ) then
             print *,"ii is ",ii
             print *,"this%y_endpoints(ii)%box_id is ",this%y_endpoints(ii)%box_id
             print *,"this%boxes(this%y_endpoints(ii)%box_id)%minendpoint_id(1)",this%boxes(-this%y_endpoints(ii)%box_id)%minendpoint_id(2)

             stop __LINE__
          endif
       else
          if (ii .ne. this%boxes(this%y_endpoints(ii)%box_id)%maxendpoint_id(2)  ) then
             stop __LINE__
          endif
       endif
    enddo

    do ii = 1, this%z_endpoints_len
       if (0 > this%z_endpoints(ii)%box_id) then
          if (ii .ne. this%boxes(-this%z_endpoints(ii)%box_id)%minendpoint_id(3)  ) then
             stop __LINE__
          endif
       else
          if (ii .ne. this%boxes(this%z_endpoints(ii)%box_id)%maxendpoint_id(3)  ) then
             stop __LINE__
          endif
       endif
    enddo

    do ii = 1, this%boxes_len
       pid = this%boxes(ii)%particle_id
       if ( 0.0001 < abs(this%x_endpoints(this%boxes(ii)%maxendpoint_id(1))%value - this%x_endpoints(this%boxes(ii)%minendpoint_id(1))%value - 2*des_radius(pid))) then

          print *,this%x_endpoints(this%boxes(ii)%maxendpoint_id(1))%value
          print *,this%x_endpoints(this%boxes(ii)%minendpoint_id(1))%value
          print *,2*des_radius(pid)

          stop __LINE__
       endif
    enddo

    do ii = 1, this%boxes_len
       pid = this%boxes(ii)%particle_id
       if ( 0.0001 < abs(this%y_endpoints(this%boxes(ii)%maxendpoint_id(2))%value - this%y_endpoints(this%boxes(ii)%minendpoint_id(2))%value - 2*des_radius(pid))) then
          stop __LINE__
       endif
    enddo

    do ii = 1, this%boxes_len
       pid = this%boxes(ii)%particle_id
       if ( 0.0001 < abs(this%z_endpoints(this%boxes(ii)%maxendpoint_id(3))%value - this%z_endpoints(this%boxes(ii)%minendpoint_id(3))%value - 2*des_radius(pid))) then
          stop __LINE__
       endif
    enddo

  end subroutine check

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: init_sap                                                !
  !                                                                      !
  !  Purpose: sap_t constructor                                          !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine init_sap(this,id)
    implicit none
    type(sap_t), intent(inout) :: this
    integer, intent(in) :: id

    this%x_endpoints_len = 0
    this%y_endpoints_len = 0
    this%z_endpoints_len = 0
    this%boxes_len = 0
    this%id = id

    allocate (this%x_endpoints(10))
    allocate (this%y_endpoints(10))
    allocate (this%z_endpoints(10))
    allocate (this%boxes(10))

    call list_init(this%deleted_boxes)

    call init_pairs(this%hashtable)

  end subroutine init_sap

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: boxes_GROW                                              !
  !                                                                      !
  !  Purpose: resize array of box_t                                      !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  SUBROUTINE boxes_GROW(box_array,new_size)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: new_size
    type(box_t), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: box_array
    type(box_t), DIMENSION(:), ALLOCATABLE :: box_tmp
    INTEGER lSIZE

    lSIZE = size(box_array,1)
    allocate(box_tmp(new_size))
    box_tmp(1:lSIZE) = box_array(1:lSIZE)
    call move_alloc(box_tmp,box_array)

  END SUBROUTINE boxes_GROW

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: endpoints_GROW                                          !
  !                                                                      !
  !  Purpose: Resizes array of endpoint_t                                !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  SUBROUTINE endpoints_GROW(endpoint_array,new_size)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: new_size
    type(endpoint_t), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: endpoint_array
    type(endpoint_t), DIMENSION(:), ALLOCATABLE :: endpoint_tmp
    INTEGER lSIZE

    lSIZE = size(endpoint_array,1)
    allocate(endpoint_tmp(new_size))
    endpoint_tmp(1:lSIZE) = endpoint_array(1:lSIZE)
    call move_alloc(endpoint_tmp,endpoint_array)

  END SUBROUTINE endpoints_GROW

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: add_box                                                 !
  !                                                                      !
  !  Purpose: Adds new box_t to sap_t that represents particle_id        !
  !           Returns id to refer to the box later                       !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine add_box(this,aabb,particle_id,id)

    implicit none
    type(sap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(in) :: particle_id
    integer, intent(out) :: id

    if (this%x_endpoints_len+2 > size(this%x_endpoints)) then
       call endpoints_GROW(this%x_endpoints,2*(this%x_endpoints_len+2))
    endif

    if (this%y_endpoints_len+2 > size(this%y_endpoints)) then
       call endpoints_GROW(this%y_endpoints,2*(this%y_endpoints_len+2))
    endif

    if (this%z_endpoints_len+2 > size(this%z_endpoints)) then
       call endpoints_GROW(this%z_endpoints,2*(this%z_endpoints_len+2))
    endif

    if (this%boxes_len+1 > size(this%boxes)) then
       call boxes_GROW(this%boxes,2*(this%boxes_len+1))
    endif

    id = list_pop(this%deleted_boxes)

    ! if there are no deleted boxes (no holes in boxes array), then append the
    ! new box at the end of the boxes array.
    if ( id < 0 ) then
       this%boxes_len = this%boxes_len + 1
       id = this%boxes_len
    endif

    ! no holes in the endpoints array; new endpoints are always appended
    this%x_endpoints_len = this%x_endpoints_len + 2
    this%y_endpoints_len = this%y_endpoints_len + 2
    this%z_endpoints_len = this%z_endpoints_len + 2

    this%boxes(id)%particle_id = particle_id
    this%boxes(id)%maxendpoint_id(1) = this%x_endpoints_len
    this%boxes(id)%maxendpoint_id(2) = this%y_endpoints_len
    this%boxes(id)%maxendpoint_id(3) = this%z_endpoints_len
    this%boxes(id)%minendpoint_id(1) = this%x_endpoints_len-1
    this%boxes(id)%minendpoint_id(2) = this%y_endpoints_len-1
    this%boxes(id)%minendpoint_id(3) = this%z_endpoints_len-1

    this%x_endpoints(this%x_endpoints_len-1)%box_id = -id
    this%x_endpoints(this%x_endpoints_len-1)%value = aabb%minendpoint(1)
    this%y_endpoints(this%y_endpoints_len-1)%box_id = -id
    this%y_endpoints(this%y_endpoints_len-1)%value = aabb%minendpoint(2)
    this%z_endpoints(this%z_endpoints_len-1)%box_id = -id
    this%z_endpoints(this%z_endpoints_len-1)%value = aabb%minendpoint(3)

    this%x_endpoints(this%x_endpoints_len)%box_id = id
    this%x_endpoints(this%x_endpoints_len)%value = aabb%maxendpoint(1)
    this%y_endpoints(this%y_endpoints_len)%box_id = id
    this%y_endpoints(this%y_endpoints_len)%value = aabb%maxendpoint(2)
    this%z_endpoints(this%z_endpoints_len)%box_id = id
    this%z_endpoints(this%z_endpoints_len)%value = aabb%maxendpoint(3)

  end subroutine add_box

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: del_box                                                 !
  !                                                                      !
  !  Purpose: Removes box from sap_t                                     !
  !           Setting all the endpoints to HUGE will remove all pairs    !
  !           associated with this box as the endpoints get sorted.      !
  !           The box is "really" deleted at the end of sort().          !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine del_box(this,id)
    implicit none
    type(sap_t), intent(inout) :: this
    integer, intent(in) :: id

    this%x_endpoints(this%boxes(id)%maxendpoint_id(1))%value = HUGE(0.0)
    this%y_endpoints(this%boxes(id)%maxendpoint_id(2))%value = HUGE(0.0)
    this%z_endpoints(this%boxes(id)%maxendpoint_id(3))%value = HUGE(0.0)
    this%x_endpoints(this%boxes(id)%minendpoint_id(1))%value = HUGE(0.0)
    this%y_endpoints(this%boxes(id)%minendpoint_id(2))%value = HUGE(0.0)
    this%z_endpoints(this%boxes(id)%minendpoint_id(3))%value = HUGE(0.0)

  end subroutine del_box

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: update_box                                              !
  !                                                                      !
  !  Purpose: Update box represented by id with values in aabb           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine update_box(this,id,aabb)
    implicit none
    type(sap_t), intent(inout) :: this
    integer, intent(in) :: id
    type(aabb_t), intent(in) :: aabb

    ! update end points
    this%x_endpoints(this%boxes(id)%maxendpoint_id(1))%value = aabb%maxendpoint(1)
    this%y_endpoints(this%boxes(id)%maxendpoint_id(2))%value = aabb%maxendpoint(2)
    this%z_endpoints(this%boxes(id)%maxendpoint_id(3))%value = aabb%maxendpoint(3)
    this%x_endpoints(this%boxes(id)%minendpoint_id(1))%value = aabb%minendpoint(1)
    this%y_endpoints(this%boxes(id)%minendpoint_id(2))%value = aabb%minendpoint(2)
    this%z_endpoints(this%boxes(id)%minendpoint_id(3))%value = aabb%minendpoint(3)

  end subroutine update_box

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: sort                                                    !
  !                                                                      !
  !  Purpose: (Re)sorts particles after each timestep,                   !
  !           which keeps the pair manager updated.                      !
  !           Calls sort_endpoints, which does an insertion sort,        !
  !           which is fast for nearly sorted arrays but slow for        !
  !           unsorted arrays.                                           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine sort(this)
    use geometry, only: do_k
    implicit none
    type(sap_t), intent(inout) :: this
    integer :: ii

    ! sort end points
    call sort_endpoints(this,this%x_endpoints(1:this%x_endpoints_len),1)
    call sort_endpoints(this,this%y_endpoints(1:this%y_endpoints_len),2)
    if (do_k) call sort_endpoints(this,this%z_endpoints(1:this%z_endpoints_len),3)

    ! cleanup HUGE endpoints
    !
    ! After sorting the endpoints, there may be endpoints at the end of the
    ! array corresponding to boxes that are being deleted. Now that the sorting
    ! is done, the corresponding ids in boxes array can be set to zero, and
    ! added to the deleted boxes array to be used next time we add a box.

    ! TODO: what if two overlapping boxes are deleted at the same time...could
    ! they leave a pair in the pair manager corresponding to deleted particles?
    ! Not sure how to handle that

    ii = this%x_endpoints_len
    do
       if (ii .le. 0) exit
       if (this%x_endpoints(ii)%value .ne. HUGE(0.0)) exit
       this%boxes(abs(this%x_endpoints(ii)%box_id))%minendpoint_id(:) = 0
       this%boxes(abs(this%x_endpoints(ii)%box_id))%maxendpoint_id(:) = 0
       call list_add(this%deleted_boxes,this%x_endpoints(ii)%box_id)
       ii = ii - 1
    enddo

    this%x_endpoints_len = ii
    this%y_endpoints_len = ii
    this%z_endpoints_len = ii

  end subroutine sort

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: sweep                                                   !
  !                                                                      !
  !  Purpose: Creates initial set of pairs of colliding boxes            !
  !           Called on startup after quicksort is called                !
  !           for each axis set of endpoints                             !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine sweep(this)

    use pair_manager, only: add_pair
    use geometry, only: no_k

    implicit none
    type(sap_t), intent(inout) :: this

    ! active list of box id's
    integer :: ii,aa,minmax,ai

    type(list_t) :: active

    call list_init(active)

    do ii=1, this%x_endpoints_len

       minmax = this%x_endpoints(ii)%box_id

       if ( minmax < 0) then
          ! add pairs for new box x ( active pairs )

          do ai=1, list_get_length(active)
             aa = list_get(active,ai)
             ! compare other two axes

             if (max(this%boxes(aa)%minendpoint_id(2),this%boxes(-minmax)%minendpoint_id(2)) <= min(this%boxes(-minmax)%maxendpoint_id(2),this%boxes(aa)%maxendpoint_id(2)) .and. &
                  (NO_K .or. max(this%boxes(aa)%minendpoint_id(3),this%boxes(-minmax)%minendpoint_id(3)) <= min(this%boxes(-minmax)%maxendpoint_id(3),this%boxes(aa)%maxendpoint_id(3)))) then
                call add_pair(this%hashtable,this%boxes(-minmax)%particle_id,this%boxes(aa)%particle_id)
             endif
          enddo

          ! add new box to active pair list
          call list_add(active,-minmax)

       else if ( 0 < minmax ) then
          ! remove box from active pair list
          call list_del(active,minmax)
       else
          print *,"SAP_ID=",this%id,"minmax shouldn't be zero: ",minmax
          ERROR_STOP __LINE__
       endif
    enddo

  end subroutine sweep

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: sort_endpoints                                          !
  !                                                                      !
  !  Purpose: Sorts the endpoints array                                  !
  !           axis is (1,2,3) to represent x,y,z                         !
  !                                                                      !
  !           Sorts endpoints array after each timestep.                 !
  !           Uses insertion sort.                                       !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine sort_endpoints(this, endpoints, axis)
    use pair_manager, only: add_pair, del_pair
    implicit none
    type(sap_t), intent(inout) :: this
    integer, intent(in) :: axis
    type(endpoint_t), dimension(:), intent(inout) :: endpoints
    integer :: ii, jj
    type(endpoint_t) :: sweeppoint, swappoint
    real :: sweepval

    do ii=2, size(endpoints)

       sweeppoint = endpoints(ii)
       sweepval = sweeppoint%value
       jj = ii-1

       do while ( 0 < jj )
          if ( endpoints(jj)%value <= sweepval ) then
             exit
          endif

          swappoint = endpoints(jj)

          if (sweeppoint%box_id < 0 .and. 0 < swappoint%box_id ) then
             if (fullcheck(this,-sweeppoint%box_id,swappoint%box_id,axis)) then
                call add_pair(this%hashtable,this%boxes(-sweeppoint%box_id)%particle_id,this%boxes(swappoint%box_id)%particle_id)
             endif
          endif

          if (0 < sweeppoint%box_id .and. swappoint%box_id < 0 ) then
             if (fullcheck(this,sweeppoint%box_id,-swappoint%box_id,axis)) then
                call del_pair(this%hashtable,this%boxes(sweeppoint%box_id)%particle_id,this%boxes(-swappoint%box_id)%particle_id)
             endif
          endif

          endpoints(jj+1) = swappoint

          ! if box_id is zero, that is an endpoint corresponding to a deleted box
          if (swappoint%box_id < 0) then
             this%boxes(-swappoint%box_id)%minendpoint_id(axis) = jj+1
          else if (swappoint%box_id > 0) then
             this%boxes(swappoint%box_id)%maxendpoint_id(axis) = jj+1
          endif

          jj = jj - 1
       enddo

       endpoints(jj+1) = sweeppoint

       ! if box_id is zero, that is an endpoint corresponding to a deleted box
       if (sweeppoint%box_id < 0) then
          this%boxes(-sweeppoint%box_id)%minendpoint_id(axis) = jj+1
       else if (sweeppoint%box_id > 0) then
          this%boxes(sweeppoint%box_id)%maxendpoint_id(axis) = jj+1
       endif

    enddo

  end subroutine sort_endpoints

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Function: fullcheck                                                 !
  !                                                                      !
  !  Purpose: Check whether or checking that the boxes for id and id2    !
  !           are correctly sorted                                       !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  logical function fullcheck(this,id,id2,curr_axis)
    use geometry, only: no_k
    implicit none
    type(sap_t), intent(inout) :: this
    ! box ids to compare
    integer, intent(in) :: id, id2, curr_axis

    fullcheck = ((curr_axis.eq.1 .or. max(this%boxes(id)%minendpoint_id(1),this%boxes(id2)%minendpoint_id(1)) &
                                      <= min(this%boxes(id2)%maxendpoint_id(1),this%boxes(id)%maxendpoint_id(1))) &
           .and. (curr_axis.eq.2 .or. max(this%boxes(id)%minendpoint_id(2),this%boxes(id2)%minendpoint_id(2)) &
                                       <= min(this%boxes(id2)%maxendpoint_id(2),this%boxes(id)%maxendpoint_id(2))) &
 .and. (NO_K .or. curr_axis.eq.3 .or. max(this%boxes(id)%minendpoint_id(3),this%boxes(id2)%minendpoint_id(3)) &
                                       <= min(this%boxes(id2)%maxendpoint_id(3),this%boxes(id)%maxendpoint_id(3))))

  end function fullcheck

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: quicksort_endpoints                                     !
  !                                                                      !
  !  Purpose: Sorts the endpoints arrays                                 !
  !                                                                      !
  !           Only used once at startup, just before a call to sweep.    !
  !           Quicksort is fast for unsorted arrays, but slow for        !
  !           nearly sorted arrays.  sort() is used instead after each   !
  !           timestep.                                                  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine quicksort(this)
    implicit none
    type(sap_t), intent(inout) :: this

    call quicksort_endpoints(this,this%x_endpoints(1:this%x_endpoints_len),1)
    call quicksort_endpoints(this,this%y_endpoints(1:this%y_endpoints_len),2)
    call quicksort_endpoints(this,this%z_endpoints(1:this%z_endpoints_len),3)

  end subroutine quicksort

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: quicksort_endpoints                                     !
  !                                                                      !
  !  Purpose: Sorts the endpoints array                                  !
  !           axis is (1,2,3) to represent x,y,z                         !
  !                                                                      !
  !           Only used once at startup, just before a call to sweep.    !
  !           Quicksort is fast for unsorted arrays, but slow for        !
  !           nearly sorted arrays.  sort() is used instead after each   !
  !           timestep.                                                  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  recursive subroutine quicksort_endpoints(this, endpoints, axis)
    type(sap_t), intent(inout) :: this
    integer, intent(in) :: axis
    type(endpoint_t), intent(inout), dimension(:) :: endpoints
    integer :: ii

    !call check_boxes(this)

    if(size(endpoints) > 1) then
       call partition(this, endpoints, ii, axis)
       call quicksort_endpoints(this,endpoints(:ii-1),axis)
       call quicksort_endpoints(this,endpoints(ii:),axis)
    endif

  end subroutine quicksort_endpoints

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  subroutine: partition                                               !
  !                                                                      !
  !  purpose: partition subroutine for quicksort of endpoints array      !
  !           box is represented by id                                   !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine partition(this, endpoints, marker, axis)
    type(sap_t), intent(inout) :: this
    integer, intent(in) :: axis
    type(endpoint_t), intent(in out), dimension(:) :: endpoints
    integer, intent(out) :: marker
    type(endpoint_t) :: temp
    integer :: tmp_ii, tmp_jj
    integer :: i, j
    real :: x      ! pivot point
    x = endpoints(1)%value
    i= 0
    j= size(endpoints) + 1

    do
       j = j-1
       do
          if (endpoints(j)%value <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (endpoints(i)%value >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange endpoints(i) and endpoints(j)

          if (endpoints(i)%box_id < 0) then
             tmp_ii = this%boxes(-endpoints(i)%box_id)%minendpoint_id(axis)
          else
             tmp_ii = this%boxes(endpoints(i)%box_id)%maxendpoint_id(axis)
          endif
          if (endpoints(j)%box_id < 0) then
             tmp_jj = this%boxes(-endpoints(j)%box_id)%minendpoint_id(axis)
          else
             tmp_jj = this%boxes(endpoints(j)%box_id)%maxendpoint_id(axis)
          endif

          temp = endpoints(i)
          endpoints(i) = endpoints(j)
          endpoints(j) = temp

          if (endpoints(i)%box_id < 0) then
             this%boxes(-endpoints(i)%box_id)%minendpoint_id(axis) = tmp_ii
          else
             this%boxes(endpoints(i)%box_id)%maxendpoint_id(axis) = tmp_ii
          endif
          if (endpoints(j)%box_id < 0) then
             this%boxes(-endpoints(j)%box_id)%minendpoint_id(axis) = tmp_jj
          else
             this%boxes(endpoints(j)%box_id)%maxendpoint_id(axis) = tmp_jj
          endif

       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine partition

end module sweep_and_prune
