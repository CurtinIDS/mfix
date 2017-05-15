#include "version.inc"

module list

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type active_t:                                                      !
  !                                                                      !
  !  Purpose: Represents an unsorted list of ids (positive integers).    !
  !           Implemented with an array, padded with zeroes.             !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
  type list_t
     integer, dimension(:), allocatable :: list
     integer :: list_len
  end type list_t

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: list_init                                               !
  !                                                                      !
  !  Purpose: Initialize list (constructor)                              !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine list_init(this)
    implicit none
    type(list_t), intent(inout) :: this

    allocate(this%list(10))
    this%list = 0
    this%list_len = 0

  end subroutine list_init

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: list_add                                                !
  !                                                                      !
  !  Purpose: Adds new value to the list                                 !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine list_add(this,value)
    use resize, only: integer_grow
    implicit none
    type(list_t), intent(inout) :: this
    integer, intent(in) :: value
    integer :: old_len

    this%list_len = this%list_len + 1
    if (size(this%list) < this%list_len) then
       old_len = size(this%list)
       call integer_GROW(this%list,this%list_len)
       this%list(old_len+1:size(this%list)) = 0
    endif

    if (0.eq.this%list(this%list_len)) then
       this%list(this%list_len) = value
    else
       print *,"LIST SHOULD END IN ZERO"
       print *,"value = ",value
       do old_len = 1, size(this%list)
          print *,"list: ",old_len,this%list(old_len)
       enddo
       ERROR_STOP __LINE__
    endif
  end subroutine list_add

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: list_pop                                                !
  !                                                                      !
  !  Purpose: Removes value from the end of the list                     !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  integer function list_pop(this)
    implicit none
    type(list_t), intent(inout) :: this

    if (this%list_len .eq. 0) then
       list_pop = -1
       return
    endif

    list_pop = this%list(this%list_len)
    this%list(this%list_len) = 0

    this%list_len = this%list_len - 1

  end function list_pop

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: list_del                                                !
  !                                                                      !
  !  Purpose: Removes value from array                                   !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine list_del(this,value)
    implicit none
    type(list_t), intent(inout) :: this
    integer, intent(in) :: value
    integer :: ai, aa

    do ai=1, size(this%list)
       aa = this%list(ai)
       if (value.eq.aa) then
          this%list(ai) = this%list(this%list_len)
          this%list(this%list_len) = 0
          this%list_len = this%list_len - 1
          exit
       endif
    enddo

  end subroutine list_del

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Function: list_get                                                  !
  !                                                                      !
  !  Purpose: Return element from list_t list                            !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  integer function list_get(this,index)
    implicit none
    type(list_t), intent(in) :: this
    integer :: index

    list_get = this%list(index)

  end function list_get

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Function: list_get_length                                           !
  !                                                                      !
  !  Purpose: Returns the length of list_t                               !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  integer function list_get_length(this)
    implicit none
    type(list_t), intent(in) :: this

    list_get_length = this%list_len

  end function list_get_length

end module list
