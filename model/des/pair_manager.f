!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: pair_manager                                                !
!                                                                      !
!  Purpose: maintains a hashtable of pairs of positive 32-bit          !
!           integers (meant to represent particle ids).                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
#include "version.inc"

module pair_manager

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type: pair_t                                                        !
  !                                                                      !
  !  Purpose: Represents a pair of particle ids.                         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  type pair_t
     integer(kind=4) :: ii
     integer(kind=4) :: jj
  end type pair_t

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Type: hashtable_t                                                   !
  !                                                                      !
  !  Purpose: Represents a hashtable of pair_t. Uses open addressing.    !
  !           Blank values are (0,0); deleted values are (0,1).          !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  type hashtable_t
     type(pair_t), dimension(:), allocatable :: table

     ! if table_size/size(table) > 50%, time to resize the hashtable
     integer :: table_size
     integer :: current_hash
  end type hashtable_t

  public :: init_pairs, reset_pairs, is_pair, get_pair, add_pair, del_pair
  private :: check_table

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: init_pairs                                              !
  !                                                                      !
  !  Purpose: hashtable_t constructor                                    !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine init_pairs(this)
    implicit none
    type(hashtable_t), intent(inout) :: this

    this%current_hash = 0
    if (.not.allocated(this%table)) allocate(this%table(0:10006))
    this%table(:)%ii = 0
    this%table(:)%jj = 0
    this%table_size = 0

  end subroutine init_pairs

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: reset_pairs                                             !
  !                                                                      !
  !  Purpose: Resets the iterator for traversing the hashtable.          !
  !           Used before calling get_pairs().                           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine reset_pairs(this)
    implicit none
    type(hashtable_t), intent(inout) :: this

    this%current_hash = 0

  end subroutine reset_pairs

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: check_table                                             !
  !                                                                      !
  !  Purpose: Debugging                                                  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  logical function check_table(this)
    implicit none
    type(hashtable_t), intent(in) :: this
    integer :: nn, blanks, deleted, full

    if (this%table_size > size(this%table)) then
       check_table = .false.
       return
    endif

    check_table = .true.
    ! return

    blanks = 0
    deleted = 0
    full = 0
    do nn=0, size(this%table)-1
       if (this%table(nn)%ii > 0 .and. this%table(nn)%jj > 0) then
          full = full + 1
       else if (this%table(nn)%ii .eq. 0 .and. this%table(nn)%jj .eq. 0) then
          blanks = blanks + 1
       else if (this%table(nn)%ii .eq. 0 .and. this%table(nn)%jj .eq. 1) then
          deleted = deleted + 1
       else
          print *,"SHOULD NEVER OCCUR"
          check_table = .false.
          return
       endif
    end do

    if (full .ne. this%table_size) then
       print *,"SIZE = ",size(this%table)
       print *,"blanks = ",blanks
       print *,"deleted = ",deleted
       print *,"full = ",full
       print *,"table_size = ",this%table_size
       check_table = .false.
    endif

    if (full+deleted+blanks .ne. size(this%table)) then
       print *,"SIZE = ",size(this%table)
       print *,"blanks = ",blanks
       print *,"deleted = ",deleted
       print *,"full = ",full
       print *,"table_size = ",this%table_size
       check_table = .false.
    endif

  end function check_table

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: get_pair                                                !
  !                                                                      !
  !  Purpose: Returns the next pair in the hashtable.  Calling           !
  !           reset_pairs() starts over again from the beginning.        !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine get_pair(this,pair)
    implicit none
    type(hashtable_t), intent(inout) :: this
    integer, intent(out) :: pair(2)

    pair(1) = 0
    pair(2) = 0

    do while (this%current_hash < size(this%table))
       if (0.ne.this%table(this%current_hash)%ii .and. 0.ne.this%table(this%current_hash)%jj) then
          pair(1) = this%table(this%current_hash)%ii
          pair(2) = this%table(this%current_hash)%jj
          this%current_hash = this%current_hash + 1
          return
       endif
       this%current_hash = this%current_hash + 1
    enddo
  end subroutine get_pair

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: is_pair                                                 !
  !                                                                      !
  !  Purpose: Returns .true. iff the pair (i0,j0) is in the hashtable.   !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  logical function is_pair(this,i0,j0)
    implicit none
    type(hashtable_t), intent(in) :: this
    integer, intent(in) :: i0, j0
    integer :: ii, jj, probe_count
    integer(kind=8) :: hash, init_hash

    ii = min(i0,j0)
    jj = max(i0,j0)

    if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       ERROR_STOP __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    probe_count = 1
    hash = mod(ii+jj*jj+probe_count*probe_count,size(this%table))
    if (hash < 0) hash = hash+size(this%table)
    init_hash = hash
    ! print *,"INIT HASH IS =",hash," TABLE IS ",this%table_size,"/",size(this%table)

    do
       if (this%table(hash)%ii .eq. ii .and. this%table(hash)%jj .eq. jj) then
          is_pair = .true.
          return
       endif
       if (this%table(hash)%ii .eq. 0 .and. this%table(hash)%jj .eq. 0) then
          is_pair = .false.
          return
       endif
       probe_count = probe_count + 1
       hash = mod(hash+probe_count*probe_count,size(this%table))
       if (hash < 0) hash = hash+size(this%table)
       if (hash .eq. init_hash) exit
    enddo

    print *,"loop in hash addressing, this should not occur"
    ERROR_STOP __LINE__

  end function is_pair

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: add_pair                                                !
  !                                                                      !
  !  Purpose: Adds (i0,j0) to the hashtable; resizes if necessary.       !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  recursive subroutine add_pair(this,i0,j0)
    implicit none
    type(hashtable_t), intent(inout) :: this
    integer, intent(in) :: i0,j0
    integer :: ii, jj, nn, old_size, old_tablesize, probe_count
    integer(kind=8) :: hash, init_hash
    type(pair_t), dimension(:), allocatable :: table_tmp

    if (i0 < 1 .or. j0 < 1) then
       print *,"invalid pair: ",i0,j0
       ERROR_STOP __LINE__
    endif

    if (size(this%table) < 2*this%table_size ) then
       old_size = size(this%table)
       old_tablesize = this%table_size
       allocate(table_tmp(0:old_size-1))
       if (size(table_tmp).ne.old_size) then
          print *,"size = ",size(table_tmp)
          print *,"old_size = ",old_size
          ERROR_STOP __LINE__
       endif
       table_tmp(0:old_size-1) = this%table(0:old_size-1)

       deallocate(this%table)
       allocate(this%table(0:2*old_size))
       this%table(:)%ii = 0
       this%table(:)%jj = 0
       this%table_size = 0
       do nn=0, old_size-1
          if ( table_tmp(nn)%ii .ne. 0 .and. table_tmp(nn)%jj .ne. 0) then
             call add_pair(this,table_tmp(nn)%ii,table_tmp(nn)%jj)
          endif
       enddo
       if (this%table_size.ne.old_tablesize) then
          print *,"size = ",this%table_size
          print *,"old_size = ",old_tablesize
          ERROR_STOP __LINE__
       endif
       deallocate(table_tmp)
    endif

    ii = min(i0,j0)
    jj = max(i0,j0)

    if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       ERROR_STOP __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    probe_count = 1
    hash = mod(ii+jj*jj+probe_count*probe_count,size(this%table))
    if (hash < 0) hash = hash+size(this%table)
    init_hash = hash

    do
       if (this%table(hash)%ii .eq. ii .and. this%table(hash)%jj .eq. jj) then
          ! already in table
          return
       endif
       if (this%table(hash)%ii .eq. 0 .or. this%table(hash)%jj .eq. 0) then
          this%table(hash)%ii = ii
          this%table(hash)%jj = jj
          this%table_size = this%table_size + 1
          ! if(.not. check_table(this)) ERROR_STOP __LINE__
          return
       endif
       probe_count = probe_count + 1
       hash = mod(hash+probe_count*probe_count,size(this%table))
       if (hash < 0) hash = hash+size(this%table)
       if (hash .eq. init_hash) exit
    enddo

    print *,"loop in hash addressing, this should not occur.  maybe hash table is full"
    ERROR_STOP __LINE__

  end subroutine add_pair

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: del_pair                                                !
  !                                                                      !
  !  Purpose: Removes (i0,j0) from the hashtable (if it exists).         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine del_pair(this,i0,j0)
    implicit none
    type(hashtable_t), intent(inout) :: this
    integer, intent(in) :: i0,j0
    integer :: ii, jj, probe_count
    integer(kind=8) :: hash, init_hash

    ii = min(i0,j0)
    jj = max(i0,j0)

    if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       ERROR_STOP __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    probe_count = 1
    hash = mod(ii+jj*jj+probe_count*probe_count,size(this%table))
    if (hash < 0) hash = hash+size(this%table)
    init_hash = hash

    do
       if (this%table(hash)%ii .eq. 0 .and. this%table(hash)%jj .eq. 0) then
          ! if(.not. check_table(this)) ERROR_STOP __LINE__
          return
       endif
       if (this%table(hash)%ii .eq. ii .and. this%table(hash)%jj .eq. jj) then
          ! 0,1 signifies DELETED hash entry
          this%table(hash)%ii = 0
          this%table(hash)%jj = 1
          this%table_size = this%table_size - 1
          ! if(.not. check_table(this)) ERROR_STOP __LINE__
          return
       endif
       probe_count = probe_count + 1
       hash = mod(hash+probe_count*probe_count,size(this%table))
       if (hash < 0) hash = hash+size(this%table)
       if (hash .eq. init_hash) exit
    enddo

    ! if(.not. check_table(this)) ERROR_STOP __LINE__

    print *,"loop in hash addressing. must be a lot of DELETED entries:  ",this%table_size,"/",size(this%table)

  end subroutine del_pair

end module pair_manager
