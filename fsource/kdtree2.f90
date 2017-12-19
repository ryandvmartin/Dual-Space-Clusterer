!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
module kdtree2_precision_module

  integer, parameter :: sp = kind(0.0)
  integer, parameter :: dp = kind(0.0d0)

  private :: sp, dp

  !
  ! You must comment out exactly one
  ! of the two lines.  If you comment
  ! out kdkind = sp then you get single precision
  ! and if you comment out kdkind = dp
  ! you get double precision.
  !

  !integer, parameter :: kdkind = sp
  integer, parameter :: kdkind = dp
  public :: kdkind

end module kdtree2_precision_module

module metric
    use kdtree2_precision_module, ONLY : kdkind
    implicit none

    !Rotation matrix
    real(kdkind), private :: rm2(2,2), rm3(3,3)

    !Vector between two point
    real(kdkind), private :: ba2(2), ba3(3)

    real(kdkind), private, pointer :: rmkd(:,:) => null() !rotation matrix stored in kdtree

    REAL*8, PRIVATE, PARAMETER :: PI=3.141592653589793238462643383279502884197D0
    REAL*8, PRIVATE, PARAMETER :: DEG2RAD = PI / 180.D0

    contains

    !set a pointer to a rotation matrix in a kdtree
    subroutine set_rmat_pointer ( rmp )
        real(kdkind), target, intent(in) :: rmp(:,:)
        rmkd => rmp
        return
    end subroutine set_rmat_pointer

    !Populate a rotation matrix based on angles (uncorrected!), and
    !ranges that describe a variogram or search ellipse.
    subroutine set_rmat ( angles, ranges )
        real*8, intent (in) :: angles( 3 )
        real*8, intent (in) :: ranges( 3 )
        real*8 :: sina, sinb, sint, cosa, cosb, cost
        real*8, dimension ( 3 ) :: a
        real*8 :: r1, r2
        real*8 :: rm(3,3)

        ! Fix the angles and get the required sines and cosines:
        a = angles
        call fix_angles ( a )
        sina = dsin(a(1)) ; sinb = dsin(a(2)) ; sint = dsin(a(3))
        cosa = dcos(a(1)) ; cosb = dcos(a(2)) ; cost = dcos(a(3))

        ! Construct the rotation matrix:
        r1 = ranges(1) / ranges(2)
        r2 = ranges(1) / ranges(3)
        rm(1,:) = (/ cosb * cosa , cosb * sina , -sinb /)
        rm(2,:) = (/ r1 * (-cost * sina + sint * sinb * cosa) , &
                     r1 * ( cost * cosa + sint * sinb * sina) , &
                     r1 * ( sint * cosb) /)
        rm(3,:) = (/ r2 * ( sint * sina + cost * sinb * cosa) , &
                     r2 * (-sint * cosa + cost * sinb * sina) , &
                     r2 * ( cost * cosb) /)

        rm2 = real(rm(1:2,1:2),kdkind)
        rm3 = real(rm,kdkind)

        Return
    end subroutine set_rmat

    !Convert three angles of anisotropy to those that are usable.  Note that
    !the vector of angles will be updated and returned
    subroutine fix_angles ( angles )
        real*8, intent (inout) :: angles( 3 )
        !1st - angle between the major axis of anisotropy and the
        !      E-W axis. Note: Counter clockwise is positive.
        !2nd - angle between major axis and the horizontal plane.
        !      (The dip of the ellipsoid measured positive down)
        !3rd - Angle of rotation of minor axis about the major axis
        !      of the ellipsoid.
        if(angles(1) >= 0.0D+00 .and. angles(1) < 270.0D+00) then
            angles(1) = ( 90.0D+00 - angles(1)) * DEG2RAD
        else
            angles(1) = (450.0D+00 - angles(1)) * DEG2RAD
        endif
        angles(2) = -1.0D+00 * angles(2) * DEG2RAD
        angles(3) =  1.0D+00 * angles(3) * DEG2RAD
        Return
    end subroutine fix_angles


    real(kdkind) function dsqrdp ( a, b )
        real(kdkind), intent (in) :: a(:), b(:)
        if(ubound(a,1) == 2)then
            if(associated(rmkd))then
                ba2 = b - a
                dsqrdp = (rmkd(1,1)*ba2(1) + rmkd(1,2)*ba2(2))**2 + &
                         (rmkd(2,1)*ba2(1) + rmkd(2,2)*ba2(2))**2
            else
                ba2 = b - a
                dsqrdp = ba2(1)**2 + ba2(2)**2
            endif
        else
            if(associated(rmkd))then
                ba3 = b - a
                dsqrdp = (rmkd(1,1)*ba3(1) + rmkd(1,2)*ba3(2) + rmkd(1,3)*ba3(3))**2 + &
                         (rmkd(2,1)*ba3(1) + rmkd(2,2)*ba3(2) + rmkd(2,3)*ba3(3))**2 + &
                         (rmkd(3,1)*ba3(1) + rmkd(3,2)*ba3(2) + rmkd(3,3)*ba3(3))**2
            else
                ba3 = b - a
                dsqrdp = ba3(1)**2 + ba3(2)**2 + ba3(3)**2
            endif
        endif
    end function dsqrdp

    !Anisotropic distance along single coordinate
    real(kdkind) function dsqrdpi ( a, b, i )
        real(kdkind), intent (in) :: a, b
        integer, intent (in) :: i

        if(associated(rmkd))then
            ba3 = 0
            ba3(i) = b - a
            dsqrdpi = (rmkd(1,1)*ba3(1) + rmkd(1,2)*ba3(2) + rmkd(1,3)*ba3(3))**2 + &
                        (rmkd(2,1)*ba3(1) + rmkd(2,2)*ba3(2) + rmkd(2,3)*ba3(3))**2 + &
                        (rmkd(3,1)*ba3(1) + rmkd(3,2)*ba3(2) + rmkd(3,3)*ba3(3))**2
        else
            ba3(i) = b - a
            dsqrdpi = ba3(i)**2
        endif

    end function dsqrdpi


    !
    ! Euclidean distance metrics
    !

    !Get the distance between two points in 2D
    real(kdkind) function dsqrd2 ( a, b )
        real(kdkind), intent (in) :: a(2), b(2)
        ba2 = b - a
        dsqrd2 = ba2(1)**2 + ba2(2)**2
    end function dsqrd2

    !Get the distance between two points in 3D
    real(kdkind) function dsqrd3 ( a, b )
        real(kdkind), intent (in) :: a(3), b(3)
        ba3 = b - a
        dsqrd3 = ba3(1)**2 + ba3(2)**2 + ba3(3)**2
    end function dsqrd3

    !
    ! Anisotropic distance metrics
    !

    !Get the anisotropic distance between two points in 2D
    real(kdkind) function dsqrd2a ( a, b )
        real(kdkind), intent (in) :: a(2), b(2)
        ba2 = b - a
        dsqrd2a = (rm2(1,1)*ba2(1) + rm2(1,2)*ba2(2))**2 + &
                  (rm2(2,1)*ba2(1) + rm2(2,2)*ba2(2))**2
    end function dsqrd2a

    !Get the anisotropic distance between two points in 3D
    real(kdkind) function dsqrd3a ( a, b )
        real(kdkind), intent (in) :: a(3), b(3)
        ba3 = b - a
        dsqrd3a = (rm3(1,1)*ba3(1) + rm3(1,2)*ba3(2) + rm3(1,3)*ba3(3))**2 + &
                  (rm3(2,1)*ba3(1) + rm3(2,2)*ba3(2) + rm3(2,3)*ba3(3))**2 + &
                  (rm3(3,1)*ba3(1) + rm3(3,2)*ba3(2) + rm3(3,3)*ba3(3))**2
    end function dsqrd3a

end module metric

module kdtree2_priority_queue_module
  use kdtree2_precision_module
  !
  ! maintain a priority queue (PQ) of data, pairs of 'priority/payload',
  ! implemented with a binary heap.  This is the type, and the 'dis' field
  ! is the priority.
  !
  type kdtree2_result
      real(kdkind) :: dis!=0.0   !distance
      integer      :: idx!=-1    !index
  end type kdtree2_result
  !
  ! A heap-based priority queue lets one efficiently implement the following
  ! operations, each in log(N) time, as opposed to linear time.
  !
  ! 1)  add a datum (push a datum onto the queue, increasing its length)
  ! 2)  return the priority value of the maximum priority element
  ! 3)  pop-off (and delete) the element with the maximum priority, decreasing
  !     the size of the queue.
  ! 4)  replace the datum with the maximum priority with a supplied datum
  !     (of either higher or lower priority), maintaining the size of the
  !     queue.
  !
  !
  ! In the k-d tree case, the 'priority' is the square distance of a point in
  ! the data set to a reference point.   The goal is to keep the smallest M
  ! distances to a reference point.  The tree algorithm searches terminal
  ! nodes to decide whether to add points under consideration.
  !
  ! A priority queue is useful here because it lets one quickly return the
  ! largest distance currently existing in the list.  If a new candidate
  ! distance is smaller than this, then the new candidate ought to replace
  ! the old candidate.  In priority queue terms, this means removing the
  ! highest priority element, and inserting the new one.
  !
  ! Algorithms based on Cormen, Leiserson, Rivest, _Introduction
  ! to Algorithms_, 1990, with further optimization by the author.
  !
  ! Originally informed by a C implementation by Sriranga Veeraraghavan.
  !
  ! This module is not written in the most clear way, but is implemented such
  ! for speed, as it its operations will be called many times during searches
  ! of large numbers of neighbors.
  !
  type pq
      !
      ! The priority queue consists of elements
      ! priority(1:heap_size), with associated payload(:).
      !
      ! There are heap_size active elements.
      ! Assumes the allocation is always sufficient.  Will NOT increase it
      ! to match.
      integer :: heap_size = 0
      type(kdtree2_result), pointer :: elems(:)
  end type pq

  public :: kdtree2_result

  public :: pq
  public :: pq_create
  public :: pq_delete, pq_insert
  public :: pq_extract_max, pq_max, pq_replace_max, pq_maxpri
  private

contains


  function pq_create(results_in) result(res)
    !
    ! Create a priority queue from ALREADY allocated
    ! array pointers for storage.  NOTE! It will NOT
    ! add any alements to the heap, i.e. any existing
    ! data in the input arrays will NOT be used and may
    ! be overwritten.
    !
    ! usage:
    !    real(kdkind), pointer :: x(:)
    !    integer, pointer :: k(:)
    !    allocate(x(1000),k(1000))
    !    pq => pq_create(x,k)
    !
    type(kdtree2_result), target:: results_in(:)
    type(pq) :: res
    !
    !
    integer :: nalloc

    nalloc = size(results_in,1)
    if (nalloc .lt. 1) then
       write (*,*) 'PQ_CREATE: error, input arrays must be allocated.'
    end if
    res%elems => results_in
    res%heap_size = 0
    return
  end function pq_create

  !
  ! operations for getting parents and left + right children
  ! of elements in a binary heap.
  !

!
! These are written inline for speed.
!
!  integer function parent(i)
!    integer, intent(in) :: i
!    parent = (i/2)
!    return
!  end function parent

!  integer function left(i)
!    integer, intent(in) ::i
!    left = (2*i)
!    return
!  end function left

!  integer function right(i)
!    integer, intent(in) :: i
!    right = (2*i)+1
!    return
!  end function right

!  logical function compare_priority(p1,p2)
!    real(kdkind), intent(in) :: p1, p2
!
!    compare_priority = (p1 .gt. p2)
!    return
!  end function compare_priority

  subroutine heapify(a,i_in)
    !
    ! take a heap rooted at 'i' and force it to be in the
    ! heap canonical form.   This is performance critical
    ! and has been tweaked a little to reflect this.
    !
    type(pq),pointer   :: a
    integer, intent(in) :: i_in
    !
    integer :: i, l, r, largest

    real(kdkind)    :: pri_i, pri_l, pri_r, pri_largest


    type(kdtree2_result) :: temp

    i = i_in

bigloop:  do
       l = 2*i ! left(i)
       r = l+1 ! right(i)
       !
       ! set 'largest' to the index of either i, l, r
       ! depending on whose priority is largest.
       !
       ! note that l or r can be larger than the heap size
       ! in which case they do not count.


       ! does left child have higher priority?
       if (l .gt. a%heap_size) then
          ! we know that i is the largest as both l and r are invalid.
          exit
       else
          pri_i = a%elems(i)%dis
          pri_l = a%elems(l)%dis
          if (pri_l .gt. pri_i) then
             largest = l
             pri_largest = pri_l
          else
             largest = i
             pri_largest = pri_i
          endif

          !
          ! between i and l we have a winner
          ! now choose between that and r.
          !
          if (r .le. a%heap_size) then
             pri_r = a%elems(r)%dis
             if (pri_r .gt. pri_largest) then
                largest = r
             endif
          endif
       endif

       if (largest .ne. i) then
          ! swap data in nodes largest and i, then heapify

          temp = a%elems(i)
          a%elems(i) = a%elems(largest)
          a%elems(largest) = temp
          !
          ! Canonical heapify() algorithm has tail-ecursive call:
          !
          !        call heapify(a,largest)
          ! we will simulate with cycle
          !
          i = largest
          cycle bigloop ! continue the loop
       else
          return   ! break from the loop
       end if
    enddo bigloop
    return
  end subroutine heapify

  subroutine pq_max(a,e)
    !
    ! return the priority and its payload of the maximum priority element
    ! on the queue, which should be the first one, if it is
    ! in heapified form.
    !
    type(pq),pointer :: a
    type(kdtree2_result),intent(out)  :: e

    if (a%heap_size .gt. 0) then
       e = a%elems(1)
    else
       write (*,*) 'PQ_MAX: ERROR, heap_size < 1'
       stop
    endif
    return
  end subroutine pq_max

  real(kdkind) function pq_maxpri(a)
    type(pq), pointer :: a

    if (a%heap_size .gt. 0) then
       pq_maxpri = a%elems(1)%dis
    else
       write (*,*) 'PQ_MAX_PRI: ERROR, heapsize < 1'
       stop
    endif
    return
  end function pq_maxpri

  subroutine pq_extract_max(a,e)
    !
    ! return the priority and payload of maximum priority
    ! element, and remove it from the queue.
    ! (equivalent to 'pop()' on a stack)
    !
    type(pq),pointer :: a
    type(kdtree2_result), intent(out) :: e

    if (a%heap_size .ge. 1) then
       !
       ! return max as first element
       !
       e = a%elems(1)

       !
       ! move last element to first
       !
       a%elems(1) = a%elems(a%heap_size)
       a%heap_size = a%heap_size-1
       call heapify(a,1)
       return
    else
       write (*,*) 'PQ_EXTRACT_MAX: error, attempted to pop non-positive PQ'
       stop
    end if

  end subroutine pq_extract_max


  real(kdkind) function pq_insert(a,dis,idx)
    !
    ! Insert a new element and return the new maximum priority,
    ! which may or may not be the same as the old maximum priority.
    !
    type(pq),pointer  :: a
    real(kdkind), intent(in) :: dis
    integer, intent(in) :: idx
    !    type(kdtree2_result), intent(in) :: e
    !
    integer :: i, isparent
    real(kdkind)    :: parentdis
    !

    !    if (a%heap_size .ge. a%max_elems) then
    !       write (*,*) 'PQ_INSERT: error, attempt made to insert element on full PQ'
    !       stop
    !    else
    a%heap_size = a%heap_size + 1
    i = a%heap_size

    do while (i .gt. 1)
       isparent = int(i/2)
       parentdis = a%elems(isparent)%dis
       if (dis .gt. parentdis) then
          ! move what was in i's parent into i.
          a%elems(i)%dis = parentdis
          a%elems(i)%idx = a%elems(isparent)%idx
          i = isparent
       else
          exit
       endif
    end do

    ! insert the element at the determined position
    a%elems(i)%dis = dis
    a%elems(i)%idx = idx

    pq_insert = a%elems(1)%dis
    return
    !    end if

  end function pq_insert

  subroutine pq_adjust_heap(a,i)
    type(pq),pointer  :: a
    integer, intent(in) :: i
    !
    ! nominally arguments (a,i), but specialize for a=1
    !
    ! This routine assumes that the trees with roots 2 and 3 are already heaps, i.e.
    ! the children of '1' are heaps.  When the procedure is completed, the
    ! tree rooted at 1 is a heap.
    real(kdkind) :: prichild
    integer :: parent, child, N

    type(kdtree2_result) :: e

    e = a%elems(i)

    parent = i
    child = 2*i
    N = a%heap_size

    do while (child .le. N)
       if (child .lt. N) then
          if (a%elems(child)%dis .lt. a%elems(child+1)%dis) then
             child = child+1
          endif
       endif
       prichild = a%elems(child)%dis
       if (e%dis .ge. prichild) then
          exit
       else
          ! move child into parent.
          a%elems(parent) = a%elems(child)
          parent = child
          child = 2*parent
       end if
    end do
    a%elems(parent) = e
    return
  end subroutine pq_adjust_heap


  real(kdkind) function pq_replace_max(a,dis,idx)
    !
    ! Replace the extant maximum priority element
    ! in the PQ with (dis,idx).  Return
    ! the new maximum priority, which may be larger
    ! or smaller than the old one.
    !
    type(pq),pointer         :: a
    real(kdkind), intent(in) :: dis
    integer, intent(in) :: idx
!    type(kdtree2_result), intent(in) :: e
    ! not tested as well!

    integer :: parent, child, N
    real(kdkind)    :: prichild, prichildp1

    type(kdtree2_result) :: etmp

    if (.true.) then
       N=a%heap_size
       if (N .ge. 1) then
          parent =1
          child=2

          loop: do while (child .le. N)
             prichild = a%elems(child)%dis

             !
             ! posibly child+1 has higher priority, and if
             ! so, get it, and increment child.
             !

             if (child .lt. N) then
                prichildp1 = a%elems(child+1)%dis
                if (prichild .lt. prichildp1) then
                   child = child+1
                   prichild = prichildp1
                endif
             endif

             if (dis .ge. prichild) then
                exit loop
                ! we have a proper place for our new element,
                ! bigger than either children's priority.
             else
                ! move child into parent.
                a%elems(parent) = a%elems(child)
                parent = child
                child = 2*parent
             end if
          end do loop
          a%elems(parent)%dis = dis
          a%elems(parent)%idx = idx
          pq_replace_max = a%elems(1)%dis
       else
          a%elems(1)%dis = dis
          a%elems(1)%idx = idx
          pq_replace_max = dis
       endif
    else
       !
       ! slower version using elementary pop and push operations.
       !
       call pq_extract_max(a,etmp)
       etmp%dis = dis
       etmp%idx = idx
       pq_replace_max = pq_insert(a,dis,idx)
    endif
    return
  end function pq_replace_max

  subroutine pq_delete(a,i)
    !
    ! delete item with index 'i'
    !
    type(pq),pointer :: a
    integer           :: i

    if ((i .lt. 1) .or. (i .gt. a%heap_size)) then
       write (*,*) 'PQ_DELETE: error, attempt to remove out of bounds element.'
       stop
    endif

    ! swap the item to be deleted with the last element
    ! and shorten heap by one.
    a%elems(i) = a%elems(a%heap_size)
    a%heap_size = a%heap_size - 1

    call heapify(a,i)

  end subroutine pq_delete

end module kdtree2_priority_queue_module


module kdtree2_module
  use metric
  use kdtree2_precision_module
  use kdtree2_priority_queue_module
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported.
  !
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module.
  !
  !-------------DATA TYPE, CREATION, DELETION---------------------
  public :: kdkind
  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create, kdtree2_destroy

  public :: kdtree2_init !initialize a kd tree for dynamic addition of points.  This
                         !will assume the total number of points is known in advance
  public :: kdtree2_insert_point !insert a new point into a kdtree, only split when the number
                                 !exceeds the bucket_size
  public :: kdtree2_reset !reset everything so a tree can get rebuild

  public :: kdtree2_ortho_nearest !orthogonal range query

  !---------------------------------------------------------------
  !-------------------SEARCH ROUTINES-----------------------------
  public :: kdtree2_n_nearest,kdtree2_n_nearest_around_point
  ! Return fixed number of nearest neighbors around arbitrary vector,
  ! or extant point in dataset, with decorrelation window.
  !
  public :: kdtree2_r_nearest, kdtree2_r_nearest_around_point
  ! Return points within a fixed ball of arb vector/extant point
  !
  public :: kdtree2_sort_results
  ! Sort, in order of increasing distance, rseults from above.
  !
  public :: kdtree2_r_count, kdtree2_r_count_around_point
  ! Count points within a fixed ball of arb vector/extant point
  !
  public :: kdtree2_n_nearest_brute_force, kdtree2_r_nearest_brute_force
  ! brute force of kdtree2_[n|r]_nearest
  !----------------------------------------------------------------

  interface dmet
    module procedure dsqrdp, dsqrdpi
  end interface dmet

  integer, parameter :: bucket_size = 9
  ! The maximum number of points to keep in a terminal node.

  type interval
      real(kdkind) :: lower,upper
  end type interval

  type :: tree_node
      ! an internal tree node
      private
      integer :: cut_dim
      ! the dimension to cut
      real(kdkind) :: cut_val
      ! where to cut the dimension
      real(kdkind) :: cut_val_left, cut_val_right
      ! improved cutoffs knowing the spread in child boxes.
      integer :: l, u
      type (tree_node), pointer :: left, right
      type(interval), pointer :: box(:) => null()
      ! child pointers
      ! Points included in this node are indexes[k] with k \in [l,u]

      integer, pointer :: ind(:) !indexes of data this node contains

  end type tree_node

  type :: kdtree2
      ! Global information about the tree, one per tree
      integer :: dimen=0, n=0
      ! dimensionality and total # of points
      real(kdkind), pointer :: the_data(:,:) => null()
      ! pointer to the actual data array
      !
      !  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
      !  which may be opposite of what may be conventional.
      !  This is, because in Fortran, the memory layout is such that
      !  the first dimension is in sequential order.  Hence, with
      !  (1:d,1:N), all components of the vector will be in consecutive
      !  memory locations.  The search time is dominated by the
      !  evaluation of distances in the terminal nodes.  Putting all
      !  vector components in consecutive memory location improves
      !  memory cache locality, and hence search speed, and may enable
      !  vectorization on some processors and compilers.

      integer, pointer :: ind(:) => null()
      ! permuted index into the data, so that indexes[l..u] of some
      ! bucket represent the indexes of the actual points in that
      ! bucket.
      logical       :: sort = .false.
      ! do we always sort output results?
      logical       :: rearrange = .false.
      real(kdkind), pointer :: rearranged_data(:,:) => null()
      ! if (rearrange .eqv. .true.) then rearranged_data has been
      ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
      ! permitting search to use more cache-friendly rearranged_data, at
      ! some initial computation and storage cost.
      type (tree_node), pointer :: root => null()
      ! root pointer of the tree

      real(kdkind) :: aniso(3,3)

  end type kdtree2


  type :: tree_search_record
      !
      ! One of these is created for each search.
      !
      private
      !
      ! Many fields are copied from the tree structure, in order to
      ! speed up the search.
      !
      integer           :: dimen
      integer           :: nn, nfound
      real(kdkind)      :: ballsize
      integer           :: centeridx=999, correltime=9999
      ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
      integer           :: nalloc  ! how much allocated for results(:)?
      logical           :: rearrange  ! are the data rearranged or original?
      ! did the # of points found overflow the storage provided?
      logical           :: overflow
      real(kdkind), pointer :: qv(:)  ! query vector
      real(kdkind), pointer :: rmin(:), rmax(:) !othogonal range query data
      type(kdtree2_result), pointer :: results(:) ! results
      type(pq) :: pq
      real(kdkind), pointer :: data(:,:)  ! temp pointer to data
      integer, pointer      :: ind(:)     ! temp pointer to indexes
  end type tree_search_record

  integer, public :: KDTREE2_NFOUND

  private
  ! everything else is private.

  type(tree_search_record), save, target :: sr   ! A GLOBAL VARIABLE for search

contains

  function kdtree2_create(input_data,dim,sort,rearrange,aniso) result (mr)
    !
    ! create the actual tree structure, given an input array of data.
    !
    ! Note, input data is input_data(1:d,1:N), NOT the other way around.
    ! THIS IS THE REVERSE OF THE PREVIOUS VERSION OF THIS MODULE.
    ! The reason for it is cache friendliness, improving performance.
    !
    ! Optional arguments:  If 'dim' is specified, then the tree
    !                      will only search the first 'dim' components
    !                      of input_data, otherwise, dim is inferred
    !                      from SIZE(input_data,1).
    !
    !                      if sort .eqv. .true. then output results
    !                      will be sorted by increasing distance.
    !                      default=.false., as it is faster to not sort.
    !
    !                      if rearrange .eqv. .true. then an internal
    !                      copy of the data, rearranged by terminal node,
    !                      will be made for cache friendliness.
    !                      default=.true., as it speeds searches, but
    !                      building takes longer, and extra memory is used.
    !
    ! .. Function Return Cut_value ..
    type (kdtree2), pointer :: mr
    integer, intent(in), optional      :: dim
    logical, intent(in), optional      :: sort
    logical, intent(in), optional      :: rearrange
    real(kdkind), intent(in), optional :: aniso(3,3)
    ! ..
    ! .. Array Arguments ..
    real(kdkind), target :: input_data(:,:)
    real(kdkind) :: rm(3,3)

    !
    integer :: i
    ! ..

    allocate (mr)
    mr%the_data => input_data
    ! pointer assignment

    if (present(dim)) then
       mr%dimen = dim
    else
       mr%dimen = size(input_data,1)
    end if
    mr%n = size(input_data,2)

    if (mr%dimen > mr%n) then
       !  unlikely to be correct
       return
       write (*,*) 'KD_TREE_TRANS: likely user error.'
       write (*,*) 'KD_TREE_TRANS: You passed in matrix with D=',mr%dimen
       write (*,*) 'KD_TREE_TRANS: and N=',mr%n
       write (*,*) 'KD_TREE_TRANS: note, that new format is data(1:D,1:N)'
       write (*,*) 'KD_TREE_TRANS: with usually N >> D.   If N =approx= D, then a k-d tree'
       write (*,*) 'KD_TREE_TRANS: is not an appropriate data structure.'
       stop
    end if

    if(present(aniso))then
        rm = aniso
    else
        rm = reshape((/1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0/),(/3,3/))
    endif
    mr%aniso = rm
    call set_rmat_pointer ( mr%aniso )

    call build_tree(mr)

    if (present(sort)) then
       mr%sort = sort
    else
       mr%sort = .false.
    endif

    if (present(rearrange)) then
       mr%rearrange = rearrange
    else
       mr%rearrange = .true.
    endif

    if (mr%rearrange) then
       allocate(mr%rearranged_data(mr%dimen,mr%n))
       do i=1,mr%n
          mr%rearranged_data(:,i) = mr%the_data(:,mr%ind(i))
       enddo
    else
       nullify(mr%rearranged_data)
    endif

  end function kdtree2_create

!
! New function added by John Manchuk, Dec 2009
!
    function kdtree2_init(input_data,dim,sort,rearrange,aniso) result (mr)
    !
    ! Initialize a kdtree with a known amount of data, but unknown insertion order
    ! input_data is assumed the correct size (dim by n)
    !
        type (kdtree2), pointer :: mr
        integer, intent(in), optional      :: dim
        logical, intent(in), optional      :: sort
        logical, intent(in), optional      :: rearrange
        real(kdkind), intent(in), optional :: aniso(3,3)
        ! ..
        ! .. Array Arguments ..
        real(kdkind), target :: input_data(:,:)
        real(kdkind) :: rm(3,3)

        !
        integer :: i
        ! ..

        allocate (mr)
        mr%the_data => input_data
        ! pointer assignment

        if (present(dim)) then
           mr%dimen = dim
        else
           mr%dimen = size(input_data,1)
        end if
        mr%n = size(input_data,2)

        if (mr%dimen > mr%n) then
           !  unlikely to be correct
           write (*,*) 'KD_TREE_TRANS: likely user error.'
           write (*,*) 'KD_TREE_TRANS: You passed in matrix with D=',mr%dimen
           write (*,*) 'KD_TREE_TRANS: and N=',mr%n
           write (*,*) 'KD_TREE_TRANS: note, that new format is data(1:D,1:N)'
           write (*,*) 'KD_TREE_TRANS: with usually N >> D.   If N =approx= D, then a k-d tree'
           write (*,*) 'KD_TREE_TRANS: is not an appropriate data structure.'
           stop
        end if

        if(present(aniso))then
            rm = aniso
        else
            rm = reshape((/1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0/),(/3,3/))
        endif
        mr%aniso = rm
        call set_rmat_pointer ( mr%aniso )


        allocate (mr%ind(mr%n))
        forall (j = 1:mr%n)
            mr%ind(j) = j
        end forall

        if (present(sort)) then
            mr%sort = sort
        else
            mr%sort = .false.
        endif

        mr%rearrange = .false.
        nullify(mr%rearranged_data)

    end function kdtree2_init

    recursive subroutine kdtree2_insert_point( tp, node, p, id )
    !
    ! Add a point to a kdtree node.
    !
    ! This routine first finds the containing node, but another implementation
    ! could use the node from a prior query instead of performing this twice
    !
        type (kdtree2), pointer :: tp
        type(tree_node), pointer :: node !the current node of the kdtree

        real(kdkind), intent (in) :: p(:) !the point being inserted
        integer, intent (in) :: id !index of the point being inserted: tp%ind(k) = id
        type(tree_node), pointer :: thisnode, parent
        !
        integer :: c, i, l, u, m
        integer :: dimen
        ! ..
        real(kdkind) :: qval, average
        integer, pointer :: ind(:)

        dimen = tp%dimen

        if(.not.associated(node))then
            allocate(node)
            node%cut_dim = 0
            node%cut_val = 0.0
            node%l = 1
            node%u = 0
            node%left => null()
            node%right => null()
            allocate(node%box(dimen))
            node%box(1:dimen)%lower = huge(real(1.0,kdkind))
            node%box(1:dimen)%upper = tiny(real(1.0,kdkind))
            allocate(node%ind(bucket_size + 1))
        endif

        thisnode => node

        if((associated(thisnode%left) .and. &
            associated(thisnode%right)) .eqv. .false.) then ! we are on a terminal node

            ind => thisnode%ind(1:)

            !Add a datum to this node
            thisnode%u = thisnode%u + 1
            thisnode%ind(thisnode%u) = id

            !Expand the bounding box as data are added
            do i = 1, dimen
                if(thisnode%box(i)%lower > tp%the_data(i,id))then
                    thisnode%box(i)%lower = tp%the_data(i,id)
                endif
                if(thisnode%box(i)%upper < tp%the_data(i,id))then
                    thisnode%box(i)%upper = tp%the_data(i,id)
                endif
            enddo

            !Did the insertion of this point exceed the bucketsize
            if(thisnode%u - thisnode%l + 1 > bucket_size)then

                l = thisnode%l
                u = thisnode%u

                !Determine index to split on
                c = maxloc(thisnode%box(1:dimen)%upper - &
                           thisnode%box(1:dimen)%lower,1)

                thisnode%cut_dim = c

                !Compute the average along dimension c
                !this is the location of the split
                average = sum(tp%the_data(c,thisnode%ind(l:u))) / real(u-l+1,kdkind)

                thisnode%cut_val = average

                !Organize the data into left and right
                m = select_on_coordinate_value(tp%the_data,thisnode%ind,c,average,l,u)

                !Now build the new left and right nodes
                allocate(thisnode%left)
                thisnode%left%l = 1
                thisnode%left%u = m
                thisnode%left%cut_dim = 0
                thisnode%left%cut_val = 0.0
                thisnode%left%left => null()
                thisnode%left%right => null()
                allocate(thisnode%left%ind(bucket_size + 1))
                thisnode%left%ind(thisnode%left%l:thisnode%left%u) = thisnode%ind(1:m)
                allocate(thisnode%left%box(dimen))
                do i = 1, dimen !Update the bounding box of this node
                    call spread_in_coordinate_node(tp,thisnode%left,i,1,m,thisnode%left%box(i))
                enddo

                allocate(thisnode%right)
                thisnode%right%l = 1
                thisnode%right%u = thisnode%u - m
                thisnode%right%cut_dim = 0
                thisnode%right%cut_val = 0.0
                thisnode%right%left => null()
                thisnode%right%right => null()
                allocate(thisnode%right%ind(bucket_size + 1))
                thisnode%right%ind(thisnode%right%l:thisnode%right%u) = thisnode%ind(m+1:u)
                allocate(thisnode%right%box(dimen))
                do i = 1, dimen !Update the bounding box of this node
                    call spread_in_coordinate_node(tp,thisnode%right,i,1,thisnode%u - m,thisnode%right%box(i))
                enddo

                !the indexes of the parent node are no longer required
                deallocate(thisnode%ind)

            endif

        else ! we are not on a terminal node

            qval = p(thisnode%cut_dim)

            if (qval < thisnode%cut_val) then
                parent => thisnode
                thisnode => parent%left
            else
                parent => thisnode
                thisnode => parent%right
            endif

            !Descend further into the kdtree
            call kdtree2_insert_point ( tp, thisnode, p, id )

            !Was the node actually split?
            ! If so, update the cut_vals
            ! otherwise, update it's parents cut_vals
            if(associated(thisnode%right) .and. associated(thisnode%left))then
                !Expand the bounding box as the tree is traversed
!                do i = 1, dimen
!                    if(thisnode%box(i)%lower > tp%the_data(i,id))then
!                       thisnode%box(i)%lower = tp%the_data(i,id)
!                    endif
!                    if(thisnode%box(i)%upper < tp%the_data(i,id))then
!                       thisnode%box(i)%upper = tp%the_data(i,id)
!                    endif
!                enddo

                !Update the cut vals
                c = thisnode%cut_dim
                thisnode%cut_val_right = thisnode%right%box(c)%lower
                thisnode%cut_val_left = thisnode%left%box(c)%upper
                !parent%cut_val = (parent%cut_val_left + parent%cut_val_right) / real(2,kdkind)
            endif

            !Expand the bounding box as the tree is traversed
            do i = 1, dimen
                if(parent%box(i)%lower > tp%the_data(i,id))then
                   parent%box(i)%lower = tp%the_data(i,id)
                endif
                if(parent%box(i)%upper < tp%the_data(i,id))then
                   parent%box(i)%upper = tp%the_data(i,id)
                endif
            enddo

            !Update the cut vals
            c = parent%cut_dim
            parent%cut_val_right = parent%right%box(c)%lower
            parent%cut_val_left = parent%left%box(c)%upper
            !parent%cut_val = (parent%cut_val_left + parent%cut_

        endif

        return
    end subroutine kdtree2_insert_point



    subroutine spread_in_coordinate_node(tp,node,c,l,u,interv)
        ! the spread in coordinate 'c', between l and u.
        ! This is identical to the original version, but uses ind array from the node
        ! Return lower bound in 'smin', and upper in 'smax',
        ! ..
        ! .. Structure Arguments ..
        type (kdtree2), pointer :: tp
        type(tree_node), pointer :: node
        type(interval), intent(out) :: interv
        ! ..
        ! .. Scalar Arguments ..
        integer, intent (In) :: c, l, u
        ! ..
        ! .. Local Scalars ..
        real(kdkind) :: last, lmax, lmin, t, smin,smax
        integer :: i, ulocal
        ! ..
        ! .. Local Arrays ..
        real(kdkind), pointer :: v(:,:)
        integer, pointer :: ind(:)
        ! ..
        v => tp%the_data(1:,1:)

        ind => node%ind(1:)

        smin = v(c,ind(l))
        smax = smin

        ulocal = u

        do i = l + 2, ulocal, 2
            lmin = v(c,ind(i-1))
            lmax = v(c,ind(i))
            if (lmin>lmax) then
                t = lmin
                lmin = lmax
                lmax = t
            end if
            if (smin>lmin) smin = lmin
            if (smax<lmax) smax = lmax
        end do
        if (i==ulocal+1) then
            last = v(c,ind(ulocal))
            if (smin>last) smin = last
            if (smax<last) smax = last
        end if

        interv%lower = smin
        interv%upper = smax

    end subroutine spread_in_coordinate_node

    subroutine kdtree2_reset(tp)
        ! Deallocates all memory for the tree, except initialized components
        type (kdtree2), pointer :: tp
        ! ..
        call destroy_node(tp%root)

        return

        contains
            recursive subroutine destroy_node(np)
                ! .. Structure Arguments ..
                type (tree_node), pointer :: np
                ! ..
                ! .. Intrinsic Functions ..
                intrinsic ASSOCIATED
                ! ..
                if (associated(np%left)) then
                    call destroy_node(np%left)
                    nullify (np%left)
                end if
                if (associated(np%right)) then
                    call destroy_node(np%right)
                    nullify (np%right)
                end if
                if (associated(np%box)) deallocate(np%box)
                if (associated(np%ind)) deallocate(np%ind)
                deallocate(np)
                return
            end subroutine destroy_node
    end subroutine kdtree2_reset


  !
  ! End of new functions
  !

    subroutine build_tree(tp)
      type (kdtree2), pointer :: tp
      ! ..
      integer :: j
      type(tree_node), pointer :: dummy => null()
      ! ..
      allocate (tp%ind(tp%n))
      forall (j=1:tp%n)
         tp%ind(j) = j
      end forall
      tp%root => build_tree_for_range(tp,1,tp%n, dummy)
    end subroutine build_tree

    recursive function build_tree_for_range(tp,l,u,parent) result (res)
      ! .. Function Return Cut_value ..
      type (tree_node), pointer :: res
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type (tree_node),pointer :: parent
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      integer :: i, c, m, dimen
      logical :: recompute
      real(kdkind)    :: average

!!$      If (.False.) Then
!!$         If ((l .Lt. 1) .Or. (l .Gt. tp%n)) Then
!!$            Stop 'illegal L value in build_tree_for_range'
!!$         End If
!!$         If ((u .Lt. 1) .Or. (u .Gt. tp%n)) Then
!!$            Stop 'illegal u value in build_tree_for_range'
!!$         End If
!!$         If (u .Lt. l) Then
!!$            Stop 'U is less than L, thats illegal.'
!!$         End If
!!$      Endif
!!$
      ! first compute min and max
      dimen = tp%dimen
      allocate (res)
      allocate(res%box(dimen))

      ! First, compute an APPROXIMATE bounding box of all points associated with this node.
      if ( u < l ) then
         ! no points in this box
         nullify(res)
         return
      end if

      if ((u-l)<=bucket_size) then
         !
         ! always compute true bounding box for terminal nodes.
         !
         do i=1,dimen
            call spread_in_coordinate(tp,i,l,u,res%box(i))
         end do
         res%cut_dim = 0
         res%cut_val = 0.0
         res%l = l
         res%u = u
         res%left =>null()
         res%right => null()
      else
         !
         ! modify approximate bounding box.  This will be an
         ! overestimate of the true bounding box, as we are only recomputing
         ! the bounding box for the dimension that the parent split on.
         !
         ! Going to a true bounding box computation would significantly
         ! increase the time necessary to build the tree, and usually
         ! has only a very small difference.  This box is not used
         ! for searching but only for deciding which coordinate to split on.
         !
         do i=1,dimen
            recompute=.true.
            if (associated(parent)) then
               if (i .ne. parent%cut_dim) then
                  recompute=.false.
               end if
            endif
            if (recompute) then
               call spread_in_coordinate(tp,i,l,u,res%box(i))
            else
               res%box(i) = parent%box(i)
            endif
         end do


         c = maxloc(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
         !
         ! c is the identity of which coordinate has the greatest spread.
         !

         if (.false.) then
            ! select exact median to have fully balanced tree.
            m = (l+u)/2
            call select_on_coordinate(tp%the_data,tp%ind,c,m,l,u)
         else
            !
            ! select point halfway between min and max, as per A. Moore,
            ! who says this helps in some degenerate cases, or
            ! actual arithmetic average.
            !
            if (.true.) then
               ! actually compute average
               average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1,kdkind)
            else
               average = (res%box(c)%upper + res%box(c)%lower)/2.0
            endif

            res%cut_val = average
            m = select_on_coordinate_value(tp%the_data,tp%ind,c,average,l,u)
         endif

         ! moves indexes around
         res%cut_dim = c
         res%l = l
         res%u = u
!         res%cut_val = tp%the_data(c,tp%ind(m))

         res%left => build_tree_for_range(tp,l,m,res)
         res%right => build_tree_for_range(tp,m+1,u,res)

         if (associated(res%right) .eqv. .false.) then
            res%box = res%left%box
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = res%cut_val_left
         elseif (associated(res%left) .eqv. .false.) then
            res%box = res%right%box
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val = res%cut_val_right
         else
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = (res%cut_val_left + res%cut_val_right)/2


            ! now remake the true bounding box for self.
            ! Since we are taking unions (in effect) of a tree structure,
            ! this is much faster than doing an exhaustive
            ! search over all points
            res%box%upper = max(res%left%box%upper,res%right%box%upper)
            res%box%lower = min(res%left%box%lower,res%right%box%lower)
         endif
      end if
    end function build_tree_for_range

    integer function select_on_coordinate_value(v,ind,c,alpha,li,ui) &
     result(res)
      ! Move elts of ind around between l and u, so that all points
      ! <= than alpha (in c cooordinate) are first, and then
      ! all points > alpha are second.

      !
      ! Algorithm (matt kennel).
      !
      ! Consider the list as having three parts: on the left,
      ! the points known to be <= alpha.  On the right, the points
      ! known to be > alpha, and in the middle, the currently unknown
      ! points.   The algorithm is to scan the unknown points, starting
      ! from the left, and swapping them so that they are added to
      ! the left stack or the right stack, as appropriate.
      !
      ! The algorithm finishes when the unknown stack is empty.
      !
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, li, ui
      real(kdkind), intent(in) :: alpha
      ! ..
      real(kdkind) :: v(1:,1:)
      integer :: ind(1:)
      integer :: tmp
      ! ..
      integer :: lb, rb
      !
      ! The points known to be <= alpha are in
      ! [l,lb-1]
      !
      ! The points known to be > alpha are in
      ! [rb+1,u].
      !
      ! Therefore we add new points into lb or
      ! rb as appropriate.  When lb=rb
      ! we are done.  We return the location of the last point <= alpha.
      !
      !
      lb = li; rb = ui

      do while (lb < rb)
         if ( v(c,ind(lb)) <= alpha ) then
            ! it is good where it is.
            lb = lb+1
         else
            ! swap it with rb.
            tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
            rb = rb-1
         endif
      end do

      ! now lb .eq. ub
      if (v(c,ind(lb)) <= alpha) then
         res = lb
      else
         res = lb-1
      endif

    end function select_on_coordinate_value

    subroutine select_on_coordinate(v,ind,c,k,li,ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, k, li, ui
      ! ..
      integer :: i, l, m, s, t, u
      ! ..
      real(kdkind) :: v(:,:)
      integer :: ind(:)
      ! ..
      l = li
      u = ui
      do while (l<u)
         t = ind(l)
         m = l
         do i = l + 1, u
            if (v(c,ind(i))<v(c,t)) then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            end if
         end do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         if (m<=k) l = m + 1
         if (m>=k) u = m - 1
      end do
    end subroutine select_on_coordinate

   subroutine spread_in_coordinate(tp,c,l,u,interv)
      ! the spread in coordinate 'c', between l and u.
      !
      ! Return lower bound in 'smin', and upper in 'smax',
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type(interval), intent(out) :: interv
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      real(kdkind) :: last, lmax, lmin, t, smin,smax
      integer :: i, ulocal
      ! ..
      ! .. Local Arrays ..
      real(kdkind), pointer :: v(:,:)
      integer, pointer :: ind(:)
      ! ..
      v => tp%the_data(1:,1:)
      ind => tp%ind(1:)
      smin = v(c,ind(l))
      smax = smin

      ulocal = u

      do i = l + 2, ulocal, 2
         lmin = v(c,ind(i-1))
         lmax = v(c,ind(i))
         if (lmin>lmax) then
            t = lmin
            lmin = lmax
            lmax = t
         end if
         if (smin>lmin) smin = lmin
         if (smax<lmax) smax = lmax
      end do
      if (i==ulocal+1) then
         last = v(c,ind(ulocal))
         if (smin>last) smin = last
         if (smax<last) smax = last
      end if

      interv%lower = smin
      interv%upper = smax

    end subroutine spread_in_coordinate


  subroutine kdtree2_destroy(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    type (kdtree2), pointer :: tp
    ! ..
    call destroy_node(tp%root)

    deallocate (tp%ind)
    nullify (tp%ind)

    if (tp%rearrange) then
       deallocate(tp%rearranged_data)
       nullify(tp%rearranged_data)
    endif

    deallocate(tp)
    return

  contains
    recursive subroutine destroy_node(np)
      ! .. Structure Arguments ..
      type (tree_node), pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
      intrinsic ASSOCIATED
      ! ..
      if (associated(np%left)) then
         call destroy_node(np%left)
         nullify (np%left)
      end if
      if (associated(np%right)) then
         call destroy_node(np%right)
         nullify (np%right)
      end if
      if (associated(np%box)) deallocate(np%box)
      if (associated(np%ind)) deallocate(np%ind)
      deallocate(np)
      return

    end subroutine destroy_node

  end subroutine kdtree2_destroy

  subroutine kdtree2_n_nearest(tp,qv,nn,results)
    ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
    ! returning their indexes and distances in 'indexes' and 'distances'
    ! arrays already allocated passed to this subroutine.
    type (kdtree2), pointer      :: tp
    real(kdkind), target, intent (In)    :: qv(:)
    integer, intent (In)         :: nn
    type(kdtree2_result), target :: results(:)

    integer :: i !***NEW***
    real*8 :: sd, newpri !***NEW***
    type(pq), pointer :: pqp

    call set_rmat_pointer ( tp%aniso ) !*****NEW******

    sr%ballsize = huge(1.0)
    sr%qv => qv
    sr%nn = nn
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false.

    sr%results => results

    sr%nalloc = nn   ! will be checked

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    call validate_query_storage(nn)
    sr%pq = pq_create(results)

    if(associated(tp%root))then
        call search(tp%root)
    else
        !Get all the data
        pqp => sr.pq
        do i = 1, tp.n
            sd = dmet ( sr.Data(:,i), qv )
            sr%nfound = sr%nfound + 1
            newpri = pq_insert(pqp,sd,i)
        enddo
    endif

    if (tp%sort) then
       call kdtree2_sort_results(sr%nfound, results)
    endif
!    deallocate(sr%pqp)

    KDTREE2_NFOUND = sr%nfound

    return
  end subroutine kdtree2_n_nearest

  subroutine kdtree2_n_nearest_around_point(tp,idxin,correltime,nn,results)
    ! Find the 'nn' vectors in the tree nearest to point 'idxin',
    ! with correlation window 'correltime', returing results in
    ! results(:), which must be pre-allocated upon entry.
    type (kdtree2), pointer        :: tp
    integer, intent (In)           :: idxin, correltime, nn
    type(kdtree2_result), target   :: results(:)

    call set_rmat_pointer ( tp%aniso ) !*****NEW******

    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin) ! copy the vector
    sr%ballsize = huge(1.0)       ! the largest real(kdkind) number
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%nn = nn
    sr%nfound = 0

    sr%dimen = tp%dimen
    sr%nalloc = nn

    sr%results => results

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif

    call validate_query_storage(nn)
    sr%pq = pq_create(results)

    call search(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
    deallocate (sr%qv)
    return
  end subroutine kdtree2_n_nearest_around_point

  subroutine kdtree2_r_nearest(tp,qv,r2,nfound,nalloc,results)
    ! find the nearest neighbors to point 'idxin', within SQUARED
    ! Euclidean distance 'r2'.   Upon ENTRY, nalloc must be the
    ! size of memory allocated for results(1:nalloc).  Upon
    ! EXIT, nfound is the number actually found within the ball.
    !
    !  Note that if nfound .gt. nalloc then more neighbors were found
    !  than there were storage to store.  The resulting list is NOT
    !  the smallest ball inside norm r^2
    !
    ! Results are NOT sorted unless tree was created with sort option.
    type (kdtree2), pointer      :: tp
    real(kdkind), target, intent (In)    :: qv(:)
    real(kdkind), intent(in)             :: r2
    integer, intent(out)         :: nfound
    integer, intent (In)         :: nalloc
    type(kdtree2_result), target :: results(:)

    call set_rmat_pointer ( tp%aniso ) !*****NEW******

    !
    sr%qv => qv
    sr%ballsize = r2
    sr%nn = 0      ! flag for fixed ball search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0

    sr%results => results

    call validate_query_storage(nalloc)
    sr%nalloc = nalloc
    sr%overflow = .false.
    sr%ind => tp%ind
    sr%rearrange= tp%rearrange

    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound, results)
    endif

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    endif

    return
  end subroutine kdtree2_r_nearest

  subroutine kdtree2_r_nearest_around_point(tp,idxin,correltime,r2,&
   nfound,nalloc,results)
    !
    ! Like kdtree2_r_nearest, but around a point 'idxin' already existing
    ! in the data set.
    !
    ! Results are NOT sorted unless tree was created with sort option.
    !
    type (kdtree2), pointer      :: tp
    integer, intent (In)         :: idxin, correltime, nalloc
    real(kdkind), intent(in)             :: r2
    integer, intent(out)         :: nfound
    type(kdtree2_result), target :: results(:)
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE

    call set_rmat_pointer ( tp%aniso ) !*****NEW******

    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin) ! copy the vector
    sr%ballsize = r2
    sr%nn = 0    ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%results => results

    sr%nalloc = nalloc
    sr%overflow = .false.

    call validate_query_storage(nalloc)

    !    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    !    sr%il = -1               ! set to invalid indexes

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%rearrange = tp%rearrange
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound,results)
    endif

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    endif

    deallocate (sr%qv)
    return
  end subroutine kdtree2_r_nearest_around_point

  function kdtree2_r_count(tp,qv,r2) result(nfound)
    ! Count the number of neighbors within square distance 'r2'.
    type (kdtree2), pointer   :: tp
    real(kdkind), target, intent (In) :: qv(:)
    real(kdkind), intent(in)          :: r2
    integer                   :: nfound
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE

    call set_rmat_pointer ( tp%aniso ) !*****NEW******
    ! ..
    sr%qv => qv
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0

    nullify(sr%results) ! for some reason, FTN 95 chokes on '=> null()'

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.
    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call search(tp%root)

    nfound = sr%nfound

    return
  end function kdtree2_r_count

  function kdtree2_r_count_around_point(tp,idxin,correltime,r2) &
   result(nfound)
    ! Count the number of neighbors within square distance 'r2' around
    ! point 'idxin' with decorrelation time 'correltime'.
    !
    type (kdtree2), pointer :: tp
    integer, intent (In)    :: correltime, idxin
    real(kdkind), intent(in)        :: r2
    integer                 :: nfound
    ! ..
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE

    call set_rmat_pointer ( tp%aniso ) !*****NEW******

    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin)
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    nullify(sr%results)

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call search(tp%root)

    nfound = sr%nfound

    return
  end function kdtree2_r_count_around_point

  !-------------------------------------------------------------------------------!
  !                                                                               !
  !                                                                               !
  ! These routines were added by John Manchuk (2009) for Orthogonal range queries !
  !                                                                               !
  !                                                                               !
  !-------------------------------------------------------------------------------!

    subroutine kdtree2_ortho_nearest(tp,rmin,rmax,nfound,nalloc,results)
        ! find the nearest neighbors to point 'idxin', within an orthogonal
        ! window that has limits rmin(:) rmax(:), which must have a length
        ! equal to the number of dimensions in the tree.
        !
        ! Upon ENTRY, nalloc must be the size of memory allocated for results(1:nalloc).
        ! Upon EXIT, nfound is the number actually found within the ball.
        !
        ! Results are not sorted, regardless of the sort option for the tree.
        type (kdtree2), pointer      :: tp
        real(kdkind), target, intent (in) :: rmin(:), rmax(:)
        integer, intent(out)         :: nfound
        integer, intent (In)         :: nalloc
        type(kdtree2_result), target :: results(:)

        !
        sr%rmin => rmin
        sr%rmax => rmax
        sr%nn = 0      ! flag for fixed ball search
        sr%nfound = 0
        sr%centeridx = -1
        sr%correltime = 0

        sr%results => results

        call validate_query_storage(nalloc)
        sr%nalloc = nalloc
        sr%overflow = .false.
        sr%ind => tp%ind
        sr%rearrange= tp%rearrange

        if (tp%rearrange) then
            sr%Data => tp%rearranged_data
        else
            sr%Data => tp%the_data
        endif
        sr%dimen = tp%dimen

        call search_ortho(tp%root)

        nfound = sr%nfound

        if (sr%overflow) then
            write (*,*) 'KD_TREE_TRANS: warning! return from kdtree2_ortho_nearest found more neighbors'
            write (*,*) 'KD_TREE_TRANS: than storage was provided for.'
        endif

        return
    end subroutine kdtree2_ortho_nearest


    recursive subroutine search_ortho(node)
        !
        ! Search for points within an orthogonal range
        !
        type (Tree_node), pointer          :: node
        ! ..
        type(tree_node),pointer            :: ncloser, nfarther
        !
        integer                            :: cut_dim, i
        ! ..
        real(kdkind) :: qval_left, qval_right
        real(kdkind),   pointer :: rmin(:),rmax(:)

        if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
            ! we are on a terminal node
            call process_terminal_node_ortho(node)
        else ! we are not on a terminal node
            rmin => sr%rmin(1:) !orthogonal range limits
            rmax => sr%rmax(1:)

            cut_dim = node%cut_dim !dimension being cut

            ! The range along cut_dim falls into one of three categories
            !  - entirely on the left of cut_val
            !  - entirely on the right of cut_val
            !  - spanning both the left and right

            qval_left = rmin(cut_dim)
            qval_right = rmax(cut_dim)

            if( qval_left > node%cut_val )then       !the range is fully on the right
                call search_ortho ( node%right )
            elseif( qval_right < node%cut_val )then  !the range is fully on the left
                call search_ortho ( node%left )
            else                                     !the range spans both
                call search_ortho ( node%right )
                call search_ortho ( node%left )
            endif
        endif
        return
    end subroutine search_ortho

    subroutine process_terminal_node_ortho(node)
        !
        ! Look for data within node that is inside the range
        !
        type (tree_node), pointer :: node
        real(kdkind), pointer :: rmin(:),rmax(:)
        real(kdkind), pointer :: data(:,:)
        integer, pointer :: ind(:)
        integer :: nfound
        integer :: dimen, i, j, ir
        logical :: rearrange

        ! copy values from sr to local variables
        rmin => sr%rmin(1:)         !orthogonal range limits
        rmax => sr%rmax(1:)
        dimen = sr%dimen            !dimensionality
        rearrange = sr%rearrange    !is data rearranged
        ind => sr%ind(1:)           !data indexes
        data => sr%Data(1:,1:)      !data points
        nfound = sr%nfound          !number found so far

        ! search through terminal bucket and add points that are
        ! within the range
        mainloop: do i = node%l, node%u
            if (rearrange) then
                j = i
            else
                j = ind(i)
            endif

            !Check against all range axes
            do ir = 1, dimen
                if(data(ir,j) < rmin(ir) .or. data(ir,j) > rmax(ir)) cycle mainloop
            enddo

            !If we get to here, keep the point
            nfound = nfound + 1
            if (nfound > sr%nalloc) then
                ! oh nuts, we have to add another one to the tree but
                ! there isn't enough room.
                sr%overflow = .true.
            else
                sr%results(nfound)%idx = ind(i)
            endif
        enddo mainloop
        sr%nfound = nfound
    end subroutine process_terminal_node_ortho

  !-------------------------------------------------------------------------------!
  !                                                                               !
  !                                                                               !
  ! End of new routine                                                            !
  !                                                                               !
  !                                                                               !
  !-------------------------------------------------------------------------------!


  subroutine validate_query_storage(n)
    !
    ! make sure we have enough storage for n
    !
    integer, intent(in) :: n

    if (size(sr%results,1) .lt. n) then
       write (*,*) 'KD_TREE_TRANS:  you did not provide enough storage for results(1:n)'
       stop
       return
    endif

    return
  end subroutine validate_query_storage

  function square_distance(d, iv,qv) result (res)
    ! distance between iv[1:n] and qv[1:n]
    ! .. Function Return Value ..
    ! re-implemented to improve vectorization.
    real(kdkind) :: res
    ! ..
    ! ..
    ! .. Scalar Arguments ..
    integer :: d
    ! ..
    ! .. Array Arguments ..
    real(kdkind) :: iv(:),qv(:)
    ! ..
    ! ..
    res = sum( (iv(1:d)-qv(1:d))**2 )

    !res = dmet(iv,qv)

  end function square_distance

  recursive subroutine search(node)
    !
    ! This is the innermost core routine of the kd-tree search.  Along
    ! with "process_terminal_node", it is the performance bottleneck.
    !
    ! This version uses a logically complete secondary search of
    ! "box in bounds", whether the sear
    !
    type (Tree_node), pointer          :: node
    ! ..
    type(tree_node),pointer            :: ncloser, nfarther
    !
    integer                            :: cut_dim, i
    ! ..
    real(kdkind)                               :: qval, dis
    real(kdkind)                               :: ballsize
    real(kdkind), pointer           :: qv(:)
    type(interval), pointer :: box(:)

    if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
       ! we are on a terminal node
       if (sr%nn .eq. 0) then
          call process_terminal_node_fixedball(node)
       else
          call process_terminal_node(node)
       endif
    else
       ! we are not on a terminal node
       qv => sr%qv(1:)
       cut_dim = node%cut_dim
       qval = qv(cut_dim)

       if (qval < node%cut_val) then
          ncloser => node%left
          nfarther => node%right
          !dis = (node%cut_val_right - qval)**2
          dis = dmet ( node%cut_val_right, qval, cut_dim )
!          extra = node%cut_val - qval
       else
          ncloser => node%right
          nfarther => node%left
          !dis = (node%cut_val_left - qval)**2
          dis = dmet ( node%cut_val_left, qval, cut_dim )
!          extra = qval- node%cut_val_left
       endif

       if (associated(ncloser)) call search(ncloser)

       ! we may need to search the second node.
       if (associated(nfarther)) then
          ballsize = sr%ballsize
!          dis=extra**2
          if (dis <= ballsize) then
             !
             ! we do this separately as going on the first cut dimen is often
             ! a good idea.
             ! note that if extra**2 < sr%ballsize, then the next
             ! check will also be false.
             !
             box => node%box(1:)
             do i=1,sr%dimen
                if (i .ne. cut_dim) then
                   dis = dis + dis2_from_bnd(qv(i),box(i)%lower,box(i)%upper,i)
                   if (dis > ballsize) then
                      return
                   endif
                endif
             end do

             !
             ! if we are still here then we need to search more.
             !
             call search(nfarther)
          endif
       endif
    end if
  end subroutine search


  real(kdkind) function dis2_from_bnd(x,amin,amax,idim) result (res)
    real(kdkind), intent(in) :: x, amin,amax
    integer, intent (in) :: idim

    if (x > amax) then
       !res = (x-amax)**2;
       res = dmet ( x, amax, idim )
       return
    else
       if (x < amin) then
          !res = (amin-x)**2;
          res = dmet ( amin, x, idim )
          return
       else
          res = 0.0
          return
       endif
    endif
    return
  end function dis2_from_bnd

  logical function box_in_search_range(node, sr) result(res)
    !
    ! Return the distance from 'qv' to the CLOSEST corner of node's
    ! bounding box
    ! for all coordinates outside the box.   Coordinates inside the box
    ! contribute nothing to the distance.
    !
    type (tree_node), pointer :: node
    type (tree_search_record), pointer :: sr

    integer :: dimen, i
    real(kdkind)    :: dis, ballsize
    real(kdkind)    :: l, u

    dimen = sr%dimen
    ballsize = sr%ballsize
    dis = 0.0
    res = .true.
    do i=1,dimen
       l = node%box(i)%lower
       u = node%box(i)%upper
       dis = dis + (dis2_from_bnd(sr%qv(i),l,u,i))
       if (dis > ballsize) then
          res = .false.
          return
       endif
    end do
    res = .true.
    return
  end function box_in_search_range


  subroutine process_terminal_node(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure.
    !
    type (tree_node), pointer          :: node
    !
    real(kdkind), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kdkind), pointer          :: data(:,:)
    !
    integer                :: dimen, i, indexofi, k, centeridx, correltime
    real(kdkind)                   :: ballsize, sd, newpri
    logical                :: rearrange
    type(pq), pointer      :: pqp
    !
    ! copy values from sr to local variables
    !
    !
    ! Notice, making local pointers with an EXPLICIT lower bound
    ! seems to generate faster code.
    ! why?  I don't know.
    qv => sr%qv(1:)
    pqp => sr%pq
    dimen = sr%dimen
    ballsize = sr%ballsize
    rearrange = sr%rearrange

    if(associated(node%ind))then   !***NEW***
        ind => node%ind
    else
        ind => sr%ind(1:)
    endif

    data => sr%Data(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime

    !    doing_correl = (centeridx >= 0)  ! Do we have a decorrelation window?
    !    include_point = .true.    ! by default include all points
    ! search through terminal bucket.

    mainloop: do i = node%l, node%u
       if (rearrange) then
          sd = 0.0

          sd = dmet ( data(:,i), qv )           !***NEW***
          if (sd>ballsize) cycle mainloop

!          do k = 1,dimen
!             sd = sd + (data(k,i) - qv(k))**2
!             if (sd>ballsize) cycle mainloop
!          end do

          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0

          sd = dmet ( data(:,indexofi), qv )    !***NEW***
          if (sd>ballsize) cycle mainloop

!          do k = 1,dimen
!             sd = sd + (data(k,indexofi) - qv(k))**2
!             if (sd>ballsize) cycle mainloop
!          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx) < correltime) cycle mainloop
       endif


       !
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

       if (sr%nfound .lt. sr%nn) then
          !
          ! add this point unconditionally to fill list.
          !
          sr%nfound = sr%nfound +1
          newpri = pq_insert(pqp,sd,indexofi)
          if (sr%nfound .eq. sr%nn) ballsize = newpri
          ! we have just filled the working list.
          ! put the best square distance to the maximum value
          ! on the list, which is extractable from the PQ.

       else
          !
          ! now, if we get here,
          ! we know that the current node has a squared
          ! distance smaller than the largest one on the list, and
          ! belongs on the list.
          ! Hence we replace that with the current one.
          !
          ballsize = pq_replace_max(pqp,sd,indexofi)
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%ballsize = ballsize

  end subroutine process_terminal_node

  subroutine process_terminal_node_fixedball(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure, i.e.
    ! save all within a fixed ball.
    !
    type (tree_node), pointer          :: node
    !
    real(kdkind), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kdkind), pointer          :: data(:,:)
    !
    integer                :: nfound
    integer                :: dimen, i, indexofi, k
    integer                :: centeridx, correltime, nn
    real(kdkind)                   :: ballsize, sd
    logical                :: rearrange

    !
    ! copy values from sr to local variables
    !
    qv => sr%qv(1:)
    dimen = sr%dimen
    ballsize = sr%ballsize
    rearrange = sr%rearrange

    if(associated(node%ind))then   !***NEW***
        ind => node%ind
    else
        ind => sr%ind(1:)
    endif

    data => sr%Data(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime
    nn = sr%nn ! number to search for
    nfound = sr%nfound

    ! search through terminal bucket.
    mainloop: do i = node%l, node%u

       !
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

       ! which index to the point do we use?

       if (rearrange) then
          sd = 0.0

          sd = dmet ( data(:,i), qv )           !***NEW***
          if (sd>ballsize) cycle mainloop

!          do k = 1,dimen
!             sd = sd + (data(k,i) - qv(k))**2
!             if (sd>ballsize) cycle mainloop
!          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0

          sd = dmet ( data(:,indexofi), qv )    !***NEW***
          if (sd>ballsize) cycle mainloop

!          do k = 1,dimen
!             sd = sd + (data(k,indexofi) - qv(k))**2
!             if (sd>ballsize) cycle mainloop
!          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx)<correltime) cycle mainloop
       endif

       nfound = nfound+1
       if (nfound .gt. sr%nalloc) then
          ! oh nuts, we have to add another one to the tree but
          ! there isn't enough room.
          sr%overflow = .true.
       else
          sr%results(nfound)%dis = sd
          sr%results(nfound)%idx = indexofi
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%nfound = nfound
  end subroutine process_terminal_node_fixedball

  subroutine kdtree2_n_nearest_brute_force(tp,qv,nn,results)
    ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (kdtree2), pointer :: tp
    real(kdkind), intent (In)       :: qv(:)
    integer, intent (In)    :: nn
    type(kdtree2_result)    :: results(:)

    integer :: i, j, k
    real(kdkind), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = square_distance(tp%dimen,qv,tp%the_data(:,i))
    end do
    ! now find 'n' smallest distances
    do i = 1, nn
       results(i)%dis =  huge(1.0)
       results(i)%idx = -1
    end do
    do i = 1, tp%n
       if (all_distances(i)<results(nn)%dis) then
          ! insert it somewhere on the list
          do j = 1, nn
             if (all_distances(i)<results(j)%dis) exit
          end do
          ! now we know 'j'
          do k = nn - 1, j, -1
             results(k+1) = results(k)
          end do
          results(j)%dis = all_distances(i)
          results(j)%idx = i
       end if
    end do
    deallocate (all_distances)
  end subroutine kdtree2_n_nearest_brute_force


  subroutine kdtree2_r_nearest_brute_force(tp,qv,r2,nfound,results)
    ! find the nearest neighbors to 'qv' with distance**2 <= r2 by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (kdtree2), pointer :: tp
    real(kdkind), intent (In)       :: qv(:)
    real(kdkind), intent (In)       :: r2
    integer, intent(out)    :: nfound
    type(kdtree2_result)    :: results(:)

    integer :: i, nalloc
    real(kdkind), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = square_distance(tp%dimen,qv,tp%the_data(:,i))
    end do

    nfound = 0
    nalloc = size(results,1)

    do i = 1, tp%n
       if (all_distances(i)< r2) then
          ! insert it somewhere on the list
          if (nfound .lt. nalloc) then
             nfound = nfound+1
             results(nfound)%dis = all_distances(i)
             results(nfound)%idx = i
          endif
       end if
    enddo
    deallocate (all_distances)

    call kdtree2_sort_results(nfound,results)


  end subroutine kdtree2_r_nearest_brute_force

  subroutine kdtree2_sort_results(nfound,results)
    !  Use after search to sort results(1:nfound) in order of increasing
    !  distance.
    integer, intent(in)          :: nfound
    type(kdtree2_result), target :: results(:)
    !
    !

    !THIS IS BUGGY WITH INTEL FORTRAN
    !    If (nfound .Gt. 1) Call heapsort(results(1:nfound)%dis,results(1:nfound)%ind,nfound)
    !
    if (nfound .gt. 1) call heapsort_struct(results,nfound)

    return
  end subroutine kdtree2_sort_results

  subroutine heapsort(a,ind,n)
    !
    ! Sort a(1:n) in ascending order, permuting ind(1:n) similarly.
    !
    ! If ind(k) = k upon input, then it will give a sort index upon output.
    !
    integer,intent(in)          :: n
    real(kdkind), intent(inout)         :: a(:)
    integer, intent(inout)      :: ind(:)

    !
    !
    real(kdkind)        :: value   ! temporary for a value from a()
    integer     :: ivalue  ! temporary for a value from ind()

    integer     :: i,j
    integer     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! Generate initial idum array
    !    end do

    if(n.eq.1) return

    do
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft); ivalue=ind(ileft)
       else
          value=a(iright); ivalue=ind(iright)
          a(iright)=a(1); ind(iright)=ind(1)
          iright=iright-1
          if (iright == 1) then
             a(1)=value;ind(1)=ivalue
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright)
          if(j < iright) then
             if(a(j) < a(j+1)) j=j+1
          endif
          if(value < a(j)) then
             a(i)=a(j); ind(i)=ind(j)
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value; ind(i)=ivalue
    end do
  end subroutine heapsort

  subroutine heapsort_struct(a,n)
    !
    ! Sort a(1:n) in ascending order
    !
    !
    integer,intent(in)                 :: n
    type(kdtree2_result),intent(inout) :: a(:)

    !
    !
    type(kdtree2_result) :: value ! temporary value

    integer     :: i,j
    integer     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! Generate initial idum array
    !    end do

    if(n.eq.1) return

    do
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft)
       else
          value=a(iright)
          a(iright)=a(1)
          iright=iright-1
          if (iright == 1) then
             a(1) = value
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright)
          if(j < iright) then
             if(a(j)%dis < a(j+1)%dis) j=j+1
          endif
          if(value%dis < a(j)%dis) then
             a(i)=a(j);
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value
    end do
  end subroutine heapsort_struct

end module kdtree2_module

