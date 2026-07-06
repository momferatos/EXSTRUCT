module types
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!
  !Data type definitions
  !!!!!!!!!!!!!!!!!!!!!!
  !Integer data types
  integer, parameter     :: i1b = 1 !integer 1 bit
  integer, parameter     :: i2b = 2 !integer 2 bits
  integer, parameter     :: i4b = 4 !integer 4 bits
  integer, parameter     :: i8b = 8 !integer 8 bits
  integer, parameter     :: ik = i8b!main integer type
  !Real data types
  integer, parameter     :: sp = 4  !single precision real
  integer, parameter     :: dp = 8  !double precision real
  integer, parameter     :: qp = 16 !quadruple precision real
  integer, parameter     :: rk = sp !main real type
  integer, parameter     :: rks = sp!real type for arrays
  !Complex data types
  integer, parameter     :: csp = 4  !single precision complex
  integer, parameter     :: cdp = 8  !double precision real
  integer, parameter     :: ck = cdp !main complex type
  integer, parameter     :: cks = rks!complex type for arrays
  integer(ik), parameter :: inrk = sp
  integer(ik), parameter :: outrk = dp
  ! Structure coordinates are stored as integer(i2b), and the periodic
  ! unwrapping in translate_structures can push them up to twice the grid
  ! extent, so each axis is limited to huge(i2b)/2 = 16383 points.
  integer(ik), parameter :: max_extent = (huge(0_i2b) - 1) / 2
  integer(ik)            :: nnintrv

  type :: interval
     !!!!!!!!!!!!!!!!!!!
     !Interval data type
     !!!!!!!!!!!!!!!!!!!
     ! Default-initialise the payload: consumers treat zero limits/indices as
     ! "skip this node", so a freshly allocated node must never carry
     ! undefined values (allocate does not zero memory).
     integer(i2b), dimension(1:2) :: limits = 0
     integer(i2b)                 :: i = 0, j = 0, k = 0
     type(interval), pointer      :: next => null()
  end type interval

  type :: structure
     !!!!!!!!!!!!!!!!!!!!
     !Structure data type
     !!!!!!!!!!!!!!!!!!!!
     type(interval), pointer                   :: intrv => null()
     integer(i2b)                              :: nintrv = 0
     logical                                   :: crosses_i = .false.
     logical                                   :: crosses_j = .false.
     logical                                   :: crosses_k = .false.
     type(structure), pointer                  :: next => null()
     integer(i8b)                              :: npoints = 0
     integer(i2b), dimension(:,:), allocatable :: points
  end type structure

  type :: list
     !!!!!!!!!!!!!!!!!!!!!!
     !Linked list data type
     !!!!!!!!!!!!!!!!!!!!!!
     integer(i2b)        :: i = 0,j = 0,k = 0
     type(list), pointer :: next => null()
  end type list

end module types

module data
  use types
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!
  !Wrapper module for data
  !!!!!!!!!!!!!!!!!!!!!!!!
  integer(ik)                                   :: VOLUME
  integer(ik)                                   :: n1,n2,n3,n4,nstruct
  integer(i2b), dimension(:,:),allocatable      :: neigh
  real(rk), dimension(:,:,:), allocatable       :: field
  integer(i2b), dimension(:,:,:,:), allocatable :: intrv
  integer(ik), dimension(:,:),allocatable       :: nintrv
  logical, dimension(:,:,:), allocatable        :: cintrv
  logical, dimension(:,:,:), allocatable        :: mask,refmask
  type(structure), target                       :: struct
  type(interval), pointer                       :: intrv_ptr => null()
  type(structure), pointer                      :: struct_ptr => null()

end module data

module routines
  use types
  use data
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Wrapper module for Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  function per(n, m)
    implicit none
    integer(ik)             :: per
    integer(ik), intent(in) :: n, m
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Handles periodic boundaries: wrap index n into [1, m]
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(n==m + 1) then
       per = 1
    else if(n==0) then
       per = m
    else
       per = n
    end if

    return

  end function per

  subroutine mean_sdev(n1,n2,n3,field,mean,sdev)
    use types
    implicit none
    integer(ik), intent(IN)                         :: n1,n2,n3
    real(rk), dimension(1:n1,1:n2,1:n3), intent(IN) :: field
    real(rk), intent(OUT)                           :: mean, sdev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean and standard deviation of a scalar field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                     :: i,j,k
    real(dp)                                        :: smean, ssdev, fac

    ! Accumulate in double precision and divide once at the end. The field is
    ! single precision and has n1*n2*n3 (~1e7-1e8) points: summing it in single
    ! precision loses the running sum to roundoff once it dwarfs the individual
    ! terms, biasing mean/sdev and hence the mean+m*sdev cutoff that defines the
    ! structures.
    fac = 1.0_dp / real(n1,dp) / real(n2,dp) / real(n3,dp)

    smean = 0.0_dp

    !$omp parallel do reduction(+:smean)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             smean = smean + real(field(i,j,k),dp)
          end do
       end do
    end do
    !$omp end parallel do

    smean = smean * fac

    ssdev = 0.0_dp

    !$omp parallel do reduction(+:ssdev)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             ssdev = ssdev + (real(field(i,j,k),dp) - smean)**2
          end do
       end do
    end do
    !$omp end parallel do

    ssdev = sqrt(ssdev * fac)

    mean = real(smean,rk)
    sdev = real(ssdev,rk)

    return

  end subroutine mean_sdev

  subroutine check
    implicit none
    integer(ik) :: n,m,i,j,k
    integer(ik) :: nerr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Check consistency of structures with initial mask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nerr = 0

    !$omp workshare
    mask = .false.
    !$omp end workshare

    struct_ptr => struct
    n = 0

    ! .and. does not short-circuit in Fortran: the struct_ptr%... operand may
    ! be evaluated even when struct_ptr is disassociated, which dereferences a
    ! null pointer. Test the component separately inside the loop (everywhere
    ! the list is walked).
    do while(associated(struct_ptr))
       if(.not.allocated(struct_ptr%points)) exit
       n = n + 1
       do m=1,struct_ptr%npoints
          i = struct_ptr%points(m,1)
          j = struct_ptr%points(m,2)
          k = struct_ptr%points(m,3)
          if(mask(i,j,k)) then
             print '(a,3i5)', 'error: double point ', i,j,k
             nerr = nerr + 1
          end if
          mask(i,j,k) = .true.
       end do
       struct_ptr => struct_ptr%next
    end do

    do k=1,n3
       do j=1,n2
          do i=1,n1
             if(mask(i,j,k) .and. .not.refmask(i,j,k)) then
                print '(a,3i5)', 'error: extra point ',i,j,k
                nerr = nerr + 1
             end if

             if(.not.mask(i,j,k) .and. refmask(i,j,k)) then
                print '(a,3i5)', 'error: missing point ',i,j,k
                nerr = nerr + 1
             end if
          end do
       end do
    end do

    if(nerr /= 0) then
       print '(a,i0,a)', 'check: ', nerr, ' inconsistent points.'
       error stop 'error: structure extraction failed the self-check.'
    end if

    return

  end subroutine check

  subroutine intervals(n1,n2,n3,intrv,nintrv)
    use types
    implicit none
    integer(ik),intent(IN)                       :: n1,n2,n3
    integer(i2b), dimension(1:n1,1:n2,1:(n3+1)/2,2)  :: intrv
    integer(ik), dimension(1:n1,1:n2)            :: nintrv
!!!!!!!!!!!!!!!!!
    !Set-up intervals
!!!!!!!!!!!!!!!!!
    integer(ik)                                  :: i,j,k,k1,k2
    logical                                      :: openintrv
    integer(i2b), dimension(:,:,:,:),allocatable :: intrv1
    integer(i2b), dimension(:,:),allocatable     :: nintrv1

    allocate(intrv1(n1,n2,(n3+1)/2,2))
    allocate(nintrv1(n1,n2))
    nnintrv = 0
    nintrv1 = 0
    intrv1 = 0

    do j=1,n2
       do i=1,n1
          openintrv = .false.
          do k=1,n3
             if(refmask(i,j,k) .and. (.not.openintrv)) then
                nintrv1(i,j) = nintrv1(i,j) + 1
                intrv1(i,j,nintrv1(i,j),1) = k
                openintrv = .true.
             else if((.not.refmask(i,j,k)) .and. openintrv) then
                intrv1(i,j,nintrv1(i,j),2) = k - 1
                openintrv = .false.
                nnintrv = nnintrv + 1
             end if
          end do
          if(openintrv) then
             intrv1(i,j,nintrv1(i,j),2) = n3
             openintrv = .false.
          end if
       end do
    end do

    intrv = intrv1
    nintrv = nintrv1

    do j=1,n2
       do i=1,n1
          if(nintrv1(i,j) >= 2) then
             if(intrv1(i,j,1,1)==1 .and. intrv1(i,j,nintrv1(i,j),2)==n3) then
                nintrv(i,j) = nintrv1(i,j) - 1
                intrv(i,j,1,1) = intrv1(i,j,nintrv1(i,j),1)
                intrv(i,j,1,2) = intrv1(i,j,1,2)
                k1 = 2
                k2 = nintrv1(i,j) - 1
             else
                nintrv(i,j) = nintrv1(i,j)
                k1 = 1
                k2 = nintrv1(i,j)
             end if

             do k=k1,k2
                intrv(i,j,k,:) = intrv1(i,j,k,:)
             end do
          end if
       end do
    end do

    deallocate(intrv1)
    deallocate(nintrv1)

    return

  end subroutine intervals

  subroutine find_structures
    use data
    use types
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Top-level routine for structure identification
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                  :: i,j,k
    type(interval)                               :: intrv1

    n4 = maxval(nintrv)
    allocate(cintrv(n1,n2,n4))
    cintrv(:,:,:) = .false.
    nstruct = 0
    nullify(struct_ptr)

    do j=1,n2
       do i=1,n1
          do k=1,nintrv(i,j)
             if(.not.cintrv(i,j,k)) then
                ! Open a structure node only once a new structure is seeded.
                ! The old pre-allocated tail node showed up downstream as a
                ! phantom zero-volume structure: it inflated the structure
                ! count by one, forced the minimum volume to 0 and appended a
                ! bogus 0 record to gs.dat.
                if(.not.associated(struct_ptr)) then
                   struct_ptr => struct
                else
                   allocate(struct_ptr%next)
                   struct_ptr => struct_ptr%next
                end if
                nstruct = nstruct + 1
                intrv1%limits(1) = intrv(i,j,k,1)
                intrv1%limits(2) = intrv(i,j,k,2)
                intrv1%i = i
                intrv1%j = j
                intrv1%k = k
                call find_neighbours(intrv1)
                cintrv(i,j,k) = .true.
             end if
          end do
       end do
    end do

    return

  end subroutine find_structures

  subroutine find_neighbours(intrv1)
    use data,only:cintrv,nintrv,intrv_ptr,struct_ptr,neigh
    use types
    implicit none
    type(interval), intent(INOUT) :: intrv1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Find all neighbours of a given interval
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                   :: first = 0,last = 0
    integer(ik), dimension(8,2)   :: neighbours
    integer(ik)                   :: ll,north,south,east,west,center1,center2
    integer(ik)                   :: i,j,k,ii,jj,kk
    integer(i2b), dimension(2)    :: limits
    logical                       :: cond
    logical                       :: lexit
    logical                       :: lcycle
    
    first = 0
    last = 0
    cond = .true.
    lexit = .false.
    
    do while(.true.)
       lcycle =.false.
       
       if(cond) then
          i = intrv1%i
          j = intrv1%j
          k = intrv1%k
          limits = intrv1%limits
       end if

       if(i==0 .or. j==0 .or. k==0&
            & .or. i>n1 .or. j>n2 .or. k>n3) then

          call sub1(cond, lexit, lcycle)

          if(lexit) then
             exit
          end if

          if(lcycle) then
             cycle
          end if

       end if

       call sub2

       call sub1(cond, lexit, lcycle)

       if(lexit) then
          exit
       end if

       if(lcycle) then
          cycle
       end if

    end do

    return

  contains

    subroutine sub1(lcond, lexit, lcycle)
      implicit none
      logical, intent(out) :: lcond, lexit, lcycle

      ! All three outputs must be defined on every path: they are intent(out)
      ! (undefined on entry) and the caller tests lexit/lcycle immediately, so
      ! leaving any unset means the flood fill exits or cycles on garbage and
      ! connected structures get torn apart.
      lcond  = .false.
      lexit  = .false.
      lcycle = .false.

      if(first<last) then
         first = first + 1
         ! Queue entries are written by sub2 before they are ever popped and
         ! always hold valid 1-based indices, so they can be used directly.
         ii = neigh(first,1)
         jj = neigh(first,2)
         kk = neigh(first,3)
         intrv1%limits(1) = intrv(ii,jj,kk,1)
         intrv1%limits(2) = intrv(ii,jj,kk,2)
         intrv1%i = ii
         intrv1%j = jj
         intrv1%k = kk
         lcond = .true.
         lcycle = .true.
      else
         lexit = .true.
      end if

      return

    end subroutine sub1

    subroutine sub2
      implicit none

      if(.not.cintrv(i,j,k)) then
         call append_interval(i, j, k)
         struct_ptr%nintrv = struct_ptr%nintrv + 1
         cintrv(i,j,k) = .true.
      end if

      center1 = i
      center2 = j

      if(center1 /= n1) then
         north = center1 + 1
      else
         north = 1
      end if

      if(center1 /= 1) then
         south = center1 - 1
      else
         south = n1
      end if

      if(center2 /= n2) then
         east = center2 + 1
      else
         east = 1
      end if

      if(center2 /= 1) then
         west = center2 - 1
      else
         west = n2
      end if

      neighbours(1,1) = south
      neighbours(1,2) = west

      neighbours(2,1) = south
      neighbours(2,2) = center2

      neighbours(3,1) = south
      neighbours(3,2) = east

      neighbours(4,1) = center1
      neighbours(4,2) = west

      neighbours(5,1) = center1
      neighbours(5,2) = east

      neighbours(6,1) = north
      neighbours(6,2) = west

      neighbours(7,1) = north
      neighbours(7,2) = center2

      neighbours(8,1) = north
      neighbours(8,2) = east

      do ll=1,8
         ii = neighbours(ll,1)
         jj = neighbours(ll,2)
         do kk=1,nintrv(ii,jj)
            if(is_neighbour(limits,intrv(ii,jj,kk,:)) .and. &
                 &.not.cintrv(ii,jj,kk)) then

               if(center1==1 .and. ii==n1 .or. center1==n1 .and. ii==1) then
                  struct_ptr%crosses_i = .true.
               end if

               if(center2==1 .and. jj==n2 .or. center2==n2 .and. jj==1) then
                  struct_ptr%crosses_j = .true.
               end if

               if(intrv(ii,jj,kk,1)>intrv(ii,jj,kk,2)) then
                  struct_ptr%crosses_k = .true.
               end if

               call append_interval(ii, jj, kk)
               cintrv(ii,jj,kk) = .true.

               last = last + 1
               neigh(last,1) = ii
               neigh(last,2) = jj
               neigh(last,3) = kk

            end if
         end do
      end do

      return

    end subroutine sub2

    subroutine append_interval(ai, aj, ak)
      implicit none
      integer(ik), intent(in) :: ai, aj, ak
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Append the interval at (ai,aj,ak) to the current structure's list.
      !The node is allocated only here, together with its data: allocating
      !ahead of the data (as the old code did) left a trailing node with
      !undefined limits on every structure.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(.not.associated(struct_ptr%intrv)) then
         allocate(struct_ptr%intrv)
         intrv_ptr => struct_ptr%intrv
      else
         allocate(intrv_ptr%next)
         intrv_ptr => intrv_ptr%next
      end if

      intrv_ptr%limits(1) = intrv(ai,aj,ak,1)
      intrv_ptr%limits(2) = intrv(ai,aj,ak,2)
      intrv_ptr%i = ai
      intrv_ptr%j = aj
      intrv_ptr%k = ak

      return

    end subroutine append_interval

  end subroutine find_neighbours

  function is_neighbour(i1,i2)
    use types
    implicit none
    integer(i2b), dimension(2),intent(IN) :: i1,i2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Check whether two intervals are neighbours
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical                               :: is_neighbour
    integer                               :: a,b,x,y
    integer                               :: nn3

    nn3 = int(n3)
    a = i1(1)
    b = i1(2)
    x = i2(1)
    y = i2(2)

    ! Recover a <= b and x <= y intervals thanks to periodicity:
    if(b < a) then
       b = b + nn3
    end if

    if(y < x) then
       y = y + nn3
    end if

    ! Check whether the intervals are connected in the three possible 
    !configurations due to periodicity.
    is_neighbour =  are_connected(a,b,x,y) .or.  &
         are_connected(a + nn3,b + nn3,x,y) .or.  &
         are_connected(a,b,x + nn3,y + nn3)

    return

  contains

    logical function are_connected(a,b,x,y)
      integer :: a,b,x,y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! two integer intervals [a,b] and [x,y] intersect iff (x <= b and a <= y)
      ! two intervals [a,b] and [x,y] connect iff [a - 1,b + 1] intersects [x,y]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      are_connected = ((x <= b + 1) .and. (a <= y + 1))

      return

    end function are_connected

  end function is_neighbour

  subroutine output_structures
    use types
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!
    !Write strucures to file
!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)    :: nfile
    character(128) :: filename

    nfile = 0
    struct_ptr => struct

    do while(associated(struct_ptr))
       if(.not.allocated(struct_ptr%points)) exit
       if(struct_ptr%npoints > VOLUME) then
          write(filename,'(a,i0,a)') 'out.',nfile,'.vtk'
          call write_vtk_file(trim(filename),struct_ptr)
          nfile = nfile + 1
       end if
       struct_ptr => struct_ptr%next
    end do

    print '(a,i0,a)', 'Wrote ', nfile, ' files.'

    return

  end subroutine output_structures

  subroutine write_vtk_file(filename,struct_ptr)
    use types 
    implicit none
    character(*)             :: filename
    type(structure), pointer, intent(IN) :: struct_ptr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Write vtk file for individual structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)              :: m,vtk_file = 999
    real(rk)                 :: x,y,z

    open(vtk_file,file=filename,form='formatted',action='write')

    write(vtk_file, '(t1,a)') '# vtk DataFile Version 2.0'
    write(vtk_file, '(t1,a)')  ' '
    write(vtk_file, '(t1,a)') 'ASCII'
    write(vtk_file, '(t1,a)')  ' '
    write(vtk_file, '(t1,a)') 'DATASET POLYDATA'
    write(vtk_file, '(t1,a,i20,a)') 'POINTS ', struct_ptr%npoints, ' double'

    do m=1,struct_ptr%npoints
       x = struct_ptr%points(m,4)
       y = struct_ptr%points(m,5)
       z = struct_ptr%points(m,6)
       write(vtk_file,'(t1,3e19.8)') x,y,z
    end do

    ! POLYDATA with bare POINTS has no cells, and viewers such as ParaView
    ! render it as empty: emit one VERTEX cell per point (cell list size is
    ! 2*npoints: a count and a 0-based point index per cell).
    write(vtk_file, '(t1,a,i0,a,i0)') 'VERTICES ', struct_ptr%npoints, &
         ' ', 2 * struct_ptr%npoints
    do m=1,struct_ptr%npoints
       write(vtk_file,'(t1,a,i0)') '1 ', m - 1
    end do

    close(vtk_file,status='keep')

    return

  end subroutine write_vtk_file

  subroutine translate_structures
    use types
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Translate structures so that they are continuous and not disrupted by
    !periodic boundaries!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, dimension(:,:,:), allocatable :: box
    integer(ik)                            :: m,i,j,k,n,step_i,step_j,step_k

    allocate(box(1:n1,1:n2,1:n3))
    struct_ptr => struct

    n = 0
    do while(associated(struct_ptr))
       if(.not.allocated(struct_ptr%points)) exit
       n = n + 1
       !$omp workshare
       box(:,:,:) = .false.
       !$omp end workshare

       !$omp parallel do private(i,j,k)
       do m=1,struct_ptr%npoints
          i=struct_ptr%points(m,1)
          j=struct_ptr%points(m,2)
          k=struct_ptr%points(m,3)
          box(i,j,k) = .true.
       end do
       !$omp end parallel do

       m = 0
       step_i = 0
       do while(connected(1_ik) .and. m<n1)
          call rotate(1_ik)
          m = m + 1
       end do
       step_i = m

       m = 0
       step_j = 0
       do while(connected(2_ik) .and. m<n2)
          call rotate(2_ik)
          m = m + 1
       end do
       step_j = m

       m = 0
       step_k = 0
       do while(connected(3_ik) .and. m<n3)
          call rotate(3_ik)
          m = m + 1
       end do
       step_k = m

       m = 0
       do k=1,n3
          do j=1,n2
             do i=1,n1
                if(box(i,j,k)) then
                   m = m + 1
                   struct_ptr%points(m,4) = i + step_i
                   struct_ptr%points(m,5) = j + step_j
                   struct_ptr%points(m,6) = k + step_k
                end if
             end do
          end do
       end do

       struct_ptr => struct_ptr%next

    end do

    return

  contains

    function connected(dim)
      implicit none
      integer(ik), intent(IN) :: dim
      logical :: connected,cell
      integer(ik) :: i,j,ii,jj,iii,jjj

      ! reduction(.or.:connected): every thread used to write the shared flag
      ! unsynchronised (a data race). With the reduction each thread works on
      ! a private copy, and 'if(connected) cycle' short-circuits only that
      ! thread's own remaining iterations.
      connected = .false.

      if(dim==1_ik) then
         ! i-face: does the box(1,:,:) plane connect to box(n1,:,:)?
         ! scan the (n2,n3) plane and wrap neighbours by n2,n3.
         !$omp parallel do private(j,ii,jj,iii,jjj,cell) reduction(.or.:connected)
         do i=1,n2
            if(connected) cycle
            do j=1,n3
               if(connected) cycle
               cell=box(1,i,j)
               if(cell) then
                  do ii=i-1,i+1
                     if(connected) cycle
                     do jj=j-1,j+1
                        if(connected) cycle
                        iii = per(ii,n2)
                        jjj = per(jj,n3)
                        if(box(n1,iii,jjj)) connected = .true.
                     end do
                  end do
               end if
            end do
         end do
         !$omp end parallel do
      else if(dim==2_ik) then
         ! j-face: does box(:,1,:) connect to box(:,n2,:)?
         ! scan the (n1,n3) plane and wrap neighbours by n1,n3.
         !$omp parallel do private(j,ii,jj,iii,jjj,cell) reduction(.or.:connected)
         do i=1,n1
            if(connected) cycle
            do j=1,n3
               if(connected) cycle
               cell=box(i,1,j)
               if(cell) then
                  do ii=i-1,i+1
                     if(connected) cycle
                     do jj=j-1,j+1
                        if(connected) cycle
                        iii = per(ii,n1)
                        jjj = per(jj,n3)
                        if(box(iii,n2,jjj)) connected = .true.
                     end do
                  end do
               end if
            end do
         end do
         !$omp end parallel do
      else if(dim==3_ik) then
         ! k-face: does box(:,:,1) connect to box(:,:,n3)?
         ! scan the (n1,n2) plane and wrap neighbours by n1,n2.
         !$omp parallel do private(j,ii,jj,iii,jjj,cell) reduction(.or.:connected)
         do i=1,n1
            if(connected) cycle
            do j=1,n2
               if(connected) cycle
               cell = box(i,j,1)
               if(cell) then
                  do ii=i-1,i+1
                     if(connected) cycle
                     do jj=j-1,j+1
                        if(connected) cycle
                        iii = per(ii,n1)
                        jjj = per(jj,n2)
                        if(box(iii,jjj,n3)) connected = .true.
                     end do
                  end do
               end if
            end do
         end do
         !$omp end parallel do
      end if

      return

    end function connected

    subroutine rotate(dim)
      implicit none
      integer(ik), intent(IN) :: dim
      ! Circular shift of the whole box by one cell along axis dim. cshift
      ! handles any extents, so this is correct on non-cubic grids; the old
      ! hand-rolled version reused a tmp plane of fixed shape (n1,n2).
      box = cshift(box, 1, int(dim))

      return

    end subroutine rotate

  end subroutine translate_structures

  subroutine read_hdf5_file(n1,n2,n3,field,filename, fieldname)
    use hdf5
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik), intent(out)                             :: n1, n2, n3
    real(sp), dimension(:,:,:), allocatable, intent(out) :: field
    character(len=*), intent(in)                         :: filename, fieldname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: space_id       ! Dataspace identifier
    integer     ::   error ! Error flag
    integer        :: rank
    integer(HSIZE_T), dimension(3) :: data_dims
    integer(HSIZE_T), dimension(3) :: max_dims

    print '(3a)', 'reading file ', trim(filename), '...'

    ! Initialize FORTRAN interface.

    call h5open_f(error)
    if(error < 0) then
       print '(a)', 'error: could not initialise the HDF5 interface.'
       stop 1
    end if

    ! Open an existing file read-only (we never write to it).

    call h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, error)
    if(error < 0) then
       print '(a)', 'error: could not open the HDF5 file (check the path).'
       stop 1
    end if

    ! Open an existing dataset. trim() is essential: a blank-padded name does
    ! not match the dataset and h5dopen_f then leaves field uninitialised, so
    ! the structures would be extracted from garbage.

    call h5dopen_f(file_id, trim(fieldname), dset_id, error)
    if(error < 0) then
       print '(a)', 'error: could not open the requested field (check fieldname).'
       stop 1
    end if

    !Get dataspace ID
    call h5dget_space_f(dset_id, space_id,error)

    ! Reject anything that is not a 3D dataset rather than reading it scrambled.
    call h5sget_simple_extent_ndims_f(space_id, rank, error)
    if(error < 0 .or. rank /= 3) then
       print '(a)', 'error: the field dataset is not 3-dimensional.'
       stop 1
    end if

    !Get dataspace dims

    call h5sget_simple_extent_dims_f(space_id, data_dims, max_dims, error)

    n1 = data_dims(1)
    n2 = data_dims(2)
    n3 = data_dims(3)

    ! Refuse grids the 16-bit coordinate storage cannot represent, instead of
    ! overflowing silently (see max_extent in the types module).
    if(max(n1, n2, n3) > max_extent) then
       print '(a,i0,a)', 'error: grid dimensions above ', max_extent, &
            ' are not supported (coordinates are stored as 16-bit integers).'
       stop 1
    end if

    !Allocate dimensions to dset_data for reading
    allocate(field(n1, n2, n3))

    !Get data
    call h5dread_f(dset_id, H5T_NATIVE_REAL, field, data_dims, error)
    if(error < 0) then
       print '(a)', 'error: could not read the field data.'
       stop 1
    end if

    ! Close the dataspace, dataset and file before closing the interface,
    ! otherwise these handles leak.
    call h5sclose_f(space_id, error)
    call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)

    print '(a)', 'done.'

    return
  end subroutine read_hdf5_file

end module routines

program exstruct
  use types
  use routines
  use data
  implicit none
!!!!!!!!!!!!!
  !Main program
!!!!!!!!!!!!!
  real(rk)                               :: mean,sdev
  integer(ik)                            :: m = 3
  type(interval), pointer                :: intrv_ptr2
  integer(ik)                            :: i,j,k,vol,k1,k2,maxvol,minvol
  integer(ik)                            :: kk
  real(dp)                               :: t1,t2
  real(8), external                      :: omp_get_wtime
  character(len=1024)                    :: fname, field_name
  logical                                :: args_ok
  logical                                :: check_failed

  call parse_args(fname, field_name, m, volume, args_ok)
  ! 'stop' with a message string exits with status 0, so every error path
  ! here uses 'stop 1' (user error) or 'error stop' (failed self-check)
  ! instead: scripts must be able to detect a failed run.
  if(.not. args_ok) stop 1

  call read_hdf5_file(n1,n2,n3,field,fname,field_name)

  ! Worst case the BFS queue holds every interval; intrv() caps each column at
  ! (n3+1)/2 intervals (an alternating column holds ceil(n3/2) of them, one
  ! more than n3/2 when n3 is odd), so this bound is valid for non-cubic grids
  ! too (and is smaller than the old n1**3 for cubic ones).
  allocate(neigh(n1 * n2 * ((n3 + 1) / 2), 3))
  allocate(mask(n1,n2,n3))

  !$omp workshare
  mask = .false.
  !$omp end workshare
  allocate(refmask(n1,n2,n3))
  !$omp workshare
  refmask = .false.
  !$omp end workshare
  allocate(intrv(n1,n2,(n3+1)/2,2))
  allocate(nintrv(n1,n2))

  !$omp parallel do
  do k=1,(n3+1)/2
     do j=1,n2
        do i=1,n1
           intrv(i,j,k,:) = 0_ik
        end do
     end do
  end do
  !$omp end parallel do

  !$omp parallel do
  do j=1,n2
     do i=1,n1
        nintrv(i,j) = 0_ik
     end do
  end do
  !$omp end parallel do

  call mean_sdev(n1,n2,n3,field,mean,sdev)

  print '(a,3es13.5)', 'mean, sdev, cutoff: ', mean, sdev, mean + m * sdev

  !$omp parallel do
  do k=1,n3
     do j=1,n2
        do i=1,n1
           if(field(i,j,k)>mean + m * sdev) refmask(i,j,k) = .true.
        end do
     end do
  end do
  !$omp end parallel do

  print '(a)', 'intervals...'

  call intervals(n1,n2,n3,intrv,nintrv)

  !$omp workshare
  mask = .false.
  !$omp end workshare
  
  check_failed = .false.
  do i=1,n1
     do j=1,n2
        do kk=1,nintrv(i,j)
           k1 = intrv(i,j,kk,1)
           k2 = intrv(i,j,kk,2)
           if(k1 <= k2) then
              do k=k1,k2
                 if(mask(i,j,k)) then
                    print '(a,6i5)', 'error: double point ', i,j,k,k1,k2,nintrv(i,j)
                    check_failed = .true.
                 end if
                 mask(i,j,k) = .true.
              end do
           else 
              do k=k1,n3
                 if(mask(i,j,k)) then
                    print '(a,6i5)', 'error: double point ', i,j,k,k1,k2,nintrv(i,j)
                    check_failed = .true.
                 end if
                 mask(i,j,k) = .true.
              end do
              do k=1,k2
                 if(mask(i,j,k)) then
                    print '(a,6i5)', 'error: double point ', i,j,k,k1,k2,nintrv(i,j)
                    check_failed = .true.
                 end if
                 mask(i,j,k) = .true.
              end do
           end if
        end do
     end do
  end do

  if(check_failed .or. .not.all(mask.eqv.refmask)) then
     error stop 'error: the interval decomposition does not match the mask.'
  end if
  print '(a)', 'Intervals ok.'

  print '(a,i10)', 'intervals: ', sum(nintrv)

  print '(a)', 'ok'

  print '(a)', 'structures...'

  call timing(t1)
  call find_structures
  call timing(t2)
  print '(a,e15.3)', 'time :', t2 - t1

  deallocate(intrv)
  deallocate(cintrv)

  struct_ptr => struct
  maxvol = 0
  minvol = 99999999
  nstruct = 0
  nnintrv = 0
  do while(associated(struct_ptr))
     if(.not.associated(struct_ptr%intrv)) exit
     intrv_ptr => struct_ptr%intrv
     nnintrv = nnintrv + struct_ptr%nintrv
     vol = 0
     do while(associated(intrv_ptr))
        if(any(intrv_ptr%limits==0)) then
           intrv_ptr => intrv_ptr%next
           cycle
        end if
        k1 = intrv_ptr%limits(1)
        k2 = intrv_ptr%limits(2)

        i = intrv_ptr%i
        j = intrv_ptr%j
        if(i /= 0 .and. j /= 0 .and. k1 /= 0 .and. k2 /= 0) then
           if(k1 <= k2) then
              vol = vol + k2 - k1 + 1
           else 
              vol = vol + n3 - k1 + 1
              vol = vol + k2
           end if
        end if
        intrv_ptr => intrv_ptr%next
     end do
     struct_ptr%npoints = vol
     allocate(struct_ptr%points(struct_ptr%npoints,6))
     struct_ptr%points = 0_i2b
     maxvol = max(maxvol,vol)
     minvol = min(minvol,vol)
     struct_ptr => struct_ptr%next
     nstruct = nstruct + 1
  end do

  !print *, nnintrv, sum(nintrv)

  struct_ptr => struct
  open(966,file='gs.dat')
  do while(associated(struct_ptr))
     if(.not.allocated(struct_ptr%points)) exit
     ! i0, not a fixed width: structures beyond 99999 points would print as
     ! '*****' and corrupt gs.dat for downstream scripts.
     write(966,'(t1,i0)') struct_ptr%npoints
     struct_ptr => struct_ptr%next
  end do

  mask(:,:,:) = .false.
  nnintrv = 0
  struct_ptr => struct
  do while(associated(struct_ptr))
     if(.not.associated(struct_ptr%intrv)) exit
     intrv_ptr => struct_ptr%intrv
     kk = 1
     do while(associated(intrv_ptr))
        if(any(intrv_ptr%limits==0)) then
           intrv_ptr => intrv_ptr%next
           cycle
        end if
        k1 = intrv_ptr%limits(1)
        k2 = intrv_ptr%limits(2)
        i = intrv_ptr%i
        j = intrv_ptr%j
        if(i /= 0 .and. j /= 0 .and. k1 /= 0 .and. k2 /= 0) then
           if(k1 <= k2) then
              do k=k1,k2
                 if(.not.mask(i,j,k)) then
                    struct_ptr%points(kk,1) = i
                    struct_ptr%points(kk,2) = j
                    struct_ptr%points(kk,3) = k
                    kk = kk + 1
                    mask(i,j,k) = .true.
                 end if
              end do
           else
              do k=k1,n3
                 if(.not.mask(i,j,k)) then
                    struct_ptr%points(kk,1) = i
                    struct_ptr%points(kk,2) = j
                    struct_ptr%points(kk,3) = k
                    kk = kk + 1
                    mask(i,j,k) = .true.
                 end if
              end do
              do k=1,k2
                 if(.not.mask(i,j,k)) then
                    struct_ptr%points(kk,1) = i
                    struct_ptr%points(kk,2) = j
                    struct_ptr%points(kk,3) = k
                    kk = kk + 1
                    mask(i,j,k) = .true.
                 end if
              end do
           end if
        end if
        intrv_ptr => intrv_ptr%next
        nnintrv = nnintrv + 1
     end do
     ! Every point counted into npoints must have been recorded exactly once:
     ! a shortfall would leave zero rows in points() and break check().
     if(kk - 1 /= struct_ptr%npoints) then
        error stop 'error: point count mismatch during extraction.'
     end if
     struct_ptr => struct_ptr%next
  end do

  print '(a,i10)', 'nnintrv: ', nnintrv

  if(any(refmask.neqv.mask)) then
     error stop 'error: the extracted structures do not cover the mask.'
  end if
  print '(a)', 'extraction ok.'

  struct_ptr => struct
  do while(associated(struct_ptr))
     if(.not.associated(struct_ptr%intrv)) exit
     intrv_ptr => struct_ptr%intrv
     do while(associated(intrv_ptr))
        intrv_ptr2 => intrv_ptr
        intrv_ptr => intrv_ptr%next
        deallocate(intrv_ptr2)
     end do
     ! Drop the dangling head so nothing below dereferences freed memory.
     nullify(struct_ptr%intrv)
     struct_ptr => struct_ptr%next
  end do

  print '(a,i10)', '# structures: ', nstruct
  print '(a,i10)', 'maximum volume: ', maxvol
  print '(a,i10)', 'minimum volume: ', minvol

  print '(a)', 'checking extraction...'
  call check
  print '(a)', 'ok'

  open(985,file='struct.bin',form='unformatted',action='write')
  ! nstruct no longer includes the phantom tail node, so write it as is.
  write(985) nstruct
  struct_ptr => struct
  do while(associated(struct_ptr))
     if(.not.allocated(struct_ptr%points)) exit
     if(struct_ptr%npoints /= 0) then
        write(985) struct_ptr%npoints
        write(985) struct_ptr%points(1:struct_ptr%npoints,1:3)
     end if
     struct_ptr => struct_ptr%next
  end do

  call translate_structures

  call output_structures

  stop

contains

  subroutine timing(t)
    implicit none
    real(dp), intent(OUT) :: t

    ! omp_get_wtime() is callable from serial code: the old parallel region
    ! around it made every thread write the shared t concurrently. Double
    ! precision, because wtime is seconds since an arbitrary epoch and a
    ! single-precision difference of two large values loses the interval.
#ifdef _OPENMP_
    t = omp_get_wtime()
#else
    call cpu_time(t)
#endif

    return

  end subroutine timing

  subroutine parse_args(infile, infield, sdevs, vol, ok)
    implicit none
    character(len=*), intent(out) :: infile, infield
    integer(ik), intent(out)      :: sdevs, vol
    logical, intent(out)          :: ok
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Parse the command line into the input file, field name and thresholds.
    !Recognises -f/--file, -n/--field, -s/--sdevs, -v/--volume and -h/--help,
    !plus the --key=value form. Returns ok=.false. on help or any error.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer             :: nargs, iarg, ios, eqpos
    character(len=1024) :: arg, key, val
    logical             :: have_file, have_field

    ! Defaults.
    infile     = ''
    infield    = ''
    sdevs      = 3_ik
    vol        = 0_ik
    have_file  = .false.
    have_field = .false.
    ok         = .true.

    nargs = command_argument_count()
    if(nargs == 0) then
       call print_usage
       ok = .false.
       return
    end if

    iarg = 1
    do while(iarg <= nargs)
       call get_command_argument(iarg, arg)

       ! Split a --key=value argument; otherwise the value (if any) is the
       ! following argument.
       eqpos = index(arg, '=')
       if(arg(1:2) == '--' .and. eqpos > 0) then
          key = arg(1:eqpos-1)
          val = arg(eqpos+1:)
       else
          key = arg
          val = ''
       end if

       select case(trim(key))
       case('-h', '--help')
          ! Help is a normal exit (status 0), unlike the error paths, which
          ! must exit nonzero so scripts can detect them.
          call print_usage
          stop

       case('-f', '--file')
          call take_value(iarg, nargs, key, val, ok)
          if(.not. ok) return
          infile = val
          have_file = .true.

       case('-n', '--field')
          call take_value(iarg, nargs, key, val, ok)
          if(.not. ok) return
          infield = val
          have_field = .true.

       case('-s', '--sdevs')
          call take_value(iarg, nargs, key, val, ok)
          if(.not. ok) return
          read(val, *, iostat=ios) sdevs
          if(ios /= 0) then
             print '(3a)', 'error: ', trim(key), ' expects an integer.'
             ok = .false.
             return
          end if

       case('-v', '--volume')
          call take_value(iarg, nargs, key, val, ok)
          if(.not. ok) return
          read(val, *, iostat=ios) vol
          if(ios /= 0) then
             print '(3a)', 'error: ', trim(key), ' expects an integer.'
             ok = .false.
             return
          end if

       case default
          print '(2a)', 'error: unknown argument: ', trim(arg)
          print '(a)', "try 'exstruct.exe --help'"
          ok = .false.
          return
       end select

       iarg = iarg + 1
    end do

    ! Validate.
    if(.not. have_file) then
       print '(a)', 'error: no input file given (-f/--file).'
       ok = .false.
       return
    end if

    if(.not. have_field) then
       print '(a)', 'error: no field name given (-n/--field).'
       ok = .false.
       return
    end if

    if(sdevs <= 0_ik) then
       print '(a)', 'error: sdevs (-s/--sdevs) must be a positive integer.'
       ok = .false.
       return
    end if

    if(vol < 0_ik) then
       print '(a)', 'error: volume (-v/--volume) must be non-negative.'
       ok = .false.
       return
    end if

    return

  end subroutine parse_args

  subroutine take_value(iarg, nargs, key, val, ok)
    implicit none
    integer, intent(inout)          :: iarg
    integer, intent(in)             :: nargs
    character(len=*), intent(in)    :: key
    character(len=*), intent(inout) :: val
    logical, intent(out)            :: ok
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Resolve the value for an option: either the text after '=' (already in
    !val) or the next command-line argument, advancing iarg in that case.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ok = .true.
    if(len_trim(val) > 0) return
    if(iarg >= nargs) then
       print '(3a)', 'error: option ', trim(key), ' requires a value.'
       ok = .false.
       return
    end if
    iarg = iarg + 1
    call get_command_argument(iarg, val)

    return

  end subroutine take_value

  subroutine print_usage
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Print a short usage/help message.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print '(a)', 'exstruct - extract high-dissipation structures from an HDF5 field'
    print '(a)', ''
    print '(a)', 'usage: exstruct.exe -f FILE -n FIELD [-s SDEVS] [-v VOLUME]'
    print '(a)', ''
    print '(a)', 'required:'
    print '(a)', '  -f, --file FILE     path to the input HDF5 file'
    print '(a)', '  -n, --field FIELD   name of the 3D dataset to read'
    print '(a)', ''
    print '(a)', 'optional:'
    print '(a)', '  -s, --sdevs N       threshold at mean + N*stddev (positive int, default 3)'
    print '(a)', '  -v, --volume N      output only structures with more than N points (default 0)'
    print '(a)', '  -h, --help          show this help and exit'

    return

  end subroutine print_usage

end program exstruct
