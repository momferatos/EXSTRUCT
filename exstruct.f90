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
  real(rk), parameter    :: LBOX = 2.0_rk * 3.14159_rk
  integer(ik)            :: nnintrv, nnintrv2

  type :: interval
     !!!!!!!!!!!!!!!!!!!
     !Interval data type
     !!!!!!!!!!!!!!!!!!!
     integer(i2b), dimension(1:2) :: limits 
     integer(i2b)                 :: i, j, k
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
     integer(i8b)                              :: npoints
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
  integer(ik)                                   :: MAXNEIGHBOURS = 512,nnn
  real(rk), parameter                           :: PI = 3.1415926535897932384626&
       &4338327950288419716939937510582097494459230781640628620899862803482534&
       &2117068_rk
  integer(ik)                                   :: VOLUME
  integer(ik)                                   :: n1,n2,n3,n4,nstruct
  integer(i2b), dimension(:,:),allocatable      :: neigh
  real(rk), dimension(:,:,:), allocatable       :: field
  integer(i2b), dimension(:,:,:,:), allocatable :: intrv
  integer(ik), dimension(:,:),allocatable       :: nintrv
  logical, dimension(:,:,:), allocatable        :: cintrv,rintrv
  integer(ik)                                   :: vol1,vol2,maxintrv
  integer(ik)                                   :: nnstruct,nintrv2
  logical, dimension(:,:,:), allocatable        :: mask,refmask
  logical, dimension(:,:,:),allocatable         :: box
  real(rk), dimension(:), allocatable           :: L, Lx, Ly, Lz, R, V, A
  real(rk), dimension(:), allocatable           :: ohm,visc,ad,thick,C,totdis
  real(rk), dimension(:), allocatable           :: ux,uy,uz,sigux,siguy,siguz
  type(structure), target                       :: struct
  type(interval), pointer                       :: intrv_ptr => null()
  type(structure), pointer                      :: struct_ptr => null()
  real(rk), dimension(:), allocatable           :: dist, diff
  logical, dimension(:), allocatable            :: diffmask
 
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
    real(rk)                                        :: fac

    fac = real(n1 * n2 * n3,rk)**(-1)

    mean = 0.0_rk
    sdev = 0.0_rk

    !$omp parallel do reduction(+:mean)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             mean = mean + fac * field(i,j,k)
          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do reduction(+:sdev)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             sdev = sdev + fac * (field(i,j,k) - mean)**2
          end do
       end do
    end do
    !$omp end parallel do

    sdev = sqrt(sdev)

    return

  end subroutine mean_sdev

  subroutine check
    implicit none
    integer(ik) :: n,m,i,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Check consistency of structures with initial mask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !$omp workshare
    mask = .false.
    !$omp end workshare

    struct_ptr => struct
    n = 0

    do while(associated(struct_ptr) .and. allocated(struct_ptr%points))
       n = n + 1
       do m=1,struct_ptr%npoints
          i = struct_ptr%points(m,1)
          j = struct_ptr%points(m,2)
          k = struct_ptr%points(m,3)
          if(mask(i,j,k)) print '(a,3i5)', 'error: double point ', i,j,k
          mask(i,j,k) = .true.
       end do
       struct_ptr => struct_ptr%next
    end do

    do k=1,n3
       do j=1,n2
          do i=1,n1
             if(mask(i,j,k) .and. .not.refmask(i,j,k)) then
                print '(a,3i5)', 'error: extra point ',i,j,k
             end if

             if(.not.mask(i,j,k) .and. refmask(i,j,k)) then
                print '(a,3i5)', 'error: missing point ',i,j,k
             end if
          end do
       end do
    end do

    return         

  end subroutine check

  subroutine intervals(n1,n2,n3,field,intrv,nintrv,mean,sdev,m)
    use types
    implicit none
    integer(ik),intent(IN)                       :: n1,n2,n3
    real(rk), dimension(1:n1,1:n2,1:n3)          :: field
    integer(i2b), dimension(1:n1,1:n2,1:n3/2,2)  :: intrv
    integer(ik), dimension(1:n1,1:n2)            :: nintrv
    real(rk)                                     :: mean,sdev
    integer(ik)                                  :: m
!!!!!!!!!!!!!!!!!
    !Set-up intervals
!!!!!!!!!!!!!!!!!
    integer(ik)                                  :: i,j,k,start,end,skip1
    integer(ik)                                  :: skip2,k1,k2,kk
    logical                                      :: openintrv
    integer(i2b), dimension(:,:,:,:),allocatable :: intrv1
    integer(i2b), dimension(:,:),allocatable     :: nintrv1
    logical, dimension(:,:),allocatable          :: contintrv

    allocate(intrv1(n1,n2,n3/2,2))
    allocate(contintrv(n1,n2))
    allocate(nintrv1(n1,n2))
    contintrv(:,:) = .false.
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
    deallocate(contintrv)
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
    integer(ik)                                  :: i,j,k,start,end
    type(interval)                               :: intrv1
    type(structure), allocatable, dimension(:,:) :: structarr


    n4 = maxval(nintrv)
    allocate(cintrv(n1,n2,n4),rintrv(n1,n2,n4))
    cintrv(:,:,:) = .false.
    rintrv(:,:,:) = .false.
    nstruct = 0
    struct_ptr => struct
    allocate(struct%intrv)
    intrv_ptr => struct%intrv
    struct_ptr%nintrv = 0

    do j=1,n2
       do i=1,n1
          do k=1,nintrv(i,j)
             if(.not.cintrv(i,j,k)) then
                intrv1%limits(1) = intrv(i,j,k,1)
                intrv1%limits(2) = intrv(i,j,k,2)
                intrv1%i = i
                intrv1%j = j
                intrv1%k = k
                call find_neighbours(intrv1)
                cintrv(i,j,k) = .true.
                allocate(struct_ptr%next)
                struct_ptr => struct_ptr%next
                nullify(struct_ptr%next)
                nstruct = nstruct + 1
                allocate(struct_ptr%intrv)
                intrv_ptr => struct_ptr%intrv
             end if
          end do
       end do
    end do

    return

  end subroutine find_structures

  subroutine find_neighbours(intrv1)
    use data,only:cintrv,rintrv,nintrv,intrv_ptr,struct_ptr,neigh
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

      if(first<last) then
         first = first + 1
         ii = neigh(first,1)
         jj = neigh(first,2)
         kk = neigh(first,3)
         if(ii /= 0 .and. jj /= 0 .and. kk /= 0) then
            intrv1%limits(1) = intrv(ii,jj,kk,1)
            intrv1%limits(2) = intrv(ii,jj,kk,2)
            intrv1%i = ii
            intrv1%j = jj
            intrv1%k = kk
            lcond = .true.
            lcycle = .true.
         else
            lcond = .false.
            lcycle = .true.
         end if
      else
         lexit = .true.
      end if

      return

    end subroutine sub1

    subroutine sub2
      implicit none

      !rintrv(i,j,k) = .true.
      if(.not.cintrv(i,j,k)) then
         intrv_ptr%limits(1) = intrv(i,j,k,1)
         intrv_ptr%limits(2) = intrv(i,j,k,2)
         intrv_ptr%i = i
         intrv_ptr%j = j
         intrv_ptr%k = k
         struct_ptr%nintrv = struct_ptr%nintrv + 1
         allocate(intrv_ptr%next)
         intrv_ptr => intrv_ptr%next
         nullify(intrv_ptr%next)
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

               intrv_ptr%limits(1) = intrv(ii,jj,kk,1)
               intrv_ptr%limits(2) = intrv(ii,jj,kk,2)

               intrv_ptr%i = ii
               intrv_ptr%j = jj
               intrv_ptr%k = kk
               allocate(intrv_ptr%next)
               intrv_ptr => intrv_ptr%next
               nullify(intrv_ptr%next)
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

  end subroutine find_neighbours

  subroutine find_neighbours3(intrv1)
    use data,only:cintrv,rintrv,nintrv,intrv_ptr,struct_ptr,neigh
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
    first = 0
    last = 0
    cond = .true.
    do while(.true.)
       print *, first, last

       if(cond) then
          i = intrv1%i
          j = intrv1%j
          k = intrv1%k
          limits = intrv1%limits
       end if

       if(i==0 .or. j==0 .or. k==0&
            & .or. i>n1 .or. j>n2 .or. k>n3) then
          if(first<last) then
             first = first + 1
             ii = neigh(first,1)
             jj = neigh(first,2)
             kk = neigh(first,3)
             if(ii /= 0 .and. jj /= 0 .and. kk /= 0) then
                intrv1%limits(1) = intrv(ii,jj,kk,1)
                intrv1%limits(2) = intrv(ii,jj,kk,2)
                intrv1%i = ii
                intrv1%j = jj
                intrv1%k = kk
                cond = .true.
                cycle
             else
                cond = .false.
                cycle
             end if
          else
             exit
          end if
       end if

       !rintrv(i,j,k) = .true.
       if(.not.cintrv(i,j,k)) then
          intrv_ptr%limits(1) = intrv(i,j,k,1)
          intrv_ptr%limits(2) = intrv(i,j,k,2)
          intrv_ptr%i = i
          intrv_ptr%j = j
          intrv_ptr%k = k
          struct_ptr%nintrv = struct_ptr%nintrv + 1
          allocate(intrv_ptr%next)
          intrv_ptr => intrv_ptr%next
          nullify(intrv_ptr%next)
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

                intrv_ptr%limits(1) = intrv(ii,jj,kk,1)
                intrv_ptr%limits(2) = intrv(ii,jj,kk,2)

                intrv_ptr%i = ii
                intrv_ptr%j = jj
                intrv_ptr%k = kk
                allocate(intrv_ptr%next)
                intrv_ptr => intrv_ptr%next
                nullify(intrv_ptr%next)
                cintrv(ii,jj,kk) = .true.

                last = last + 1
                neigh(last,1) = ii
                neigh(last,2) = jj
                neigh(last,3) = kk

             end if
          end do
       end do

       if(first<last) then
          first = first + 1
          ii = neigh(first,1)
          jj = neigh(first,2)
          kk = neigh(first,3)
          if(ii /= 0 .and. jj /= 0 .and. kk /= 0) then
             intrv1%limits(1) = intrv(ii,jj,kk,1)
             intrv1%limits(2) = intrv(ii,jj,kk,2)
             intrv1%i = ii
             intrv1%j = jj
             intrv1%k = kk
             cond = .true.
             cycle
          else
             cond = .false.
             cycle
          end if
       else
          exit
       end if
    end do

    return

  end subroutine find_neighbours3

  subroutine find_neighbours2(intrv1)
    use data,only:cintrv,rintrv,nintrv,intrv_ptr,struct_ptr,neigh
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

    first = 0
    last = 0

100 continue

    i = intrv1%i
    j = intrv1%j
    k = intrv1%k
    limits = intrv1%limits

    if(i==0 .or. j==0 .or. k==0&
         & .or. i>n1 .or. j>n2 .or. k>n3) then
       !first=first + 1
       goto 200
    end if

    !rintrv(i,j,k) = .true.
    if(.not.cintrv(i,j,k)) then
       intrv_ptr%limits(1) = intrv(i,j,k,1)
       intrv_ptr%limits(2) = intrv(i,j,k,2)
       intrv_ptr%i = i
       intrv_ptr%j = j
       intrv_ptr%k = k
       struct_ptr%nintrv = struct_ptr%nintrv + 1
       allocate(intrv_ptr%next)
       intrv_ptr => intrv_ptr%next
       nullify(intrv_ptr%next)
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

             intrv_ptr%limits(1) = intrv(ii,jj,kk,1)
             intrv_ptr%limits(2) = intrv(ii,jj,kk,2)

             intrv_ptr%i = ii
             intrv_ptr%j = jj
             intrv_ptr%k = kk
             allocate(intrv_ptr%next)
             intrv_ptr => intrv_ptr%next
             nullify(intrv_ptr%next)
             cintrv(ii,jj,kk) = .true.

             last = last + 1
             neigh(last,1) = ii
             neigh(last,2) = jj
             neigh(last,3) = kk

          end if
       end do
    end do

200 if(first<last) then
       first = first + 1
       ii = neigh(first,1)
       jj = neigh(first,2)
       kk = neigh(first,3)
       if(ii /= 0 .and. jj /= 0 .and. kk /= 0) then
          intrv1%limits(1) = intrv(ii,jj,kk,1)
          intrv1%limits(2) = intrv(ii,jj,kk,2)
          intrv1%i = ii
          intrv1%j = jj
          intrv1%k = kk
          goto 100
       else
          goto 200
       end if
    end if

    return

  end subroutine find_neighbours2

  function is_neighbour(i1,i2)
    use types
    implicit none
    integer(i2b), dimension(2),intent(IN) :: i1,i2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Check whether two intervals are neighbours
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical                               :: is_neighbour
    integer(ik)                           :: start1,start2,end1,end2
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
    integer(ik)    :: n,i,j,k,nfile,m
    character(128) :: filename

    nfile = 0
    struct_ptr => struct

    do while(associated(struct_ptr) .and. allocated(struct_ptr%points))
       if(struct_ptr%npoints > VOLUME) then
          write(filename,'(a,i0,a)') 'out.',nfile,'.vtk'
          call write_vtk_file(trim(filename),struct_ptr)
          nfile = nfile + 1
       end if
       struct_ptr => struct_ptr%next
    end do

    print '(a,i4,a)', 'Wrote ', nfile, ' files.'

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
    integer(ik)              :: i,j,k,m,vtk_file = 999
    integer(ik)              :: tmp
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

    close(vtk_file,status='keep')

    return

  end subroutine write_vtk_file

  !#end pass

  subroutine surface_area(astruct,nsurf)
    use types
    use data,only:refmask
    implicit none
    type(structure), intent(IN), pointer :: astruct
    integer(ik), intent(OUT)             :: nsurf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate surface area of structrure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, dimension(3,3,3)            :: box
    integer(ik)                          :: i,j,k,m,ii,jj,kk

    nsurf = 0

    do m=1,astruct%npoints
       i = astruct%points(m,1)
       j = astruct%points(m,2)
       k = astruct%points(m,3)
       box = .false.
       do ii=1,3 ; do jj=1,3 ; do kk=1,3
          box(ii,jj,kk) = refmask(per(i+ii-2,n1),per(j+jj-2,n2),per(k+kk-2,n3))
       enddo; enddo ; enddo
       if (.not.all(box)) nsurf = nsurf + 1
    enddo

    return

  end subroutine surface_area

  function in_box(n)
    use types
    use data,only:nnn
    implicit none
    integer(i2b), intent(IN) :: n
    integer(ik)              :: in_box
!!!!!!!!!!!!!!!!!!!!!!!!
    !???????????????????????
!!!!!!!!!!!!!!!!!!!!!!!!

    if(n<1) then
       in_box = n + nnn
    else if(n>nnn) then
       in_box = n - nnn
    else
       in_box = n
    end if

    return

  end function in_box

  subroutine statistics
    use types
    use data
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate structure statistics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                           :: n,m,i,j,k,ix,jx,kx,nsurf
    real(rk)                              :: dx,dv,ds,tmpohm,tmpvisc,tmpad
    real(rk)                              :: distx,disty,distz,distl,disttot
    real(rk)                              :: tmpux,tmpuy,tmpuz,tmpsigux
    real(rk)                              :: tmpsiguy,tmpsiguz
    real(rk), dimension(:,:), allocatable :: x

    allocate(L(1:nstruct))
    allocate(Lx(1:nstruct))
    allocate(Ly(1:nstruct))
    allocate(Lz(1:nstruct))
    allocate(R(1:nstruct))
    allocate(V(1:nstruct))
    allocate(A(1:nstruct))
    allocate(thick(1:nstruct))
    allocate(C(1:nstruct))
    allocate(ohm(1:nstruct))
    allocate(visc(1:nstruct))
    allocate(ad(1:nstruct))
    allocate(ux(1:nstruct))
    allocate(uy(1:nstruct))
    allocate(uz(1:nstruct))
    allocate(sigux(1:nstruct))
    allocate(siguy(1:nstruct))
    allocate(siguz(1:nstruct))
    dx = LBOX/(n1 - 1)
    ds = dx**2
    dv = dx**3



    n = 1
    struct_ptr => struct
    do while(associated(struct_ptr) .and. allocated(struct_ptr%points))
       allocate(x(struct_ptr%npoints,3))
       x(:,:) = 0.0_rk
       nsurf = 0
       call surface_area(struct_ptr,nsurf)
       tmpohm = 0.0_rk
       tmpvisc = 0.0_rk
       tmpad = 0.0_rk
       tmpux = 0
       tmpuy = 0
       tmpuz = 0
       !$omp parallel do private(i,j,k,ix,jx,kx) reduction(+:tmpohm,tmpvisc,&
       !$omp tmpad,tmpux,tmpuy,tmpuz,tmpsigux,tmpsiguy,tmpsiguz)
       do m=1,struct_ptr%npoints
          ix = struct_ptr%points(m,4)
          jx = struct_ptr%points(m,5)
          kx = struct_ptr%points(m,6)

          i = struct_ptr%points(m,1)
          j = struct_ptr%points(m,2)
          k = struct_ptr%points(m,3)
          x(m,1) = dx * (ix - 1)
          x(m,2) = dx * (jx - 1)
          x(m,3) = dx * (kx - 1)
!!$          tmpohm = tmpohm + dv * j2(i,j,k)
!!$          tmpvisc = tmpvisc + dv * e(i,j,k)
!!$          tmpad = tmpad + dv * jxb2(i,j,k)
!!$          tmpux = tmpux + u1(i,j,k)
!!$          tmpuy = tmpuy + u2(i,j,k)
!!$          tmpuz = tmpuz + u3(i,j,k)
       end do
       !$omp end parallel do
       if(struct_ptr%npoints /= 0) then
          tmpux = tmpux/struct_ptr%npoints
          tmpuy = tmpuy/struct_ptr%npoints
          tmpuz = tmpuz/struct_ptr%npoints
       else
          tmpux = 0.
          tmpuy = 0.
          tmpuz = 0.
       end if


       tmpsigux = 0.0
       tmpsiguy = 0.0
       tmpsiguz = 0.0
       do m=1,struct_ptr%npoints
!!$          tmpsigux = tmpsigux + (u1(i,j,k) - tmpux)**2
!!$          tmpsiguy = tmpsiguy + (u2(i,j,k) - tmpuy)**2
!!$          tmpsiguz = tmpsiguz + (u3(i,j,k) - tmpuz)**2
       end do
       tmpsigux = sqrt(tmpsigux/struct_ptr%npoints)
       tmpsiguy = sqrt(tmpsiguy/struct_ptr%npoints)
       tmpsiguz = sqrt(tmpsiguz/struct_ptr%npoints)

       distx = 0.0_rk
       disty = 0.0_rk
       distz = 0.0_rk
       distl = 0.0_rk
       !$omp parallel do reduction(max:distx,disty,distz,distl) private(disttot,n)
       do i=1,struct_ptr%npoints
          do j=i + 1,struct_ptr%npoints
             distx = max(abs(x(i,1) - x(j,1)),distx)
             disty = max(abs(x(i,2) - x(j,2)),disty)
             distz = max(abs(x(i,3) - x(j,3)),distz)
             disttot = sqrt((x(i,1) - x(j,1))**2 + (x(i,2) - x(j,2))**2 + &
                  &(x(i,3) - x(j,3))**2)
             distl = max(distl,disttot)
          end do
       end do
       !$omp end parallel do
       if(distx==0.0_rk) distx = dx
       if(disty==0.0_rk) disty = dx
       if(distz==0.0_rk) distz = dx
       if(distl==0.0_rk) distl = dx
       L(n) = distl
       Lx(n) = distx
       Ly(n) = disty
       Lz(n) = distz
       R(n) = sqrt(distx**2 + disty**2 + distz**2)
       V(n) = dv * struct_ptr%npoints
       A(n) = ds * nsurf
       thick(n) = V(n)/(A(n)/2._rk)
       C(n) = PI * (L(n)/2._rk)**2/(A(n)/2._rk)
       ohm(n) = tmpohm
       visc(n) = tmpvisc
       ad(n) = tmpad
       ux(n) = tmpux
       uy(n) = tmpuy
       uz(n) = tmpuz
       sigux(n) = tmpsigux
       siguy(n) = tmpsiguy
       siguz(n) = tmpsiguz
       deallocate(x)
       struct_ptr => struct_ptr%next
       n = n + 1
    end do

  end subroutine statistics


  subroutine compute_pdf_hist(npoints,dat,nbins,fname)
    implicit none
    integer(ik), intent(IN)                  :: npoints
    real(rk), dimension(npoints), intent(IN) :: dat
    integer(ik), intent(IN)                  :: nbins
    character(*), intent(IN)                 :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate histogram of structures' statistics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                              :: i,j,nn,n,ii
    real(rk)                                 :: bw,dmax,dmin,vpmax,vpmin
    real(rk)                                 :: accelmax,xdens,acceldens
    real(rk)                                 :: vpdens,d,dx,c1,c2,c3,dir
    real(rk)                                 :: area,xmin,xmax
    real(rk), dimension(3)                   :: tmpx,tmpvp,tmpaccel
    real(rk), dimension(:,:), allocatable    :: pdf
    integer                                  :: pdf_file = 99
    print '(2a)', 'pdf: ', trim(fname)
    allocate(pdf(nbins,3))

    dmax = 0.0_rk
    dmin = 1.0e11_rk
    !$omp parallel do private(d) reduction(max:dmax) reduction(min:dmin)
    do n=1,npoints
       if(dat(n)>0.0) then
          d = dat(n)
          dmax = max(dmax,d)
          dmin = min(dmin,d)
       end if
    end do
    !$omp end parallel do

    dx = (log(dmax) - log(dmin))/real(nbins)
    xmin = log(dmin)
    pdf(:,:) = 0.0
    !$omp parallel do private(ii)
    do i=1,npoints
       if(dat(i)>0.0) then
          ii = int((log(dat(i)) - xmin)/dx) + 1
          if(ii==nbins + 1) ii = nbins
          if(ii<1 .or. ii>nbins) print '(a,2i5)', 'histogram error:', ii,xmin
          pdf(ii,2) = pdf(ii,2) + 1
       end if
    end do
    !$omp end parallel do

    dx = (log(dmax) - log(dmin))/real(nbins)
    !print '(a,3f10.3)', fname, dmin,dmax,dx
    !$omp parallel do private(xmin,xmax)
    do ii=1,nbins
       xmin = exp(log(dmin) + (ii - 1) * dx)
       xmax = exp(log(dmin) + (ii * dx))
       pdf(ii,1) = 0.5 * (xmin + xmax)
       pdf(ii,2) = pdf(ii,2)/(xmax - xmin)
    end do
    !$omp end parallel do

    area = 0.0
    !$omp parallel do private(dx) reduction(+:area)
    do i=1,nbins-1
       dx = pdf(i + 1,1) - pdf(i,1)
       area = area + (pdf(i + 1,2) + pdf(i,2)) * dx/2
    end do
    !$omp end parallel do

    if(area /= 0.0) then
       !$omp parallel do
       do i=1,nbins
          pdf(i,2) = pdf(i,2)/area
          pdf(i,3) = pdf(i,3)/area
       end do
       !$omp end parallel do
    end if

    open(pdf_file,file=fname,form='formatted',action='write')

    do i=1,nbins
       write(pdf_file,*) pdf(i,1:3)
    end do

    close(pdf_file,status='keep')

    deallocate(pdf)

    return

  end subroutine compute_pdf_hist

  subroutine output_statistics
    implicit none
    integer(ik) :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Write structure statistics to file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open(777,file='stats.dat',action='write',form='formatted')
    do i=1,nstruct
       write(777,'(9e30.15)') V(i),A(i),L(i),R(i),ad(i),ohm(i),visc(i),&
            &thick(i),C(i)
    end do

    open(778,file='stats-new.dat',action='write',form='formatted')
    do i=1,nstruct
       write(778,'(11e30.15)') V(i),A(i),L(i),R(i),ad(i),ohm(i),visc(i),&
            &thick(i),C(i),ux(i),uy(i),uz(i),sigux(i),siguy(i),siguz(i)
    end do

    close(777,status='keep')

    return

  end subroutine output_statistics

  subroutine translate_structures
    use types
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Translate structures so that they are continuous and not disrupted by
    !periodic boundaries!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, dimension(:,:,:), allocatable :: box,box2
    logical, dimension(:,:), allocatable   :: tmp
    integer(ik)                            :: m,i,j,k,n,step_i,step_j,step_k

    allocate(box2(1:n1,1:n2,1:n3))
    allocate(box(1:n1,1:n2,1:n3))
    allocate(tmp(1:n1,1:n2))
    struct_ptr => struct

    n = 0
    do while(associated(struct_ptr) .and. allocated(struct_ptr%points))
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

      connected = .false.

      if(dim==1_ik) then
         ! i-face: does the box(1,:,:) plane connect to box(n1,:,:)?
         ! scan the (n2,n3) plane and wrap neighbours by n2,n3.
         !$omp parallel do private(j,ii,jj,iii,jjj,cell) shared(connected)
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
         !$omp parallel do private(j,ii,jj,iii,jjj,cell) shared(connected)
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
         !$omp parallel do private(j,ii,jj,iii,jjj,cell) shared(connected)
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
      box2 = box
      if(dim==1) then
         !$omp workshare
         tmp(:,:) = box2(1,:,:)
         box(1:n1 - 1,:,:) = box2(2:n1,:,:)
         box(n1,:,:) = tmp(:,:)
         !$omp end workshare
      else if(dim==2) then
         !$omp workshare
         tmp(:,:) = box2(:,1,:)
         box(:,1:n2 - 1,:) = box2(:,2:n2,:)
         box(:,n2,:) = tmp(:,:)
         !$omp end workshare
      else if(dim==3) then
         !$omp workshare
         tmp(:,:) = box2(:,:,1)
         box(:,:,1:n3 - 1) = box2(:,:,2:n3)
         box(:,:,n3) = tmp(:,:)
         !$omp end workshare
      end if

      return

    end subroutine rotate

  end subroutine translate_structures

  function newpoint(points,n,i,j,k)
    implicit none
    integer(i2b), dimension(:,:),intent(IN) :: points
    integer(ik)                             :: i,j,k,n
    logical                                 :: newpoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Check whether we have a new point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                             :: ii

    newpoint = .true.
    do ii=1,n
       if(all(points(ii,:)==(/i,j,k/))) then
          newpoint = .false.
          return
       end if
    end do

    return

  end function newpoint

  subroutine read_hdf5_file(n1,n2,n3,field,filename)
    use hdf5
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik), intent(out)                             :: n1, n2, n3
    real(sp), dimension(:,:,:), allocatable, intent(out) :: field
    character(len=*), intent(in)                         :: filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=4), parameter :: dsetname = "diss"     ! Dataset name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: space_id       ! Dataspace identifier
    integer(HID_T) :: dtype_id       ! Dataspace identifier
    integer     ::   error ! Error flag    
    integer(HSIZE_T), dimension(3) :: data_dims
    integer(HSIZE_T), dimension(3) :: max_dims                  


    print '(3a)', 'reading file ', trim(filename), '...'



    ! Initialize FORTRAN interface.

    call h5open_f(error)


    ! Open an existing file.

    call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)


    ! Open an existing dataset.

    call h5dopen_f(file_id, dsetname, dset_id, error)


    !Get dataspace ID
    call h5dget_space_f(dset_id, space_id,error)


    !Get dataspace dims

    call h5sget_simple_extent_dims_f(space_id, data_dims, max_dims, error)

    n1 = data_dims(1)
    n2 = data_dims(2)
    n3 = data_dims(3)

    !Allocate dimensions to dset_data for reading
    allocate(field(n1, n2, n3))


    !Get data
    call h5dread_f(dset_id, H5T_NATIVE_REAL, field, data_dims, error)

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
  use data,only:nnn
  implicit none
!!!!!!!!!!!!!
  !Main program
!!!!!!!!!!!!!
  real(rk)                               :: mean,sdev
  integer(ik)                            :: m = 3
  type(interval), pointer                :: intrv_ptr2
  type(structure),pointer                :: struct_ptr2
  integer(ik)                            :: i,j,k,vol,k1,k2,maxvol,minvol,n
  integer(ik)                            :: kk,none,which,nppdf,nstruct2
  real(rk)                               :: t1,t2,tmp,totohm,totad,totvisc
  real(rk)                               :: sohm,svisc,sad,tottot,stot,svol
  character(256)                         :: path,fmt
  real(8), external                      :: omp_get_wtime
  integer(i2b),  dimension(:,:), pointer :: point_array
  character(len=1024)                    :: fname, strdevs, strvolume


  if(command_argument_count() /= 3) then
     print '(a)', 'usage: exstruct filename sdevs volume'
     stop
  endif

  call get_command_argument(1, fname)
  call get_command_argument(2, strdevs)
  call get_command_argument(3, strvolume)


  read(strdevs, *, err=100) m
  read(strvolume, *, err=200) volume
  if(m <= 0_ik) then
     print '(a)', 'error: argument sdev should be a positive integer greater than zero.'
     stop
  end if

  if(volume < 0_ik) then
     print '(a)', 'error: argument volume should be a positive integer.'
     stop
  end if
  
  call read_hdf5_file(n1,n2,n3,field,fname)


  nnn = n1
  allocate(neigh(nnn**3,3))
  neigh = 0
  allocate(mask(n1,n2,n3))

  !$omp workshare
  mask = .false.
  !$omp end workshare
  allocate(refmask(n1,n2,n3))
  !$omp workshare
  refmask = .false.
  !$omp end workshare
  allocate(intrv(n1,n2,n3/2,2))
  allocate(nintrv(n1,n2))

  allocate(dist(n1 * n2 * n3))
  allocate(diff(n1 * n2 * n3))
  allocate(diffmask(n1 * n2 * n3))

  !$omp parallel do
  do k=1,n3/2
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

  call intervals(n1,n2,n3,field,intrv,nintrv,mean,sdev,m)

  !$omp workshare
  mask = .false.
  !$omp end workshare
  
  do i=1,n1
     do j=1,n2
        do kk=1,nintrv(i,j)
           k1 = intrv(i,j,kk,1)
           k2 = intrv(i,j,kk,2)
           if(k1 <= k2) then
              do k=k1,k2
                 if(mask(i,j,k)) print '(a,6i5)', 'error: double point ', i,j,k,k1,k2,nintrv(i,j)
                 mask(i,j,k) = .true.
              end do
           else 
              do k=k1,n3
                 if(mask(i,j,k)) print '(a,6i5)', 'error: double point ', i,j,k,k1,k2,nintrv(i,j)
                 mask(i,j,k) = .true.
              end do
              do k=1,k2
                 if(mask(i,j,k)) print '(a,6i5)', 'error: double point ', i,j,k,k1,k2,nintrv(i,j)
                 mask(i,j,k) = .true.
              end do
           end if
        end do
     end do
  end do

  if(all(mask.eqv.refmask)) print '(a)', 'Intervals ok.'

  print '(a,i10)', 'intervals: ', sum(nintrv)

  print '(a)', 'ok'

  print '(a)', 'structures...'

  call timing(t1)
  call find_structures
  call timing(t2)
  print '(a,e15.3)', 'time :', t2 - t1

  deallocate(intrv)
  deallocate(cintrv)
  deallocate(rintrv)

  struct_ptr => struct
  maxvol = 0
  minvol = 99999999
  nstruct = 0
  n = 0
  nnintrv = 0
  do while(associated(struct_ptr) .and. associated(struct_ptr%intrv))
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
  vol = 0
  open(966,file='gs.dat')
  do while(associated(struct_ptr) .and. allocated(struct_ptr%points))
     write(966,'(t1,i5)') struct_ptr%npoints
     vol = vol + struct_ptr%npoints
     struct_ptr => struct_ptr%next
  end do

  mask(:,:,:) = .false.
  nnintrv = 0
  struct_ptr => struct
  n = 0
  do while(associated(struct_ptr) .and. associated(struct_ptr%intrv))
     n = n + 1
     intrv_ptr => struct_ptr%intrv
     kk = 1
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
              do k=k1,k2
                 if(i /= 0 .and. j /= 0 .and. k /= 0 .and. .not.mask(i,j,k)) then
                    struct_ptr%points(kk,1) = i
                    struct_ptr%points(kk,2) = j
                    struct_ptr%points(kk,3) = k
                    kk = kk + 1
                    vol = vol + 1
                    if(mask(i,j,k)) print '(a,4i5)', 'error: double point ',i,j,k,n
                    mask(i,j,k) = .true.
                 end if
              end do
           else 
              do k=k1,n3
                 if(i /= 0 .and. j /= 0 .and. k /= 0 .and. .not.mask(i,j,k)) then
                    struct_ptr%points(kk,1) = i
                    struct_ptr%points(kk,2) = j
                    struct_ptr%points(kk,3) = k
                    kk = kk + 1
                    vol = vol + 1
                    if(mask(i,j,k)) print '(a,4i5)','error: double point ', i,j,k,n
                    mask(i,j,k) = .true.
                 end if
              end do
              do k=1,k2
                 if(i /= 0 .and. j /= 0 .and. k /= 0 .and. .not.mask(i,j,k)) then
                    struct_ptr%points(kk,1) = i
                    struct_ptr%points(kk,2) = j
                    struct_ptr%points(kk,3) = k
                    kk = kk + 1
                    vol = vol + 1
                    if(mask(i,j,k)) print '(a,4i5)','error: double point ', i,j,k,n
                    mask(i,j,k) = .true.
                 end if
              end do
           end if
        end if
        intrv_ptr => intrv_ptr%next
        nnintrv = nnintrv + 1
        n = n + 1
     end do
     struct_ptr => struct_ptr%next
  end do

  print '(a,i10)', 'nnintrv: ', nnintrv

  if(any(refmask.neqv.mask)) then
     print '(a)', 'extraction error.'
  else
     print '(a)', 'extraction ok.'
  end if

  struct_ptr => struct
  do while(associated(struct_ptr) .and. associated(struct_ptr%intrv))
     intrv_ptr => struct_ptr%intrv
     do while(associated(intrv_ptr))
        intrv_ptr2 => intrv_ptr
        intrv_ptr => intrv_ptr%next
        deallocate(intrv_ptr2)
     end do
     struct_ptr => struct_ptr%next
  end do

  print '(a,i10)', '# structures: ', nstruct
  print '(a,i10)', 'maximum volume: ', maxvol
  print '(a,i10)', 'minimum volume: ', minvol

  print '(a)', 'checking extraction...'
  call check
  print '(a)', 'ok'

  open(985,file='struct.bin',form='unformatted',action='write')
  write(985) nstruct - 1
  struct_ptr => struct
  do while(associated(struct_ptr) .and. associated(struct_ptr%intrv))
     if(struct_ptr%npoints /= 0) then
        write(985) struct_ptr%npoints
        write(985) struct_ptr%points(1:struct_ptr%npoints,1:3)
     end if
     struct_ptr => struct_ptr%next
  end do

  call translate_structures

  call output_structures

  stop

100 print '(a)', 'error: argument sdev should be a positive integer greater than zero.'
  stop
200 print '(a)', 'error: argument volume should be a positive integer.'
  stop

contains

  subroutine timing(t)
    implicit none
    real(rk), intent(OUT) :: t

#ifdef _OPENMP_
    !$omp parallel
    t = omp_get_wtime()
    !$omp end parallel
#else
    call cpu_time(t)
#endif

    return

  end subroutine timing

end program exstruct
