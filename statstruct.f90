module types
  implicit none
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer, parameter     :: i1b=1 !integer 1 bit
  integer, parameter     :: i2b=2 !integer 2 bits
  integer, parameter     :: i4b=4 !integer 4 bits
  integer, parameter     :: i8b=8 !integer 8 bits
  integer, parameter     :: ik=i8b!main integer type
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer, parameter     :: sp=4  !single precision real
  integer, parameter     :: dp=8  !double precision real
  integer, parameter     :: qp=16 !quadruple precision real
  integer, parameter     :: rk=sp !main real type
  integer, parameter     :: rks=sp!real type for arrays
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer, parameter     :: csp=4  !single precision complex
  integer, parameter     :: cdp=8  !double precision real
  integer, parameter     :: ck=cdp !main complex type
  integer, parameter     :: cks=rks!complex type for arrays
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer(ik), parameter :: inrk=sp
  integer(ik), parameter :: outrk=dp
  
  real(rk), parameter :: LBOX=2.0_rk*3.14159_rk
  integer(ik) :: nnintrv,nnintrv2


  type :: struct
     integer(i2b), dimension(:,:),allocatable :: points
  end type struct
  
  type :: arr_pointer
     type(struct), pointer :: struct_ptr
     integer(i8b) :: npoints
  end type arr_pointer


  
end module types

module data
  use types
  implicit none
  integer(ik) :: MAXNEIGHBOURS=512,nnn
  real(rk), parameter    :: PI=3.14159265358979323846264338327950&
       &2884197169399375105820974944592307816406286208998628034825&
       &342117068_rk
  integer(ik) :: VOLUME
  integer(ik) :: n1,n2,n3,n4,nstruct
  real(rk) :: ADCOEFF,ETA
  integer(i2b), dimension(:,:),allocatable :: neigh
  real(rk), dimension(:,:,:), allocatable :: j2,e,jxb2,total,u1,u2,u3
  integer(i2b), dimension(:,:,:,:), allocatable :: intrv
  integer(ik), dimension(:,:),allocatable :: nintrv
  logical, dimension(:,:,:), allocatable :: cintrv,rintrv
  integer(ik) :: vol1,vol2,maxintrv,nnstruct,nintrv2
  logical, dimension(:,:,:), allocatable :: mask,refmask
  logical, dimension(:,:,:),allocatable :: box
  real(rk), dimension(:), allocatable :: L, Lx, Ly, Lz, R, V, A, ohm,&
       &visc,ad,thick,C,totdis,ux,uy,uz,sigux,siguy,siguz
  real(rk), dimension(:), allocatable :: dist, diff
  logical, dimension(:), allocatable :: diffmask
  type(arr_pointer), dimension(:), allocatable :: structs
  integer(i2b), dimension(:,:,:), allocatable :: structures

end module data

module sort
  use types
  implicit none
  ! Recursive Fortran 95 quicksort routine
  ! sorts real numbers into ascending numerical order
  ! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
  ! Based on algorithm from Cormen et al., Introduction to Algorithms,
  ! 1997 printing

  ! Made F conformant by Walt Brainerd
  public :: Qsort
  private :: Partition

contains

  recursive subroutine Qsort(A,n)
    integer(ik), intent(in) :: n
    real(rk), intent(in out), dimension(1:) :: A
    
    integer(ik) :: iq

    if(n > 1) then
       call Partition(A, iq,n)
       call Qsort(A(:iq-1),iq-1)
       call Qsort(A(iq:),n-iq)
    endif
  end subroutine Qsort

  subroutine Partition(A, marker,n)
    integer(ik), intent(in) :: n
    real(rk), intent(in out), dimension(1:) :: A
    integer(ik), intent(out) :: marker
    integer(ik) :: i, j
    real(rk) :: temp
    real(rk) :: x      ! pivot point
    x = A(1)
    i= 0
    j= n + 1

    do
       j = j-1
       do
          if (A(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine Partition

end module sort

module routines
  use types
  use data
  implicit none

contains

  function per(n)
    implicit none
    integer(ik) :: per,n
    if(n==n1+1) then
       per=1
    else if(n==0) then
       per=n1
    else
       per=n
    end if
    return
  end function per

  subroutine read_field(path,n1,n2,n3,j2,e,jxb2,u1,u2,u3)
    use types
    implicit none
    character(*), intent(IN) :: path
    integer(ik), intent(IN) :: n1,n2,n3
    real(rk), dimension(1:n1,1:n2,1:n3) :: j2,e,jxb2,u1,u2,u3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k,nnewl
    real(sp) :: jj1,jj2,jj3
    character(512) :: infile
    character :: c,newline

    j=0.0
    e=0.0
    jxb2=0.0

    newline=char(10)

    print *, 'j00000.vtk...'
    write(infile,'(2a)') trim(path), '/j00000.vtk'
    open(777,file=infile,access='stream',&
         &action='read',status='old')
    nnewl=0

    do while(nnewl<9)
       read(777) c
       if(c==newline) nnewl=nnewl+1
    end do
    do k=1,n3
       do j=1,n2
          do i=1,n1
             read(777) jj1,jj2,jj3
             j2(i,j,k)=ETA*(jj1**2+jj2**2+jj3**2)
          end do
       end do
    end do
    close(777,status='keep')



    print *, 'e00000.vtk...'
    write(infile,'(2a)') trim(path), '/e00000.vtk'
    open(777,file=infile,access='stream',&
         &action='read',status='old')
    nnewl=0
    do while(nnewl<10)
       read(777) c
       if(c==newline) nnewl=nnewl+1
    end do
    do k=1,n3
       do j=1,n2
          do i=1,n1
             read(777) jj1
             e(i,j,k)=jj1
          end do
       end do
    end do
    close(777,status='keep')



    print *, 'jxb00000.vtk...'
    write(infile,'(2a)') trim(path), '/jxb00000.vtk'
    open(777,file=infile,access='stream',&
         &action='read',err=100,status='old')
    nnewl=0
    do while(nnewl<9)
       read(777) c
       if(c==newline) nnewl=nnewl+1
    end do
    do k=1,n3
       do j=1,n2
          do i=1,n1
             read(777) jj1,jj2,jj3
             jxb2(i,j,k)=ADCOEFF*(jj1**2+jj2**2+jj3**2)
          end do
       end do
    end do
    close(777,status='keep')

!!$    print *, 'u00000.vtk...'
!!$    write(infile,'(2a)') trim(path), '/u00000.vtk'
!!$    open(777,file=infile,access='stream',&
!!$         &action='read',err=100,status='old')
!!$    nnewl=0
!!$    do while(nnewl<9)
!!$       read(777) c
!!$       if(c==newline) nnewl=nnewl+1
!!$    end do
!!$    do k=1,n3
!!$       do j=1,n2
!!$          do i=1,n1
!!$             read(777) jj1,jj2,jj3
!!$             !u1(i,j,k)=jj1
!!$             !u2(i,j,k)=jj2
!!$             !u3(i,j,k)=jj3
!!$          end do
!!$       end do
!!$    end do
!!$    close(777,status='keep')

100 print *, 'unable to open jxb file'

    return

  end subroutine read_field

  subroutine mean_sdev(n1,n2,n3,field,mean,sdev)
    use types
    implicit none
    integer(ik), intent(IN) :: n1,n2,n3
    real(rk), dimension(1:n1,1:n2,1:n3), intent(IN) :: field
    real(rk), intent(OUT) :: mean, sdev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k
    real(rk) :: fac

    fac=real(n1*n2*n3,rk)**(-1)

    mean=0.0_rk
    sdev=0.0_rk
    !$omp parallel do reduction(+:mean)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             mean=mean+fac*field(i,j,k)
          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do reduction(+:sdev)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             sdev=sdev+fac*(field(i,j,k)-mean)**2
          end do
       end do
    end do
    !$omp end parallel do

    sdev=sqrt(sdev)

    return
  end subroutine mean_sdev


  subroutine surface_area(n,nsurf)
    use types
    use data
    implicit none
    integer(ik), intent(IN) :: n
    integer(ik), intent(OUT) :: nsurf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, dimension(3,3,3) :: box
    integer(ik) :: i,j,k,m,ii,jj,kk
    nsurf=0
    do m=1,structs(n)%npoints
       i=structs(n)%struct_ptr%points(m,1)
       j=structs(n)%struct_ptr%points(m,2)
       k=structs(n)%struct_ptr%points(m,3)
       box=.false.
       do ii=1,3 ; do jj=1,3 ; do kk=1,3
          box(ii,jj,kk)=refmask(per(i+ii-2),per(j+jj-2),per(k+kk-2))
       enddo; enddo ; enddo
       if (.not.all(box)) nsurf=nsurf+1
    enddo
    return
  end subroutine surface_area



  function in_box(n)
    use types
    use data,only:nnn
    implicit none
    integer(i2b), intent(IN) :: n
    integer(ik) :: in_box

    if(n<1) then
       in_box=n+nnn
    else if(n>nnn) then
       in_box=n-nnn
    else
       in_box=n
    end if

    return

  end function in_box

  subroutine statistics
    use types
    use data
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: n,m,i,j,k,ix,jx,kx,nsurf,mmm
    real(rk) :: dx,dv,ds,tmpohm,tmpvisc,tmpad,distx,disty,distz,distl,disttot
    real(rk) :: tmpux,tmpuy,tmpuz,tmpsigux,tmpsiguy,tmpsiguz
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
    dx=LBOX/(n1-1)
    ds=dx**2
    dv=dx**3



    mmm=1
    n=1
    do n=1,nstruct
       print *, n, structs(n)%npoints
       allocate(x(structs(n)%npoints,3))
       x(:,:)=0.0_rk
       nsurf=0
       call surface_area(n,nsurf)
       tmpohm=0.0_rk
       tmpvisc=0.0_rk
       tmpad=0.0_rk
       tmpux=0
       tmpuy=0
       tmpuz=0
       !$omp parallel do private(i,j,k,ix,jx,kx) reduction(+:tmpohm,tmpvisc,&
       !$omp tmpad,tmpux,tmpuy,tmpuz,tmpsigux,tmpsiguy,tmpsiguz)
       do m=1,structs(n)%npoints
          ix=structs(n)%struct_ptr%points(m,4)
          jx=structs(n)%struct_ptr%points(m,5)
          kx=structs(n)%struct_ptr%points(m,6)

          i=structs(n)%struct_ptr%points(m,1)
          j=structs(n)%struct_ptr%points(m,2)
          k=structs(n)%struct_ptr%points(m,3)
          x(m,1)=dx*(ix-1)
          x(m,2)=dx*(jx-1)
          x(m,3)=dx*(kx-1)
          tmpohm=tmpohm+dv*j2(i,j,k)
          tmpvisc=tmpvisc+dv*e(i,j,k)
          tmpad=tmpad+dv*jxb2(i,j,k)
!!$          tmpux=tmpux+u1(i,j,k)
!!$          tmpuy=tmpuy+u2(i,j,k)
!!$          tmpuz=tmpuz+u3(i,j,k)
       end do
       !$omp end parallel do
       if(structs(n)%npoints/=0) then
          tmpux=tmpux/structs(n)%npoints
          tmpuy=tmpuy/structs(n)%npoints
          tmpuz=tmpuz/structs(n)%npoints
       else
          tmpux=0.
          tmpuy=0.
          tmpuz=0.
       end if


       tmpsigux=0.0
       tmpsiguy=0.0
       tmpsiguz=0.0
       !do m=1,npoints(n)
       !tmpsigux=tmpsigux+(u1(i,j,k)-tmpux)**2
       !tmpsiguy=tmpsiguy+(u2(i,j,k)-tmpuy)**2
       !tmpsiguz=tmpsiguz+(u3(i,j,k)-tmpuz)**2
       !end do
       !tmpsigux=sqrt(tmpsigux/npoints(n))
       !tmpsiguy=sqrt(tmpsiguy/npoints(n))
       !tmpsiguz=sqrt(tmpsiguz/npoints(n))

       distx=0.0_rk
       disty=0.0_rk
       distz=0.0_rk
       distl=0.0_rk
       !!$omp parallel do reduction(max:distx,disty,distz,distl) private(disttot,n)
       do i=1,structs(n)%npoints
          do j=i,structs(n)%npoints
             distx=max(abs(x(i,1)-x(j,1)),distx)
             disty=max(abs(x(i,2)-x(j,2)),disty)
             distz=max(abs(x(i,3)-x(j,3)),distz)
             disttot=sqrt((x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+&
                  &(x(i,3)-x(j,3))**2)
             distl=max(distl,disttot)
          end do
       end do
       !!$omp end parallel do
       if(distx==0.0_rk) distx=dx
       if(disty==0.0_rk) disty=dx
       if(distz==0.0_rk) distz=dx
       if(distl==0.0_rk) distl=dx
       L(n)=distl
       Lx(n)=distx
       Ly(n)=disty
       Lz(n)=distz
       R(n)=sqrt(distx**2+disty**2+distz**2)
       V(n)=dv*structs(n)%npoints
       A(n)=ds*nsurf
       thick(n)=V(n)/(A(n)/2._rk)
       C(n)=PI*(L(n)/2._rk)**2/(A(n)/2._rk)
       ohm(n)=tmpohm
       visc(n)=tmpvisc
       ad(n)=tmpad
       ux(n)=tmpux
       uy(n)=tmpuy
       uz(n)=tmpuz
       sigux(n)=tmpsigux
       siguy(n)=tmpsiguy
       siguz(n)=tmpsiguz
       deallocate(x)
    end do

  end subroutine statistics


  subroutine compute_pdf_hist(npoints,dat,nbins,fname)
    implicit none
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik), intent(IN) :: npoints
    real(rk), dimension(npoints), intent(IN) :: dat
    integer(ik), intent(IN) :: nbins
    character(*), intent(IN) :: fname
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik) :: i,j,nn,n,ii
    real(rk) :: bw,dmax,dmin,vpmax,vpmin,accelmax,xdens,&
         &acceldens,vpdens,d,dx,c1,c2,c3,dir,area,xmin,xmax
    real(rk), dimension(3) :: tmpx,tmpvp,tmpaccel
    real(rk), dimension(:,:), allocatable :: pdf
    integer :: pdf_file=99
    print *, 'pdf: ', fname
    allocate(pdf(nbins,3))

    dmax=0.0_rk
    dmin=1.0e11_rk
    !$omp parallel do private(d) reduction(max:dmax) reduction(min:dmin)
    do n=1,npoints
       if(dat(n)>0.0) then
          d=dat(n)
          dmax=max(dmax,d)
          dmin=min(dmin,d)
       end if
    end do
    !$omp end parallel do

    dx=(log(dmax)-log(dmin))/real(nbins)
    xmin=log(dmin)
    pdf(:,:)=0.0
    !print *, '!!!!', dmax,dmin,xmin
    !$omp parallel do private(ii)
    do i=1,npoints
       if(dat(i)>0.0) then
          ii=int((log(dat(i))-xmin)/dx)+1
          if(ii==nbins+1) ii=nbins
          if(ii<1.or.ii>nbins) print *, 'histogram error:', ii,xmin
          pdf(ii,2)=pdf(ii,2)+1
       end if
    end do
    !$omp end parallel do

    dx=(log(dmax)-log(dmin))/real(nbins)
    !print *, fname, dmin,dmax,dx
    !$omp parallel do private(xmin,xmax)
    do ii=1,nbins
       xmin=exp(log(dmin)+(ii-1)*dx)
       xmax=exp(log(dmin)+(ii*dx))
       pdf(ii,1)=0.5*(xmin+xmax)
       pdf(ii,2)=pdf(ii,2)/(xmax-xmin)
    end do
    !$omp end parallel do

    area=0.0
    !$omp parallel do private(dx) reduction(+:area)
    do i=1,nbins-1
       dx=pdf(i+1,1)-pdf(i,1)
       area=area+(pdf(i+1,2)+pdf(i,2))*dx/2
    end do
    !$omp end parallel do

    if(area/=0.0) then
       !$omp parallel do
       do i=1,nbins
          pdf(i,2)=pdf(i,2)/area
          pdf(i,3)=pdf(i,3)/area
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

  subroutine compute_pdf_hist_uz(npoints,dat,nbins,fname)
    implicit none
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik), intent(IN) :: npoints
    real(rk), dimension(npoints), intent(IN) :: dat
    integer(ik), intent(IN) :: nbins
    character(*), intent(IN) :: fname
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik) :: i,j,nn,n,ii
    real(rk) :: bw,dmax,dmin,vpmax,vpmin,accelmax,xdens,&
         &acceldens,vpdens,d,dx,c1,c2,c3,dir,area,xmin,xmax
    real(rk), dimension(3) :: tmpx,tmpvp,tmpaccel
    real(rk), dimension(:,:), allocatable :: pdf
    integer :: pdf_file=99

    allocate(pdf(nbins,3))

    dmax=0.0_rk
    dmin=1.0e11_rk
    !$omp parallel do private(d) reduction(max:dmax) reduction(min:dmin)
    do n=1,npoints
       if(.true.) then
          d=dat(n)
          dmax=max(dmax,d)
          dmin=min(dmin,d)
       end if
    end do
    !$omp end parallel do
    print *, dmax,dmin
    dx=((dmax)-(dmin))/real(nbins)
    xmin=(dmin)
    pdf(:,:)=0.0
    !$omp parallel do private(ii)
    do i=1,npoints
       if(.true.) then
          ii=int(((dat(i))-xmin)/dx)+1
          if(ii==nbins+1) ii=nbins
          if(ii<1.or.ii>nbins) print *, 'histogram error:', ii,xmin
          pdf(ii,2)=pdf(ii,2)+1
       end if
    end do
    !$omp end parallel do

    dx=((dmax)-(dmin))/real(nbins)
    !print *, fname, dmin,dmax,dx
    !$omp parallel do private(xmin,xmax)
    do ii=1,nbins
       xmin=((dmin)+(ii-1)*dx)
       xmax=((dmin)+(ii*dx))
       pdf(ii,1)=0.5*(xmin+xmax)
       pdf(ii,2)=pdf(ii,2)/(xmax-xmin)
    end do
    !$omp end parallel do

    area=0.0
    !$omp parallel do private(dx) reduction(+:area)
    do i=1,nbins-1
       dx=pdf(i+1,1)-pdf(i,1)
       area=area+(pdf(i+1,2)+pdf(i,2))*dx/2
    end do
    !$omp end parallel do

    if(area/=0.0) then
       !$omp parallel do
       do i=1,nbins
          pdf(i,2)=pdf(i,2)/area
          pdf(i,3)=pdf(i,3)/area
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

  end subroutine compute_pdf_hist_uz


  subroutine output_statistics
    implicit none
    integer(ik) :: i

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
    logical, dimension(:,:,:), allocatable :: box,box2
    logical, dimension(:,:), allocatable :: tmp
    integer(ik) :: m,i,j,k,n,step_i,step_j,step_k
    allocate(box2(1:n1,1:n2,1:n3))
    allocate(box(1:n1,1:n2,1:n3))
    allocate(tmp(1:n1,1:n2))


    n=0
    do n=1,nstruct

       !$omp workshare
       box(:,:,:)=.false.
       !$omp end workshare
       
       !$omp parallel do private(i,j,k)
       do m=1,structs(n)%npoints
          i=structs(n)%struct_ptr%points(m,1)       
          j=structs(n)%struct_ptr%points(m,2)       
          k=structs(n)%struct_ptr%points(m,3)       
          box(i,j,k)=.true.
       end do
       !$omp end parallel do

       m=0
       step_i=0
       do while(connected(1_ik).and.m<n1)
          call rotate(1_ik)
          m=m+1
       end do
       step_i=m
       
       m=0
       step_j=0
       do while(connected(2_ik).and.m<n2)
          call rotate(2_ik)
          m=m+1
       end do
       step_j=m
       
       m=0
       step_k=0
       do while(connected(3_ik).and.m<n3)
          call rotate(3_ik)
          m=m+1
       end do
       step_k=m
       
       m=0
       do k=1,n3
          do j=1,n2
             do i=1,n1
                if(box(i,j,k)) then
                   m=m+1
                   structs(n)%struct_ptr%points(m,4)=i+step_i
                   structs(n)%struct_ptr%points(m,5)=j+step_j
                   structs(n)%struct_ptr%points(m,6)=k+step_k
                end if
             end do
          end do
       end do
       

!       struct_ptr=>struct_ptr%next

    end do


    return

  contains

    function connected(dim)
      implicit none
      integer(ik), intent(IN) :: dim
      logical :: connected,cell
      integer(ik) :: i,j,ii,jj,iii,jjj

      if(dim==1_ik) then
         connected=.false.
         !$omp parallel do private(i,j,k,ii,jj,iii,jjj,cell) shared(connected)
         do i=1,n1
            if(connected) cycle
            do j=1,n2
               if(connected) cycle
               cell=box(1,i,j)
               if(cell) then
                  do ii=i-1,i+1
                     if(connected) cycle
                     do jj=j-1,j+1
                        if(connected) cycle
                        iii=per(ii)
                        jjj=per(jj)
                        if(cell.and.box(n1,iii,jjj)) then
                           !$omp critical
                           connected=.true.
                           !$omp end critical
                        end if
                     end do
                  end do
               end if
            end do
         end do
         !$omp end parallel do
      else if(dim==2_ik) then
         connected=.false.
         !$omp parallel do private(i,j,k,ii,jj,iii,jjj,cell) shared(connected)
         do i=1,n1
            if(connected) cycle
            do j=1,n2
               if(connected) cycle
               cell=box(i,1,j)
               if(cell) then
                  do ii=i-1,i+1
                     if(connected) cycle
                     do jj=j-1,j+1
                        if(connected) cycle
                        iii=per(ii)
                        jjj=per(jj)
                        if(cell.and.box(iii,n2,jjj)) then
                           !$omp critical
                           connected=.true.
                           !$omp end critical
                        end if
                     end do
                  end do
               end if
            end do
         end do
         !$omp end parallel do
      else if(dim==3_ik) then
         connected=.false.
         !$omp parallel do private(i,j,k,ii,jj,iii,jjj,cell) shared(connected)
         do i=1,n1
            if(connected) cycle
            do j=1,n2
               if(connected) cycle
               cell=box(i,j,1)
               if(cell) then
                  do ii=i-1,i+1
                     if(connected) cycle
                     do jj=j-1,j+1
                        if(connected) cycle
                        iii=per(ii)
                        jjj=per(jj)
                        if(cell.and.box(iii,jjj,n3)) then
                           !$omp critical
                           connected=.true.
                           !$omp end critical
                        end if
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
      box2=box
      if(dim==1) then
         !$omp workshare
         tmp(:,:)=box2(1,:,:)
         box(1:n1-1,:,:)=box2(2:n1,:,:)
         box(n1,:,:)=tmp(:,:)
         !$omp end workshare
      else if(dim==2) then
         !$omp workshare
         tmp(:,:)=box2(:,1,:)
         box(:,1:n2-1,:)=box2(:,2:n2,:)
         box(:,n2,:)=tmp(:,:)
         !$omp end workshare
      else if(dim==3) then
         !$omp workshare
         tmp(:,:)=box2(:,:,1)
         box(:,:,1:n3-1)=box2(:,:,2:n3)
         box(:,:,n3)=tmp(:,:)
         !$omp end workshare
      end if

      return

    end subroutine rotate

  end subroutine translate_structures


  function newpoint(points,n,i,j,k)
    implicit none
    integer(i2b), dimension(:,:),intent(IN) :: points
    integer(ik) :: i,j,k,n
    logical :: newpoint
    integer(ik) :: ii

    newpoint=.true.
    do ii=1,n
       if(all(points(ii,:)==(/i,j,k/))) then
          newpoint=.false.
          return
       end if
    end do

    return

  end function newpoint

  subroutine output_structures
    use types
    implicit none
    integer(ik) :: n,i,j,k,nfile,m
    character(128) :: filename

    
    nfile=0
    do n=1,nstruct
       if(structs(n)%npoints>VOLUME) then
          write(filename,'(a,i0,a)') 'out.',nfile,'.vtk'
          call write_vtk_file(filename,n)
          !print *, n,nfile
          nfile=nfile+1
       end if
    end do
    print '(a,i4,a)', 'Wrote ', nfile, ' files.'

    !call write_vtk_file_2

    return

  end subroutine output_structures

   subroutine write_vtk_file(filename,n)
    use types 
    implicit none
    character(*) :: filename
    integer(ik) :: n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k,m,vtk_file=999
    integer(ik) :: tmp
    real(rk) :: x,y,z

    open(vtk_file,file=filename,form='formatted',action='write')

    write(vtk_file, '(t1,a)') '# vtk DataFile Version 2.0'
    write(vtk_file, '(t1,a)')  ' '
    write(vtk_file, '(t1,a)') 'ASCII'
    write(vtk_file, '(t1,a)')  ' '
    write(vtk_file, '(t1,a)') 'DATASET POLYDATA'
    write(vtk_file, '(t1,a,i20,a)') 'POINTS ', structs(n)%npoints, ' double'


    do m=1,structs(n)%npoints
       x=LBOX*structs(n)%struct_ptr%points(m,1)/(n1-1)
       y=LBOX*structs(n)%struct_ptr%points(m,2)/(n2-1)
       z=LBOX*structs(n)%struct_ptr%points(m,3)/(n3-1)
       write(vtk_file,'(t1,3e19.8)') x,y,z
    end do


    close(vtk_file,status='keep')

    return

  end subroutine write_vtk_file

end module routines

program statstruct
  use types
  use routines
  use data,only:nnn
  implicit none
  real(rk) :: mean,sdev
  integer(ik) :: m=3
  integer(ik) :: i,j,k,vol,k1,k2,maxvol,minvol,n,kk,none,which,nppdf,&
       &nstruct2
  real(rk) :: t1,t2,tmp,totohm,totad,totvisc,sohm,svisc,sad,tottot,stot,svol
  character(256) :: path,fmt
  real(8), external :: omp_get_wtime
  integer(i2b),  dimension(:,:), pointer :: point_array
  integer(i8b) :: maxpoints, itmp
  integer(i2b), dimension(:,:), allocatable :: ptmp
  integer(i2b) :: i2btmp
  integer(i8b) :: i8btmp
  VOLUME=0

  print *, 'Resolution: '
  nnn=128
  read(*,*) nnn
  print *, 'Path: '
  path='~/runs/mhd128uf'
  read(*,'(a)') path
  print *, 'Which?'
  which=0
  read(*,*) which
  print *, 'Sdevs: '
  m=3
  read(*,*) m
  print *, 'AD coefficient: '
  ADCOEFF=1.!0.01
  read(*,*) ADCOEFF
  print *, 'eta: '
  ETA=1.!0.181e-2
  read(*,*) ETA
  n1=nnn
  n2=nnn
  n3=nnn

 
  allocate(j2(1:n1,1:n2,1:n3))
  allocate(e(1:n1,1:n2,1:n3))
  allocate(jxb2(1:n1,1:n2,1:n3))
  allocate(total(1:n1,1:n2,1:n3))
  allocate(mask(n1,n2,n3))
  !$omp workshare
  mask=.false.
  !$omp end workshare
  allocate(refmask(n1,n2,n3))
  !$omp workshare
  refmask=.false.
  !$omp end workshare
  allocate(intrv(n1,n2,n3/2,2))
  allocate(nintrv(n1,n2))

  allocate(dist(n1*n2*n3))
  allocate(diff(n1*n2*n3))
  allocate(diffmask(n1*n2*n3))

  !$omp parallel do
  do k=1,n3
     do j=1,n2
        do i=1,n1
           j2(i,j,k)=0.0_rk
        end do
     end do
  end do
  !$omp end parallel do

  !$omp parallel do
  do k=1,n3/2
     do j=1,n2
        do i=1,n1
           intrv(i,j,k,:)=0_ik
        end do
     end do
  end do
  !$omp end parallel do


  !$omp parallel do
  do j=1,n2
     do i=1,n1
        nintrv(i,j)=0_ik
     end do
  end do
  !$omp end parallel do

  print '(3a)', 'reading files in ', trim(path), '...'
  call read_field(path,n1,n2,n3,j2,e,jxb2,u1,u2,u3)
  call cpu_time(t1)

  !$omp parallel do 
  do k=1,n3
     do j=1,n2
        do i=1,n1
           total(i,j,k)=j2(i,j,k)+e(i,j,k)+jxb2(i,j,k)
        end do
     end do
  end do
  !$omp end parallel do

  print *, 'ok'

  if(which==0) then
     call mean_sdev(n1,n2,n3,j2,mean,sdev)
  else if(which==1) then
     call mean_sdev(n1,n2,n3,e,mean,sdev)
  else if(which==2) then
     call mean_sdev(n1,n2,n3,jxb2,mean,sdev)
  else if(which==3) then
     call mean_sdev(n1,n2,n3,total,mean,sdev)
  end if



  !$omp parallel do private(tmp)
  do k=1,n3
     do j=1,n2
        do i=1,n1
           if(which==0) then
              tmp=j2(i,j,k)
           else if(which==1) then
              tmp=e(i,j,k)
           else if(which==2) then
              tmp=jxb2(i,j,k)
           else if(which==3) then
              tmp=total(i,j,k)
           end if
           if(tmp>mean+m*sdev) refmask(i,j,k)=.true.
        end do
     end do
  end do
  !$omp end parallel do


  open(986,file='struct.bin',form='unformatted',action='read')
  read(986) nstruct
  print *, nstruct
  allocate(structs(nstruct))
  print *, 'nstruct: ', nstruct
  do n=1,nstruct
     read(986) structs(n)%npoints
     allocate(structs(n)%struct_ptr)
     allocate(structs(n)%struct_ptr%points(structs(n)%npoints,6))
     read(986) structs(n)%struct_ptr%points(1:structs(n)%npoints,1:3)
  end do
  
  close(986)


  print *, 'translating structures...'

  call translate_structures
  print *, 'ok'

  print *, 'statistics...'
  call statistics

  print *, 'ok'

  print *, 'pdfs...'
  nppdf=int(sqrt(real(nstruct)))
  allocate(totdis(1:nstruct))
  totdis=visc+ohm+ad

  call compute_pdf_hist(nstruct,totdis,nppdf,'tot-pdf.dat')
  call compute_pdf_hist(nstruct,ad,nppdf,'ad-pdf.dat')
  call compute_pdf_hist(nstruct,ohm,nppdf,'ohm-pdf.dat')
  call compute_pdf_hist(nstruct,visc,nppdf,'visc-pdf.dat')
  call compute_pdf_hist(nstruct,V,nppdf,'vol-pdf.dat')
  call compute_pdf_hist(nstruct,L,nppdf,'L-pdf.dat')
  call compute_pdf_hist(nstruct,R,nppdf,'R-pdf.dat')
  call compute_pdf_hist(nstruct,A,nppdf,'A-pdf.dat')
  call compute_pdf_hist(nstruct,thick,nppdf,'H-pdf.dat')
  !call random_number(C)
  !call compute_pdf_hist(nstruct,C,nppdf,'C-pdf.dat')
!!$  call compute_pdf_hist_uz(nstruct,ux,nppdf,'ux-pdf.dat')
!!$  call compute_pdf_hist_uz(nstruct,uy,nppdf,'uy-pdf.dat')
!!$  call compute_pdf_hist_uz(nstruct,uz,nppdf,'uz-pdf.dat')
!!$  call compute_pdf_hist(nstruct,sigux,nppdf,'sigux-pdf.dat')
!!$  call compute_pdf_hist(nstruct,siguy,nppdf,'siguy-pdf.dat')
!!$  call compute_pdf_hist(nstruct,siguz,nppdf,'siguz-pdf.dat')

  print *, 'ok'
  call cpu_time(t2)

  print *, 'time: ', t2-t1

  print *, 'output...'
  call output_statistics
  call output_structures
  print *, 'ok.'

  tmp=0.0
  totohm=0.0
  totvisc=0.0
  totad=0.0
  tottot=0.0
  !$omp parallel do reduction(+:totohm,totvisc,totad,tottot)
  do k=1,n3
     do j=1,n2
        do i=1,n1
           totohm=totohm+j2(i,j,k)
           totvisc=totvisc+e(i,j,k)
           totad=totad+jxb2(i,j,k)
           tottot=tottot+j2(i,j,k)+e(i,j,k)+jxb2(i,j,k)
        end do
     end do
  end do
  !$omp end parallel do

  sohm=sum(j2,mask=refmask)
  svisc=sum(e,mask=refmask)
  sad=sum(jxb2,mask=refmask)
  stot=sohm+svisc+sad
  svol=count(refmask)

  totohm=sum(j2)
  totvisc=sum(e)
  totad=sum(jxb2)

  tmp=(LBOX/(nnn-1))**3
  !print '(a,e25.7)', 'global sigux: ', sqrt(sum(u1(:,:,:)**2)/(n1*n2*n3))
  !print '(a,e25.7)', 'global siguy: ', sqrt(sum(u2(:,:,:)**2)/(n1*n2*n3))
  !print '(a,e25.7)', 'global siguz: ', sqrt(sum(u3(:,:,:)**2)/(n1*n2*n3))

  print '(a,f7.2,a)','Ohm: ', sohm/totohm*100.,'%'
  print '(a,f7.2,a)','Visc: ', svisc/totvisc*100.,'%'
  print '(a,f7.2,a)','AD: ', sad/totad*100.,'%'
  print '(a,f7.2,a)','Total: ', stot/tottot*100.,'%'
  print '(a,f7.2,a)','Volume: ', svol/(nnn**3)*100.,'%'
  open(888,file='percent.dat',form='formatted',action='write')
  write(888,'(a,i5)') 'n= ', nnn
  write(888,'(a,i5)') 'field= ', which
  write(888,'(3(a,e15.3))') 'tohm= ', tmp*totohm, ' totvisc= ', tmp*totvisc,&
       & ' totad=', tmp*totad
  write(888,'(i6,a)') nstruct, ' structures'
  write(888,'(a,i15)') 'maxvol: ', maxvol
  write(888,'(a,f7.2,a)') 'Ohm: ', sohm/totohm*100.,'%'
  write(888,'(a,f7.2,a)') 'Visc: ', svisc/totvisc*100.,'%'
  write(888,'(a,f7.2,a)') 'AD: ', sad/totad*100.,'%'
  write(888,'(a,f7.2,a)') 'Total: ', stot/tottot*100.,'%'
  write(888,'(a,f7.2,a)') 'Volume: ', svol/(nnn**3)*100.,'%'

  write(888,'(a,f7.2,a)') 'Total: ', (sohm+svisc+sad)/(totohm+totvisc+&
       &totad)*100.,'%'


  stop

contains

  subroutine timing(t)
    implicit none
    real(rk), intent(OUT) :: t

#ifdef OPENMP
    !$omp parallel
    t=omp_get_wtime()
    !$omp end parallel
#else
    call cpu_time(t)
#endif

    return

  end subroutine timing

end program statstruct
