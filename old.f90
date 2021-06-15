
subroutine read_field(path,n1,n2,n3,j2,e,jxb2,u1,u2,u3)
  use types
  implicit none
  character(*), intent(IN)            :: path
  integer(ik), intent(IN)             :: n1,n2,n3
  real(rk), dimension(1:n1,1:n2,1:n3) :: j2,e,jxb2,u1,u2,u3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Reads an input scalar field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(ik)                         :: i,j,k,nnewl
  real(dp)                            :: jj1,jj2,jj3
  character(512)                      :: infile
  character                           :: c,newline

  !?????????????????????????????????????????????????

  j=0.0
  !e=0.0
  !jxb2=0.0

  newline=char(10)
!!$ 
!!$    print *, 'j00000.vtk...'
!!$    write(infile,'(2a)') trim(path), '/j00000.vtk'
!!$    open(777,file=infile,access='stream',&
!!$         &action='read',status='old')
!!$    nnewl=0
!!$
!!$    do while(nnewl<9)
!!$       read(777) c
!!$       if(c==newline) nnewl=nnewl+1
!!$    end do
!!$    do k=1,n3
!!$       print *,  k
!!$       do j=1,n2
!!$          do i=1,n1
!!$             read(777) jj1,jj2,jj3
!!$             j2(i,j,k)=ETA*(jj1**2+jj2**2+jj3**2)
!!$          end do
!!$       end do
!!$    end do
!!$    close(777,status='keep')
!!$
!!$

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
!!$
!!$
!!$
!!$    print *, 'jxb00000.vtk...'
!!$    write(infile,'(2a)') trim(path), '/jxb00000.vtk'
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
!!$             jxb2(i,j,k)=ADCOEFF*(jj1**2+jj2**2+jj3**2)
!!$          end do
!!$       end do
!!$    end do
!!$    close(777,status='keep')
!!$
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
!!$             u1(i,j,k)=jj1
!!$             u2(i,j,k)=jj2
!!$             u3(i,j,k)=jj3
!!$          end do
!!$       end do
!!$    end do
!!$    close(777,status='keep')

100 print *, 'unable to open jxb file'

  return

end subroutine read_field

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

subroutine translate_structures_old
  use types
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(ik) :: i,j,k,m,step_i,step_j,step_k,si,sj,sk
  integer(ik), dimension(1) :: loc


  struct_ptr=>struct
  do while(associated(struct_ptr).and.allocated(struct_ptr%points))

     step_i=0
     step_j=0
     step_k=0

     dist(1:struct_ptr%npoints)=0.0_rk
     diff(1:struct_ptr%npoints)=0.0_rk

     diffmask(1:struct_ptr%npoints-1)=.true.
     diffmask(struct_ptr%npoints:)=.false.


     if(struct_ptr%crosses_i) then
        !$omp parallel do
        do m=1,struct_ptr%npoints
           dist(m)=real(struct_ptr%points(m,1),rk)
        end do
        !$omp end parallel do
        call qsort(dist,struct_ptr%npoints)
        !$omp parallel do
        do m=1,struct_ptr%npoints-1
           diff(m)=dist(m+1)-dist(m)
        end do
        !$omp end parallel do
        loc=maxloc(diff,diffmask)
        step_i=int(dist(loc(1)),ik)


        do m=1,struct_ptr%npoints
           if(struct_ptr%points(m,1)<=step_i) then
              struct_ptr%points(m,1)=struct_ptr%points(m,1)+nnn
           end if
        end do
     end if


     if(struct_ptr%crosses_j) then
        !$omp parallel do
        do m=1,struct_ptr%npoints
           dist(m)=real(struct_ptr%points(m,2),rk)
        end do
        !$omp end parallel do
        call qsort(dist,struct_ptr%npoints)
        !$omp parallel do
        do m=1,struct_ptr%npoints-1
           diff(m)=dist(m+1)-dist(m)
        end do
        !$omp end parallel do
        loc=maxloc(diff,diffmask)
        step_j=int(dist(loc(1)),ik)

        do m=1,struct_ptr%npoints
           if(struct_ptr%points(m,2)<=step_j) then
              struct_ptr%points(m,2)=struct_ptr%points(m,2)+nnn
           end if
        end do

     end if


     if(struct_ptr%crosses_k) then
        !$omp parallel do
        do m=1,struct_ptr%npoints
           dist(m)=real(struct_ptr%points(m,3),rk)
        end do
        !$omp end parallel do
        call qsort(dist,struct_ptr%npoints)
        !$omp parallel do
        do m=1,struct_ptr%npoints-1
           diff(m)=dist(m+1)-dist(m)
        end do
        !$omp end parallel do
        loc=maxloc(diff,diffmask)
        step_k=int(dist(loc(1)),ik)

        do m=1,struct_ptr%npoints
           if(struct_ptr%points(m,3)<=step_k) then
              struct_ptr%points(m,3)=struct_ptr%points(m,3)+nnn
           end if
        end do

     end if

     struct_ptr=>struct_ptr%next

  end do



  return

end subroutine translate_structures_old


function is_neighbour_old(i1,i2)
  use types
  implicit none
  integer(i2b), dimension(2),intent(IN) :: i1,i2
  logical                               :: is_neighbour_old
  integer(ik)                           :: start1,start2,end1,end2

  start1=i1(1)
  end1=i1(2)
  start2=i2(1)
  end2=i2(2)

  if(start1<=end1.and.start2<=end2) then
     if(end1<start2-1.or.end2<start1-1) then
        is_neighbour_old=.false.
     else
        is_neighbour_old=.true.
     end if
  else if(start1>end1.and.start2<=end2) then
     if(start2>end1+1.and.end2<start1-1) then
        is_neighbour_old=.false.
     else
        is_neighbour_old=.true.
     end if
  else if(start1<=end1.and.start2>end2) then
     if(start1>end2+1.and.end1<start2-1) then
        is_neighbour_old=.false.
     else
        is_neighbour_old=.true.
     end if
  else if(start1>end1.and.start2>end2) then
     is_neighbour_old=.true.
  endif

  return

end function is_neighbour_old

subroutine intervals_old(n1,n2,n3,field,intrv,nintrv,mean,sdev,m)
  use types
  implicit none
  integer(ik),intent(IN)                       :: n1,n2,n3
  real(rk), dimension(1:n1,1:n2,1:n3)          :: field
  integer(i2b), dimension(1:n1,1:n2,1:n3/2,2)  :: intrv
  integer(ik), dimension(1:n1,1:n2)            :: nintrv
  real(rk)                                     :: mean,sdev
  integer(ik)                                  :: m
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(ik)                                  :: i,j,k,start,end,skip1
  integer(ik)                                  :: skip2,k1,k2,kk
  logical                                      :: openintrv
  integer(i2b), dimension(:,:,:,:),allocatable :: intrv1
  integer(i2b), dimension(:,:),allocatable     :: nintrv1
  logical, dimension(:,:),allocatable          :: contintrv

  allocate(intrv1(n1,n2,n3/2,2))
  allocate(contintrv(n1,n2))
  allocate(nintrv1(n1,n2))
  contintrv(:,:)=.false.
  nnintrv=0
  nintrv1=0
  intrv1=0


  do j=1,n2
     do i=1,n1
        k=1
        openintrv=.false.
        do while (k<=n3+1)
           kk=k
           if(kk==n3+1) kk=1
           if(refmask(i,j,kk)) then
              if(.not.(openintrv)) then
                 nintrv1(i,j)=nintrv1(i,j)+1
                 intrv1(i,j,nintrv1(i,j),1)=k
                 openintrv=.true.
              end if
           else
              if(openintrv) then
                 intrv1(i,j,nintrv1(i,j),2)=k-1
                 openintrv=.false.
                 nnintrv=nnintrv+1
              end if
           end if
           k=k+1
        end do
        if(openintrv) contintrv(i,j)=.true.
     end do
  end do



  intrv=0
  nintrv=0


  do j=1,n2
     do i=1,n1
        if(intrv1(i,j,1,1)==1.and.intrv1(i,j,nintrv1(i,j),2)==n3) then
           nintrv(i,j)=nintrv(i,j)+1
           intrv(i,j,1,1)=intrv1(i,j,nintrv1(i,j),1)
           intrv(i,j,1,2)=intrv1(i,j,1,2)
           k1=2
           k2=nintrv1(i,j)-1
        else
           k1=1
           k2=nintrv1(i,j)
        end if

        do k=k1,k2
           nintrv(i,j)=nintrv(i,j)+1
           intrv(i,j,nintrv(i,j),:)=intrv1(i,j,k,:)
        end do
     end do
  end do


  deallocate(intrv1)
  deallocate(contintrv)
  deallocate(nintrv1)

  return

end subroutine intervals_old

subroutine write_vtk_file_2
  use types 
  implicit none
  integer(ik) :: n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(ik) :: i,j,k,m,ivtk_file=999,dvtk_file=888,itot,dtot
  integer(ik) :: tmp
  real(rk)    :: x,y,z,idown=0.1,iup=1.3,ddown=0.01,dup=0.08

  open(ivtk_file,file='inertial.vtk',form='formatted',action='write')
  open(dvtk_file,file='dissipative.vtk',form='formatted',action='write')

  write(ivtk_file, '(t1,a)') '# vtk DataFile Version 2.0'
  write(ivtk_file, '(t1,a)')  ' '
  write(ivtk_file, '(t1,a)') 'ASCII'
  write(ivtk_file, '(t1,a)')  ' '
  write(ivtk_file, '(t1,a)') 'DATASET POLYDATA'

  write(dvtk_file, '(t1,a)') '# vtk DataFile Version 2.0'
  write(dvtk_file, '(t1,a)')  ' '
  write(dvtk_file, '(t1,a)') 'ASCII'
  write(dvtk_file, '(t1,a)')  ' '
  write(dvtk_file, '(t1,a)') 'DATASET POLYDATA'
  itot=0
  dtot=0
  n=1
  struct_ptr=>struct
  n=1
  do while(associated(struct_ptr).and.allocated(struct_ptr%points))
     if(idown<L(n).and.L(n)<iup) itot=itot+struct_ptr%npoints
     if(ddown<L(n).and.L(n)<dup) dtot=dtot+struct_ptr%npoints
     struct_ptr=>struct_ptr%next
     n=n+1
  end do

  write(ivtk_file, '(t1,a,i20,a)') 'POINTS ', itot, ' double'

  write(dvtk_file, '(t1,a,i20,a)') 'POINTS ', dtot, ' double'

  struct_ptr=>struct
  n=1
  do while(associated(struct_ptr).and.allocated(struct_ptr%points))
     if(idown<L(n).and.L(n)<iup) then
        do m=1,struct_ptr%npoints
           x=LBOX*struct_ptr%points(m,1)/(n1-1)
           y=LBOX*struct_ptr%points(m,2)/(n2-1)
           z=LBOX*struct_ptr%points(m,3)/(n3-1)
           write(ivtk_file,'(t1,3e19.8)') x,y,z
        end do
     else if(ddown<L(n).and.L(n)<dup) then
        do m=1,struct_ptr%npoints
           x=LBOX*struct_ptr%points(m,1)/(n1-1)
           y=LBOX*struct_ptr%points(m,2)/(n2-1)
           z=LBOX*struct_ptr%points(m,3)/(n3-1)
           write(dvtk_file,'(t1,3e19.8)') x,y,z
        end do
     end if
     struct_ptr=>struct_ptr%next
     n=n+1
  end do

  close(ivtk_file,status='keep')

  close(dvtk_file,status='keep')

  return
end subroutine write_vtk_file_2

subroutine surface_area_old(astruct,nsurf)
  use types
  implicit none
  type(structure), intent(IN), pointer   :: astruct
  integer(ik), intent(OUT)               :: nsurf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, dimension(:,:,:), allocatable :: box
  integer(ik)                            :: is,js,ks,ie,je,ke,ir,jr
  integer(ik)                            :: kr,i,j,k,m,ii,jj,kk,ll
  real(rk)                               :: gradi,gradj,gradk,grad
  integer(ik), dimension(3,3,3)          :: neighbours
  logical                                :: boundary
  integer(i2b), dimension(6,3)           :: neigh
  if(astruct%npoints<27) then
     nsurf=astruct%npoints
  else
     is=minval(astruct%points(:,4))
     js=minval(astruct%points(:,5))
     ks=minval(astruct%points(:,6))
     ie=maxval(astruct%points(:,4))
     je=maxval(astruct%points(:,5))
     ke=maxval(astruct%points(:,6))
     ir=ie-is+10
     jr=je-js+10
     kr=ke-ks+10
     allocate(box(1:ir,1:jr,1:kr))
     box=.false.
     !$omp parallel do private(i,j,k)
     do m=1,astruct%npoints
        i=astruct%points(m,4)-is+3
        j=astruct%points(m,5)-js+3
        k=astruct%points(m,6)-ks+3
        box(i,j,k)=.true.
     end do
     !$omp end parallel do

     nsurf=0
     !$omp parallel do private(boundary) reduction(+:nsurf)
     do k=1,kr
        do j=1,jr
           do i=1,ir
              if(box(i,j,k)) then
                 boundary=.false.
                 neigh(1,:)=(/i+1,j,k/)
                 neigh(2,:)=(/i-1,j,k/)
                 neigh(3,:)=(/i,j+1,k/)
                 neigh(4,:)=(/i,j-1,k/)
                 neigh(5,:)=(/i,j,k+1/)
                 neigh(6,:)=(/i,j,k-1/)
                 do ll=1,6
                    ii=neigh(ll,1)
                    jj=neigh(ll,2)
                    kk=neigh(ll,3)
                    if(.not.box(ii,jj,kk)) then
                       boundary=.true.
                       exit
                    end if
                 end do
                 if(boundary) nsurf=nsurf+1
              end if
           end do
        end do
     end do
     !$omp end parallel do
  end if

  return 
end subroutine surface_area_old


subroutine compute_pdf_hist_uz(npoints,dat,nbins,fname)
  implicit none
  integer(ik), intent(IN)                  :: npoints
  real(rk), dimension(npoints), intent(IN) :: dat
  integer(ik), intent(IN)                  :: nbins
  character(*), intent(IN)                 :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(ik) :: i,j,nn,n,ii
  real(rk)                                 :: bw,dmax,dmin,vpmax,vpmin
  real(rk)                                 :: accelmax,xdens,acceldens
  real(rk)                                 :: vpdens,d,dx,c1,c2,c3,dir
  real(rk)                                 :: area,xmin,xmax
  real(rk), dimension(3)                   :: tmpx,tmpvp,tmpaccel
  real(rk), dimension(:,:), allocatable    :: pdf
  integer                                  :: pdf_file=99

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

!!$  open(986,file='struct-old.dat',form='formatted',action='read')
!!$  read(986,*) nstruct
!!$  struct_ptr=>struct
!!$  do n=1,nstruct
!!$     if(n/=1) allocate(struct_ptr)
!!$     read(986,*) struct_ptr%npoints
!!$     if(struct_ptr%npoints/=0) then
!!$        write(fmt,'(i0)') 3*(struct_ptr%npoints)
!!$        fmt='('//trim(fmt)//'i7'//')'
!!$        point_array=>struct_ptr%points
!!$        allocate(point_array(1:struct_ptr%npoints,1:3))
!!$        read(986,trim(fmt)) point_array(1:struct_ptr%npoints,1:3)
!!$     end if
!!$     struct_ptr=>struct_ptr%next
!!$  end do
!!$
!!$
!!$  open(985,file='struct.dat',form='formatted',action='write')
!!$  write(985,*) nstruct-1
!!$  struct_ptr=>struct
!!$  do while(associated(struct_ptr).and.associated(struct_ptr%intrv))
!!$     if(struct_ptr%npoints/=0) then
!!$        write(985,*) struct_ptr%npoints
!!$        write(fmt,'(i0)') 3*(struct_ptr%npoints)
!!$        fmt='('//trim(fmt)//'i7'//')'
!!$        write(985,fmt) struct_ptr%points(1:struct_ptr%npoints,1:3)
!!$     end if
!!$     struct_ptr=>struct_ptr%next
!!$  end do

  
!!$  stop
!!$  
!!$  print *, 'Resolution: '
!!$  nnn=128
!!$  read(*,*) nnn
!!$  print *, 'Path: '
!!$  path='~/runs/mhd128uf'
!!$  read(*,'(a)') path
!!$  print *, 'Which?'
!!$  which=0
!!$  read(*,*) which
!!$  print *, 'Sdevs: '
!!$  m=3
!!$  read(*,*) m
!!$  print *, 'AD coefficient: '
!!$  ADCOEFF=1.!0.01
!!$  read(*,*) ADCOEFF
!!$  print *, 'eta: '
!!$  ETA=1.!0.181e-2
!!$  read(*,*) ETA
!!$  n1=nnn
!!$  n2=nnn
!!$  n3=nnn

!!$  print '(3a)', 'reading files in ', trim(path), '...'
!!$  call read_field(path,n1,n2,n3,j2,e,jxb2,u1,u2,u3)
!!$  call cpu_time(t1)

!!$omp parallel do 
!!$  do k=1,n3
!!$     do j=1,n2
!!$        do i=1,n1
!!$           !total(i,j,k)=j2(i,j,k)+e(i,j,k)+jxb2(i,j,k)
!!$        end do
!!$     end do
!!$  end do
!!$omp end parallel do

!!$  print *, 'ok'
!!$
!!$  if(which==0) then
!!$     call mean_sdev(n1,n2,n3,j2,mean,sdev)
!!$  else if(which==1) then
!!$     call mean_sdev(n1,n2,n3,e,mean,sdev)
!!$  else if(which==2) then
!!$     call mean_sdev(n1,n2,n3,jxb2,mean,sdev)
!!$  else if(which==3) then
!!$     call mean_sdev(n1,n2,n3,total,mean,sdev)
!!$  end if
