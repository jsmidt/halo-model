program test
use halo
use interp
implicit none

real(dl) :: m,rombint,kg,mm,mmm,kk,mmmm,zi
real(dl), allocatable,dimension(:) :: xx,yy
real(dl), dimension(1) :: bix,biy
real(dl), dimension(1,1) :: biz
integer :: i,ier

allocate(hm%k(mpts),hm%Pk(mpts))
allocate(xx(10),yy(10))
call init_params()
call linear_pk(hm%k,hm%Pk)
call init_halo()


kk  = -4
m = 1e11
mm = 4.64e12
mmm = 2.15e14
mmmm = 1e16
open(unit=10,file='output/' // trim(hm%run_name) // '_ukm.dat',form='formatted',status='unknown')
do i=1,48
    kg = 10**kk
    !write(*,*) kg
    write(10,'(6Es12.3)') kg,ukm(kg,m),ukm(kg,mm),ukm(kg,mmm),ukm(kg,mmmm)
    kk = kk+0.15
end do

kk  = -4
open(unit=10,file='output/' // trim(hm%run_name) // '_powerspec.dat',form='formatted',status='unknown')
do i=1,48
    kg = 10**kk
    !write(*,*) kg
    write(10,'(4Es12.3)') kg,interpf(hm%k,dble(hm%Pk),kg),P1h(kg),P2h(kg)
    kk = kk+0.15
end do

end program test

