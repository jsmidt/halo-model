program test
use halo
use interp
use bicubic
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
write(*,*) hm%k(200), hm%Pk(200)

call init_halo()

!write(*,*) size(hm%ukm,1)
!write(*,*) exp(hm%lnk2d)

!kk  = -4
!m = 1e11
!mm = 4.64e12
!mmm = 2.15e14
!mmmm = 1e16
!open(unit=10,file='output/' // trim(hm%run_name) // '_ukm.dat',form='formatted',status='unknown')
!do i=1,50
!    kg = 10**kk
!    bix(1) = log(kg)
!    biy(1) = log(mm)
!    call rgbi3p(1, size(hm%lnk2d),size(hm%lnm2d), hm%lnk2d, hm%lnm2d, hm%ukm, 1, bix, biy, biz, ier)
!    write(10,'(6Es12.3)') kg,ukm(kg,m),ukm(kg,mm),ukm(kg,mmm),ukm(kg,mmmm),biz(1,1)
!    kk = kk+0.15
!end do

!do i = 1,160
!    write(*,*) exp(hm%lnk2d(i)), exp(hm%lnm2d(180)),hm%ukm(i,180)
!enddo

kk  = -4
open(unit=10,file='output/' // trim(hm%run_name) // '_powerspec.dat',form='formatted',status='unknown')
do i=1,48
    kg = 10**kk
    write(*,*) kg
    write(10,'(4Es12.3)') kg,interpf(hm%k,dble(hm%Pk),kg),P1h(kg),P2h(kg)
    kk = kk+0.15
end do

end program test

