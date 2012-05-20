program test
use halo
implicit none

real(dp) :: m,rombint,kk
real, allocatable,dimension(:) :: k,Pk
integer :: i

z = 0.1
m = 1.0

write(*,*) w_top(m)
write(*,*) sig_2(m)

allocate(k(mpts),Pk(mpts))
!z = 0.0
!call linear_pk(k,Pk)
!open(unit=10,file='matterpower0.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(10,*) k(i),Pk(i)
!end do
!close(10)

!z = 0.1
!call linear_pk(k,Pk)
!open(unit=10,file='matterpower01.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(10,*) k(i),Pk(i)
!end do
!close(10)


!z = 1.0
!call linear_pk(k,Pk)
!open(unit=10,file='matterpower1.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(10,*) k(i),Pk(i)
!end do
!close(10)

!z = 10.0
!call linear_pk(k,Pk)
!open(unit=10,file='matterpower10.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(10,*) k(i),Pk(i)
!end do
!close(10)

!z = 0.1
!write(*,*) nu(m)
!z = 10
!write(*,*) nu(m)
!z = 0.1
!z = 0.1
!m = 1e6
!write(*,*) fnu(m)
!z = 10
!write(*,*) fnu(m)
write(*,*) tol


!z = 0.1
!m = 1e6
!write(*,*) ukm(z,m)
!z = 10
!write(*,*) ukm(z,m)


z = 1.0
kk = 0.1
write(*,*) P1h(kk)


end program test
