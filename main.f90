program test
use halo
implicit none

real(dp) :: m,rombint
real, allocatable,dimension(:) :: k,Pk
integer :: i

z = 0.1
m = 1.0

write(*,*) w_top(m)
write(*,*) sig_2(m)
allocate(k(mpts),Pk(mpts))
call linear_pk(k,Pk)

   open(unit=10,file='matterpower01.dat',form='formatted',status='unknown')
   do i=1,mpts
     write(10,*) k(i),Pk(i)
   end do
   close(10)

z = 0.1
write(*,*) nu(m)
z = 10
write(*,*) nu(m)
z = 0.1
write(*,*) nu_fnu(m)
z = 10
write(*,*) nu_fnu(m)


z = 0.1
write(*,*) ukm(z,m)
z = 10
write(*,*) ukm(z,m)


end program test
