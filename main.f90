program test
use halo
use interp
implicit none

real(dl) :: m,rombint,kg,mm,mmm,kk
real(dl), allocatable,dimension(:) :: xx,yy
integer :: i

!! Define cosmology (WMAP 7)
!omegab  = 0.045
!omegac  = 0.227
!omegam  = omegab+omegac
!omegal  = 0.728
!omegan  = 0
!H0      = 70.2
!YHe     = 0.24
!Num_Nu_massless =3.04
!Num_Nu_massive  =0


! To get the critical density rho_c we need to calculate 3H^2/(8\piG).  
! We want units of M_\sun/Mpc^3. We know H0 has units of km/s/Mpc and 
! G = 4.302e-9 Mpc(km/s)^2/M_\sun So that 3H^2/(8\piG) has the right units.
!rho_c = 3.0d0*H0**2/8.0d0/pi/4.302e-9

allocate(hm%k(mpts),hm%Pk(mpts))
allocate(xx(10),yy(10))
call init_params()
call linear_pk(hm%k,hm%Pk)

do i = 1,10
   xx(i) = i
   yy(i) = xx(i)**2
enddo

!stop
!write(*,*) xx
!write(*,*) interpf(xx,yy,dble(2.5))
!m = 1e9
!write(*,*) sig_2(m)



open(unit=10,file='matterpower0.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(10,*) k(i),Pk(i)
!end do
!close(10)

!kg = 2
!do i = 1,150
!m = 10**kg
!write(10,*) m,nu(m),sig_2(m),nu_fnu(m)
!kg = kg+0.1
!end do
!call linear_pk(k,Pk)
!open(unit=10,file='matterpower01.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(10,*) k(i),Pk(i)
!end do
!close(10)


!type(HMdata) :: hm    
!open(unit=10,file='matterpower1.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(*,*) hm%k(i),hm%Pk(i)
!end do
!close(10)

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


!z = 0.1
!m = 1e6
!write(*,*) ukm(z,m)
!z = 10
!write(*,*) ukm(z,m)


!cz = 1.0
!kk = 0.1
!write(*,*) P1h(kk)

!z = 1.0
!kk = 0.1
!write(*,*) P2h(kk)

!z = 0.0
!call linear_pk(k,Pk)
!open(unit=10,file='matterpower0.dat',form='formatted',status='unknown')
!do i=1,mpts
!    write(10,*) k(i),Pk(i)
!end do
!close(10)

kk  = -2
open(unit=10,file='hey.dat',form='formatted',status='unknown')
do i=1,40
    kg = 10**kk
    write(10,*) kg,P1h(kg)
    kk = kk+0.1
end do
!close(10)
!kk  = -2
!open(unit=10,file='hey.dat',form='formatted',status='unknown')
!do i=1,56
!    kg = 10**kk
!    m = 1e11
!    mm = 1e13
!    mmm = 1e16
!    write(*,*) kg,ukm(kg,m),ukm(kg,mm),ukm(kg,mmm)
!    kk = kk+0.1
!end do
!close(10)






end program test
