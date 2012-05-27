module halo
use interp
use CAMB
use hm_init
implicit none

! Set types to be consistant with c double
!integer,parameter :: idp = selected_int_kind(13)
!integer,parameter :: dp = selected_real_kind(p=15,r=307)
!
real(dl), parameter:: kmin = 2.0e-4
real(dl), parameter:: kmax = 2000.0d0
real(dl), parameter:: mmin = 2e9
real(dl), parameter:: mmax = 1e12
real(dl), parameter:: dlnk = 0.105d0
integer, parameter:: mpts = 160
real(dl), parameter:: rho_bar = 1.0d0
!real(dl) :: z,kk,crv
!real(dl) :: omegab, omegac, omegal, omegan, H0, YHe, Num_Nu_massless, Num_Nu_massive, omegam,rho_c,r_s,kkkk
!real(dl) :: rho_c

contains

subroutine init_halo()
    real(dl) :: lnm,m
    integer :: i,N
    N = 190
    allocate(hm%m(N),hm%sig_2(N),hm%nu(N),hm%nu_fnu(N),hm%bias_1(N),hm%bias_2(N))
    lnm = 2
    open(unit=10,file='nu_fnu.dat',form='formatted',status='unknown')
    do i = 1,N
        m = 10**lnm
        hm%m(i) = m
        hm%sig_2(i) = sig_2(m)
        hm%nu(i) = nu(m)
        hm%nu_fnu(i) = nu_fnu(m)
        hm%bias_1(i) = bias_1(m)
        hm%bias_1(i) = bias_2(m)
        write(10,'(4Es12.4)') hm%m(i), hm%sig_2(i),hm%nu(i), hm%nu_fnu(i),hm%bias_1(i),hm%bias_2(i)
        lnm = lnm+0.08
    end do
    close(10)
end subroutine init_halo

! Tophat window function. See eq. 58 of astro-ph/0206508.
elemental real(dl) function win_top(x)
    real(dl),intent(in) :: x
    win_top = (3.0d0/x**3)*(sin(x)-x*cos(x))
end function win_top


! Variance in the initial density fluctuation (\sigma^2(m))
! Equation 58 of astro-ph/0206508.
real(dl) function sig_2(m)
    real(dl),intent(in) :: m
    real(dl) :: rombint,R
    R = (3.0d0*m/4.0d0/pi/hm%rho_c/0.3d0)**0.3333333333d0
    CALL qromb(sig_2int,log(kmin),dlog(kmax),sig_2,R)
end function sig_2
real(dl) function sig_2int(lnk,R)
    real(dl),intent(in) :: lnk,R
    real(dl) :: k,P_lin
    k = exp(lnk)
    P_lin = interpf(log(hm%k),dble(hm%Pk),lnk)
    sig_2int = (k**3.0d0/2.0d0/pi**2)*P_lin*win_top(k*R)**2
end function sig_2int

! The growth Function normalized s.t. D(z=0) = 1. Eq. 25 of astro-ph/0206508.
elemental real(dl) function growth(zz)
    real(dl),intent(in) :: zz
    real(dl) :: D,D0, omegamz, omegalz,omegam,omegal,z
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegav

    omegamz = omegam*(1.0+hm%z)**3/(omegam*(1.0+hm%z)**3+omegal)
    omegalz = 1.0-omegamz
    D = 5.0d0/2.0d0*omegamz/(1.0d0+z) &
        /(omegamz**(4.0/7.0d0)-omegalz+(1.0+0.5d0*omegamz)*(1.0+omegalz/70.0d0))
    D0 = 5.0d0/2.0d0*omegam &
        /(omegam**(4.0/7.0d0)-omegal+(1.0+0.5d0*omegam)*(1.0+omegal/70.0d0))
    growth = D/D0
end function growth

! Calculate nu for mass m and redshift z.  Equation 57 of astro-ph/0206508.
real(dl) function nu(m)
    real(dl), intent(in) :: m
    real(dl) :: z
    nu = (1.686470199841145d0/growth(hm%z))**2/sig_2(m)
end function nu

! Calculate f_nu for mass m and redshift z. Equation 59 of astro-ph/0206508.
real(dl) function nu_fnu(m)
    real(dl), intent(in) :: m
    real(dl):: A,q
    A = 0.3222d0
    q = 0.75d0
    nu_fnu=A*(1.0d0+(q*nu(m))**0.3)*(q*nu(m)/2.0d0/pi)**0.5*exp(-q*nu(m)/2.0d0)
end function nu_fnu

! The bias parameters.  Eq. 68 of astro-ph/0206508
real(dl) function bias_1(m)
real(dl), intent(in) :: m
real(dl) :: q,p,ep1,E1
q = 0.75d0
p = 0.3
ep1 = (q*nu(m) - 1.0)/(1.686470199841145d0/growth(hm%z))
E1  = 2.0*p/(1.686470199841145d0/growth(hm%z))/(1.0+(q*nu(m))**p)
bias_1 = 1.0+ep1+E1
end function bias_1
real(dl) function bias_2(m)
real(dl), intent(in) :: m
real(dl) :: q,p,ep1,E1,ep2,E2
q = 0.75d0
p = 0.3
ep1 = (q*nu(m)-1.0)/(1.686470199841145d0/growth(hm%z))
ep2 = (q*nu(m))*(q*nu(m)-3.0)/(1.686470199841145d0/growth(hm%z))**2
E1  = 2.0*p/(1.686470199841145d0/growth(hm%z))/(1.0+(q*nu(m))**p)
E2  = E1*((1.0+2.0*p)/(1.686470199841145d0/growth(hm%z))+2*ep1)
bias_2 = 2.0*(1.0-17.0/21.0)*(ep1+E1)+ep2+E2
end function bias_2



! Calculate the Fourier transform of the dark matter distribution u(k|m)
! Equation 80 & 74 of astro-ph/0206508
real(dl) function ukm(k,m)
    real(dl), intent(in) :: k
    real(dl) :: x,delta_c,r_vir,p_s,c,m,gg,rombint,Ez2,omegam,omegal,ms
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegan

    ! Consentration parameter. Eq. 78 of astro-ph/0206508
    ms = interpf(log(hm%nu),hm%m,log(1.0d0))
    c = 9.0/(1.0+hm%z)*(m/ms)**(-0.13d0)

    ! \Delta_c & E(z)^2. Eq. 5-6 of arxiv:0907.4387
    x = omegam*(1.0+hm%z)**3/(omegam*(1.0+hm%z)**3+omegal)-1.0d0
    delta_c = 18.0*pi**2+82.0*x-39.0*x**2
    Ez2 = omegam*(1.0+hm%z)**3+omegal

    ! Virial radius. Eq. 4 of arxiv:0907.4387
    r_vir = (3.0*m/(hm%rho_c*Ez2*4.0*pi*delta_c))**(1.0/3.0)

    ! r_s = r_vir/c.  Above eq. 77 of astro-ph/0206508
    hm%r_s = r_vir/c

    ! p_s  Eq. 3 of arxiv:0907.4387
    p_s = c**3*m/(4.0*pi*r_vir**3)/(log(1.0+c)-c/(1.0+c))

    gg = 1.0e-4
    CALL qromb(ukmi,log(gg),log(r_vir),ukm,k)
    ukm = p_s*ukm/m
end function ukm
real(dl) function ukmi(lnr,k)
    real(dl), intent(in) :: lnr,k
    real(dl) :: p,r
    r = exp(lnr)
    p = 1.0/(r/hm%r_s)/(1.0+(r/hm%r_s))**2
    ukmi = 4.0*pi*r**3*sin(k*r)/(k*r)*p
end function ukmi


! Get linear power spectrum Pk at k from CAMB. 
subroutine linear_pk(k,Pk)
   real(dl), allocatable,dimension(:),intent(inout) :: k
   real(sp), allocatable,dimension(:),intent(inout) :: Pk
   integer :: i,error
   type(CAMBdata) :: P

   ! Get k.
   do i=1,mpts
     k(i)=exp(log(kmin)+dlnk*(i-1))
   end do

   ! Get P(k)
   call CAMB_GetTransfers(hm%Params, P, error)
   call Transfer_GetMatterPower(P%MTrans,Pk,1,1,real(kmin),real(dlnk),mpts)

end subroutine linear_pk

! Get 1-Halo Term
real(dl) function P1h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint,tol2
    !P1h = rombint(P1hi,log(mmin),log(mmax),tol)
    CALL qromb(P1hi,log(mmin),dlog(mmax),P1h,kg)
end function P1h
real(dl) function P1hi(lnm,k)
    real(dl),intent(in) :: lnm
    real(dl) :: m,k,inu_fnu
    m = exp(lnm)
    inu_fnu = interpf(log(hm%m),hm%nu_fnu,lnm)
    P1hi = m/hm%rho_c*inu_fnu*ukm(k,m)**2
    !P1hi = m/hm%rho_c*nu_fnu(m)*ukm(k,m)**2
end function P1hi

! Get 2-Halo Term
real(dl) function P2h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint
    CALL qromb(P2h_int1,log(mmin),dlog(mmax),P2h,kg)
end function P2h
real(dl) function P2h_int1(lnm,kg)
    real(dl), intent(in) :: lnm,kg
    real(dl) :: m,inu_fnu,ibias_1
    m = exp(lnm)
    inu_fnu = interpf(log(hm%m),hm%nu_fnu,lnm)
    ibias_1 = interpf(log(hm%m),hm%bias_1,lnm)
    CALL qromb(P2h_int2,log(mmin),dlog(mmax),P2h_int1,kg)
    P2h_int1 = ibias_1*inu_fnu*ukm(kg,m)*P2h_int1
end function P2h_int1
real(dl) function P2h_int2(lnm,kg)
    real(dl),intent(in) :: lnm,kg
    real(dl) :: m,inu_fnu,ibias_2
    m = exp(lnm)
    inu_fnu = interpf(log(hm%m),hm%nu_fnu,lnm)
    ibias_2 = interpf(log(hm%m),hm%bias_2,lnm)
    P2h_int2 = inu_fnu*ukm(kg,m)*interpf(log(hm%k),dble(hm%Pk),log(kg))
    !P2h_int2 = bias_2(m)*nu_fnu(m)*ukm(kg,m)*interpf(log(hm%k),dble(hm%Pk),log(kg))
end function P2h_int2


end module halo
