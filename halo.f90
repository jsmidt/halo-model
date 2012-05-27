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
real(dl), parameter:: kmax = 9.610024d0
real(dl), parameter:: mmin = 2e9
real(dl), parameter:: mmax = 1e12
real(dl), parameter:: dlnk = 0.07d0
integer, parameter:: mpts = 155
real(dl), parameter:: rho_bar = 1.0d0
!real(dl) :: z,kk,crv
!real(dl) :: omegab, omegac, omegal, omegan, H0, YHe, Num_Nu_massless, Num_Nu_massive, omegam,rho_c,r_s,kkkk
!real(dl) :: rho_c

contains


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
    R = (3.0d0*m/4.0d0/pi/hm%rho_c)**0.3333333333d0
    CALL qromb(sig_2int,log(kmin),dlog(kmax),sig_2,R)
end function sig_2
real(dl) function sig_2int(lnk,R)
    real(dl),intent(in) :: lnk,R
    real(dl) :: k,P_lin
    k = exp(lnk)
    P_lin = interpf(log(hm%k),dble(hm%Pk),lnk)
    !P_lin = 2.077e4*exp(-((lnk+4.072)/1.205)**2)+1.682e4*exp(-((lnk+4.589)/2.23)**2)
    !P_lin = 2.077e4*exp(-((lnk+4.072)/1.205)**2)+1.682e4*exp(-((lnk+4.589)/2.23)**2)
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

!! Calculate the Fourier transform of the dark matter distribution u(k|m)
!! Equation 81 & 82 of astro-ph/0206508
!real(dl) function ukm(k,m)
!    real(dl), intent(in) :: k,m
!    real(dl) :: x,delta_c,r_vir,p_s,gg,c
!
!    ! Consentration parameter. Eq. 78 of astro-ph/0206508
!    c = 9.0/(1.0+z)*(m/4.822e13)**(-0.13)
!    r_s = (m/4.0/pi/rho_c/(log(1.0+c)-c/(1.0+c)))**0.3333333d0
!
!    !ukm = 4.0d0*pi*rho_c*r_s**3/m*(sin(k*r_s)*(Si((1.0d0+c)*k*r_s)-Si(k*r_s)) &
!    !- sin(c*k*r_s)/((1.0d0+c)*k*r_s)+cos(k*r_s)*(Ci((1.0d0+c)*k*r_s)-Ci(k*r_s)))
!    ukm = Ci(k*r_s)
!end function ukm
!real(dl) function Ci(x)
!   real(dl), intent(in) :: x
!   real(dl) :: rombint,maxa
!   maxa = 30000.0
!   Ci = rombint(Cii,log(x),log(maxa),tol)
!end function Ci
!real(dl) function Si(x)
!   real(dl), intent(in) :: x
!   real(dl) :: rombint,maxa
!   maxa = 1e-10
!   Si = rombint(Sii,log(maxa),log(x),tol)
!end function Si
!real(dl) function Cii(t)
!   real(dl), intent(in) :: t
!   Cii = -cos(exp(t))
!end function Cii
!real(dl) function Sii(t)
!   real(dl), intent(in) :: t
!   Sii = sin(exp(t))
!end function Sii


! Calculate the Fourier transform of the dark matter distribution u(k|m)
! Equation 80 & 74 of astro-ph/0206508
real(dl) function ukm(k,m)
    real(dl), intent(in) :: k
    real(dl) :: x,delta_c,r_vir,p_s,c,m,gg,rombint,Ez2,omegam,omegal,z
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegan

    ! Consentration parameter. Eq. 78 of astro-ph/0206508
    c = 9.0/(1.0+z)*(m/1.259e13)**(-0.13d0)

    ! \Delta_c & E(z)^2. Eq. 5-6 of arxiv:0907.4387
    x = omegam*(1.0+hm%z)**3/(omegam*(1.0+hm%z)**3+omegal)-1.0d0
    delta_c = 18.0*pi**2+82.0*x-39.0*x**2
    Ez2 = omegam*(1.0+z)**3+omegal

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
    real(dl) :: m,k
    m = exp(lnm)
    P1hi = m/hm%rho_c*nu_fnu(m)*ukm(k,m)**2
end function P1hi

! Get 2-Halo Term
real(dl) function P2h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint
    P2h = rombint(P2hi1,log(mmin),log(mmax),tol) + rombint(P2hi2,log(mmin),log(mmax),tol)
end function P2h
real(dl) function P2hi1(lnm)
    real(dl), intent(in) :: lnm
    real(dl) :: m,kk
    m = exp(lnm)
    P2hi1 = nu_fnu(m)*ukm(kk,m)
end function P2hi1
real(dl) function P2hi2(lnm)
    real(dl),intent(in) :: lnm
    real(dl) :: Phh,m,kk
    Phh = 1.0
    m = exp(lnm)
    P2hi2 = nu_fnu(m)*ukm(kk,m)*Phh
end function P2hi2


end module halo
