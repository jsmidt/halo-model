module halo
use CAMB
implicit none

integer, parameter::dp = kind(1.0d0)
real(dp), parameter:: kmin = 2.0e-4
real(dp), parameter:: kmax = 9.610024d0
real(dp), parameter:: mmin = 2e9
real(dp), parameter:: mmax = 1e13
real(dp), parameter:: dlnk = 0.07d0
integer, parameter:: mpts = 155
real(dp), parameter:: rho_bar = 1.0d0
real(dp) :: z,kk,crv
real(dp) :: omegab, omegac, omegal, omegan, H0, YHe, Num_Nu_massless, Num_Nu_massive, omegam,rho_c,r_s

contains


! Tophat window function. See eq. 58 of astro-ph/0206508.
elemental real(dp) function win_top(x)
    real(dp),intent(in) :: x
    win_top = (3.0d0/x**3)*(sin(x)-x*cos(x))
end function win_top


! Variance in the initial density fluctuation (\sigma^2(m))
! Equation 58 of astro-ph/0206508.
real(dp) function sig_2(m)
    real(dp),intent(in) :: m
    real(dp) :: rombint,R
    R = (3.0d0*m/4.0d0/pi/rho_c)**0.3333333333d0
    CALL qromb(sig_2int,log(kmin),dlog(kmax),sig_2,R)
end function sig_2
real(dp) function sig_2int(lnk,R)
    real(dp),intent(in) :: lnk,R
    real(dp) :: k,P_lin
    k = exp(lnk)
    !P_lin = 1e9*k**(3.0)*exp(-sqrt(sqrt(k))*20.0)
    P_lin = 2.077e4*exp(-((lnk+4.072)/1.205)**2)+1.682e4*exp(-((lnk+4.589)/2.23)**2)
    sig_2int = (k**3.0d0/2.0d0/pi**2)*P_lin*win_top(k*R)**2
end function sig_2int

! The growth Function normalized s.t. D(z=0) = 1. Eq. 25 of astro-ph/0206508.
elemental real(dp) function growth(zz)
    real(dp),intent(in) :: zz
    real(dp) :: D,D0, omegamz, omegalz
    omegamz = omegam*(1.0+z)**3/(omegam*(1.0+z)**3+omegal)
    omegalz = 1.0-omegamz
    D = 5.0d0/2.0d0*omegamz/(1.0d0+z) &
        /(omegamz**(4.0/7.0d0)-omegalz+(1.0+0.5d0*omegamz)*(1.0+omegalz/70.0d0))
    D0 = 5.0d0/2.0d0*omegam &
        /(omegam**(4.0/7.0d0)-omegal+(1.0+0.5d0*omegam)*(1.0+omegal/70.0d0))
    growth = D/D0
end function growth

! Calculate nu for mass m and redshift z.  Equation 57 of astro-ph/0206508.
real(dp) function nu(m)
    real(dp), intent(in) :: m
    nu = (1.686470199841145d0/growth(z))**2/sig_2(m)
end function nu

! Calculate f_nu for mass m and redshift z. Equation 59 of astro-ph/0206508.
real(dp) function nu_fnu(m)
    real(dp), intent(in) :: m
    real(dp):: A,q
    A = 0.3222d0
    q = 0.75d0
    nu_fnu=A*(1.0d0+(q*nu(m))**0.3)*(q*nu(m)/2.0d0/pi)**0.5*exp(-q*nu(m)/2.0d0)
end function nu_fnu

!! Calculate the Fourier transform of the dark matter distribution u(k|m)
!! Equation 81 & 82 of astro-ph/0206508
!real(dp) function ukm(k,m)
!    real(dp), intent(in) :: k
!    real(dp) :: m
!    real(dp) :: ps, rs,c,zz,cp
!    real(dp) :: si_z,ci_z,si_cz,ci_cz
!    ps = 1.0
!    rs = 2e-5*(m)**0.3333
!    c = 9.0/(1.0+z)*m**(-0.13)*90
!    cp = c+1.0
!    zz  = k*rs
!    !ukm = 4.0d0*pi*ps*rs**3/m*(sin(zz)*(Si(cp*zz)-Si(zz)) )
!    !       - sin(c*zz)/(cp*zz)+cos(zz)*(Ci(cp*zz)-Ci(zz)) )
!    ukm = 4.0d0*pi*ps*rs**3/m*(sin(k*rs)*(Si((1.0d0+C)*k*rs)-Si(k*rs)) &
!        - sin(c*k*rs)/((1.0d0+c)*k*rs)+cos(k*rs)*(Ci((1.0d0+C)*k*rs)-Ci(k*rs)))
!end function ukm
!real(dp) function Ci(x)
!   real(dp), intent(in) :: x
!   real(dp) :: rombint,maxa
!   maxa = 1000.0
!   Ci = rombint(Cii,log(x),log(maxa),tol)
!end function Ci
!real(dp) function Si(x)
!   real(dp), intent(in) :: x
!   real(dp) :: rombint,maxa
!   maxa = 10000
!   Si = rombint(Sii,log(x),log(maxa),tol)
!end function Si
!real(dp) function Cii(t)
!   real(dp), intent(in) :: t
!   Cii = cos(exp(t))
!end function Cii
!real(dp) function Sii(t)
!   real(dp), intent(in) :: t
!   Sii = sin(exp(t))
!end function Sii


! Calculate the Fourier transform of the dark matter distribution u(k|m)
! Equation 80 & 74 of astro-ph/0206508
real(dp) function ukm(k,m)
    real(dp), intent(in) :: k
    real(dp) :: x,delta_c,r_vir,p_s,c,m,gg

    ! Consentration parameter. Eq. 78 of astro-ph/0206508
    c = 9.0/(1.0+z)*(m/4.822e13)**(-0.13)
    r_s = (m/4.0/pi/rho_c/(log(1.0+c)-c/(1.0+c)))**0.33333

    ! \Delta_c. Eq. 6 of arxiv:0907.4387
    x = omegam*(1.0+z)**3/(omegam*(1.0+z)**3+omegal)-1.0d0
    delta_c = 18.8*pi**2+82.0*x+39.0*x**2

    ! Virial radius. Eq. 4 of arxiv:0907.4387
    r_vir = c*r_s

    ! p_s.  Eq. 3 of arxiv:0907.4387
    !p_s = c**3*m/(4.0*pi*r_vir**3)/(log(1.0+c)-c/(1.0+c))
    !crv =  c/r_vir

    ! Want change of variables from r -> c*r/r_vir so that Eq. 1
    ! of arxiv:0907.4387 becomes simple p_s/r/(1+r)^2.  This means we
    ! ar integrating to c and integral must be multiplied by an additional
    ! p_s*c**3/r_vir**3
    !r_vir = 0.1
    gg = 1e-6
    CALL qromb(ukmi,log(gg),log(r_vir),ukm,k)
    ukm = ukm/m
end function ukm
real(dp) function ukmi(lnr,k)
    real(dp), intent(in) :: lnr,k
    real(dp) :: p,r
    r = exp(lnr)
    p = rho_c/(r/r_s)/(1.0+(r/r_s))**2
    ukmi = 4.0*pi*r**3*sin(k*r)/(k*r)*p
    !ukmi = 4.0*pi*r**3*p
end function ukmi




! Get linear power spectrum Pk at k from CAMB. 
subroutine linear_pk(k,Pk)
   real, allocatable,dimension(:),intent(inout) :: k,Pk
   integer :: i,error
   type(CAMBdata) :: P

   !allocate(k(mpts),Pk(mpts))
   do i=1,mpts
     k(i)=exp(log(kmin)+dlnk*(i-1))
   end do

   ! Get Default Cosmology Parameters.
   call CAMB_SetDefParams(P%Params)
   P%Params%WantTransfer= .true.
   P%Params%Transfer%redshifts=z
   P%Params%InitPower%ScalarPowerAmp(1) = 2.50e-9
   P%Params%omegab  = omegab
   P%Params%omegac  = omegac
   P%Params%omegav  = omegal
   P%Params%omegan  = omegan
   P%Params%H0      = H0
   P%Params%YHe     = 0.24
   P%Params%Num_Nu_massless = Num_Nu_massless
   P%Params%Num_Nu_massive  = Num_Nu_massive



   ! Get transfer functions and power spectra Pk
   call CAMB_GetTransfers(P%Params, P, error)
   call Transfer_GetMatterPower(P%MTrans,Pk,1,1,real(kmin),real(dlnk),mpts)
   !Pk = Pk/3.0e8

end subroutine linear_pk

! Get 1-Halo Term
real(dp) function P1h(kg)
    real(dp), intent(in) :: kg
    real(dp) :: rombint,tol2
    !P1h = rombint(P1hi,log(mmin),log(mmax),tol)
    CALL qromb(P1hi,log(mmin),dlog(mmax),P1h,kg)
end function P1h
real(dp) function P1hi(lnm,k)
    real(dp),intent(in) :: lnm
    real(dp) :: m,k
    m = exp(lnm)
    P1hi = m/rho_c*nu_fnu(m)*ukm(k,m)**2
end function P1hi

! Get 2-Halo Term
real(dp) function P2h(kg)
    real(dp), intent(in) :: kg
    real(dp) :: rombint
    P2h = rombint(P2hi1,log(mmin),log(mmax),tol) + rombint(P2hi2,log(mmin),log(mmax),tol)
end function P2h
real(dp) function P2hi1(lnm)
    real(dp), intent(in) :: lnm
    real(dp) :: m
    m = exp(lnm)
    P2hi1 = nu_fnu(m)*ukm(kk,m)
end function P2hi1
real(dp) function P2hi2(lnm)
    real(dp),intent(in) :: lnm
    real(dp) :: Phh,m
    Phh = 1.0
    m = exp(lnm)
    P2hi2 = nu_fnu(m)*ukm(kk,m)*Phh
end function P2hi2


end module halo
