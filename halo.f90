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
real(dp) :: z,kk

contains

! Tophat window function. See eq. 58 of astro-ph/0206508.
elemental real(dp) function w_top(x)
    real(dp),intent(in) :: x
    w_top = (3.0d0/x**3)*(sin(x)-x*cos(x))
end function w_top


! Variance in the initial density fluctuation (\sigma^2(m))
! Equation 58 of astro-ph/0206508.
real(dp) function sig_2(m)
    real(dp),intent(in) :: m
    real(dp) :: rombint,tol2,R
    R = (3.0d0*m/4.0d0/pi/rho_bar)**0.3333333333d0
    tol2 = 1e-6
    !sig_2 = qromb(sig_2int,log(kmin),log(kmax),tol2)
    CALL qromb(sig_2int,log(kmin),dlog(kmax),sig_2,R)
end function sig_2
real(dp) function sig_2int(lnk,R)
    real(dp),intent(in) :: lnk,R
    real(dp) :: k,P_lin
    k = exp(lnk)
    P_lin = 1e9*k**(3.0)*exp(-sqrt(sqrt(k))*20.0)
    sig_2int = (k**3.0d0/2.0d0/pi**2)*P_lin*w_top(k*R)**2
end function sig_2int

! Calculate nu for mass m and redshift z.  Equation 57 of astro-ph/0206508.
real(dp) function nu(m)
    real(dp), intent(in) :: m
    nu = 1.686470199841145d0*(1.0d0+z)*sig_2(m)
end function nu

! Calculate f_nu for mass m and redshift z. Equation 59 of astro-ph/0206508.
real(dp) function fnu(m)
    real(dp), intent(in) :: m
    real(dp):: A,q
    A = 0.3222d0
    q = 0.75d0
    fnu=A*(1.0d0+(q*nu(m))**0.3)*(q*nu(m)/2.0d0/pi)**0.5*exp(-q*nu(m)/2.0d0)/nu(m)
end function fnu

! Calculate the Fourier transform of the dark matter distribution u(k|m)
! Equation 81 & 82 of astro-ph/0206508
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
!    !ukm = 4.0d0*pi*ps*rs**3/m*(cos(zz)*Ci(zz))
!    ukm = Ci(zz)
!    !ukm = 4.0d0*pi*ps*rs**3/m*(sin(zz)*(Si(cp*zz)-Si(zz)) )
!    !       - sin(c*zz)/(cp*zz)+cos(zz)*(Ci(cp*zz)-Ci(zz)) )
!    !ukm = 4.0d0*pi*ps*rs**3/m*(sin(k*rs)*(Si((1.0d0+C)*k*rs)-Si(k*rs)) &
!    !    - sin(c*k*rs)/((1.0d0+c)*k*rs)+cos(k*rs)*(Ci((1.0d0+C)*k*rs)-Ci(k*rs)))
!end function ukm
!real(dp) function Ci(x)
!   real(dp), intent(in) :: x
!   real(dp) :: rombint,maxa
!   maxa = 1000.0
!   !Ci = rombint(Cii,x,maxa,tol)
!    CALL qromb(Cii,x,maxa,Ci)
!end function Ci
!real(dp) function Si(x)
!   real(dp), intent(in) :: x
!   real(dp) :: rombint,maxa
!   maxa = 10000
!   Si = rombint(Sii,x,maxa,tol)
!end function Si
!real(dp) function Cii(t)
!   real(dp), intent(in) :: t
!   Cii = -sin(t+3.14159/2.0)/t
!end function Cii
!real(dp) function Sii(t)
!   real(dp), intent(in) :: t
!   Sii = sin(t)/t
!end function Sii


! Calculate the Fourier transform of the dark matter distribution u(k|m)
! Equation 81 & 82 of astro-ph/0206508
real(dp) function ukm(k,m)
    real(dp), intent(in) :: k
    real(dp) :: m
    real(dp) :: ps, rs,c,zz,cp
    real(dp) :: si_z,ci_z,si_cz,ci_cz,x
    ps = 1.0
    rs = 2e-5*(m)**0.3333
    c = 9.0/(1.0+z)*m**(-0.13)*90
    cp = c+1.0
    zz  = k*rs
    x = 1e-2
    CALL qromb(ukmi,x,rs,ukm,k)
    ukm = 2.25*ukm/m
end function ukm

real(dp) function ukmi(r,k)
    real(dp), intent(in) :: r,k
    real(dp) :: p
    p = 1.0/r/(1.0+r)**3
    ukmi = 4.0*pi*r**2*sin(k*r)/(k*r)*p
    !ukmi = 4.0*pi*r**2*sin(k*r)
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

   ! Get transfer functions and power spectra Pk
   call CAMB_GetTransfers(P%Params, P, error)
   call Transfer_GetMatterPower(P%MTrans,Pk,1,1,real(kmin),real(dlnk),mpts)
   Pk = Pk/3.0e8

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
    P1hi = nu(m)*fnu(m)*(m/rho_bar)*ukm(k,m)**2
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
    P2hi1 = fnu(m)*ukm(kk,m)
end function P2hi1
real(dp) function P2hi2(lnm)
    real(dp),intent(in) :: lnm
    real(dp) :: Phh,m
    Phh = 1.0
    m = exp(lnm)
    P2hi2 = fnu(m)*ukm(kk,m)*Phh
end function P2hi2


end module halo