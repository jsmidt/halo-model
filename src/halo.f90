module halo
use interp
use CAMB
use hm_init
implicit none


contains

subroutine init_halo()
    real(dl) :: lnm,m,lnk,k
    integer :: i,j
    allocate(hm%m(hm%mpts),hm%sig_2(hm%mpts),hm%nu_m(hm%mpts),hm%nu_fnum(hm%mpts),hm%bias_1(hm%mpts),hm%bias_2(hm%mpts))
    allocate(hm%sp_nu_fnum(hm%mpts),hm%sp_bias_1(hm%mpts))
    allocate(hm%ukm(hm%mpts,hm%mpts),hm%sp_ukm(hm%mpts,hm%mpts))
    open(unit=10,file='output/' // trim(hm%run_name) // '_nu_fnu.dat',form='formatted',status='unknown')
    do i = 1,hm%mpts
        m = 10**(dlog10(hm%mmin)+(dlog10(hm%mmax)-dlog10(hm%mmin))/(hm%mpts-1.0)*(i-1.0))
        hm%m(i) = m
        hm%sig_2(i) = sig_2(m)
        hm%nu_m(i) = nu_m(m)
        hm%nu_fnum(i) = nu_fnu(hm%nu_m(i))
        hm%bias_1(i) = bias_1(hm%nu_m(i))
        write(10,'(6Es13.3)') hm%m(i), hm%sig_2(i),hm%nu_m(i), hm%nu_fnum(i), hm%bias_1(i)
    end do
    close(10)

    do j=1,size(hm%k)
        do i=1,size(hm%nu_m)
            if (hm%k(j) .lt. 3000) then
            hm%ukm(j,i) = ukm(hm%k(j),hm%m(i))
            else
                hm%ukm(j,i) = 0
                endif
        enddo
    enddo


    call spline(dlog(hm%nu_m),hm%nu_fnum,size(hm%nu_m),1d40,1d40,hm%sp_nu_fnum)
    call spline(dlog(hm%nu_m),hm%bias_1,size(hm%nu_m),1d40,1d40,hm%sp_bias_1)
    call splie2(dlog(hm%k), dlog(hm%m),hm%ukm, hm%mpts,hm%mpts, hm%sp_ukm)

    ! Write out redshift and cosmology being used
    write(*,*) ' ' 
    write(*,'("   Halo model initialized for redshift z = ", F4.2, ". With cosmology:")') hm%z
    write(*,'("   \Omega_b = ", F5.3, ". \Omega_c = ", F5.3, ". \Omega_l = ", &
    F5.3, ". \Omega_nu = ", F5.3)') hm%Params%omegab,hm%Params%omegac,hm%Params%omegav,hm%Params%omegan
    write(*,'("    H_0 = ", F4.1, ". Num_Nu_massless = ", F4.2, ". Num_Nu_massive = ", &
    F4.2, ". YHe = ", F4.2)') hm%Params%H0,hm%Params%Num_Nu_massless,real(hm%Params%Num_Nu_massive),hm%Params%YHe
    write(*,*) ' ' 
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
    real(dl) :: rombint,R,omegam,omegal,omegamz
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegav
    omegamz = omegam*(1.0+hm%z)**3/(omegam*(1.0+hm%z)**3+omegal)
    R = (3.0d0*m/(4.0d0*pi*hm%rho_c*omegamz))**(1.0d0/3.0d0)
    CALL qromb(sig_2int,dlog(hm%kmin),dlog(hm%kmax),sig_2,R)
end function sig_2
real(dl) function sig_2int(lnk,R)
    real(dl),intent(in) :: lnk,R
    real(dl) :: k,P_lin
    k = exp(lnk)
    P_lin = interpf(dlog(hm%k),dble(hm%Pk),lnk)
    sig_2int = (k**3/2.0d0/pi**2)*P_lin*win_top(k*R)**2
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
real(dl) function nu_m(m)
    real(dl), intent(in) :: m
    real(dl) :: z
    !nu = (1.686470199841145d0/growth(hm%z))**2/sig_2(m)
    nu_m = (1.686470199841145d0)**2/sig_2(m)
end function nu_m

! Calculate f_nu for mass m and redshift z. Equation 59 of astro-ph/0206508.
real(dl) function nu_fnu(nu)
    real(dl), intent(in) :: nu
    real(dl):: A,q,p,nup
    !A = 0.3222d0
    A = 0.413078984
    q = 0.75d0
    p = 0.3d0
    nup = q*nu
    nu_fnu=A*(1.0d0+nup**(-p))*(nup/2.0d0/pi)**0.5*exp(-nup/2.0d0)
end function nu_fnu


! Calculate f_nu for mass m and redshift z. Equation 59 of astro-ph/0206508.
real(dl) function nfnu(m)
    real(dl), intent(in) :: m
    real(dl) :: rombint
    nfnu = rombint(nfnui,dlog(4.434d-019),dlog(2.386d06),tol)
end function nfnu


! Calculate f_nu for mass m and redshift z. Equation 59 of astro-ph/0206508.
real(dl) function nfnui(lnm)
    real(dl), intent(in) :: lnm
    real(dl):: q,p,num
    q = 0.75d0
    p = 0.3d0
    num = exp(lnm)
    nfnui=(1.0d0+(q*num)**(-p))*(q*num/2.0d0/pi)**0.5*exp(-q*num/2.0d0)
end function nfnui


! The bias parameters.  Eq. 68 of astro-ph/0206508
real(dl) function bias_1(nu)
real(dl), intent(in) :: nu
real(dl) :: q,p,ep1,E1,delc,nup
q = 0.75d0
p = 0.3
nup = nu*q
!delc = (1.686470199841145d0/growth(hm%z))
delc = (1.686470199841145d0)
ep1 = (nup-1.0)/delc
E1  = 2.0*p/delc/(1.0+nup**p)
bias_1 = 1.0+ep1+E1
end function bias_1
!real(dl) function bias_2(m)
!real(dl), intent(in) :: m
!real(dl) :: q,p,ep1,E1,ep2,E2
!q = 0.75d0
!p = 0.3
!ep1 = (q*nu(m)-1.0)/(1.686470199841145d0/growth(hm%z))
!ep2 = (q*nu(m))*(q*nu(m)-3.0)/(1.686470199841145d0/growth(hm%z))**2
!E1  = 2.0*p/(1.686470199841145d0/growth(hm%z))/(1.0+(q*nu(m))**p)
!E2  = E1*((1.0+2.0*p)/(1.686470199841145d0/growth(hm%z))+2*ep1)
!bias_2 = 2.0*(1.0-17.0/21.0)*(ep1+E1)+ep2+E2
!end function bias_2



! Calculate the Fourier transform of the dark matter distribution u(k|m)
! Equation 80 & 74 of astro-ph/0206508
real(dl) function ukm(k,m)
    real(dl), intent(in) :: k
    real(dl) :: x,delta_c,r_vir,p_s,c,m,gg,rombint,Ez2,omegam,omegal,ms
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegan

    ! Consentration parameter. Eq. 78 of astro-ph/0206508
    ms = interpf(log(hm%nu_m),hm%m,log(1.0d0))
    !ms = 3.6d12
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
    p_s = c**3*m/(4.0*pi*r_vir**3)/(dlog(1.0+c)-c/(1.0+c))

    gg = 1.0e-4
    CALL qromb(ukmi,dlog(gg),dlog(r_vir),ukm,k)
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
   real(dl) :: dlnk
   type(CAMBdata) :: P

   ! Get k.
   do i=1,hm%mpts
     k(i)=exp(dlog(hm%kmin)+(dlog(hm%kmax)-dlog(hm%kmin))/(hm%mpts-1.0)*(i-1.0))
   end do
   dlnk = (dlog(hm%kmax)-dlog(hm%kmin))/(hm%mpts-1.0)

   ! Get P(k)
   call CAMB_GetTransfers(hm%Params, P, error)
   call Transfer_GetMatterPower(P%MTrans,Pk,1,1,real(hm%kmin),real(dlnk),hm%mpts)

end subroutine linear_pk

! Get 1-Halo Term
real(dl) function P1h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint,tol2
    CALL qromb(P1hi,dlog(minval(hm%nu_m)),dlog(maxval(hm%nu_m)),P1h,kg)
end function P1h
real(dl) function P1hi(lnnu,k)
    real(dl),intent(in) :: lnnu
    real(dl) :: m,k,inu_fnu,omegam,omegal,omegamz,nu,iukm
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegav
    omegamz = omegam*(1.0+hm%z)**3/(omegam*(1.0+hm%z)**3+omegal)
    nu = exp(lnnu)
    m = interpf(dlog(hm%nu_m),hm%m,lnnu)
  call splint(dlog(hm%nu_m),hm%nu_fnum,hm%sp_nu_fnum,size(hm%nu_m),lnnu,inu_fnu)
call splin2(dlog(hm%k),dlog(hm%m),hm%ukm,hm%sp_ukm,hm%mpts,hm%mpts,dlog(k),dlog(m),iukm)
    P1hi = m/hm%rho_c/omegamz*inu_fnu*iukm**2
end function P1hi

! Get 2-Halo Term
real(dl) function P2h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint
    CALL qromb(P2hi,log(minval(hm%nu_m)),dlog(maxval(hm%nu_m)),P2h,kg)
    P2h = P2h**2*interpf(dlog(hm%k),dble(hm%Pk),log(kg))
end function P2h
real(dl) function P2hi(lnnu,k)
    real(dl), intent(in) :: lnnu,k
    real(dl) :: m,inu_fnu,ibias_1,nu,iukm
    nu = exp(lnnu)
    m = interpf(dlog(hm%nu_m),hm%m,lnnu)
  call splint(dlog(hm%nu_m),hm%nu_fnum,hm%sp_nu_fnum,size(hm%nu_m),lnnu,inu_fnu)
  call splint(dlog(hm%nu_m),hm%bias_1,hm%sp_bias_1,size(hm%nu_m),lnnu,ibias_1)
call splin2(dlog(hm%k),dlog(hm%m),hm%ukm,hm%sp_ukm,hm%mpts,hm%mpts,dlog(k),dlog(m),iukm)
    P2hi = ibias_1*inu_fnu*iukm
end function P2hi


end module halo
