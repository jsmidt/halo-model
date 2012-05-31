module hm_init
use Precision
use CAMB
implicit none

! Use camb derived type for halo-model parameters HP.
!integer,parameter :: dp = selected_real_kind(p=15,r=307)

! The global halo-model parameter structure. Consistant with CAMB's. 
type, extends(CAMBdata) :: HMdata
   real(dl) :: rho_c
   real(dl) :: r_s
   real(dl) :: z
   real(dl) :: kmin,kmax
   real(dl) :: mmin,mmax
   integer :: mpts
   real(dl), allocatable,dimension(:) :: k,m,nu_m,sig_2,bias_1,bias_2,nu_fnum
   real(dl), allocatable,dimension(:) :: sp_bias_1,sp_nu_fnum
   real(sp), allocatable,dimension(:) :: Pk
   real(dl), allocatable,dimension(:) :: lnk2d,lnm2d
   real(dl), allocatable,dimension(:,:) :: ukm,sp_ukm
   character*8 :: run_name 
end type HMdata

type(HMdata) :: hm

contains

subroutine init_params()
    hm%kmin = 8.4d-5
    hm%kmax = 4.0d6
    hm%mmin = 1.0d-8
    hm%mmax = 1.0d16
    hm%mpts = 100


    ! Get Default Cosmology Parameters. Match to WMAP 7.
    call CAMB_SetDefParams(hm%Params)

    hm%run_name = 'test'
    hm%z = 0.0
    hm%Params%WantTransfer= .true.
    hm%Params%WantCls= .false.              
    hm%Params%Transfer%redshifts=hm%z       ! Redshift
    hm%Params%InitPower%ScalarPowerAmp(1) = 2.45e-9 ! Scalar amplitude
    hm%Params%omegab  = 0.0455              ! Baryon density
    hm%Params%omegac  = 0.227               ! CDM density
    hm%Params%omegav  = 0.728               ! Dark energy density  
    hm%Params%omegan  = 0                   ! Neutrino density
    hm%Params%H0      = 70.2                ! Hubble Constant
    hm%Params%YHe     = 0.24                ! Helium fraction
    hm%Params%Num_Nu_massless = 3.04        ! Massless neutrinos
    hm%Params%Num_Nu_massive  = 0           ! # Massive neutrinos
    hm%Params%InitPower%an(1)  = 0.961      ! Spectral Index


    ! To get the critical density rho_c we need to calculate 3H^2/(8\piG).  
    ! We want units of M_\sun/Mpc^3. We know H0 has units of km/s/Mpc and 
    ! G = 4.302e-9 Mpc(km/s)^2/M_\sun So that 3H^2/(8\piG) has the right units.
    hm%rho_c =3.0d0*hm%Params%H0**2/8.0d0/pi/4.302e-9/(hm%Params%H0/100.0)**2* &
     ((hm%Params%omegab+hm%Params%omegac)*(1.0+hm%z)**3.0+hm%Params%omegav) 
end subroutine init_params

end module hm_init
