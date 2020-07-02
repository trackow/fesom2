! Modules of cell-vertex ocean model
! S. Danilov, 2012 (sergey.danilov@awi.de)
! SI units are used
  
!==========================================================
MODULE o_PARAM
integer, parameter            :: WP=8        ! Working precision
integer		                  :: mstep
real(kind=WP), parameter      :: pi=3.14159265358979
real(kind=WP), parameter      :: rad=pi/180.0_WP
real(kind=WP), parameter      :: density_0=1030.0_WP
real(kind=WP), parameter      :: density_0_r=1.0_WP/density_0 ! [m^3/kg]         
real(kind=WP), parameter      :: g=9.81_WP
real(kind=WP), parameter      :: r_earth=6367500.0_WP
real(kind=WP), parameter      :: omega=2*pi/(3600.0_WP*24.0_WP)
real(kind=WP), parameter      :: vcpw=4.2e6   ![J/m^3/K] water heat cap
real(kind=WP), parameter      :: inv_vcpw = 1._WP / vcpw  ! inverse, to replace divide by multiply
real(kind=WP), parameter      :: small=1.0e-8 !small value

real(kind=WP)                 :: C_d= 0.0025_WP ! Bottom drag coefficient
real(kind=WP)	              :: kappa=0.4      !von Karman's constant
real(kind=WP)                 :: mix_coeff_PP=0.01_WP   ! mixing coef for PP scheme
real(kind=WP)                 :: A_hor=100.0_WP		! Horizontal harm. visc.    
real(kind=WP)                 :: A_hor_max=1500.0_WP	! Maximum viscosity allowed (to limit Smag and Leith
							! contributions when they are too large
real(kind=WP)                 :: Leith_c=0.		! Leith viscosity. It needs vorticity, which is only computed for 
							! the vector invariant form of momentum advection (mom_adv=4)
real(kind=WP)                 :: tau_c=0.4		! Controls the strength of filters. Should be about 0.4
real(kind=WP)                 :: A_ver=0.001_WP ! Vertical harm. visc.
real(kind=WP)                 :: Div_c=0.5_WP
real(kind=WP)                 :: Smag_c=0.0_WP  ! 0.2   ! (C/pi)^2
real(kind=WP)                 :: Abh0=8.0e12    ! 
logical                       :: laplacian=.false.
logical                       :: biharmonic=.true.
real(kind=WP)                 :: K_hor=10._WP
real(kind=WP)                 :: K_ver=0.00001_WP
real(kind=WP)                 :: scale_area=2.0e8
real(kind=WP)                 :: surf_relax_T= 0.0_WP
real(kind=WP)                 :: surf_relax_S= 10.0_WP/(60*3600.0_WP*24)
logical                       :: balance_salt_water =.true.
real(kind=WP)                 :: clim_relax= 1.0_WP/(10*3600.0_WP*24)
real(kind=WP)                 :: clim_decay, clim_growth
                                 ! set to 0.0 if no relaxation
logical                       :: ref_sss_local=.false.
real(kind=WP)                 :: ref_sss=34.7
logical                       :: Fer_GM =.false.  !flag for Ferrari et al. (2010) GM scheme
real(kind=WP)                 :: K_GM=1000.
logical			      :: scaling_Ferreira   =.true.
logical			      :: scaling_Rossby     =.false.
logical			      :: scaling_resolution =.true.
logical			      :: scaling_FESOM14    =.false.
logical			      :: Redi               =.false.  !flag for Redi scheme

real(kind=WP)                 :: visc_sh_limit=5.0e-3      !for KPP, max visc due to shear instability
real(kind=WP)                 :: diff_sh_limit=5.0e-3      !for KPP, max diff due to shear instability
logical                       :: Kv0_const=.true.		    !use Kv0 varying with depth and latitude 
logical                       :: double_diffusion=.false.  !for KPP,dd switch
                                 ! KPP parametrization
 character(5)                 :: mix_scheme='KPP'	   !'KPP','PP'
real(KIND=WP)                 :: Ricr   = 0.3_WP  ! critical bulk Richardson Number
real(KIND=WP)                 :: concv  = 1.6_WP  ! constant for pure convection (eqn. 23) (Large 1.5-1.6; MOM default 1.8)

logical                       :: hbl_diag =.false.   ! writen boundary layer depth

! Time stepping                               
! real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
                                                 ! elevation and divergence
real(kind=WP)                 :: epsilon=0.1_WP  ! AB2 offset 
! Tracers
logical                       :: i_vert_diff= .true.
logical                       :: i_vert_visc= .true.
integer                       :: tracer_adv=1

logical                       :: w_split  =.false.
real(kind=WP)                 :: w_exp_max=1.e-5_WP

logical                       :: SPP=.false.

				! 1 MUSCL
				! 2 MUSCL-FCT
integer	                       :: num_tracers=2
integer, dimension(100)        :: tracer_ID  = RESHAPE((/0, 1/), (/100/), (/0/)) ! ID for each tracer for treating the initialization and surface boundary condition
                                                                                 ! 0=temp, 1=salt etc.

! Momentum
logical                       :: free_slip=.false.
                                ! false=no slip 
integer                       :: mom_adv=2
                                ! 1 vector control volumes, p1 velocities
				! 2 scalar control volumes  
				! 3 vector invariant 

logical                       :: open_b=.false.   ! Reserved    

logical                       :: mo_on=.true. !.false. !Monin-Obukhov
real(kind=WP) :: modiff=0.01                   !for PP, mixing coefficient within MO length

  ! *** active tracer cutoff
logical          :: limit_salinity=.true.         !set an allowed range for salinity
real(kind=WP)    :: salinity_min=5.0              !minimal salinity 
real(kind=WP)    :: coeff_limit_salinity=0.0023   !m/s, coefficient to restore s to s_min

  namelist /tracer_cutoff/ limit_salinity, salinity_min, coeff_limit_salinity

! *** others ***
 real(kind=WP)                        :: time_sum=0.0 ! for runtime estimate


 NAMELIST /oce_dyn/ C_d, A_ver, laplacian, A_hor, A_hor_max, Leith_c, tau_c, Div_c, Smag_c, &
                    biharmonic, Abh0, scale_area, mom_adv, free_slip, i_vert_visc, w_split, w_exp_max, SPP,&
                    Fer_GM, K_GM, scaling_Ferreira, scaling_Rossby, scaling_resolution, scaling_FESOM14, Redi, visc_sh_limit, mix_scheme, Ricr, concv
 NAMELIST /oce_tra/ diff_sh_limit, Kv0_const, double_diffusion, K_ver, K_hor, surf_relax_T, surf_relax_S, balance_salt_water, clim_relax, &
		    ref_sss_local, ref_sss, i_vert_diff, tracer_adv, num_tracers, tracer_ID
END MODULE o_PARAM  
!==========================================================

!==========================================================
MODULE o_MESH
USE o_PARAM
USE, intrinsic :: ISO_FORTRAN_ENV
! All variables used to keep the mesh structure +
! auxiliary variables involved in implementation 
! of open boundaries and advection schemes
! 

integer, parameter                         :: MAX_ADJACENT=32 ! Max allowed number of adjacent nodes
integer                                    ::   nod2D      ! the number of 2D nodes
real(kind=WP)                              ::   ocean_area
real(kind=WP), allocatable, dimension(:,:) ::   coord_nod2D, geo_coord_nod2D
integer                                    ::   edge2D     ! the number of 2D edges
integer                                    ::   edge2D_in  
                                              ! the number of internal 2D edges
integer                                    ::   elem2D     ! the number of 2D elements
integer, allocatable, dimension(:,:)       ::   elem2D_nodes
                                              ! elem2D_nodes(:,n) lists
				              ! 3 nodes of element n   
integer, allocatable, dimension(:,:)       ::   edges
                                              ! edge(:,n) lists 2 nodes
				              ! edge n
integer, allocatable, dimension(:,:)       ::   edge_tri
                                              ! edge_tri(:,n) lists 2 
				              ! elements containing edge n
				              ! The first one is to left 
				              ! of the line directed
				              ! to the second node
integer, allocatable, dimension(:,:)       ::   elem_edges
                                              ! elem_edges(:,n) are edges of 
                                              ! element n.  
real(kind=WP), allocatable, dimension(:)   ::   elem_area
real(kind=WP), allocatable, dimension(:,:) ::   edge_dxdy, edge_cross_dxdy
real(kind=WP), allocatable, dimension(:)   ::   elem_cos, metric_factor
integer,allocatable,dimension(:,:)         ::   elem_neighbors
integer,allocatable,dimension(:,:)         ::   nod_in_elem2D
real(kind=WP),allocatable,dimension(:,:)   ::   x_corners, y_corners ! cornes for the scalar points
integer,allocatable,dimension(:)           ::   nod_in_elem2D_num
real(kind=WP),allocatable,dimension(:)     ::   depth
                                              ! depth(n) is the depths at 
				              ! node n 
real(kind=WP),allocatable,dimension(:,:)    ::   gradient_vec 
                                              ! Coefficients of linear reconstruction
					      ! of velocities on elements
real(kind=WP),allocatable,dimension(:,:)    ::   gradient_sca
                                              ! Coefficients to compute
					      ! gradient of scalars on elements
INTEGER,       ALLOCATABLE, DIMENSION(:)    :: bc_index_nod2D(:)
! Vertical structure             
integer                                    :: nl
real(kind=WP), allocatable, dimension(:)    :: zbar, Z,elem_depth
integer, allocatable, dimension(:)         :: nlevels, nlevels_nod2D
real(kind=WP), allocatable, dimension(:,:)  :: area, area_inv
real(kind=WP), allocatable, dimension(:)   :: mesh_resolution


  type sparse_matrix 
     integer :: nza
     integer :: dim
     real(kind=WP), allocatable, dimension(:)      :: values
     integer(int32), allocatable,   dimension(:) :: colind
     integer(int32), allocatable,   dimension(:) :: rowptr
     integer(int32), allocatable,   dimension(:) :: colind_loc
     integer(int32), allocatable,   dimension(:) :: rowptr_loc
  end type sparse_matrix
! Elevation stiffness matrix
type(sparse_matrix)                           :: ssh_stiff

! Auxiliary arrays. They are not related to mesh structure, but are 
! kept here because they are just used for temporary storage in computations

! Open boundary:
integer                                       :: ob_num  ! number of OB fragments

TYPE ob_type
    integer      :: len
    integer, allocatable, dimension(:)       :: list
END TYPE ob_type

TYPE ob_rhs_type
    integer      :: len
    real(kind=WP), allocatable, dimension(:) :: list
END TYPE ob_rhs_type

type(ob_type), allocatable                    ::  ob_info(:)
type(ob_rhs_type), allocatable                ::  ob_2rhs(:)
!
! The fct part
integer                                       :: fct_iter=1
real(kind=WP),allocatable,dimension(:,:)      :: fct_LO, fct_HO           ! Low-order solution
real(kind=WP),allocatable,dimension(:,:)      :: fct_adf_h, fct_adf_h2    ! Antidif. horiz. contrib. from edges / backup for iterafive fct scheme
real(kind=WP),allocatable,dimension(:,:)      :: fct_adf_v, fct_adf_v2    ! Antidif. vert. fluxes from nodes    / backup for iterafive fct scheme

real(kind=WP),allocatable,dimension(:,:)      :: fct_ttf_max,fct_ttf_min
real(kind=WP),allocatable,dimension(:,:)      :: fct_plus,fct_minus
! Quadratic reconstruction part
integer,allocatable,dimension(:)              :: nlevels_nod2D_min, nn_num, nboundary_lay
real(kind=WP),allocatable,dimension(:,:,:)    :: quad_int_mat, quad_int_coef
integer,allocatable,dimension(:,:)            :: nn_pos
! MUSCL type reconstruction
integer,allocatable,dimension(:,:)            :: edge_up_dn_tri
real(kind=WP),allocatable,dimension(:,:,:)    :: edge_up_dn_grad

#if defined (__oasis)
  real(kind=WP), allocatable, dimension(:)      :: lump2d_south, lump2d_north  
  integer, allocatable, dimension(:)           :: ind_south, ind_north    
#endif  
end module o_MESH
!==========================================================

!==========================================================
MODULE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
! Arrays are described in subroutine array_setup  
real(kind=WP), allocatable    :: UV(:,:,:)
real(kind=WP), allocatable    :: UV_rhs(:,:,:), UV_rhsAB(:,:,:)
real(kind=WP), allocatable    :: eta_n(:), d_eta(:)
real(kind=WP), allocatable    :: ssh_rhs(:), Wvel(:,:), hpressure(:,:)
real(kind=WP), allocatable    :: Wvel_e(:,:), Wvel_i(:,:)
real(kind=WP), allocatable    :: CFL_z(:,:)
real(kind=WP), allocatable    :: stress_surf(:,:)
real(kind=WP), allocatable    :: T_rhs(:,:) 
real(kind=WP), allocatable    :: heat_flux(:), Tsurf(:) 
real(kind=WP), allocatable    :: heat_flux_old(:), Tsurf_old(:)  !PS
real(kind=WP), allocatable    :: S_rhs(:,:)
real(kind=WP), allocatable    :: tr_arr(:,:,:),tr_arr_old(:,:,:)
real(kind=WP), allocatable    :: del_ttf(:,:)

real(kind=WP), allocatable    :: water_flux(:), Ssurf(:)
real(kind=WP), allocatable    :: virtual_salt(:), relax_salt(:)
real(kind=WP), allocatable    :: water_flux_old(:), Ssurf_old(:) !PS
real(kind=WP), allocatable    :: Tclim(:,:), Sclim(:,:)
real(kind=WP), allocatable    :: Visc(:,:)
real(kind=WP), allocatable    :: Tsurf_t(:,:), Ssurf_t(:,:)
real(kind=WP), allocatable    :: tau_x_t(:,:), tau_y_t(:,:) 
real(kind=WP), allocatable    :: heat_flux_t(:,:), heat_rel_t(:,:), heat_rel(:) 
real(kind=WP), allocatable    :: coriolis(:), coriolis_node(:)
real(kind=WP), allocatable    :: relax2clim(:)
real(kind=WP), allocatable    :: MLD1(:), MLD2(:)
integer,       allocatable    :: MLD1_ind(:), MLD2_ind(:)

! Passive and age tracers
real(kind=WP), allocatable    :: tracer(:,:,:), tracer_rhs(:,:,:)   
!Tracer gradients&RHS      
real(kind=WP), allocatable :: ttrhs(:,:)
real(kind=WP), allocatable :: tr_xy(:,:,:)
real(kind=WP), allocatable :: tr_z(:,:)

! Auxiliary arrays for vector-invariant form of momentum advection
real(kind=WP), allocatable,dimension(:,:)   :: vorticity

!Viscosity and diff coefs
real(kind=WP), allocatable,dimension(:,:)   :: Av,Kv
real(kind=WP), allocatable,dimension(:,:,:) :: Kv_double
real(kind=WP), allocatable,dimension(:)     :: Kv0
!Velocities interpolated to nodes
real(kind=WP), allocatable,dimension(:,:,:)   :: Unode

! Auxiliary arrays to store Redi-GM fields
real(kind=WP), allocatable,dimension(:,:,:) :: neutral_slope
real(kind=WP), allocatable,dimension(:,:,:) :: slope_tapered
real(kind=WP), allocatable,dimension(:,:,:) :: sigma_xy
real(kind=WP), allocatable,dimension(:,:)   :: sw_beta, sw_alpha
!real(kind=WP), allocatable,dimension(:,:,:) :: tsh, tsv, tsh_nodes
!real(kind=WP), allocatable,dimension(:,:)   :: hd_flux,vd_flux
!Isoneutral diffusivities (or xy diffusivities if Redi=.false)
real(kind=WP), allocatable :: Ki(:,:)

!_______________________________________________________________________________
! Arrays added for ALE implementation:
! --> layer thinkness at node and depthlayer for t=n and t=n+1
real(kind=WP), allocatable,dimension(:,:)   :: hnode, hnode_new, zbar_3d_n, Z_3d_n

! --> layer thinkness at elements, interpolated from hnode
real(kind=WP), allocatable,dimension(:,:)   :: helem

! --> thinkness of bottom elem (important for partial cells)
real(kind=WP), allocatable,dimension(:)     :: bottom_elem_thickness 
real(kind=WP), allocatable,dimension(:)     :: bottom_node_thickness 

! --> The increment of total fluid depth on elements. It is used to update the matrix
real(kind=WP), allocatable,dimension(:)     :: dhe

! --> hbar, hbar_old: correspond to the elevation, but on semi-integer time steps.
real(kind=WP), allocatable,dimension(:)     :: hbar, hbar_old

! --> auxiliary array to store an intermediate part of the rhs computations.
real(kind=WP), allocatable,dimension(:)     :: ssh_rhs_old !, ssh_rhs_old2 !PS

! --> auxiliary array to store depth of layers and depth of mid level due to changing 
!     layer thinkness at every node
real(kind=WP), allocatable,dimension(:)     :: zbar_n, Z_n

! new bottom depth at node and element due to partial cells
real(kind=WP), allocatable,dimension(:)     :: zbar_n_bot
real(kind=WP), allocatable,dimension(:)     :: zbar_e_bot

! --> multiplication factor for surface boundary condition in 
!     diff_ver_part_impl_ale(tr_num) between linfs -->=0.0 and noninfs 
!     (zlevel,zstar...) --> = 1.0
real(kind=WP)                               :: is_nonlinfs

!_______________________________________________________________________________
! Arrays added for pressure gradient force calculation
real(kind=WP), allocatable,dimension(:,:)   :: density_m_rho0
real(kind=WP), allocatable,dimension(:,:)   :: pgf_x, pgf_y
!_______________________________________________________________________________
!Monin-Obukhov correction
real(kind=WP),allocatable :: mo(:,:),mixlength(:)
!GM_stuff
real(kind=WP),allocatable :: bvfreq(:,:),mixlay_dep(:),bv_ref(:)

real(kind=WP),         allocatable    :: fer_UV(:,:,:), fer_wvel(:,:)
real(kind=WP), target, allocatable    :: fer_c(:), fer_K(:,:), fer_gamma(:,:,:)

real(kind=WP),         allocatable    :: ice_rejected_salt(:)
END MODULE o_ARRAYS
!==========================================================


!==========================================================
MODULE bgc
! Parameters, variables and functions for transient tracer simulations.

  implicit none
  save

! Normalized and fractionation-corrected atmospheric 14CO2 / 12CO2 ratios
  real(kind=8) :: r14c_a  = 1.0, & ! Global average and/or value in air-sea flux calculation
                  r14c_nh = 1.0, & ! Northern Hemisphere
                  r14c_tz = 1.0, & ! Tropics
                  r14c_sh = 1.0    ! Southern Hemisphere
! Normalized and fractionation-corrected atmospheric 39Ar/40Ar ratio
  real(kind=8) :: r39ar_a  = 1.0   ! Global average and/or value in air-sea flux calculation
! Atmospheric CO2 concentration
! CMIP6 & OMIP-BGC: xCO2_a = 284.32 ppm for 1700-1850 CE
! PMIP4:            xCO2_a = 190.00 ppm for 21 kcal BP
  real(kind=8) :: xCO2_a = 284.23e-6  ! mole fraction in dry air
! Atmospheric CFC-12 concentration (mole fraction in dry air)
  real(kind=8) :: xf12_a  = 0.00e-12, &  ! value passed in air-sea flux calculation
                  xf12_nh = 0.00e-12, &  ! Northern Hemisphere
                  xf12_sh = 0.00e-12     ! Southern Hemisphere
! Atmospheric SF6 concentration (mole fraction in dry air)
  real(kind=8) :: xsf6_a  = 0.00e-12, &  ! value passed in air-sea flux calculation
                  xsf6_nh = 0.00e-12, &  ! Northern Hemisphere
                  xsf6_sh = 0.00e-12     ! Southern Hemisphere
! Global-mean DIC concentration in the mixed layer (mol / m**3)
  real(kind=8) :: dic_0 = 2.00        ! GLODAPv2, 0-50 m: TCO2 ~ 2050 umol / kg
! Radioactive decay constants (1 / s; default values assume that 1 year = 365.00 days)
  real(kind=8) :: decay14 = 3.8561e-12 , & ! 14C; t1/2 = 5700 a following OMIP-BGC
                  decay39 = 8.1708e-11     ! 39Ar; t1/2 = 269 a
! Latitude of/for atmospheric boundary conditions
  real(kind=8) :: y_abc
! Switches for off-line simulations
  logical ::  offline = .false., online = .true. ! on-line simulations (default setup)
! logical :: offline = .true., online = .true.  ! diagnose dynamic fields to be used in off-line simulations
! logical :: offline = .true., online = .false. ! enable off-line simulations
  
! Namelist to modify default parameter settings
  namelist / bgc_param / r14c_a, r14c_nh, r14c_tz, r14c_sh, xco2_a, &  ! atmospheric values for Delta14C
                         xf12_a, xf12_nh, xf12_sh, &                   ! atmospheric values for CFC-12
                         xsf6_a, xsf6_nh, xsf6_sh, &                   ! atmospheric values for SF6
                         dic_0, decay14, &                             ! both are needed for Delta14C
                         offline, online


  contains


    function iso_flux(temp_c, sal, u_10, v_10, f_ice, p_air, x_gas, iso_a, iso_w, c_surf)
!     Calculate isotopic air-sea exchange fluxes in m / s, positive values mean oceanic uptake.
!     This function should become obsolete as it could be replaced with function gas_flux().
      implicit none

      real(kind=8) :: iso_flux
!     Input parameters
      real(kind=8), intent(in) :: temp_c, sal, u_10, v_10, f_ice, p_air, x_gas, iso_a, iso_w, c_surf
!     temp_c = temperature (SST) in deg C
!     sal = salinity in PSU or permil
!     u_10, v_10 = zonal and meridional wind speed (m / s) at 10 m height
!     f_ice = sea-ice fractional coverage
!     p_air = total atmospheric pressure (Pa)
!     x_gas = mole fraction of atmospheric trace gas (for F14C: 12CO2)
!     iso_a, iso_w = normalized isotopic ratios (14C/12C or 39/40Ar) in atmosphere and surface water
!     c_surf = tracer concentration (e.g. of DIC) in surface water

!     Computation of 14CO2 fluxes following Wanninkhof (2014, equation 5)
!     assuming local air-sea flux equilibrium for CO2 and using SI units:
!     0.251 (cm / h) / (m / s)**2 by Wanninkhof (2014) -> 6.9722e-7 s / m
!!      flux_r14co2 = 6.9722e-7 * solub("co2", temp_c, sal) * sc_660("co2", temp_c)**(-0.5) * & 
!!                    (u_10**2 + v_10**2) * partial_press(xco2_a, p_air, temp_c, sal) * & 
!!                    (r14c_a - r14c_w) * (1. - f_ice) / dic_0
      iso_flux = transfer_vel("co2", temp_c, u_10, v_10) * &
                 solub("co2", temp_c, sal) * partial_press(x_gas, p_air, temp_c, sal) * &
                 (iso_a - iso_w) * (1. - f_ice) / c_surf
      return
    end function iso_flux

    function gas_flux(which_gas, temp_c, sal, u_10, v_10, f_ice, p_air, x_gas, c_surf)
!     Computes air-sea exchange gas fluxes in m / s, positive values mean oceanic uptake.
      implicit none

      real(kind=8) :: gas_flux
!     Input parameters
      character(len=3), intent(in) :: which_gas  ! trace gas name
      real(kind=8), intent(in) :: temp_c, sal, & ! SST (deg C) and SSS ("PSU" or permil)
                                  u_10, v_10,  & ! wind speed at 10 m height (m / s)
                                  f_ice, &       ! sea-ice fractional coverage
                                  p_air, &       ! total atmospheric pressure (Pa)
                                  x_gas, &       ! atmospheric mole fraction 
                                  c_surf         ! marine surface water concentration (mol / m**3)
!     Internal variables
      real(kind=8) :: c_sat                      ! marine saturation concentration (mol / m**3)
      c_sat = solub(which_gas, temp_c, sal) * partial_press(x_gas, p_air, temp_c, sal) * x_gas
      gas_flux = transfer_vel(which_gas, temp_c, u_10, v_10) * (c_sat - c_surf) * (1. - f_ice)

      return
    end function gas_flux


    function partial_press(x_gas, p_air, temp_c, sal)
!     Converts the mole fraction of a gas to its partial pressure in marine air.
      implicit none

      real(kind=8) :: partial_press
!     Input parameters
      real(kind=8), intent(in) :: x_gas, &       ! mole fraction of the gas to be converted
                                  p_air, &       ! total atmospheric pressure in Pa
                                  temp_c, sal    ! temperature (deg C) and salinity ("PSU" or permil)
!     Internal variables
      real(kind=8) :: p_h2o, &                   ! vapor pressure of sea water in atm
                      temp_k100                  ! water temperature in K / 100

      temp_k100 = (temp_c + 273.15) * 0.01

!     Vapor pressure over seawater in atm (Weiss & Price 1980, eq. 10)
      p_h2o = exp(24.4543 - 67.4509 / temp_k100 - 4.8489 * log(temp_k100) - 0.000544 * sal)
!     Partial pressure in moist air in atm, 1 atm = 101325 Pa
      partial_press = x_gas  * (1. - p_h2o) * p_air / 1.01325e5 

      return 
    end function partial_press


    function solub(which_gas, temp_c, sal)
!     Computes the solubility of trace gases in seawater.
      implicit none
      real(kind=8) :: solub                      ! solubility (mol / (m**3 * atm))
!     Input parameters
      character(len=3), intent(in) :: which_gas  ! tracer name
      real(kind=8), intent(in) :: temp_c, &      ! temperature (deg C)
                                  sal            ! salinity ("PSU" or permil)
      real(kind=8) :: a1, a2, a3, b1, b2, b3, &  ! polynomial coefficients of the solubility function 
                      temp_k100, &               ! water temperature in K / 100
                      per_m3 = 1000.             ! factor to convert from 1 / L to 1 / m**3

      temp_k100 = (temp_c + 273.15) * 0.01

      select case (which_gas)
      case ("co2")
!       CO2 solubility according to Weiss, 1974 (eq. 12 and tab. I). 
!       We follow the OMPIC-BCG recommendations for "abiotic" carbon tracers 
!       (Orr et al. 2017, eq. 17-18 and tab. 3) to use this solubility function
!       instead of the one by Weiss & Price (1985). 
        a1 = -58.0931 ; a2 =  90.5069;  a3 = 22.2940
        b1 =  0.027766; b2 = -0.025888; b3 = 0.0050578
!!        solub = exp(-58.0931 + 90.5069 / temp_k100 + 22.2940 * log(temp_k100) + & 
!!                    sal * (0.027766 - 0.025888 * temp_k100 + 0.0050578 * temp_k100**2)) * 1000
      case ("f12") 
!       CFC-12 !! CHECK, ADD REF
        a1 = -122.3246; a2 = 182.5306; a3 = 50.5898
        b1 = -0.145633; b2 = 0.092509; b3 = -0.0156627 
        per_m3 = per_m3 * 1e12 ! -> solubility in pmol / m**3
      case ("sf6") 
!       SF6 !! CHECK, ADD REF
        a1 = -96.5975; a2 = 139.883; a3 = 37.8193
        b1 =  0.0310693; b2 = -0.0356385; b3 = 0.00743254
        per_m3 = per_m3 * 1e12 ! -> solubility in pmol / m**3
      case("arg")
!       Ar-39, Bunsen coefficient -> convert to K' <- UNDER CONSTRUCTION
        a1 = -55.6578; a2 = 82.0262; a3 = 22.5929
        b1 = -0.036267; b2 = 0.016241; b3 = -0.0020114
!       per_m3 =  
      end select

      solub = exp(       a1 + a2 / temp_k100 + a3 * log(temp_k100) + & 
                  sal * (b1 - b2 * temp_k100 + b3 * temp_k100**2))  ! in mol / L
      solub = solub * per_m3


      return
    end function solub


    function sc_660(which_gas, temp_c)
!     Schmidt numbers of trace gases in sea water with S = 35 
!     normalized to 20 degC (Sc(CO2) ~660; Wanninkhof 2014, tab. 1)).
      implicit none

!     Result
      real(kind=8) :: sc_660                       ! Schmidt number
!     Input parameters
      character(len=3), intent(in) :: which_gas    ! tracer name
      real(kind=8), intent(in) :: temp_c           ! temperature (deg C)
!     Internal parameters and/or variables
      real(kind=8) :: as, bs, cs, ds, es           ! polynomial coefficients

      select case (which_gas)
      case ("co2") ! CO2
        as = 2116.8; bs = -136.25; cs = 4.7353; ds = -0.092307; es = 0.0007555
!!        sc_660 = (2116.8 - 136.25 *temp_c + 4.7353 * temp_c**2 - 0.092307 * temp_c**3 + 0.0007555 * temp_c**4) / 660
      case ("f12") ! CFC-12
        as = 3828.1; bs = -249.86; cs = 8.7603; ds = -0.171600; es = 0.0014080
      case ("sf6") ! SF6
        as = 3177.5; bs = -200.57; cs = 6.8865; ds = -0.133350; es = 0.0010877
      case ("arg") ! Ar-39
        as = 2078.1; bs = -146.74; cs = 5.6403; ds = -0.118380; es = 0.0010148
      end select
      
      sc_660 = (as + bs *temp_c + cs * temp_c**2 + ds * temp_c**3 + es * temp_c**4) / 660.
      
      return
    end function sc_660


    function transfer_vel(which_gas, temp_c, u_10, v_10)
!     Compute gas transfer velocities of / for tracers
!     Result
      real(kind=8) :: transfer_vel                 ! transfer velocity (m / s)
!     Input parameters
      character(len=3), intent(in) :: which_gas    ! tracer name
      real(kind=8), intent(in) :: temp_c, &        ! temperature (deg C)
                                  u_10, v_10       ! wind speed at 10 m height (m / s)

!     Wanninkhof (2014), eq. (4) with a = 0.251 (cm / h) / (m / s)**2 -> 6.9722e-7 s / m 
!     to obtain the gas transfer velocity in m / s
      transfer_vel = 6.9722e-7 * sc_660(which_gas, temp_c)**(-0.5) * (u_10**2 + v_10**2) 

      return
    end function transfer_vel

END MODULE bgc
