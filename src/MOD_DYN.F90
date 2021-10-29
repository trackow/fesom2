!==========================================================
MODULE MOD_DYN
USE O_PARAM
USE, intrinsic :: ISO_FORTRAN_ENV
USE MOD_WRITE_BINARY_ARRAYS
USE MOD_READ_BINARY_ARRAYS
IMPLICIT NONE
SAVE

!
!
!_______________________________________________________________________________
! set data array strucutre for instant zonal, meridionasl and vertical velocities + 
! Adams-Bashfort rhs
TYPE T_DYN_DATA
    ! instant zonal merdional velocity & Adams-Bashfort rhs
    real(kind=WP), allocatable, dimension(:,:,:):: uv, uv_rhsAB  

    ! instant vertical velm explicite+implicite part
    real(kind=WP), allocatable, dimension(:,:)  :: w, w_e, w_i   

    !___________________________________________________________________________
    contains
        procedure WRITE_T_DYN_DATA
        procedure READ_T_DYN_DATA
        generic :: write(unformatted) => WRITE_T_DYN_DATA
        generic :: read(unformatted)  => READ_T_DYN_DATA
END TYPE T_DYN_DATA

!
!
!_______________________________________________________________________________
! set working arrays structure 
TYPE T_DYN_WORK
    real(kind=WP), allocatable, dimension(:,:,:):: uv_rhs, uvnode, uvnode_rhs
    
    !___________________________________________________________________________
    contains
        procedure WRITE_T_DYN_WORK
        procedure READ_T_DYN_WORK
        generic :: write(unformatted) => WRITE_T_DYN_WORK
        generic :: read(unformatted)  => READ_T_DYN_WORK
END TYPE T_DYN_WORK

!
!
!_______________________________________________________________________________
! set main structure for dynamicss, contains viscosity options and parameters + 
! option for momentum advection 
TYPE T_DYN
    
    ! contains: %uv, %uv_rhsAB, %w, %w_e, %w_i 
    type(t_dyn_data)                         :: data
    
    ! contains: %uv_rhs, %uvnode, %uvnode_rhsAB
    type(t_dyn_work)                         :: work
    
    ! visc_option=...
    ! 1=Harmonic Leith parameterization;
    ! 2=Laplacian+Leith+biharmonic background
    ! 3=Biharmonic Leith parameterization
    ! 4=Biharmonic flow aware
    ! 5=Kinematic (easy) Backscatter
    ! 6=Biharmonic flow aware (viscosity depends on velocity Laplacian)
    ! 7=Biharmonic flow aware (viscosity depends on velocity differences)
    ! 8=Dynamic Backscatter
    integer                                     :: visc_opt=5      

    ! gamma0 [m/s],   backgroung viscosity= gamma0*len, it should be as small 
    !                 as possible (keep it < 0.01 m/s).
    ! gamma1 [nodim], for computation of the flow aware viscosity
    ! gamma2 [s/m],   is only used in easy backscatter option
    real(kind=WP)                               :: gamma0_visc  = 0.03
    real(kind=WP)                               :: gamma1_visc  = 0.1
    real(kind=WP)                               :: gamma2_visc  = 0.285

    ! div_c the strength of the modified Leith viscosity, nondimensional, 0.3 -- 1.0
    ! leith the strength of the Leith viscosity
    real(kind=WP)                               :: div_c_visc   = 0.5
    real(kind=WP)                               :: leith_c_visc = 0.05

    logical                                     :: use_ivertvisc = .true.
    integer                                     :: momadv_opt   = 2
    
    ! Switch on free slip
    logical                                     :: use_freeslip  = .false. 
    
    ! do implicite, explicite spliting of vertical velocity
    logical                                     :: use_wsplit    = .false.
    ! maximum allowed CFL criteria in vertical (0.5 < w_max_cfl < 1.) 
    ! in older FESOM it used to be w_exp_max=1.e-3
    real(kind=WP)                               :: wsplit_maxcfl= 1.0     

    !___________________________________________________________________________
    contains
        procedure WRITE_T_DYN
        procedure READ_T_DYN
        generic :: write(unformatted) => WRITE_T_DYN
        generic :: read(unformatted)  => READ_T_DYN
END TYPE T_DYN

contains
!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN_DATA
subroutine WRITE_T_DYN_DATA(dyndata, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN_DATA),    intent(in)     :: dyndata
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    call write_bin_array(dyndata%uv      , unit, iostat, iomsg)
    call write_bin_array(dyndata%uv_rhsAB, unit, iostat, iomsg)
    call write_bin_array(dyndata%w       , unit, iostat, iomsg)
    call write_bin_array(dyndata%w_e     , unit, iostat, iomsg)
    call write_bin_array(dyndata%w_i     , unit, iostat, iomsg)
end subroutine WRITE_T_DYN_DATA

subroutine READ_T_DYN_DATA(dyndata, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN_DATA),    intent(inout)  :: dyndata
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    call read_bin_array(dyndata%uv      , unit, iostat, iomsg)
    call read_bin_array(dyndata%uv_rhsAB, unit, iostat, iomsg)
    call read_bin_array(dyndata%w       , unit, iostat, iomsg)
    call read_bin_array(dyndata%w_e     , unit, iostat, iomsg)
    call read_bin_array(dyndata%w_i     , unit, iostat, iomsg)
end subroutine READ_T_DYN_DATA

!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN_WORK
subroutine WRITE_T_DYN_WORK(dynwork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN_WORK),    intent(in)     :: dynwork
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    call write_bin_array(dynwork%uv_rhs    , unit, iostat, iomsg)
    call write_bin_array(dynwork%uvnode    , unit, iostat, iomsg)
    call write_bin_array(dynwork%uvnode_rhs, unit, iostat, iomsg)
end subroutine WRITE_T_DYN_WORK

subroutine READ_T_DYN_WORK(dynwork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN_WORK),    intent(inout)  :: dynwork
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    call read_bin_array(dynwork%uv_rhs    , unit, iostat, iomsg)
    call read_bin_array(dynwork%uvnode    , unit, iostat, iomsg)
    call read_bin_array(dynwork%uvnode_rhs, unit, iostat, iomsg)
end subroutine READ_T_DYN_WORK

!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN
subroutine WRITE_T_DYN(dynamics, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN),         intent(in)     :: dynamics
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%data
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%work
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_opt
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma0_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma1_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma2_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%div_c_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%leith_c_visc
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%use_ivertvisc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%momadv_opt
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%use_freeslip
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%use_wsplit
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%wsplit_maxcfl
end subroutine WRITE_T_DYN

subroutine READ_T_DYN(dynamics, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN),         intent(inout)  :: dynamics
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%data
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%work
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_opt
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma0_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma1_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma2_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%div_c_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%leith_c_visc
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_ivertvisc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%momadv_opt
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_freeslip
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_wsplit
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%wsplit_maxcfl
end subroutine READ_T_DYN

END MODULE MOD_DYN