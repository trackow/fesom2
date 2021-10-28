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
    real(kind=WP), allocatable, dimension(:,:,:):: w, w_e, w_i   

    !___________________________________________________________________________
    contains
        procedure WRITE_T_DYN_DATA
        procedure READ_T_DYN_DATA
        generic :: write(unformatted) => WRITE_T_DYN_DATA
        generic :: read(unformatted)  => READ_T_DYN_DATA
END TYPE T_TRACER_DATA

!
!
!_______________________________________________________________________________
! set working arrays structure 
TYPE T_DYN_WORK
    real(kind=WP), allocatable, dimension(:,:,:):: uv_rhs, uvnode, uvnode_rhsAB
    
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
! set main structure for dynamics, contains viscosity options and parameters + 
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
    integer                                     :: visc_opt      

    ! gamma0 [m/s],   backgroung viscosity= gamma0*len, it should be as small 
    !                 as possible (keep it < 0.01 m/s).
    ! gamma1 [nodim], for computation of the flow aware viscosity
    ! gamma2 [s/m],   is only used in easy backscatter option
    real(kind=WP)                               :: gamma0_visc, gamma1_visc, gamma2_visc

    ! div_c the strength of the modified Leith viscosity, nondimensional, 0.3 -- 1.0
    ! leith the strength of the Leith viscosity
    real(kind=WP)                               :: div_c_visc, leith_c_visc

    logical                                     :: i_vert_visc =.true.
    integer                                     :: mom_adv=2
    
    !___________________________________________________________________________
    contains
        procedure WRITE_T_DYN
        procedure READ_T_DYN
        generic :: write(unformatted) => WRITE_T_DYN
        generic :: read(unformatted)  => READ_T_DYN
END TYPE T_DYN

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
    call write_bin_array(dynwork%uv_rhs      , unit, iostat, iomsg)
    call write_bin_array(dynwork%uvnode      , unit, iostat, iomsg)
    call write_bin_array(dynwork%uvnode_rhsAB, unit, iostat, iomsg)
end subroutine WRITE_T_DYN_WORK

subroutine READ_T_DYN_WORK(dynwork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN_WORK),    intent(inout)  :: dynwork
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    call read_bin_array(dynwork%uv_rhs      , unit, iostat, iomsg)
    call read_bin_array(dynwork%uvnode      , unit, iostat, iomsg)
    call read_bin_array(dynwork%uvnode_rhsAB, unit, iostat, iomsg)
end subroutine READ_T_DYN_WORK

!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN
subroutine WRITE_T_DYN(dynamic, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN),         intent(in)     :: dynamic
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%data
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%work
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%visc_opt
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%gamma0_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%gamma1_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%gamma2_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%div_c_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%leith_c_visc
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%i_vert_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamic%mom_adv
end subroutine WRITE_T_DYN

subroutine READ_T_DYN(dynamic, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN),         intent(in)     :: dynamic
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%data
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%work
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%visc_opt
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%gamma0_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%gamma1_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%gamma2_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%div_c_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%leith_c_visc
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%i_vert_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamic%mom_adv
end subroutine READ_T_DYN

END MODULE MOD_DYN