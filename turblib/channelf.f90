    !	Copyright 2011 Johns Hopkins University
    !
    !  Licensed under the Apache License, Version 2.0 (the "License");
    !  you may not use this file except in compliance with the License.
    !  You may obtain a copy of the License at
    !
    !      http://www.apache.org/licenses/LICENSE-2.0
    !
    !  Unless required by applicable law or agreed to in writing, software
    !  distributed under the License is distributed on an "AS IS" BASIS,
    !  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    !  See the License for the specific language governing permissions and
    !  limitations under the License.


    program TurbTest

    use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_float, c_f_pointer
    implicit none

    type, bind(c) :: thresholdinfo
        integer(c_int) :: x, y, z
        real(c_float) :: value
    end type thresholdinfo

    integer, parameter :: RP=4 ! Number of bytes for reals (single precision)

    ! ---- Temporal Interpolation Options ----
    integer, parameter :: NoTInt = 0   ! No temporal interpolation
    integer, parameter :: PCHIPInt = 1 ! Piecewise cubic Hermit interpolation in time

    ! ---- Spatial Interpolation Flags for GetVelocity & GetVelocityAndPressure ----
    integer, parameter :: NoSInt = 0 ! No spatial interpolation
    integer, parameter :: Lag4 = 4   ! 4th order Lagrangian interpolation in space
    integer, parameter :: Lag6 = 6   ! 6th order Lagrangian interpolation in space
    integer, parameter :: Lag8 = 8   ! 8th order Lagrangian interpolation in space

    ! ---- Spatial Differentiation & Interpolation Flags for GetVelocityGradient & GetPressureGradient ----
    integer, parameter :: FD4NoInt = 40 ! 4th order finite differential scheme for grid values, no spatial interpolation
    integer, parameter :: FD6NoInt = 60 ! 6th order finite differential scheme for grid values, no spatial interpolation
    integer, parameter :: FD8NoInt = 80 ! 8th order finite differential scheme for grid values, no spatial interpolation
    integer, parameter :: FD4Lag4 = 44  ! 4th order finite differential scheme for grid values, 4th order Lagrangian interpolation in space

    ! ---- Spline interpolation and differentiation Flags for getVelocity,
    !      getPressure, getVelocityGradient, getPressureGradient,
    !      getVelocityHessian, getPressureHessian
    integer, parameter :: M1Q4 = 104   ! Splines with smoothness 1 (3rd order) over 4 data points. Not applicable for Hessian.
    integer, parameter :: M2Q8 = 208   ! Splines with smoothness 2 (5th order) over 8 data points.
    integer, parameter :: M2Q14 = 214  ! Splines with smoothness 2 (5th order) over 14 data points.

    !
    ! Choose which dataset to use in this query
    ! Currently, only valid datasets are:
    !   'isotropic1024coarse', 'isotropic1024fine', 'mhd1024', 'channel', 'mixing' and 'isotropic4096'
    !
    character(*), parameter :: dataset = 'channel' // CHAR(0)

    !
    ! Specify your access key here.
    ! If you need one, please visit http://turbulence.pha.jhu.edu/
    ! (We just want to know a bit about our users!)
    !
    character(*), parameter :: authkey = 'edu.jhu.pha.turbulence.testing-201406' // CHAR(0)

    real(RP), parameter :: time = 0.364_RP;

    real(RP) :: points(3, 10)    ! input
    real(RP) :: dataout1(10)     ! p
    real(RP) :: dataout3(3, 10)  ! x,y,z
    real(RP) :: dataout4(4, 10)  ! x,y,z,p
    real(RP) :: dataout6(6, 10)  ! results from Pressure Hessian
    real(RP) :: dataout9(9, 10)  ! results from Velocity Gradient
    real(RP) :: dataout18(18, 10) ! results from Velocity Hessian

    integer,parameter :: x_start=1, y_start=1, z_start=1, x_end=4, y_end=4, z_end=4

    real(RP) :: threshold = 0.5_RP; ! threshold level
    character(*), parameter :: threshold_field = 'vorticity' // CHAR(0) ! field used for threshold query
    integer(c_int) :: result_size; ! the size of the array returned from the getThreshold function
    type(c_ptr) :: cptr_to_array   ! C pointer to the dynamically allocated array from the getThreshold function
    type(thresholdinfo), pointer :: points_above_threshold(:) => NULL() ! Fortran pointer used to extract the array
    ! of points from the C pointer above

    ! Declare the return type of the turblib functions as integer.
    ! This is required for custom error handling (see the README).
    integer :: getvelocity, getvelocityandpressure, getvelocitygradient
    integer :: getvelocitylaplacian, getvelocityhessian
    integer :: getpressure, getpressuregradient, getpressurehessian
    integer :: getthreshold
    ! return code
    integer :: rc

    ! loop iterator
    integer :: i

    ! Formatting rules
    character(*), parameter :: format1='(i3,1(a,e13.6))'
    character(*), parameter :: rawformat1='(i4,1(a,e13.6))'
    character(*), parameter :: format3='(i3,3(a,e13.6))'
    character(*), parameter :: rawformat3='(i4,3(a,e13.6))'
    character(*), parameter :: format4='(i3,4(a,e13.6))'
    character(*), parameter :: format6='(i3,6(a,e13.6))'
    character(*), parameter :: format9='(i3,9(a,e13.6))'
    character(*), parameter :: format18='(i3,18(a,e13.6))'
    character(*), parameter :: formatT='(3(a,i4),1(a,e13.6))'
    !
    ! Intialize the gSOAP runtime.
    ! This is required before any WebService routines are called.
    !
    CALL soapinit()

    ! Enable exit on error.  See README for details.
    CALL turblibSetExitOnError(1)

    do i = 1, 10
        points(1, i) = 0.20 * i
        points(2, i) = 0.09 * i
        points(3, i) = 0.15 * i
    end do

    write(*,*)
    write(*,'(a)') "Coordinates of 10 points where variables are requested:"
    do i = 1, 10
        write(*,format3) i, ': ', points(1,i), ', ', points(2,i), ', ', points(3,i)
    end do


    write(*,*)
    write(*,'(a)') 'Requesting velocity at 10 points...'
    rc = getvelocity(authkey, dataset, time, Lag6, NoTInt, 10, points, dataout3)
    do i = 1, 10
        write(*,format3) i, ': ', dataout3(1,i), ', ', dataout3(2,i), ', ', dataout3(3,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting velocity and pressure at 10 points...'
    rc = getvelocityandpressure(authkey, dataset, time, Lag6, NoTInt, 10, points, dataout4)
    do i = 1, 10
        write(*,format4) i, ': ', dataout4(1,i), ', ', dataout4(2,i), ', ', dataout4(3,i), ', ', dataout4(4,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting velocity gradient at 10 points...'
    rc = getvelocitygradient(authkey, dataset, time, FD4Lag4, NoTInt, 10, points, dataout9)
    do i = 1, 10
        write(*,format9) i, ': duxdx=', dataout9(1,i), ', duxdy=', dataout9(2,i), &
            ', duxdz=', dataout9(3,i), ', duydx=', dataout9(4,i),  &
            ', duydy=', dataout9(5,i), ', duydz=', dataout9(6,i),  &
            ', duzdx=', dataout9(7,i), ', duzdy=', dataout9(8,i),  &
            ', duzdz=', dataout9(9,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting velocity laplacian at 10 points...'
    rc = getvelocitylaplacian(authkey, dataset, time, FD4Lag4, NoTInt, 10, points, dataout3)
    do i = 1, 10
        write(*,format3) i, ': grad2ux=', dataout3(1,i), ', grad2uy=', dataout3(2,i), ', grad2uz=', dataout3(3,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting velocity hessian at 10 points...'
    rc = getvelocityhessian(authkey, dataset, time, FD4Lag4, NoTInt, 10, points, dataout18)
    do i = 1, 10
        write(*,format18) i, ': d2uxdxdx=', dataout18(1,i), &
            ', d2uxdxdy=', dataout18(2,i), &
            ', d2uxdxdz=', dataout18(3,i), &
            ', d2uxdydy=', dataout18(4,i), &
            ', d2uxdydz=', dataout18(5,i), &
            ', d2uxdzdz=', dataout18(6,i), &
            ', d2uydxdx=', dataout18(7,i), &
            ', d2uydxdy=', dataout18(8,i), &
            ', d2uydxdz=', dataout18(9,i), &
            ', d2uydydy=', dataout18(10,i), &
            ', d2uydydz=', dataout18(11,i), &
            ', d2uydzdz=', dataout18(12,i), &
            ', d2uzdxdx=', dataout18(13,i), &
            ', d2uzdxdy=', dataout18(14,i), &
            ', d2uzdxdz=', dataout18(15,i), &
            ', d2uzdydy=', dataout18(16,i), &
            ', d2uzdydz=', dataout18(18,i), &
            ', d2uzdzdz=', dataout18(18,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting pressure at 10 points...'
    rc = getpressure(authkey, dataset, time, Lag6, NoTInt, 10, points, dataout1)
    do i = 1, 10
        write(*,format1) i, ': ', dataout1(i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting pressure gradient at 10 points...'
    rc = getpressuregradient(authkey, dataset, time, FD4Lag4, NoTInt, 10, points, dataout3)
    do i = 1, 10
        write(*,format3) i, ': dpdx=', dataout3(1,i), ', dpdy=', dataout3(2,i), ', dpdz=', dataout3(3,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting pressure hessian at 10 points...'
    rc = getpressurehessian(authkey, dataset, time, FD4Lag4, NoTInt, 10, points, dataout6)
    do i = 1, 10
        write(*,format6) i, ': d2pdxdx=', dataout6(1,i), ', d2pdxdy=', dataout6(2,i), &
            ', d2pdxdz=', dataout6(3,i), ', d2pdydy=', dataout6(4,i),  &
            ', d2pdydz=', dataout6(5,i), ', d2pdzdz', dataout6(6,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting points above threshold ...'
    rc = getthreshold(authkey, dataset, threshold_field, time, threshold, FD4NoInt, x_start, y_start, z_start, &
        x_end, y_end, z_end, cptr_to_array, result_size)
    ! convert the C pointer to Fortran pointer
    call c_f_pointer(cptr_to_array, points_above_threshold, [result_size])
    do i = 1, result_size
        write(*,formatT) '(', points_above_threshold(i)%x, ', ', points_above_threshold(i)%y, ', ', &
            points_above_threshold(i)%z, ') value = ', points_above_threshold(i)%value
    end do
    ! The dynamically allocated array for the getThreshold function needs to be deallocated.
    ! Call this function to deallocate the array, otherwise a memory leak can occur.
    call deallocate_array(cptr_to_array)

    !
    ! Destroy the gSOAP runtime.
    ! No more WebService routines may be called.
    !
    CALL soapdestroy()

    end program TurbTest

