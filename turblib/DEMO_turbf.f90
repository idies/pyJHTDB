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
    character(*), parameter :: dataset = 'isotropic1024coarse' // CHAR(0)

    !
    ! Specify your access key here.
    ! If you need one, please visit http://turbulence.pha.jhu.edu/
    ! (We just want to know a bit about our users!)
    !
    character(*), parameter :: authkey = 'edu.jhu.pha.turbulence.testing-201406' // CHAR(0)

    character(*), parameter :: field = 'velocity' // CHAR(0) ! field used for box filter
    character(*), parameter :: scalar_fields = 'pp' // CHAR(0) ! two scalar fields ("p" and "p") used
    ! for getBoxFilterSGSscalar
    character(*), parameter :: vector_scalar_fields = 'up' // CHAR(0) ! a vector and a scalar field ("u" and "p")
    ! used for getBoxFilterSGSvector
    real(RP), parameter :: filterwidth = 0.055223308_RP ! 9 * dx, where dx = 2*PI/1024

    ! dx = 2*PI/4096 if using the isotropic4096 dataset.

    real(RP), parameter :: spacing = 0.030680_RP ! 5 * dx

    real(RP), parameter :: time = 0.364_RP;
    real(RP), parameter :: startTime = 0.364_RP;
    real(RP), parameter :: endTime = 0.376_RP;
    real(RP), parameter :: lag_dt = 0.0004_RP     ! fraction of database timestep to use for
    ! getposition method (Lagrangian integration time step)
    real(RP) :: points(3, 10)    ! input
    real(RP) :: dataout1(10)     ! p
    real(RP) :: dataout2(2, 10)  ! results from invariant
    real(RP) :: dataout3(3, 10)  ! x,y,z
    real(RP) :: dataout4(4, 10)  ! x,y,z,p
    real(RP) :: dataout6(6, 10)  ! results from Pressure Hessian
    real(RP) :: dataout9(9, 10)  ! results from Velocity Gradient
    real(RP) :: dataout18(18, 10) ! results from Velocity Hessian

    integer,parameter :: x_start=1, y_start=1, z_start=1, x_end=4, y_end=4, z_end=4

    real(RP) :: threshold = 30.0_RP; ! threshold level
    character(*), parameter :: threshold_field = 'vorticity' // CHAR(0) ! field used for threshold query
    integer(c_int) :: result_size; ! the size of the array returned from the getThreshold function
    type(c_ptr) :: cptr_to_array   ! C pointer to the dynamically allocated array from the getThreshold function
    type(thresholdinfo), pointer :: points_above_threshold(:) => NULL() ! Fortran pointer used to extract the array
    ! of points from the C pointer above

    ! Declare the return type of the turblib functions as integer.
    ! This is required for custom error handling (see the README).
    integer :: getvelocity, getforce, getvelocityandpressure, getvelocitygradient
    integer :: getboxfilter, getboxfiltergradient
    integer :: getboxfiltersgssymtensor, getboxfiltersgsscalar, getboxfiltersgsvector, getboxfiltersgstensor
    integer :: getvelocitylaplacian, getvelocityhessian
    integer :: getpressure, getpressuregradient, getpressurehessian
    integer :: getinvariant
    integer :: getposition
    integer :: getthreshold
    ! If working with hdf5 cutout files, uncomment the lines below
    ! to load the file and use the functions locally:
    ! character(*), parameter :: filename = 'isotropic1024coarse.h5' // CHAR(0)
    ! integer :: turblibaddlocalsource

    ! If selecting a dataset without time evolution (e.g. isotropic4096), certain getFunctions such as getPosition do not work.

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

    ! If working with cutout files, CUTOUT_SUPPORT should be defined during compilation.
    ! Make sure to run make with the "CUTOUT_SUPPORT=1" option, e.g.:
    ! $ make turbf CUTOUT_SUPPORT=1
    ! Change the filename to the name of the downloaded cutout file
    ! and uncomment the line below:
    ! rc = turblibaddlocalsource(filename);


    !   The client library implements all of the server-side functionality "locally" (except
    !   for particle tracking and filtering). Therefore, if an hdf5 file with cutout data is
    !   available and loaded as above all queries for data that are within the region defined
    !   in the file will be evaluated locally (without being sent to the server). An example
    !   is provided below.

    !   Please note that the use of this feature of the client library requires an hdf5
    !   installation. The standard approach of simply executing queries through the server
    !   does not require hdf5 or downloading any cutout data.

    !   For users that make frequent calls for data in a particular region and would like to
    !   download this data to their local machine the only change in their existing code
    !   will be to use the turblibAddLocalSource function to load the hdf5 file that they
    !   have downloaded. Each function call will determine whether the data is available
    !   locally and will evaluate the query locally if it is and will make a call to the
    !   Web-services on the server if it is not.

    !   The steps below can be followed to download data cutouts in hdf5 format and use the
    !   client library locally:
    !   1) Download an hdf5 file containing a cubic region of the velocity data at timestep 0
    !      with the following download link:
    !	http://turbulence.pha.jhu.edu/download.aspx/[authorization Token]/isotropic1024coarse/u/0,1/0,16/0,16/0,16/
    !   2) Compile this sample code with the "CUTOUT_SUPPORT=1" option, e.g.:
    !      $ make turbf CUTOUT_SUPPORT=1
    !   3) Uncomment the line below to Load the hdf5 cutout file:

    ! rc = turblibAddLocalSource(filename);

    !   4) Uncomment the code below to restrict target locations to within the data region downloaded:

    ! do i = 1, 10
    !   points(1, i) = 2 * 3.141592 / 1024.0 * i
    !   points(2, i) = 2 * 3.141592 / 1024.0 * i
    !   points(3, i) = 2 * 3.141592 / 1024.0 * i
    ! end do

    !  5) Uncomment the code below to call getVelocity, which will be evaluated locally as
    !     the time chosen is 0.0, which corresponds to timestep 0 and no spatial or temporal
    !     interpolation is requested:

    ! write(*,*)
    ! write(*,'(a)') 'Requesting velocity at 10 points...'
    ! rc = getvelocity(authkey, dataset,  0.0, NoSInt, NoTInt, 10, points, dataout3)
    ! do i = 1, 10
    !   write(*,format3) i, ': ', dataout3(1,i), ', ', dataout3(2,i), ', ', dataout3(3,i)
    ! end do

    ! In this sample code, the default time chosen is 0.364, so all of the function calls in the
    !   remainder of the sample code will be evaluated at the server.

    do i = 1, 10
        points(1, i) = 0.20 * i
        points(2, i) = 0.50 * i
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
    write(*,'(a)') 'Requesting forcing at 10 points...'
    rc = getforce(authkey, dataset, time, Lag6, NoTInt, 10, points, dataout3)
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
    write(*,'(2(a,f8.6),a)') 'Requesting position at 10 points, starting at time ', &
        startTime, ' and ending at time ', endTime, '...'
    rc = getposition(authkey, dataset, startTime, endTime, lag_dt, Lag6, 10, points, dataout3)

    write(*,*)
    write(*,'(a)') 'Coordinates of 10 points at startTime:'
    do i = 1, 10
        write(*,format3) i, ': ', points(1,i), ', ', points(2,i), ', ', points(3,i)
    end do

    write(*,*)
    write(*,'(a)') 'Coordinates of 10 points at endTime:'
    do i = 1, 10
        write(*,format3) i, ': ', dataout3(1,i), ', ', dataout3(2,i), ', ', dataout3(3,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting box filter of velocity at 10 points...'
    rc = getboxfilter(authkey, dataset, field, time, filterwidth, 10, points, dataout3)
    do i = 1, 10
        write(*,format3) i, ': ', dataout3(1,i), ', ', dataout3(2,i), ', ', dataout3(3,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting sub-grid stress symmetric tensor at 10 points...'
    rc = getboxfiltersgssymtensor(authkey, dataset, field, time, filterwidth, 10, points, dataout6)
    do i = 1, 10
        write(*,format6) i, ': xx=', dataout6(1,i), ', yy=', dataout6(2,i), &
            ', zz=', dataout6(3,i), ', xy=', dataout6(4,i),  &
            ', xz=', dataout6(5,i), ', yz', dataout6(6,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting sub-grid stress of two scalar fields at 10 points...'
    rc = getboxfiltersgsscalar(authkey, dataset, scalar_fields, time, filterwidth, 10, points, dataout1)
    do i = 1, 10
        write(*,format1) i, ': ', dataout1(i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting sub-grid stress of a vector and a scalar field at 10 points...'
    rc = getboxfiltersgsvector(authkey, dataset, vector_scalar_fields, time, filterwidth, 10, points, dataout3)
    do i = 1, 10
        write(*,format3) i, ': ', dataout3(1,i), ', ', dataout3(2,i), ', ', dataout3(3,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting box filter of velocity gradient at 10 points...'
    rc = getboxfiltergradient(authkey, dataset, field, time, filterwidth, spacing, 10, points, dataout9)
    do i = 1, 10
        write(*,format9) i, ': duxdx=', dataout9(1,i), ', duxdy=', dataout9(2,i), &
            ', duxdz=', dataout9(3,i), ', duydx=', dataout9(4,i),  &
            ', duydy=', dataout9(5,i), ', duydz=', dataout9(6,i),  &
            ', duzdx=', dataout9(7,i), ', duzdy=', dataout9(8,i),  &
            ', duzdz=', dataout9(9,i)
    end do

    write(*,*)
    write(*,'(a)') 'Requesting invariant at 10 points...'
    rc = getinvariant(authkey, dataset, time, FD4Lag4, NoTInt, 10, points, dataout2)
    do i = 1, 10
        write(*,format3) i, ': S2=', dataout2(1,i), ', O2=', dataout2(2,i)
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

