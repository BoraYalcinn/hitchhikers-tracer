program bench_rays
    use geodesic, only: trace_ray
    use photon, only: init_ray
    implicit none

    double precision, parameter :: pi = 3.14159265358979323846d0
    double precision, parameter :: r0 = 30.0d0         ! start radius
    double precision, parameter :: r_capture = 2.05d0  ! just outside the horizon
    double precision, parameter :: r_escape = 40.0d0
    double precision, parameter :: h = 0.02d0          ! RK4 step
    integer, parameter :: max_steps = 200000           ! limit per ray

    integer :: nrays
    double precision, allocatable :: b(:), phi_serial(:), phi_par(:)
    integer, allocatable :: status_serial(:), status_par(:)

    double precision :: y(8), f, nphi
    integer :: status, nsteps

    integer :: k
    character(32) :: arg
    integer(8) :: t0, t1, rate
    double precision :: ms_serial, ms_par

    ! Number of rays: default 2000, optional first command-line argument.
    nrays = 2000
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) nrays
    end if

    allocate(b(nrays), phi_serial(nrays), phi_par(nrays))
    allocate(status_serial(nrays), status_par(nrays))

    do k = 1, nrays
        b(k) = 4.0d0 + 8.0d0 * dble(k - 1) / dble(nrays - 1)
    end do

    call system_clock(count_rate=rate)

    ! ---- serial reference ----------------------------------------------
    call system_clock(t0)
    do k = 1, nrays
        f = 1.0d0 - 2.0d0 / r0
        nphi = b(k) * sqrt(f) / r0
        call init_ray(r0, pi/2.0d0, 0.0d0, &
                      [-sqrt(1.0d0 - nphi**2), 0.0d0, nphi], y)
        call trace_ray(y, h, max_steps, r_capture, r_escape, status, nsteps)
        status_serial(k) = status
        phi_serial(k) = y(4)
    end do
    call system_clock(t1)
    ms_serial = 1000.0d0 * dble(t1 - t0) / dble(rate)

    ! ---- parallel version ----------------------------------------------
    call system_clock(t0)
    do concurrent (k = 1:nrays) local(y, f, nphi, status, nsteps)
        f = 1.0d0 - 2.0d0 / r0
        nphi = b(k) * sqrt(f) / r0
        call init_ray(r0, pi/2.0d0, 0.0d0, &
                      [-sqrt(1.0d0 - nphi**2), 0.0d0, nphi], y)
        call trace_ray(y, h, max_steps, r_capture, r_escape, status, nsteps)
        status_par(k) = status
        phi_par(k) = y(4)
    end do
    call system_clock(t1)
    ms_par = 1000.0d0 * dble(t1 - t0) / dble(rate)

    ! ---- compare -------------------------------------------------------
    print '(a,i0)', "rays:                 ", nrays
    print '(a,f10.1,a)', "serial:               ", ms_serial, " ms"
    print '(a,f10.1,a)', "do concurrent:        ", ms_par, " ms"
    print '(a,f10.2)',   "speedup:              ", ms_serial / ms_par
    print '(a,i0)',      "status mismatches:    ", count(status_serial /= status_par)
    print '(a,es10.2)',  "max |phi difference|: ", maxval(abs(phi_serial - phi_par))
    print '(a,i0,a,i0,a,i0)', "captured / escaped / limit: ", &
          count(status_par == 1), " / ", count(status_par == 2), " / ", count(status_par == 0)

end program bench_rays