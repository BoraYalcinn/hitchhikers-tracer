program test_trace
    ! End-to-end physics checks for Schwarzschild null geodesics, using the
    ! full chain: photon::init_ray -> geodesic::trace_ray (RK4 + rhs + metric).
    !
    !   1. photon sphere: a tangential ray at r = 3 stays on its circular orbit
    !   2. capture and escape: b below / above the critical value
    !   3. critical impact parameter b_c = 3*sqrt(3) found by bisection
    !   4. weak-field deflection angle against the series in 1/b
    use geodesic, only: trace_ray
    use photon, only: init_ray, hamiltonian
    implicit none

    double precision, parameter :: pi = 3.14159265358979323846d0
    double precision, parameter :: b_crit = 5.196152422706632d0   ! 3*sqrt(3)
    double precision, parameter :: r_capture = 2.05d0             ! just outside r = 2
    integer :: nfail

    nfail = 0

    call test_photon_sphere()
    call test_capture_and_escape()
    call test_critical_b()
    call test_deflection(1.0d4, 100.0d0, 2.0d0)
    call test_deflection(2.0d4, 200.0d0, 4.0d0)

    if (nfail > 0) then
        print '(a,i0,a)', "FAILED: ", nfail, " check(s) failed."
        stop 1
    end if
    print '(a)', "test_trace: all checks passed."

contains

    subroutine launch(r0, b, y)
        ! Ingoing ray in the equatorial plane from radius r0 with impact
        ! parameter b = L/E. Uses b = r sin(theta) n_phi / sqrt(f).
        double precision, intent(in) :: r0, b
        double precision, intent(out) :: y(8)
        double precision :: f, nphi

        f = 1.0d0 - 2.0d0 / r0
        nphi = b * sqrt(f) / r0
        call init_ray(r0, pi/2.0d0, 0.0d0, &
                      [-sqrt(1.0d0 - nphi**2), 0.0d0, nphi], y)
    end subroutine launch

    subroutine test_photon_sphere()
        ! Tangential ray at r = 3 has b = 3*sqrt(3), p_r = 0: a circular
        ! (unstable) orbit. With E = sqrt(f), L = r, the orbit has
        !   dphi/dlambda = L/r^2 = 1/3,
        ! so after lambda = 20 we expect phi = 20/3 and r still 3.
        double precision :: y(8), h
        integer :: status, nsteps, nsteps_wanted

        h = 0.01d0
        nsteps_wanted = 2000
        call init_ray(3.0d0, pi/2.0d0, 0.0d0, [0.0d0, 0.0d0, 1.0d0], y)
        call trace_ray(y, h, nsteps_wanted, r_capture, 1.0d9, status, nsteps)

        call check_int("photon sphere status (step limit)", status, 0)
        call check_int("photon sphere steps", nsteps, nsteps_wanted)
        call check("photon sphere r", y(2), 3.0d0, 1.0d-9)
        call check("photon sphere phi", y(4), 20.0d0/3.0d0, 1.0d-9)
        call check("photon sphere H", hamiltonian(y), 0.0d0, 1.0d-10)
    end subroutine test_photon_sphere

    subroutine test_capture_and_escape()
        double precision :: y(8)
        integer :: status, nsteps

        call launch(30.0d0, 4.0d0, y)
        call trace_ray(y, 0.02d0, 100000, r_capture, 40.0d0, status, nsteps)
        call check_int("b = 4 is captured", status, 1)

        call launch(30.0d0, 7.0d0, y)
        call trace_ray(y, 0.02d0, 100000, r_capture, 40.0d0, status, nsteps)
        call check_int("b = 7 escapes", status, 2)
    end subroutine test_capture_and_escape

    subroutine test_critical_b()
        ! Bisection between a captured and an escaping impact parameter.
        double precision :: y(8), lo, hi, mid
        integer :: status, nsteps, it
        logical :: inconclusive

        lo = 4.0d0
        hi = 7.0d0
        inconclusive = .false.
        do it = 1, 40
            mid = 0.5d0 * (lo + hi)
            call launch(30.0d0, mid, y)
            call trace_ray(y, 0.02d0, 100000, r_capture, 40.0d0, status, nsteps)
            if (status == 1) then
                lo = mid
            else if (status == 2) then
                hi = mid
            else
                inconclusive = .true.
                exit
            end if
        end do

        if (inconclusive) then
            print '(a)', "FAIL critical b: step limit reached during bisection"
            nfail = nfail + 1
        else
            print '(a,f12.9,a,es10.2)', "  critical b = ", 0.5d0*(lo + hi), &
                  "  error = ", 0.5d0*(lo + hi) - b_crit
            call check("critical b", 0.5d0*(lo + hi), b_crit, 1.0d-6)
        end if
    end subroutine test_critical_b

    subroutine test_deflection(r0, b, h)
        ! Total change of phi minus the straight-line (flat space) value for
        ! the same start and end radii gives the deflection angle alpha.
        ! Reference series (M = 1):
        !   alpha = 4/b + 15*pi/(4 b^2) + 128/(3 b^3)
        double precision, intent(in) :: r0, b, h
        double precision :: y(8), flat, alpha, series
        integer :: status, nsteps
        character(len=24) :: label

        write(label, '(a,f6.1)') "deflection b = ", b

        call launch(r0, b, y)
        call trace_ray(y, h, 1000000, r_capture, r0 * 1.0001d0, status, nsteps)
        call check_int(trim(label) // " escapes", status, 2)

        flat = acos(b / r0) + acos(b / y(2))
        alpha = y(4) - flat
        series = 4.0d0/b + 15.0d0*pi/(4.0d0*b**2) + 128.0d0/(3.0d0*b**3)

        print '(a,f6.1,a,f10.7,a,f10.7)', "  b = ", b, "  alpha = ", alpha, &
              "  series = ", series
        call check_rel(trim(label) // " vs series", alpha, series, 5.0d-4)
        call check_rel(trim(label) // " vs 4/b", alpha, 4.0d0/b, 5.0d-2)
    end subroutine test_deflection

    subroutine check(name, got, want, tol_abs)
        character(len=*), intent(in) :: name
        double precision, intent(in) :: got, want, tol_abs
        double precision :: tol

        tol = tol_abs * max(1.0d0, abs(want))
        if (abs(got - want) > tol) then
            print '(a,a,a,es23.15,a,es23.15)', "FAIL ", name, &
                  ": got ", got, " expected ", want
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine check_rel(name, got, want, tol_rel)
        character(len=*), intent(in) :: name
        double precision, intent(in) :: got, want, tol_rel

        if (abs(got - want) > tol_rel * abs(want)) then
            print '(a,a,a,es23.15,a,es23.15)', "FAIL ", name, &
                  ": got ", got, " expected ", want
            nfail = nfail + 1
        end if
    end subroutine check_rel

    subroutine check_int(name, got, want)
        character(len=*), intent(in) :: name
        integer, intent(in) :: got, want

        if (got /= want) then
            print '(a,a,a,i0,a,i0)', "FAIL ", name, ": got ", got, &
                  " expected ", want
            nfail = nfail + 1
        end if
    end subroutine check_int

end program test_trace
