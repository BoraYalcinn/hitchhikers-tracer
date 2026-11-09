program test_rk4
    ! Integrates an ingoing null ray in the equatorial plane and checks:
    !   1. E = -p_t and L = p_phi stay constant (they have zero derivative),
    !   2. the null condition H = 0 holds up to the RK4 truncation error,
    !   3. halving the step size reduces |H| by about 2^4 = 16 (4th order).
    use geodesic, only: rk4_step
    use metric, only: metric_inv
    implicit none

    double precision, parameter :: pi = 3.14159265358979323846d0
    double precision, parameter :: lambda_end = 2.0d0
    double precision :: y0(8), y_coarse(8), y_fine(8)
    double precision :: h_coarse, h_fine
    double precision :: H_start, H_coarse, H_fine, ratio
    integer :: n_coarse, n_fine, nfail

    nfail = 0

    call initial_ray(y0)
    H_start = hamiltonian(y0)

    ! Start must be null (up to rounding).
    call check_small("H at start", H_start, 1.0d-14)

    ! Same path with step h and h/2.
    n_coarse = 20
    n_fine = 40
    h_coarse = lambda_end / n_coarse
    h_fine = lambda_end / n_fine

    y_coarse = integrate(y0, n_coarse, h_coarse)
    y_fine = integrate(y0, n_fine, h_fine)

    ! 1. Conserved quantities: p_t and p_phi must not change at all.
    call check_small("E drift (coarse)", abs(y_coarse(5) - y0(5)), 1.0d-14)
    call check_small("L drift (coarse)", abs(y_coarse(8) - y0(8)), 1.0d-14)
    call check_small("E drift (fine)", abs(y_fine(5) - y0(5)), 1.0d-14)
    call check_small("L drift (fine)", abs(y_fine(8) - y0(8)), 1.0d-14)

    ! 2. Null condition drift stays small.
    H_coarse = abs(hamiltonian(y_coarse))
    H_fine = abs(hamiltonian(y_fine))
    call check_small("|H| coarse", H_coarse, 1.0d-6)
    call check_small("|H| fine", H_fine, 1.0d-7)

    ! 3. Convergence order: expect a ratio near 16.
    ratio = H_coarse / H_fine
    print '(a,es12.4,a,es12.4,a,f8.3)', "  |H| coarse = ", H_coarse, &
          "  |H| fine = ", H_fine, "  ratio = ", ratio
    if (ratio < 8.0d0 .or. ratio > 24.0d0) then
        print '(a,f8.3,a)', "FAIL convergence ratio ", ratio, &
              " is outside [8, 24] (expected about 16)."
        nfail = nfail + 1
    end if

    ! The ray must still be outside the photon sphere region we assumed.
    if (y_fine(2) < 3.0d0) then
        print '(a)', "FAIL ray reached r < 3, test setup is no longer valid."
        nfail = nfail + 1
    end if

    if (nfail > 0) then
        print '(a,i0,a)', "FAILED: ", nfail, " check(s) failed."
        stop 1
    end if
    print '(a)', "test_rk4: all checks passed."

contains

    subroutine initial_ray(y)
        ! Ingoing photon at r = 6 in the equatorial plane.
        ! p_t = -1 (E = 1), p_phi = 4 (b = 4 < b_c, so it falls inward),
        ! p_r from the null condition H = 0, negative sign = ingoing.
        double precision, intent(out) :: y(8)
        double precision :: r0, f, pt, pphi, pr

        r0 = 6.0d0
        f = 1.0d0 - 2.0d0 / r0
        pt = -1.0d0
        pphi = 4.0d0
        pr = -sqrt((pt**2 / f - pphi**2 / r0**2) / f)

        y = [0.0d0, r0, pi/2.0d0, 0.0d0, pt, pr, 0.0d0, pphi]
    end subroutine initial_ray

    function hamiltonian(y) result(H)
        ! H = 1/2 g^{mu nu} p_mu p_nu
        double precision, intent(in) :: y(8)
        double precision :: H
        double precision :: ginv(4,4), dginv(4,4,4), p(4)
        integer :: mu, nu

        call metric_inv(y(1:4), ginv, dginv)
        p = y(5:8)
        H = 0.0d0
        do mu = 1, 4
            do nu = 1, 4
                H = H + 0.5d0 * ginv(mu,nu) * p(mu) * p(nu)
            end do
        end do
    end function hamiltonian

    function integrate(ystart, nsteps, h) result(y)
        double precision, intent(in) :: ystart(8), h
        integer, intent(in) :: nsteps
        double precision :: y(8)
        integer :: i

        y = ystart
        do i = 1, nsteps
            call rk4_step(y, h)
        end do
    end function integrate

    subroutine check_small(name, value, limit)
        character(len=*), intent(in) :: name
        double precision, intent(in) :: value, limit

        if (abs(value) > limit) then
            print '(a,a,a,es12.4,a,es12.4)', "FAIL ", name, &
                  ": ", value, " exceeds ", limit
            nfail = nfail + 1
        end if
    end subroutine check_small

end program test_rk4