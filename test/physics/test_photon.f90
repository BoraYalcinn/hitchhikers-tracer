program test_photon
    ! Checks photon::init_ray: initial state, null condition, normalization
    ! of the direction, and the conserved quantities E and L.
    use photon, only: init_ray, hamiltonian
    implicit none

    double precision, parameter :: pi = 3.14159265358979323846d0
    double precision, parameter :: b_crit = 5.196152422706632d0   ! 3*sqrt(3)
    double precision :: y(8), y2(8), expected(8)
    double precision :: n(3), nhat(3), norm
    double precision :: r, theta, phi, f, E, L, b
    character(len=12) :: names(8)
    integer :: i, nfail

    nfail = 0
    names = [character(len=12) :: "t", "r", "theta", "phi", &
             "p_t", "p_r", "p_theta", "p_phi"]

    ! ------------------------------------------------------------------
    ! Case 1: equatorial ray, direction (-0.6, 0, 0.8) at r = 6.
    ! f = 2/3, so p_t = -sqrt(f), p_r = n_r/sqrt(f), p_phi = r*n_phi.
    ! ------------------------------------------------------------------
    r = 6.0d0
    call init_ray(r, pi/2.0d0, 0.0d0, [-0.6d0, 0.0d0, 0.8d0], y)

    expected = [0.0d0, 6.0d0, pi/2.0d0, 0.0d0, &
                -sqrt(2.0d0/3.0d0), &
                -0.6d0/sqrt(2.0d0/3.0d0), &
                0.0d0, &
                4.8d0]
    do i = 1, 8
        call check(trim(names(i)), y(i), expected(i), 1.0d-12)
    end do

    ! The ray is light: H must vanish.
    call check("H (case 1)", hamiltonian(y), 0.0d0, 1.0d-14)

    ! Conserved quantities and impact parameter.
    E = -y(5)
    L = y(8)
    b = L / E
    call check("E", E, sqrt(2.0d0/3.0d0), 1.0d-12)
    call check("L", L, 4.8d0, 1.0d-12)
    call check("b", b, 4.8d0/sqrt(2.0d0/3.0d0), 1.0d-12)
    if (b <= b_crit) then
        print '(a,f9.5)', "FAIL b is not above the critical value: ", b
        nfail = nfail + 1
    end if

    ! ------------------------------------------------------------------
    ! Case 2: a direction of length 5 must give the same state as the
    ! unit direction (init_ray normalizes it).
    ! ------------------------------------------------------------------
    call init_ray(r, pi/2.0d0, 0.0d0, [-3.0d0, 0.0d0, 4.0d0], y2)
    do i = 1, 8
        call check("normalized " // trim(names(i)), y2(i), y(i), 1.0d-12)
    end do

    ! ------------------------------------------------------------------
    ! Case 3: off the equator, all three directions non-zero.
    ! H = 0, E = sqrt(f), L = r sin(theta) n_phi (unit n).
    ! ------------------------------------------------------------------
    r = 10.0d0
    theta = pi/3.0d0
    phi = 1.0d0
    n = [0.3d0, -0.5d0, 0.2d0]
    norm = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
    nhat = n / norm
    f = 1.0d0 - 2.0d0/r

    call init_ray(r, theta, phi, n, y)
    call check("H (case 3)", hamiltonian(y), 0.0d0, 1.0d-14)
    call check("E (case 3)", -y(5), sqrt(f), 1.0d-12)
    call check("L (case 3)", y(8), r*sin(theta)*nhat(3), 1.0d-12)
    call check("p_theta (case 3)", y(7), r*nhat(2), 1.0d-12)
    call check("p_r (case 3)", y(6), nhat(1)/sqrt(f), 1.0d-12)

    if (nfail > 0) then
        print '(a,i0,a)', "FAILED: ", nfail, " check(s) failed."
        stop 1
    end if
    print '(a)', "test_photon: all checks passed."

contains

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

end program test_photon