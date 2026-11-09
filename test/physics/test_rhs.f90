program test_rhs
    ! Checks geodesic::rhs at one point against values worked out by hand.
    ! State y = (t, r, theta, phi, p_t, p_r, p_theta, p_phi).
    use geodesic, only: rhs
    implicit none

    double precision, parameter :: pi = 3.14159265358979323846d0
    double precision :: y(8), dy(8), expected(8)
    character(len=12) :: names(8)
    integer :: i, nfail

    nfail = 0

    ! Point: r = 6, theta = pi/3, momentum p = (-1, 0.5, 0, 2)
    y = [0.0d0, 6.0d0, pi/3.0d0, 0.0d0, &
         -1.0d0, 0.5d0, 0.0d0, 2.0d0]

    ! Hand-derived (Schwarzschild, M = 1):
    !   dx^mu/dlambda = g^{mu nu} p_nu
    !   dp_mu/dlambda = -1/2 (d_mu g^{ab}) p_a p_b
    expected = [ 1.5d0, &
                 1.0d0/3.0d0, &
                 0.0d0, &
                 2.0d0/27.0d0, &
                 0.0d0, &
                 -29.0d0/648.0d0, &
                 4.0d0/(27.0d0*sqrt(3.0d0)), &
                 0.0d0 ]

    names = [character(len=12) :: "dt", "dr", "dtheta", "dphi", &
             "dp_t", "dp_r", "dp_theta", "dp_phi"]

    call rhs(y, dy)

    do i = 1, 8
        call check(trim(names(i)), dy(i), expected(i))
    end do

    if (nfail > 0) then
        print '(a,i0,a)', "FAILED: ", nfail, " check(s) failed."
        stop 1
    end if
    print '(a)', "test_rhs: all checks passed."

contains

    subroutine check(name, got, want)
        character(len=*), intent(in) :: name
        double precision, intent(in) :: got, want
        double precision :: tol

        tol = 1.0d-12 * max(1.0d0, abs(want))
        if (abs(got - want) > tol) then
            print '(a,a,a,es23.15,a,es23.15)', "FAIL ", name, &
                  ": got ", got, " expected ", want
            nfail = nfail + 1
        end if
    end subroutine check

end program test_rhs