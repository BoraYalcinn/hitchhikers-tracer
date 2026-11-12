program test_observer
    ! Checks observer::static_tetrad: the four leg vectors of a static
    ! observer must be orthonormal in the metric, g(e_a, e_b) = eta_ab with
    ! eta = diag(-1,+1,+1,+1), and match the hand-computed values.
    use observer, only: static_tetrad
    use metric, only: metric_cov
    implicit none

    double precision, parameter :: pi = 3.14159265358979323846d0
    double precision :: x(4), e(4,4), gcov(4,4)
    double precision :: g_ab, want
    integer :: a, b, nfail

    nfail = 0

    ! ------------------------------------------------------------------
    ! Case 1: r = 20, theta = pi/2 (hand-computed values).
    ! f = 0.9, so e_t^t = 1/sqrt(0.9), e_r^r = sqrt(0.9),
    ! e_theta^theta = 1/20, e_phi^phi = 1/(20*sin(pi/2)) = 1/20.
    ! ------------------------------------------------------------------
    x = [0.0d0, 20.0d0, pi/2.0d0, 0.0d0]
    call static_tetrad(x, e)

    call check("e_t^t",         e(1,1), 1.0d0/sqrt(0.9d0), 1.0d-12)
    call check("e_r^r",         e(2,2), sqrt(0.9d0),       1.0d-12)
    call check("e_theta^theta", e(3,3), 0.05d0,            1.0d-12)
    call check("e_phi^phi",     e(4,4), 0.05d0,            1.0d-12)

    ! Each leg lives on a single coordinate axis: everything off the
    ! diagonal must be zero.
    do b = 1, 4
        do a = 1, 4
            if (a /= b) call check("off-diagonal", e(a,b), 0.0d0, 1.0d-14)
        end do
    end do

    ! ------------------------------------------------------------------
    ! Case 2: orthonormality g_{mu nu} e_a^mu e_b^nu = eta_ab, away from
    ! the equator (theta = 1.0) so that sin(theta) /= 1 is tested too.
    ! ------------------------------------------------------------------
    x = [0.0d0, 6.0d0, 1.0d0, 0.5d0]
    call static_tetrad(x, e)
    call metric_cov(x, gcov)

    do b = 1, 4
        do a = 1, 4
            g_ab = dot_product(e(:,a), matmul(gcov, e(:,b)))
            if (a /= b) then
                want = 0.0d0
            else if (a == 1) then
                want = -1.0d0
            else
                want = 1.0d0
            end if
            call check("g(e_a,e_b)", g_ab, want, 1.0d-12)
        end do
    end do

    if (nfail > 0) then
        print '(a,i0,a)', "FAILED: ", nfail, " check(s) failed."
        stop 1
    end if
    print '(a)', "test_observer: all checks passed."

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

end program test_observer
