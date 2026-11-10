! Tests for ray_plane.
!
! Why these tests exist: a wrong rotation does not crash, it only makes the
! lensed sky appear mirrored or turned, which is hard to spot on screen. With
! no gravity a ray is a straight line, so the round trip
! to_equatorial -> (straight-line phi_end) -> from_equatorial must reproduce the
! camera direction exactly. That gives an exact reference without any physics.
program test_ray_plane
    use ray_plane
    implicit none

    integer :: failures
    double precision, parameter :: tol = 1.0d-12

    failures = 0

    call test_worked_example()
    call test_straight_radial_ray()
    call test_unit_length_kept()
    call test_flat_space_round_trip()

    if (failures > 0) then
        print '(a,i0,a)', 'FAILED: ', failures, ' check(s)'
        stop 1
    end if
    print '(a)', 'test_ray_plane: all checks passed'

contains

    subroutine check(name, got, expected, tolerance)
        character(*), intent(in) :: name
        double precision, intent(in) :: got, expected, tolerance

        if (abs(got - expected) > tolerance) then
            print '(a,a,es12.4,a,es12.4)', 'FAIL ', name, got, ' expected ', expected
            failures = failures + 1
        end if
    end subroutine check

    ! Why: the numbers worked out by hand in the module comment, so a reader can
    ! follow the code with a concrete case: n = (-0.6, 0.48, 0.64).
    subroutine test_worked_example()
        double precision :: n(3), n_plane(3), tangent(2)

        n = [-0.6d0, 0.48d0, 0.64d0]
        call to_equatorial(n, n_plane, tangent)
        call check('example n_r', n_plane(1), -0.6d0, tol)
        call check('example n_theta', n_plane(2), 0.0d0, tol)
        call check('example n_t', n_plane(3), 0.8d0, tol)
        call check('example tangent theta', tangent(1), 0.6d0, tol)
        call check('example tangent phi', tangent(2), 0.8d0, tol)
    end subroutine test_worked_example

    ! Why: the centre pixel has no tangential part, so the plane is undefined.
    ! A division by zero here would put NaN into the very middle of every image.
    subroutine test_straight_radial_ray()
        double precision :: n(3), n_plane(3), tangent(2)

        n = [-1.0d0, 0.0d0, 0.0d0]
        call to_equatorial(n, n_plane, tangent)
        call check('radial n_r', n_plane(1), -1.0d0, tol)
        call check('radial n_t', n_plane(3), 0.0d0, tol)
        call check('radial tangent length', sqrt(sum(tangent**2)), 1.0d0, tol)
    end subroutine test_straight_radial_ray

    ! Why: init_ray assumes a unit direction; the rotation must not change the
    ! length, otherwise every ray would start with a wrong energy.
    subroutine test_unit_length_kept()
        double precision :: n(3), n_plane(3), tangent(2)

        n = [-0.3d0, -0.5d0, sqrt(1.0d0 - 0.09d0 - 0.25d0)]
        call to_equatorial(n, n_plane, tangent)
        call check('length kept', sqrt(sum(n_plane**2)), 1.0d0, tol)
        call check('tangent unit', sqrt(sum(tangent**2)), 1.0d0, tol)
    end subroutine test_unit_length_kept

    ! Why: the exact end-to-end check. Without bending the final angle is
    ! phi_end = atan2(n_t, n_r), and the sky direction must equal the camera
    ! direction in Cartesian axes (n_r, n_phi, -n_theta). Cases cover looking
    ! up, down, left, right, outward and straight at the hole.
    subroutine test_flat_space_round_trip()
        double precision :: n(3, 6), n_plane(3), tangent(2), d(3), phi_end
        integer :: k

        n(:, 1) = [-0.6d0, 0.48d0, 0.64d0]
        n(:, 2) = [-0.6d0, -0.48d0, 0.64d0]
        n(:, 3) = [-0.6d0, 0.48d0, -0.64d0]
        n(:, 4) = [0.6d0, -0.48d0, 0.64d0]
        n(:, 5) = [-1.0d0, 0.0d0, 0.0d0]
        n(:, 6) = [0.0d0, 0.0d0, 1.0d0]

        do k = 1, 6
            call to_equatorial(n(:, k), n_plane, tangent)
            phi_end = atan2(n_plane(3), n_plane(1))
            call from_equatorial(phi_end, tangent, d)
            call check('round trip x', d(1), n(1, k), tol)
            call check('round trip y', d(2), n(3, k), tol)
            call check('round trip z', d(3), -n(2, k), tol)
        end do
    end subroutine test_flat_space_round_trip

end program test_ray_plane
