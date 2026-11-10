! Tests for trace_camera_ray (the CPU side of geodesic_kernel).
!
! Why these tests exist: the GPU kernel itself cannot be unit tested, but the
! ray loop inside it is host,device, so we run it on the CPU. If the loop is
! wrong (swapped status codes, missing escape condition, a bad step size) the
! picture just turns black or noisy, and finding the cause on screen is hard.
! These cases have answers we know independently of the code.
program test_geodesic_kernel
    use geodesic_kernel
    implicit none

    integer :: failures
    double precision, parameter :: tol = 1.0d-6

    ! Observer 20 units from the hole, capture just outside the horizon (r = 2).
    double precision, parameter :: r_obs = 20.0d0
    double precision, parameter :: r_capture = 2.05d0
    double precision, parameter :: r_escape = 100.0d0
    double precision, parameter :: step_scale = 0.005d0
    integer, parameter :: max_steps = 200000

    failures = 0

    call test_captured_below_critical()
    call test_escapes_above_critical()
    call test_radial_inward_and_outward()
    call test_step_limit()
    call test_critical_impact_parameter()
    call test_single_precision_matches()

    if (failures > 0) then
        print '(a,i0,a)', 'FAILED: ', failures, ' check(s)'
        stop 1
    end if
    print '(a)', 'test_geodesic_kernel: all checks passed'

contains

    subroutine check(name, got, expected, tolerance)
        character(*), intent(in) :: name
        double precision, intent(in) :: got, expected, tolerance

        if (abs(got - expected) > tolerance) then
            print '(a,a,es12.4,a,es12.4)', 'FAIL ', name, got, ' expected ', expected
            failures = failures + 1
        end if
    end subroutine check

    subroutine check_int(name, got, expected)
        character(*), intent(in) :: name
        integer, intent(in) :: got, expected

        if (got /= expected) then
            print '(a,a,i0,a,i0)', 'FAIL ', name, got, ' expected ', expected
            failures = failures + 1
        end if
    end subroutine check_int

    ! Why: the tests think in impact parameter b (a known physical number), but
    ! trace_camera_ray takes a direction. For a static observer at r_obs the
    ! sideways part of the direction is n_t = b*sqrt(f)/r_obs with
    ! f = 1 - 2/r_obs, and the radial part points inward so the ray flies at the hole.
    subroutine make_n_plane(b, n_plane)
        double precision, intent(in) :: b
        double precision, intent(out) :: n_plane(3)

        double precision :: f, n_t

        f = 1.0d0 - 2.0d0/r_obs
        n_t = b*sqrt(f)/r_obs
        n_plane(1) = -sqrt(1.0d0 - n_t**2)
        n_plane(2) = 0.0d0
        n_plane(3) = n_t
    end subroutine make_n_plane

    ! Why: b = 3 is well below the critical value 3*sqrt(3) = 5.196, so the
    ! ray must fall in. This fails if the capture test sets the wrong status.
    subroutine test_captured_below_critical()
        double precision :: n_plane(3), phi_end
        integer :: status

        call make_n_plane(3.0d0, n_plane)
        call trace_camera_ray(n_plane, r_obs, r_capture, r_escape, step_scale, max_steps, status, phi_end)
        call check_int('b=3 status', status, 1)
    end subroutine test_captured_below_critical

    ! Why: b = 8 is above the critical value, so the ray swings round the hole
    ! and escapes. phi_end is compared with 3.50968, obtained independently by
    ! integrating the same ray with a 10x smaller step in a separate Python
    ! script. The tolerance 1e-3 covers the coarse end-of-ray step.
    subroutine test_escapes_above_critical()
        double precision :: n_plane(3), phi_end
        integer :: status

        call make_n_plane(8.0d0, n_plane)
        call trace_camera_ray(n_plane, r_obs, r_capture, r_escape, step_scale, max_steps, status, phi_end)
        call check_int('b=8 status', status, 2)
        call check('b=8 phi_end', phi_end, 3.58974d0, 1.0d-3)
    end subroutine test_escapes_above_critical

    ! Why: the centre pixel of every image is a ray straight at the hole and the
    ! opposite ray straight away from it. Inward must fall in; outward must
    ! escape with phi unchanged (nothing bends a radial ray sideways). A missing
    ! p_r > 0 test or a wrong sign would show up here.
    subroutine test_radial_inward_and_outward()
        double precision :: n_plane(3), phi_end
        integer :: status

        n_plane = [-1.0d0, 0.0d0, 0.0d0]
        call trace_camera_ray(n_plane, r_obs, r_capture, r_escape, step_scale, max_steps, status, phi_end)
        call check_int('radial in status', status, 1)

        n_plane = [1.0d0, 0.0d0, 0.0d0]
        call trace_camera_ray(n_plane, r_obs, r_capture, r_escape, step_scale, max_steps, status, phi_end)
        call check_int('radial out status', status, 2)
        call check('radial out phi_end', phi_end, 0.0d0, tol)
    end subroutine test_radial_inward_and_outward

    ! Why: the loop must never run forever. With only 5 steps the ray cannot
    ! reach either boundary, so status must stay 0 (the "step limit" code the
    ! shader treats as black).
    subroutine test_step_limit()
        double precision :: n_plane(3), phi_end
        integer :: status

        call make_n_plane(8.0d0, n_plane)
        call trace_camera_ray(n_plane, r_obs, r_capture, r_escape, step_scale, 5, status, phi_end)
        call check_int('step limit status', status, 0)
    end subroutine test_step_limit

    ! Why: the strongest end-to-end check. Bisecting on b between "falls in" and
    ! "escapes" must land on the exact critical value 3*sqrt(3). It only works
    ! if the equations, the initial conditions, the step rule and both
    ! termination tests are all right at once.
    subroutine test_critical_impact_parameter()
        double precision :: lo, hi, mid, n_plane(3), phi_end
        integer :: status, k

        lo = 4.0d0
        hi = 7.0d0
        do k = 1, 40
            mid = 0.5d0*(lo + hi)
            call make_n_plane(mid, n_plane)
            call trace_camera_ray(n_plane, r_obs, r_capture, r_escape, step_scale, max_steps, status, phi_end)
            if (status == 1) then
                lo = mid
            else
                hi = mid
            end if
        end do
        call check('critical b', 0.5d0*(lo + hi), 3.0d0*sqrt(3.0d0), tol)
    end subroutine test_critical_impact_parameter

    ! Why: the GPU runs trace_camera_ray_sp, the single precision copy of the
    ! loop. It must give the same fate and the same phi_end to float accuracy
    ! as the double precision reference, or the picture silently drifts from
    ! the physics the tests above pin down. 1e-3 covers float round-off over
    ! the ~1000 steps of a b = 8 ray.
    subroutine test_single_precision_matches()
        double precision :: n_plane(3), phi_end
        real(4) :: n_plane_sp(3), phi_end_sp
        integer :: status, status_sp

        call make_n_plane(8.0d0, n_plane)
        call trace_camera_ray(n_plane, r_obs, r_capture, r_escape, step_scale, max_steps, status, phi_end)
        n_plane_sp = real(n_plane, 4)
        call trace_camera_ray_sp(n_plane_sp, real(r_obs,4), real(r_capture,4), real(r_escape,4), &
                                 real(step_scale,4), max_steps, status_sp, phi_end_sp)
        call check_int('sp b=8 status', status_sp, status)
        call check('sp b=8 phi_end', dble(phi_end_sp), phi_end, 1.0d-3)

        call make_n_plane(3.0d0, n_plane)
        n_plane_sp = real(n_plane, 4)
        call trace_camera_ray_sp(n_plane_sp, real(r_obs,4), real(r_capture,4), real(r_escape,4), &
                                 real(step_scale,4), max_steps, status_sp, phi_end_sp)
        call check_int('sp b=3 status', status_sp, 1)
    end subroutine test_single_precision_matches

end program test_geodesic_kernel