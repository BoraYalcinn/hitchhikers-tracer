! Tests for camera_rays.
!
! Why these tests exist: pixel_direction is pure geometry, and the mistakes it
! can hide (a flipped sign, a swapped axis, integer division in the aspect
! ratio) do not crash anything. They only show up later as a mirrored or
! stretched black hole image, which is very hard to debug on screen. Checking
! a few pixels whose direction we can compute by hand catches them now.
program test_camera_rays
    use camera_rays
    implicit none

    integer :: failures
    double precision, parameter :: tol = 1.0d-12

    failures = 0

    call test_degrees_to_radians()
    call test_centre_pixel()
    call test_unit_length()
    call test_edge_directions()
    call test_aspect_ratio()
    call test_fov_scaling()

    if (failures > 0) then
        print '(a,i0,a)', 'FAILED: ', failures, ' check(s)'
        stop 1
    end if
    print '(a)', 'test_camera_rays: all checks passed'

contains

    subroutine check(name, got, expected, tolerance)
        character(*), intent(in) :: name
        double precision, intent(in) :: got, expected, tolerance

        if (abs(got - expected) > tolerance) then
            print '(a,a,es12.4,a,es12.4)', 'FAIL ', name, got, ' expected ', expected
            failures = failures + 1
        end if
    end subroutine check

    ! Why: the FOV conversion is the only place degrees enter the program; a
    ! wrong factor here would silently change the field of view of every image.
    subroutine test_degrees_to_radians()
        call check('deg 180', degrees_to_radians(180.0d0), pi, tol)
        call check('deg 60', degrees_to_radians(60.0d0), pi/3.0d0, tol)
        call check('deg 0', degrees_to_radians(0.0d0), 0.0d0, tol)
    end subroutine test_degrees_to_radians

    ! Why: the middle of the screen must look straight at the black hole,
    ! n = (-1, 0, 0). If this fails, everything else is shifted or rotated.
    ! A 4x4 image has no centre pixel, so the centre point is the corner of
    ! pixel (2,2), reached with offset 0; a 5x5 image has a true centre pixel.
    subroutine test_centre_pixel()
        double precision :: n(3)

        call pixel_direction(2, 2, 0.0d0, 0.0d0, 4, 4, degrees_to_radians(60.0d0), n)
        call check('centre4 n_r', n(1), -1.0d0, tol)
        call check('centre4 n_theta', n(2), 0.0d0, tol)
        call check('centre4 n_phi', n(3), 0.0d0, tol)

        call pixel_direction(2, 2, 0.5d0, 0.5d0, 5, 5, degrees_to_radians(60.0d0), n)
        call check('centre5 n_r', n(1), -1.0d0, tol)
        call check('centre5 n_theta', n(2), 0.0d0, tol)
        call check('centre5 n_phi', n(3), 0.0d0, tol)
    end subroutine test_centre_pixel

    ! Why: init_ray and the later rotation step assume a unit direction. A
    ! missing normalisation would not fail there, it would only bend rays
    ! slightly wrong, so we check it at the source for many pixels and offsets.
    subroutine test_unit_length()
        double precision :: n(3)
        integer :: i, j

        do j = 0, 8, 4
            do i = 0, 16, 4
                call pixel_direction(i, j, 0.25d0, 0.75d0, 17, 9, degrees_to_radians(75.0d0), n)
                call check('unit length', sqrt(sum(n**2)), 1.0d0, tol)
            end do
        end do
    end subroutine test_unit_length

    ! Why: this fixes the orientation convention (up = -theta, right = +phi).
    ! With a 90 degree FOV on a square image the image plane edge is at
    ! distance 1, so the edge directions are known exactly: (-1, 0, 1)/sqrt(2)
    ! at the right edge and (-1, -1, 0)/sqrt(2) at the top edge. A flipped sign
    ! would give a mirrored sky and is invisible to the length test above.
    subroutine test_edge_directions()
        double precision :: n(3)
        double precision :: s
        s = 1.0d0/sqrt(2.0d0)

        ! right edge, middle row (x_ndc = +1, y_ndc = 0)
        call pixel_direction(4, 2, 0.0d0, 0.0d0, 4, 4, degrees_to_radians(90.0d0), n)
        call check('right n_r', n(1), -s, tol)
        call check('right n_theta', n(2), 0.0d0, tol)
        call check('right n_phi', n(3), s, tol)

        ! top edge, middle column (x_ndc = 0, y_ndc = +1)
        call pixel_direction(2, 0, 0.0d0, 0.0d0, 4, 4, degrees_to_radians(90.0d0), n)
        call check('top n_r', n(1), -s, tol)
        call check('top n_theta', n(2), -s, tol)
        call check('top n_phi', n(3), 0.0d0, tol)
    end subroutine test_edge_directions

    ! Why: integer division in w/h would give 1 for a 16:9 image and stretch
    ! the picture horizontally without any error message. At the right edge
    ! of a 90 degree 16:9 image, n_phi / (-n_r) must equal 16/9.
    subroutine test_aspect_ratio()
        double precision :: n(3)

        call pixel_direction(16, 4, 0.0d0, 0.5d0, 16, 9, degrees_to_radians(90.0d0), n)
        call check('aspect ratio', n(3)/(-n(1)), 16.0d0/9.0d0, tol)
        call check('aspect n_theta', n(2), 0.0d0, tol)
    end subroutine test_aspect_ratio

    ! Why: fov_y is the vertical angle, so at the top edge the tilt must be
    ! tan(fov/2), not tan(fov). Checking 60 degrees catches a missing 0.5.
    subroutine test_fov_scaling()
        double precision :: n(3)

        call pixel_direction(2, 0, 0.0d0, 0.0d0, 4, 4, degrees_to_radians(60.0d0), n)
        call check('fov top tilt', n(2)/n(1), tan(degrees_to_radians(30.0d0)), tol)
    end subroutine test_fov_scaling

end program test_camera_rays
