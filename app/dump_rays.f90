program dump_rays
    ! Traces a handful of rays past a Schwarzschild black hole and writes
    ! their paths to output/rays.csv, so tools/plot_rays.py can draw them.
    ! The output/ folder is for generated files and is not committed.
    !
    ! Why this exists: the tests prove the numbers are right; this lets us
    ! see the paths. It uses the same chain as the tests (init_ray, then
    ! trace_ray), so the picture shows exactly what the validated code does.
    !
    ! trace_ray is called in short chunks. After each chunk the current
    ! position is written to the file; the status tells us when the ray's
    ! fate is known. This way trace_ray itself needs no sampling code.
    !
    ! Output columns: ray_id, b, r, phi   (equatorial plane, M = 1)
    ! Run from the repository root (paths are relative to it):
    !   fpm run dump_rays --compiler nvfortran --flag "-cuda" --link-flag "-cuda"
    use geodesic, only: trace_ray
    use photon, only: init_ray
    implicit none

    double precision, parameter :: pi = 3.14159265358979323846d0
    double precision, parameter :: r0 = 30.0d0         ! start radius
    double precision, parameter :: r_capture = 2.05d0  ! just outside the horizon
    double precision, parameter :: r_escape = 40.0d0
    double precision, parameter :: h = 0.02d0          ! RK4 step
    integer, parameter :: chunk = 10                   ! steps between samples
    integer, parameter :: max_chunks = 20000           ! limit per ray
    integer, parameter :: nrays = 9

    double precision :: bs(nrays)
    double precision :: y(8), f, nphi
    integer :: unit, i, c, status, nsteps

    ! Impact parameters: below, near and above b_c = 3*sqrt(3) = 5.196.
    bs = [2.0d0, 3.0d0, 4.0d0, 5.0d0, 5.15d0, 5.25d0, 6.0d0, 8.0d0, 12.0d0]

    ! Fortran cannot create folders, so ask the shell (Linux/macOS).
    call execute_command_line("mkdir -p output")

    open(newunit=unit, file="output/rays.csv", status="replace", action="write")
    write(unit, '(a)') "ray_id,b,r,phi"

    do i = 1, nrays
        ! Ingoing equatorial ray from r0 with impact parameter bs(i);
        ! same conversion as in test_trace: b = r n_phi / sqrt(f).
        f = 1.0d0 - 2.0d0 / r0
        nphi = bs(i) * sqrt(f) / r0
        call init_ray(r0, pi/2.0d0, 0.0d0, &
                      [-sqrt(1.0d0 - nphi**2), 0.0d0, nphi], y)

        call write_point(unit, i, bs(i), y)
        do c = 1, max_chunks
            call trace_ray(y, h, chunk, r_capture, r_escape, status, nsteps)
            call write_point(unit, i, bs(i), y)
            if (status /= 0) exit      ! captured or escaped
        end do
    end do

    close(unit)
    print '(a)', "Wrote output/rays.csv"

contains

    subroutine write_point(unit, id, b, y)
        integer, intent(in) :: unit, id
        double precision, intent(in) :: b, y(8)

        write(unit, '(i0,a,f0.4,a,es16.8,a,es16.8)') &
              id, ",", b, ",", y(2), ",", y(4)
    end subroutine write_point

end program dump_rays