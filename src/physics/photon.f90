module photon
    use cudafor
    use metric, only: metric_cov, metric_inv, metric_cov_sp
    implicit none

    contains
    attributes(host,device) pure subroutine init_ray(r,theta,phi,n,y)
    
        double precision, intent(in) :: r
        double precision, intent(in) :: theta
        double precision, intent(in) :: phi
        double precision, intent(in) :: n(3)
        double precision, intent(out) :: y(8)

        double precision :: gcov(4,4)
        double precision :: x(4)
        double precision :: normalized_n(3)

        double precision :: k(4)
        double precision :: p(4)

        integer :: nu
        integer :: mu

        x = [0.0d0,r,theta,phi]
        call metric_cov(x,gcov)
        normalized_n = n / (sqrt(n(1)**2 + n(2)**2 +n(3)**2))

        k(2) = normalized_n(1)/sqrt(gcov(2,2))
        k(3) = normalized_n(2)/sqrt(gcov(3,3))
        k(4) = normalized_n(3)/sqrt(gcov(4,4))
        k(1) = 1.0d0/sqrt(-gcov(1,1))

        p = 0.0d0
        do mu = 1, 4
            do nu = 1, 4
                p(mu) = p(mu) + gcov(mu,nu) * k(nu)
            end do
        end do

        y(1:4) = x 
        y(5:8) = p

    end subroutine
    
    ! Single precision copy of init_ray for the GPU renderer (see rhs_sp).
    attributes(host,device) pure subroutine init_ray_sp(r,theta,phi,n,y)

        real(4), intent(in) :: r
        real(4), intent(in) :: theta
        real(4), intent(in) :: phi
        real(4), intent(in) :: n(3)
        real(4), intent(out) :: y(8)

        real(4) :: gcov(4,4)
        real(4) :: x(4)
        real(4) :: normalized_n(3)

        real(4) :: k(4)
        real(4) :: p(4)

        integer :: nu
        integer :: mu

        x = [0.0e0,r,theta,phi]
        call metric_cov_sp(x,gcov)
        normalized_n = n / (sqrt(n(1)**2 + n(2)**2 +n(3)**2))

        k(2) = normalized_n(1)/sqrt(gcov(2,2))
        k(3) = normalized_n(2)/sqrt(gcov(3,3))
        k(4) = normalized_n(3)/sqrt(gcov(4,4))
        k(1) = 1.0e0/sqrt(-gcov(1,1))

        p = 0.0e0
        do mu = 1, 4
            do nu = 1, 4
                p(mu) = p(mu) + gcov(mu,nu) * k(nu)
            end do
        end do

        y(1:4) = x
        y(5:8) = p

    end subroutine init_ray_sp

    attributes(host,device) function hamiltonian(y) result(h)
        
        double precision, intent(in) :: y(8)
        double precision :: h 
        double precision :: p(4)
        double precision :: ginv(4,4)
        double precision :: dginv(4,4,4)

        integer :: nu
        integer :: mu

        call metric_inv(y(1:4), ginv, dginv)
        p = y(5:8)
        h = 0.0d0

        do mu = 1, 4
            do nu = 1, 4
                h = h + 0.5d0 * ginv(mu,nu) * p(mu) * p(nu)
            end do
        end do


    end function
    

end module photon

