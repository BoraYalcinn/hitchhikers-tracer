module geodesic
    use cudafor
    use metric, only: metric_inv
    implicit none

    contains
    attributes(host,device) subroutine rhs(y,dy) 
    ! right hand side of the equation
    ! tells us the speed at that time  
    ! y(1)=t, y(2)=r, y(3)=θ, y(4)=φ, y(5)=p_t, y(6)=p_r, y(7)=p_θ, y(8)=p_φ
    double precision, intent(in)  :: y(8)
    double precision, intent(out) :: dy(8)

    double precision :: x(1:4)  
    double precision :: p(1:4)

    integer :: mu,nu,alpha,beta
    double precision :: ginv(4,4)
    double precision :: dginv(4,4,4)

    x(1:4) = y(1:4)
    p(1:4) = y(5:8)
    
    call metric_inv(x,ginv,dginv)
    dy = 0.0d0

    ! dx^mu/dlambda = sum_nu g^{mu nu} p_nu
    do mu = 1, 4
        do nu = 1, 4
            dy(mu) = dy(mu) + ginv(mu,nu) * p(nu)
        end do
    end do

    ! dp_mu/dlambda = -1/2 sum_{alpha,beta} (d_mu g^{alpha beta}) p_alpha p_beta
    do mu = 1, 4
        do alpha = 1, 4
            do beta = 1, 4
                dy(4+mu) = dy(4+mu) - 0.5d0 * dginv(mu,alpha,beta) * p(alpha) * p(beta)
            end do
        end do
    end do
    
    end subroutine rhs
 
    attributes(host,device) subroutine rk4_step(y,h)
    !
    double precision, intent(in)  :: h 
    double precision, intent(inout) :: y(8)

    double precision :: k1(8)
    double precision :: k2(8)
    double precision :: k3(8)
    double precision :: k4(8)
    double precision :: ytemp(8)

    end subroutine rk4_step


end module geodesic