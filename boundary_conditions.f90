module boundary_conditions
    
    use parameters

    implicit none
    
contains

    subroutine Neumann(Rho,u,v)

        real(PR), intent(inout) :: Rho(1:imax+1,1:jmax+1), u(1:imax+1,1:jmax+1), v(1:imax+1,1:jmax+1)
        integer :: i, j

        do j = 2, jmax

            ! left
            Rho(1,j) = Rho(2,j)
            u(1,j) = u(2,j)
            v(1,j) = v(2,j)
   
            ! right
            Rho(imax+1,j) = Rho(imax,j)
            u(imax+1,j) = u(imax,j)
            v(imax+1,j) = v(imax,j)
   
        end do
   
        do i = 2, imax
   
            ! top
            Rho(i,1) = Rho(i,1)
            u(i,1) = u(i,1)
            v(i,1) = v(i,1)
   
            ! buttom
            Rho(i,jmax+1) = Rho(i,jmax)
            u(i,jmax+1) = u(i,jmax)
            v(i,jmax+1) = v(i,jmax)
   
        end do
        
    end subroutine Neumann
    
end module boundary_conditions