module boundary_conditions
    
    use parameters

    implicit none
    
contains

    subroutine Neumann(Rho,u,v)

        real(PR), intent(inout) :: Rho(0:imax+1,0:jmax+1), u(0:imax+1,0:jmax+1), v(0:imax+1,0:jmax+1)
        integer :: i, j

        do j = 1, jmax

            ! left
            Rho(0,j) = Rho(1,j)
            u(0,j) = u(1,j)
            v(0,j) = v(1,j)
   
            ! right
            Rho(imax+1,j) = Rho(imax,j)
            u(imax+1,j) = u(imax,j)
            v(imax+1,j) = v(imax,j)
   
        end do
   
        do i = 1, imax
   
            ! bottom
            Rho(i,0) = Rho(i,1)
            u(i,0) = u(i,1)
            v(i,0) = v(i,1)
   
            ! top
            Rho(i,jmax+1) = Rho(i,jmax)
            u(i,jmax+1) = u(i,jmax)
            v(i,jmax+1) = v(i,jmax)
   
        end do
        
    end subroutine Neumann

    subroutine Dirichlet(Rho,u,v,x,y,t,k)

        integer, intent(in) :: k
        real(PR), intent(in) :: x(0:imax), y(0:jmax), t
        real(PR), intent(inout) :: Rho(0:imax+1,0:jmax+1), u(0:imax+1,0:jmax+1), v(0:imax+1,0:jmax+1)
        integer :: i, j

        if (case == 'one') then ! The only case in which the exact solution is known

            Rho(i,j) = exp(-(x(i)+y(j)))
            u(i,j) = exp(-t)/(k*exp(-(x(i)+y(j))))
            v(i,j) = exp(-t)/(k*exp(-(x(i)+y(j))))

            do j = 1, jmax

                ! left
                Rho(0,j) = exp(-(x(0)+y(j)))
                u(0,j) = exp(-t)/(k*exp(-(x(0)+y(j))))
                v(0,j) = exp(-t)/(k*exp(-(x(0)+y(j))))
       
                ! right
                Rho(imax+1,j) = exp(-(x(imax)+y(j)))
                u(imax+1,j) = exp(-t)/(k*exp(-(x(imax)+y(j))))
                v(imax+1,j) = exp(-t)/(k*exp(-(x(imax)+y(j))))
       
            end do

            do i = 1, imax
   
                ! bottom
                Rho(i,0) = exp(-(x(i)+y(0)))
                u(i,0) = exp(-t)/(k*exp(-(x(i)+y(0))))
                v(i,0) = exp(-t)/(k*exp(-(x(i)+y(0))))
       
                ! top
                Rho(i,jmax+1) = exp(-(x(i)+y(jmax)))
                u(i,jmax+1) = exp(-t)/(k*exp(-(x(i)+y(jmax))))
                v(i,jmax+1) = exp(-t)/(k*exp(-(x(i)+y(jmax))))
       
            end do

        else

            print*, "Incompatible type of boundary condition for this test case"

        end if
        
    end subroutine Dirichlet
    
end module boundary_conditions