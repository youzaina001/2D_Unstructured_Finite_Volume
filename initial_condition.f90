module initial_condition

    use parameters
    use flux

    implicit none
    
contains

    subroutine exact_solution(Rho,u,v,p,sigma,x,y,k,t,UU)

        integer, intent(in) :: k
        real(PR), intent(in) :: t
        real(PR), intent(in) :: x(1:imax+1), y(1:jmax+1)
        real(PR), intent(inout) :: Rho(1:imax+1,1:jmax+1), u(1:imax+1,1:jmax+1), v(1:imax+1,1:jmax+1)
        real(PR), intent(inout) :: p(1:imax+1,1:jmax+1), sigma(1:imax+1,1:jmax+1)
        real(PR), intent(out) ::  UU(1:imax+1,1:jmax+1,1:3)
        integer :: i, j

        if (case == 'one') then

            do j = 1, jmax + 1

                do i = 1, imax + 1

                    Rho(i,j) = exp(-(x(i)+y(j)))
                    u(i,j) = exp(-t)/(k*exp(-(x(i)+y(j))))
                    v(i,j) = exp(-t)/(k*exp(-(x(i)+y(j))))
                    p(i,j) = pressure(Rho(i,j))
                    sigma(i,j) = compute_sigma(x(i),y(j),t,k)
                    
                    call non_conservative_to_conservative(Rho(i,j),u(i,j),v(i,j),UU(i,j,1:3))
                   
                end do
          
            end do
            
        end if

    end subroutine exact_solution

    subroutine init(Rho,u,v,p,sigma,x,y,k,t,UU)

        integer, intent(in) :: k
        real(PR), intent(in) :: t
        real(PR), intent(in) :: x(1:imax+1), y(1:jmax+1)
        real(PR), intent(inout) :: Rho(1:imax+1,1:jmax+1), u(1:imax+1,1:jmax+1), v(1:imax+1,1:jmax+1)
        real(PR), intent(inout) :: p(1:imax+1,1:jmax+1), sigma(1:imax+1,1:jmax+1)
        real(PR), intent(out) ::  UU(1:imax+1,1:jmax+1,1:3)
        integer :: i, j

        call exact_solution(Rho,u,v,p,sigma,x,y,k,t,UU)
 
    end subroutine init

end module initial_condition
