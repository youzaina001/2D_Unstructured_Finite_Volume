module initial_condition

    use parameters
    use flux

    implicit none
    
contains

    subroutine exact_solution(Rho,u,v,x,y,k,t)

        integer, intent(in) :: k
        real(PR), intent(in) :: t
        real(PR), intent(in) :: x(1:imax), y(1:jmax)
        real(PR), intent(inout) :: Rho(0:imax+1,0:jmax+1), u(0:imax+1,0:jmax+1), v(0:imax+1,0:jmax+1)
        integer :: i, j

        if (case == 'one') then

            do j = 1, jmax
                do i = 1, imax

                    Rho(i,j) = exp(-(x(i)+y(j)))
                    u(i,j) = exp(-t)/(k*exp(-(x(i)+y(j))))
                    v(i,j) = exp(-t)/(k*exp(-(x(i)+y(j))))                    
                   
                end do          
            end do

        else if (case == 'two') then

            print*, "There is no exact solution for this test case"

        else if (case == 'thr') then

            print*, "There is no exact solution for this test case"
            
        end if

    end subroutine exact_solution

    subroutine init(Rho,u,v,x,y,k,t)

        integer, intent(in) :: k
        real(PR), intent(in) :: t
        real(PR), intent(in) :: x(1:imax), y(1:jmax)
        real(PR), intent(inout) :: Rho(0:imax+1,0:jmax+1), u(0:imax+1,0:jmax+1), v(0:imax+1,0:jmax+1)
        integer :: i, j

        if (case == 'one') then

            call exact_solution(Rho,u,v,x,y,k,t)

        else if (case == 'two') then

            do j = 1, jmax
                do i = 1, imax

                    Rho(i,j) = exp(-(x(i)-50._PR)**2-(y(i)-50._PR)**2) + 1._PR
                    u(i,j) = 0._PR
                    v(i,j) = 0._PR
                                   
                end do
            end do

        else if (case == 'thr') then

            do j = 1, jmax
                do i = 1, imax

                    if (((x(i)-0.5_PR)**2 + (y(j)-0.4_PR)**4) <= 0.3_PR**2) then

                        Rho(i,j) = 5._PR
                        u(i,j) = 0._PR
                        v(i,j) = 0._PR

                    else

                        Rho(i,j) = 1._PR
                        u(i,j) = 0._PR
                        v(i,j) = 0._PR
                        
                    end if
                   
                end do
            end do
            
        end if
 
    end subroutine init

end module initial_condition
