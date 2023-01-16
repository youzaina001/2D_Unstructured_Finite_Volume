module time_scheme

    use parameters
    use flux
    use solver
    use initial_condition

    implicit none
    
contains

    subroutine time_step(Rho,u,v,xm,ym,t,k,dt)

        real(PR), dimension(0:imax+1,0:jmax+1), intent(in) :: Rho, u, v
        real(PR), intent(in) :: xm(1:imax), ym(1:jmax)
        real(PR), intent(in) :: t
        integer, intent(in) :: k
        real(PR), intent(out) :: dt
        real(PR), dimension(1:imax,1:jmax) :: sgm
        real(PR) :: be, dK, sigmaK, PP
        integer :: i, j

        ! Initialisation
        PP = pressure_prime(maxval(Rho))
        be = max(maxval(abs(u)+sqrt(PP)),maxval(abs(v)+sqrt(PP)))
        dK = dx*dy/2._PR*(dx+dy)

        do j = 1, jmax

            do i = 1, imax

                sgm = compute_sigma(Rho(i,j),xm(i),ym(j),t,k)

            end do

        end do

        sigmaK = maxval(sgm)

        dt = CFL * dK / (4._PR * be &
           & + sigmaK*dK) 

    end subroutine time_step

    subroutine explicit_euler(Un,t,dt,kms,xm,ym,Unp1)

        integer, intent(in) :: kms 
        real(PR), intent(in) :: Un(0:imax+1,0:jmax+1,1:3), xm(1:imax), ym(1:jmax), t, dt
        real(PR), intent(out) ::  Unp1(0:imax+1,0:jmax+1,1:3)
        real(PR) :: Sn(0:imax+1,0:jmax+1,1:3), phi(0:imax+1,0:jmax+1,1:3), sigma(0:imax+1,0:jmax+1)
        integer :: i, j

        do j = 1, jmax

            do i = 1, imax
   
                call spatial_discretization(Un(i,j,1:3),Un(i-1,j,1:3),Un(i+1,j,1:3),Un(i,j+1,1:3),Un(i,j-1,1:3),phi(i,j,1:3))
                call source_term(Un(i,j,1:3),Sn(i,j,1:3))
                sigma(i,j) = compute_sigma(Un(i,j,1),xm(i),ym(j),t,kms)
                Unp1(i,j,1:3) = Un(i,j,1:3) + dt * phi(i,j,1:3) + dt * sigma(i,j) * Sn(i,j,1:3)
               
            end do
            
        end do
        
    end subroutine explicit_euler

    subroutine RK3(Un,t,dt,kms,xm,ym,Unp1)

        integer, intent(in) :: kms 
        real(PR), intent(in) :: Un(0:imax+1,0:jmax+1,1:3), xm(1:imax), ym(1:jmax), t, dt
        real(PR), intent(out) ::  Unp1(0:imax+1,0:jmax+1,1:3)
        real(PR) :: Sn(0:imax+1,0:jmax+1,1:3), phi(0:imax+1,1:jmax+1,1:3), sigma(0:imax+1,0:jmax+1)
        real(PR) :: phi1(0:imax+1,0:jmax+1,1:3), phi2(0:imax+1,1:jmax+1,1:3)
        real(PR) :: Sn1(0:imax+1,0:jmax+1,1:3), Sn2(0:imax+1,1:jmax+1,1:3)
        real(PR) :: Un1(0:imax+1,0:jmax+1,1:3), Un2(0:imax+1,1:jmax+1,1:3), Un3(0:imax+1,0:jmax+1,1:3)
        real(PR) :: Unn(0:imax+1,0:jmax+1,1:3)
        integer :: i, j

        do j = 1, jmax

            do i = 1, imax
   
                ! First intermediate solution
                call spatial_discretization(Un(i,j,1:3),Un(i-1,j,1:3),Un(i+1,j,1:3),Un(i,j+1,1:3),Un(i,j-1,1:3),phi(i,j,1:3))
                call source_term(Un(i,j,1:3),Sn(i,j,1:3))
                sigma(i,j) = compute_sigma(Un(i,j,1),xm(i),ym(j),t,kms)
                Un1(i,j,1:3) = Un(i,j,1:3) + dt * phi(i,j,1:3) + dt * sigma(i,j) * Sn(i,j,1:3)

                ! Second intermediate solution
                call spatial_discretization(Un1(i,j,1:3),Un1(i-1,j,1:3),Un1(i+1,j,1:3),Un1(i,j+1,1:3),Un1(i,j-1,1:3),phi1(i,j,1:3))
                call source_term(Un1(i,j,1:3),Sn1(i,j,1:3))
                Un2(i,j,1:3) = Un1(i,j,1:3) + dt * phi1(i,j,1:3) + dt * sigma(i,j) * Sn1(i,j,1:3)

                ! Third intermediate solution
                Unn(i,j,1:3) = (3._PR*Un(i,j,1:3) + Un2(i,j,1:3))/4._PR
                call spatial_discretization(Unn(i,j,1:3),Unn(i-1,j,1:3),Unn(i+1,j,1:3),Unn(i,j+1,1:3),Unn(i,j-1,1:3),phi2(i,j,1:3))
                call source_term(Unn(i,j,1:3),Sn2(i,j,1:3))
                Un3(i,j,1:3) = Unn(i,j,1:3) + dt * phi2(i,j,1:3) + dt * sigma(i,j) * Sn2(i,j,1:3)

                ! Final step : solution at t + dt
                Unp1(i,j,1:3) = (Un(i,j,1:3) + 2._PR*Un3(i,j,1:3))/3._PR
               
            end do
            
        end do


    
        
    end subroutine RK3
    
end module time_scheme
