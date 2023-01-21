module time_scheme

    use parameters
    use flux
    use solver
    use initial_condition

    implicit none
    
contains

    subroutine time_step(Rho,u,v,xm,ym,PRS,sgm,k,dt)

        real(PR), dimension(0:imax+1,0:jmax+1), intent(in) :: Rho, u, v
        real(PR), dimension(1:imax,1:jmax), intent(in) :: sgm
        real(PR), dimension(0:imax+1,0:jmax+1), intent(in) :: PRS
        real(PR), intent(in) :: xm(1:imax), ym(1:jmax)
        integer, intent(in) :: k
        real(PR), intent(out) :: dt
        real(PR) :: be, dK, sigmaK
        integer :: i, j

        ! Initialisation
        be = max(maxval(abs(u)+sqrt(abs(PRS))),maxval(abs(v)+sqrt(abs(PRS))))
        dK = (dx*dy)/(2._PR*(dx+dy))

        sigmaK = minval(sgm)

        dt = CFL * dK / (4._PR * be &
           & + sigmaK*dK)

    end subroutine time_step

    subroutine explicit_euler(Un,PRS,PPRS,t,dt,kms,sigma,Unp1)

        integer, intent(in) :: kms 
        real(PR), intent(in) :: Un(0:imax+1,0:jmax+1,1:3), sigma(1:imax,1:jmax)
        real(PR), intent(in) :: PRS(0:imax+1,0:jmax+1), PPRS(0:imax+1,0:jmax+1)
        real(PR), intent(in) :: t, dt
        real(PR), intent(out) ::  Unp1(1:imax,1:jmax,1:3)
        real(PR) :: Sn(1:imax,1:jmax,1:3), phi(1:imax,1:jmax,1:3)
        integer :: i, j

        !$omp parallel do private(i, j)
        do j = 1, jmax
            
            do i = 1, imax
   
                theta_e = 1._PR!/(1._PR + sigma(i,j)*t)
                call spatial_discretization(Un(i,j,1:3),Un(i-1,j,1:3),Un(i+1,j,1:3),Un(i,j+1,1:3),Un(i,j-1,1:3), &
                & PRS(i,j),PRS(i-1,j),PRS(i+1,j),PRS(i,j+1),PRS(i,j-1),PPRS(i,j),PPRS(i-1,j), &
                & PPRS(i+1,j),PPRS(i,j+1),PPRS(i,j-1),phi(i,j,1:3))
                call source_term(Un(i,j,1:3),Sn(i,j,1:3))
                Unp1(i,j,1:3) = Un(i,j,1:3) + dt * phi(i,j,1:3) + dt * sigma(i,j) * Sn(i,j,1:3)
               
            end do
            
        end do
        
    end subroutine explicit_euler

!    subroutine RK3(Un,t,dt,kms,sigma,Unp1)
!
!        integer, intent(in) :: kms 
!        real(PR), intent(in) :: Un(0:imax+1,0:jmax+1,1:3), sigma(1:imax,1:jmax), t, dt
!        real(PR), intent(out) ::  Unp1(1:imax,1:jmax,1:3)
!        real(PR) :: Sn(0:imax+1,0:jmax+1,1:3), phi(0:imax+1,0:jmax+1,1:3)
!        real(PR) :: phi1(1:imax,1:jmax,1:3), phi2(0:imax+1,0:jmax+1,1:3)
!        real(PR) :: Sn1(0:imax+1,0:jmax+1,1:3), Sn2(0:imax+1,0:jmax+1,1:3)
!        real(PR) :: Un1(0:imax+1,0:jmax+1,1:3), Un2(0:imax+1,0:jmax+1,1:3), Un3(0:imax+1,0:jmax+1,1:3)
!        real(PR) :: Unn(0:imax+1,0:jmax+1,1:3)
!        integer :: i, j
!
!        do j = 1, jmax
!
!            do i = 1, imax
!   
!                ! First intermediate solution
!                call spatial_discretization(Un(i,j,1:3),Un(i-1,j,1:3),Un(i+1,j,1:3),Un(i,j+1,1:3),Un(i,j-1,1:3),phi(i,j,1:3))
!                call source_term(Un(i,j,1:3),Sn(i,j,1:3))
!                Un1(i,j,1:3) = Un(i,j,1:3) + dt * phi(i,j,1:3) + dt * sigma(i,j) * Sn(i,j,1:3)
!
!                ! Second intermediate solution
!                call spatial_discretization(Un1(i,j,1:3),Un1(i-1,j,1:3),Un1(i+1,j,1:3),Un1(i,j+1,1:3),Un1(i,j-1,1:3),phi1(i,j,1:3))
!                call source_term(Un1(i,j,1:3),Sn1(i,j,1:3))
!                Un2(i,j,1:3) = Un1(i,j,1:3) + dt * phi1(i,j,1:3) + dt * sigma(i,j) * Sn1(i,j,1:3)
!
!                ! Third intermediate solution
!                Unn(i,j,1:3) = (3._PR*Un(i,j,1:3) + Un2(i,j,1:3))/4._PR
!                call spatial_discretization(Unn(i,j,1:3),Unn(i-1,j,1:3),Unn(i+1,j,1:3),Unn(i,j+1,1:3),Unn(i,j-1,1:3),phi2(i,j,1:3))
!               call source_term(Unn(i,j,1:3),Sn2(i,j,1:3))
!                Un3(i,j,1:3) = Unn(i,j,1:3) + dt * phi2(i,j,1:3) + dt * sigma(i,j) * Sn2(i,j,1:3)
!
!                ! Final step : solution at t + dt
!                Unp1(i,j,1:3) = (Un(i,j,1:3) + 2._PR*Un3(i,j,1:3))/3._PR
!               
!            end do
!            
!        end do
!
!
!   
!      
!    end subroutine RK3
    
end module time_scheme
