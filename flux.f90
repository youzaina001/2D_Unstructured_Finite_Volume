module flux

    use parameters

    implicit none
    
contains

    subroutine conservative_to_non_conservative(UU,Rho,u,v)

        real(PR), intent(in) :: UU(1:3)
        real(PR), intent(out) :: Rho, u, v

        Rho = UU(1)
        if (UU(1) /= 0) then

            u = UU(2)/UU(1)
            v = UU(3)/UU(1)

        else

            u = 0._PR
            v = 0._PR
            
        end if

    end subroutine conservative_to_non_conservative

    subroutine non_conservative_to_conservative(Rho,u,v,UU)

        real(PR), intent(in) :: Rho, u, v
        real(PR), intent(out) :: UU(1:3)

        UU = (/Rho,Rho*u,Rho*v/)

    end subroutine non_conservative_to_conservative

    subroutine swapping(Rhop1,up1,vp1,Rho,u,v)

        real(PR), intent(in) :: Rhop1, up1, vp1
        real(PR), intent(out) :: Rho, u, v

        Rho = Rhop1
        u = up1
        v = vp1

    end subroutine swapping

    !function compute_kappa(Rho) result(kappa)

    !    real(PR), intent(in) :: Rho(1:imax+1,1:jmax+1)
    !    real(PR) :: kappa(1:imax+1,1:jmax+1)

    !    kappa = 0.04_PR * matmul(matmul(Rho,Rho),Rho)
        
    !end function compute_kappa

    function compute_sigma(x,y,t,k) result(sigma)

        integer, intent(in) :: k
        real(PR), intent(in) :: x, y, t
        real(PR) :: sigma

        if (case == 'one') then

            sigma = 1._PR + k * (-1._PR)**gamma * (gamma - 1) * exp(-(x+y))**(gamma - 1)/exp(-t) &
              & -2._PR * exp(-t) / (k * exp(-(x+y)))

        else if (case == 'two') then

            sigma = 100._PR
            
        end if
        
    end function compute_sigma

    function pressure(Rho) result(P)

        real(PR), intent(in) :: Rho
        real(PR) :: P

        if (case == 'one' .or. case == 'two') then

            P = Rho**gamma
            
        end if
        
    end function pressure

    function pressure_prime(Rho) result(Pd)

        real(PR), intent(in) :: Rho
        real(PR) :: Pd

        if (case == 'one' .or. case == 'two') then

            Pd = (gamma - 1)*Rho**(gamma - 1)
            
        end if
        
    end function pressure_prime

    function Exact_flux(U) result(FU)

        real(PR), intent(in) :: U(1:3)
        real(PR) :: FU(1:3,1:2)
        real(PR) :: Rho, ux, vy, P

        call conservative_to_non_conservative(U,Rho,ux,vy)

        P = pressure(Rho)

        FU(1,1) = Rho * ux
        FU(1,2) = Rho * vy
        FU(2,1) = Rho * ux**2 + P 
        FU(2,2) = Rho * ux * vy
        FU(3,1) = Rho * ux * vy
        FU(3,2) = Rho * vy**2 + P
        
    end function Exact_flux

    function max_celerity(Uk,Ul,Pd) result(c)

        real(PR), intent(in) :: Uk(1:3), Ul(1:3)
        real(PR), intent(in) :: Pd
        real(PR) :: Rhok, Rhol, uxk, uxl, vyk, vyl
        real(PR) :: cuk, cul, cvk, cvl
        real(PR) :: c

        call conservative_to_non_conservative(Uk,Rhok,uxk,vyk)
        call conservative_to_non_conservative(Ul,Rhol,uxl,vyl)

        cuk = uxk + sqrt(Pd)
        cul = uxl + sqrt(Pd)
        cvk = vyk + sqrt(Pd)
        cvl = vyl + sqrt(Pd)

        c = max(cuk,cul,cvk,cvl)
        
    end function max_celerity

    function Rusanov(Uk,Ul,normal) result(Fe)

        real(PR), intent(in) :: Uk(1:3), Ul(1:3)
        real(PR), intent(in) :: normal(1:2)
        real(PR) :: FUk(1:3,1:2), FUl(1:3,1:2)
        real(PR) :: Fe(1:3)
        real(PR) :: Rho1, u1, v1, Rho2, u2, v2
        real(PR) :: Pd
        real(PR) :: be

        FUk = Exact_flux(Uk)
        Ful = Exact_flux(Ul)

        call conservative_to_non_conservative(Uk,Rho1,u1,v1)
        call conservative_to_non_conservative(Ul,Rho2,u2,v2)

        Pd = pressure_prime(max(Rho1,Rho2))

        be = max_celerity(Uk,Ul,Pd)

        Fe = 0.5_PR * (matmul(FUk,normal) + matmul(FUl,normal)) - 0.5_PR * be * theta_e * (Ul - Uk)

    end function Rusanov

    
end module flux
