module flux

    use parameters

    implicit none
    
contains

    subroutine conservative_to_non_conservative(UU,Rho,u,v)

        real(PR), intent(in) :: UU(1:3)
        real(PR), intent(out) :: Rho, u, v

        Rho = UU(1)
        if (UU(1) > 0._PR) then

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

    function compute_sigma(Rho,x,y,t,k) result(sigma)

        integer, intent(in) :: k
        real(PR), intent(in) :: Rho
        real(PR), intent(in) :: x, y, t
        real(PR) :: sigma

        if (case == 'one') then

            sigma = -k*exp(t)*((-exp(-t)/k) - gamma*exp(-gamma*(x+y)) - (2*exp(-2*t)/(k**2*exp(-(x+y)))))

        else if (case == 'two') then

            sigma = 100._PR

        else if (case == 'thr') then

            sigma = (Rho**3)/25._PR
            
        end if
        
    end function compute_sigma

    function pressure(Rho) result(P)

        real(PR), intent(in) :: Rho
        real(PR) :: P

        P = Rho**gamma
        
    end function pressure

    function pressure_prime(Rho) result(PP)

        real(PR), intent(in) :: Rho
        real(PR) :: PP

        PP = (gamma - 1)*Rho**(gamma - 1)
        
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

    function max_celerity(Uk,Ul,PPk,PPl) result(c)

        real(PR), intent(in) :: Uk(1:3), Ul(1:3)
        real(PR), intent(in) :: PPk, PPl
        real(PR) :: Rhok, Rhol, uxk, uxl, vyk, vyl
        real(PR) :: cuk, cul, cvk, cvl
        real(PR) :: c

        call conservative_to_non_conservative(Uk,Rhok,uxk,vyk)
        call conservative_to_non_conservative(Ul,Rhol,uxl,vyl)

        cuk = uxk + sqrt(PPk)
        cul = uxl + sqrt(PPl)
        cvk = vyk + sqrt(PPk)
        cvl = vyl + sqrt(PPl)

        c = max(cuk,cul,cvk,cvl)
        
    end function max_celerity

    function Rusanov(Uk,Ul,nKe) result(Fe) ! nKe : normale

        real(PR), intent(in) :: Uk(1:3), Ul(1:3)
        real(PR), intent(in) :: nKe(1:2)
        real(PR) :: FUk(1:3,1:2), FUl(1:3,1:2)
        real(PR) :: Fe(1:3)
        real(PR) :: Rho1, u1, v1, Rho2, u2, v2
        real(PR) :: PPk, PPl
        real(PR) :: be

        FUk = Exact_flux(Uk)
        Ful = Exact_flux(Ul)
        PPk = pressure_prime(Uk(1))
        PPl = pressure_prime(Ul(1))

        call conservative_to_non_conservative(Uk,Rho1,u1,v1)
        call conservative_to_non_conservative(Ul,Rho2,u2,v2)

        be = max_celerity(Uk,Ul,PPk,PPl)

        Fe(1) = 0.5_PR * (FUk(1,1)*nKe(1) + FUk(1,2)*nKe(2) + FUl(1,1)*nKe(1) + FUl(1,2)*nKe(2)) &
              & - 0.5_PR * be * theta_e * (Ul(1) - Uk(1))

        Fe(2) = 0.5_PR * (FUk(2,1)*nKe(1) + FUk(2,2)*nKe(2) + FUl(2,1)*nKe(1) + FUl(2,2)*nKe(2)) &
              & - 0.5_PR * be * (Ul(2) - Uk(2))

        Fe(3) = 0.5_PR * (FUk(3,1)*nKe(1) + FUk(3,2)*nKe(2) + FUl(3,1)*nKe(1) + FUl(3,2)*nKe(2)) &
              & - 0.5_PR * be * (Ul(3) - Uk(3))

    end function Rusanov

    
end module flux
