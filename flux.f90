module flux

    use parameters

    implicit none
    
contains

    subroutine conservative_to_non_conservative(UU,Rho,u,v)

        real(PR), intent(in) :: UU(1:3)
        real(PR), intent(out) :: Rho, u, v

        Rho = UU(1)
        u = UU(2)/UU(1)
        v = UU(3)/UU(1)

    end subroutine conservative_to_non_conservative

    function pressure(Rho) result(P)

        real(PR), intent(in) :: Rho
        real(PR) :: P
    
        P = kappa * Rho
        
    end function pressure

    function pressure_prime() result(Pd)

        real(PR) :: Pd
    
        Pd = kappa
        
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
        real(PR) :: Pd
        real(PR) :: be

        FUk = Exact_flux(Uk)
        Ful = Exact_flux(Ul)

        Pd = pressure_prime()

        be = max_celerity(Uk,Ul,Pd)

        Fe = 0.5_PR * (matmul(FUk,normal) + matmul(FUl,normal)) - 0.5_PR * be * theta_e * (Ul - Uk)

    end function Rusanov

    
end module flux
