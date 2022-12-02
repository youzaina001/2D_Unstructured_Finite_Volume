module Numerical_flux

    use read_parameters

    implicit none
    
contains

!***** Fonction permettant l'obtention des composantes de *****!
!***************** la première colonne du flux ****************!
    function flux_u(rho,u,v,p) result(Fu)

                        !** Variables **!

        !** rho : densité volumique de masse
        !** u : Composante de la vitesse dans la direction x 
        !** v : Composante de la vitesse dans la direction y
        !** p : La pression

        real(PR), intent(in) :: rho, u, v, p
        real(PR), dimension(1:3) :: Fu

        Fu(1) = rho*u
        Fu(2) = rho*u*u + p
        Fu(3) = rho*u*v

    end function flux_u

!***** Fonction permettant l'obtention du maximum des *****!
!************************ célérités ***********************!
    function max_ce(ug,ud,vg,vd,p_derive) result(be)

                        !** Variables **!

        !** u : Composante de la vitesse dans la direction x 
        !** v : Composante de la vitesse dans la direction y
        !** p_derive :: La dérivée de p par rapport à rho
        !** be : le max des valeurs absolues des célérités
        !** g : gauche
        !** d : droite

        real(PR), intent(in) :: ug, ud, vg, vd, p_derive
        real(PR) :: c_ug, c_ud, c_vg, c_vd

        ! On va pas prendre en compte les vitesses
        ! V.n - c et V.n on va prendre le max
        c_ug = ug + sqrt(p_derive)
        c_ud = ud + sqrt(p_derive)
        c_vg = vg + sqrt(p_derive)
        c_vd = vd + sqrt(p_derive)

        be = max(c_ug,c_ud,c_vg,c_vd)

    end function max_ce

!***** Fonction permettant l'obtention des composantes de *****!
!***************** la deuxième colonne du flux ****************!
    function flux_v(rho,u,v,p) result(Fv)

                        !** Variables **!

        !** rho : densité volumique de masse
        !** u : Composante de la vitesse dans la direction x 
        !** v : Composante de la vitesse dans la direction y
        !** p : La pression

        real(PR), intent(in) :: rho, u, v, p
        real(PR), dimension(1:3) :: Fv

        Fv(1) = rho*v
        Fv(2) = rho*u*v
        Fv(3) = rho*v*v + p

    end function flux_v

!***** Fonction permettant l'obtention des composantes de *****!
!************ la première colonne du flux de Rusanov **********!
    function Rusanov_u(rhog,rhod,ug,ud,vg,vd,pg,pd) result(Rus_u)

                        !** Variables **!

        !** rho : densité volumique de masse
        !** u : Composante de la vitesse dans la direction x 
        !** v : Composante de la vitesse dans la direction y
        !** p : La pression
        !** be : le max des valeurs absolues des célérités
        !** g : gauche
        !** d : droite

        real(PR), intent(in) :: rhog, rhod, ug, ud, vg, vd, pg, pd
        real(PR), dimension(1:3) :: Fu_g, Fu_d, Rus_u
        real(PR) :: be

        Fu_g = flux_u(rhog,ug,vg,pg)
        Fu_d = flux_u(rhod,ud,vd,pd)

        be = max_ce(ug,ud,vg,vd,p_derive)

        Rus_u(1) = 0.5_PR*(Fu_g(1) + Fu_d(1)) - 0.5_PR*be*(Fu_d(1) - Fu_g(1))
        Rus_u(2) = 0.5_PR*(Fu_g(2) + Fu_d(2)) - 0.5_PR*be*(Fu_d(2) - Fu_g(2))
        Rus_u(3) = 0.5_PR*(Fu_g(3) + Fu_d(3)) - 0.5_PR*be*(Fu_d(3) - Fu_g(3))

    end function Rusanov_u

!***** Fonction permettant l'obtention des composantes de *****!
!************ la deuxième colonne du flux de Rusanov **********!
    function Rusanov_v(rhog,rhod,ug,ud,vg,vd,pg,pd) result(Rus_v)

                        !** Variables **!

        !** rho : densité volumique de masse
        !** u : Composante de la vitesse dans la direction x 
        !** v : Composante de la vitesse dans la direction y
        !** p : La pression
        !** be : le max des valeurs absolues des célérités
        !** g : gauche
        !** d : droite

        real(PR), intent(in) :: rhog, rhod, ug, ud, vg, vd, pg, pd
        real(PR), dimension(1:3) :: Fv_g, Fv_d, Rus_v
        real(PR) :: be

        Fv_g = flux_v(rhog,ug,vg,pg)
        Fv_d = flux_v(rhod,ud,vd,pd)

        be = max_ce(ug,ud,vg,vd,p_derive)

        Rus_v(1) = 0.5_PR*(Fv_g(1) + Fv_d(1)) - 0.5_PR*be*(Fv_d(1) - Fv_g(1))
        Rus_v(2) = 0.5_PR*(Fv_g(2) + Fv_d(2)) - 0.5_PR*be*(Fv_d(2) - Fv_g(2))
        Rus_v(3) = 0.5_PR*(Fv_g(3) + Fv_d(3)) - 0.5_PR*be*(Fv_d(3) - Fv_g(3))

    end function Rusanov_v

    
end module Numerical_flux
