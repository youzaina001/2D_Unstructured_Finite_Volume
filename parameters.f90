module parameters

    implicit none
    
    integer, parameter :: PR = 8 ! Double precision
    real(PR), parameter :: PI = 4._PR*atan(1._PR), kappa = 1._PR, theta_e(1:3) = (/ 1, 1, 1/)
    integer, parameter :: imax = 1000, jmax = 1000, nmax = 100, pas_affichage = 10
    real(PR), parameter :: Lx = 1.0, Ly = 1.0, dx = Lx/(imax+1), dy = Ly/(jmax+1), CFL = 0.9_PR
    
end module parameters