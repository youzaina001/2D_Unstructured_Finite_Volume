module parameters

    implicit none
    
    character(3), parameter :: case = 'one' ! one := First test case : convergence to a manufactured solution
                                            ! two := Second test case : convergence to the diffusive limit
    integer, parameter :: PR = 4 ! Simple precision (4) or Double precision (8)
    real(PR), parameter :: PI = 4._PR*atan(1._PR), tf = 1._PR
    real(PR), public :: theta_e(1:3)
    integer, parameter :: imax = 300, jmax = 300, pas_affichage = 100, gamma = 2
    real(PR), parameter :: Lx = 1.0, Ly = 1.0, dx = Lx/(imax+1), dy = Ly/(jmax+1), CFL = 0.5_PR
    
end module parameters