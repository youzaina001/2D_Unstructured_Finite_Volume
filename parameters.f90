module parameters

    implicit none
    
    character(3), parameter :: case = 'one' ! one := First test case : convergence to a manufactured solution
                                            ! two := Second test case : convergence to the diffusive limit
                                            ! thr := Blast with friction
    integer, parameter :: PR = 4 ! Simple precision (4) or Double precision (8)
    real(PR), parameter :: PI = 4._PR*atan(1._PR), tf = 1._PR
    real(PR), public :: gamma = 1.4_PR
    real(PR), public :: theta_e
    integer, parameter :: imax = 100, jmax = 100, pas_affichage = 5000
    real(PR), parameter :: Lx = 1._PR, Ly = 1._PR, dx = Lx/imax, dy = Ly/jmax, CFL = 0.6_PR
    
end module parameters